/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clement Jamin
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - 2019/08 Vincent Rouvreau: Fix issue #10 for CGAL and Eigen3
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef TANGENTIAL_COMPLEX_H_
#define TANGENTIAL_COMPLEX_H_

#include <gudhi/Tangential_complex/config.h>
#include <gudhi/Tangential_complex/Simplicial_complex.h>
#include <gudhi/Tangential_complex/utilities.h>
#include <gudhi/Kd_tree_search.h>
#include <gudhi/console_color.h>
#include <gudhi/Clock.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Debug_utils.h>

#include <CGAL/Default.h>
#include <CGAL/Dimension.h>
#include <CGAL/function_objects.h>  // for CGAL::Identity
#include <CGAL/Epick_d.h>
#include <CGAL/Regular_triangulation_traits_adapter.h>
#include <CGAL/Regular_triangulation.h>
#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Combination_enumerator.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/version.h>  // for CGAL_VERSION_NR

#include <Eigen/Core>
#include <Eigen/Eigen>
#include <Eigen/src/Core/util/Macros.h>  // for EIGEN_VERSION_AT_LEAST

#include <boost/optional.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/counting_range.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/container/flat_set.hpp>

#include <tuple>
#include <vector>
#include <set>
#include <utility>
#include <sstream>
#include <iostream>
#include <limits>
#include <algorithm>
#include <functional>
#include <iterator>
#include <cmath>  // for std::sqrt
#include <string>
#include <cstddef>  // for std::size_t

#ifdef GUDHI_USE_TBB
#include <tbb/parallel_for.h>
#include <tbb/combinable.h>
#include <tbb/mutex.h>
#endif

// #define GUDHI_TC_EXPORT_NORMALS // Only for 3D surfaces (k=2, d=3)

// Make compilation fail - required for external projects - https://github.com/GUDHI/gudhi-devel/issues/10
#if CGAL_VERSION_NR < 1041101000
# error Alpha_complex_3d is only available for CGAL >= 4.11
#endif

#if !EIGEN_VERSION_AT_LEAST(3,1,0)
# error Alpha_complex_3d is only available for Eigen3 >= 3.1.0 installed with CGAL
#endif

namespace sps = Gudhi::spatial_searching;

namespace Gudhi {

namespace tangential_complex {

using namespace internal;

class Vertex_data {
 public:
  Vertex_data(std::size_t data = (std::numeric_limits<std::size_t>::max)()) : m_data(data) {}

  operator std::size_t() { return m_data; }

  operator std::size_t() const { return m_data; }

 private:
  std::size_t m_data;
};

/**
 * \class Tangential_complex Tangential_complex.h gudhi/Tangential_complex.h
 * \brief Tangential complex data structure.
 *
 * \ingroup tangential_complex
 *
 * \details
 *  The class Tangential_complex represents a tangential complex.
 *  After the computation of the complex, an optional post-processing called perturbation can
 *  be run to attempt to remove inconsistencies.
 *
 * \tparam Kernel_ requires a <a target="_blank"
 * href="http://doc.cgal.org/latest/Kernel_d/classCGAL_1_1Epick__d.html">CGAL::Epick_d</a> class, which
 * can be static if you know the ambiant dimension at compile-time, or dynamic if you don't.
 * \tparam DimensionTag can be either <a target="_blank"
 * href="http://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Dimension__tag.html">Dimension_tag<d></a>
 * if you know the intrinsic dimension at compile-time,
 * or <a target="_blank"
 * href="http://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Dynamic__dimension__tag.html">CGAL::Dynamic_dimension_tag</a>
 * if you don't.
 * \tparam Concurrency_tag enables sequential versus parallel computation. Possible values are `CGAL::Parallel_tag` (the
 * default) and `CGAL::Sequential_tag`. \tparam Triangulation_ is the type used for storing the local regular
 * triangulations. We highly recommend to use the default value (`CGAL::Regular_triangulation`).
 *
 */
template <typename Kernel_,       // ambiant kernel
          typename DimensionTag,  // intrinsic dimension
          typename Concurrency_tag = CGAL::Parallel_tag, typename Triangulation_ = CGAL::Default>
class Tangential_complex {
  typedef Kernel_ K;
  typedef typename K::FT FT;
  typedef typename K::Point_d Point;
  typedef typename K::Weighted_point_d Weighted_point;
  typedef typename K::Vector_d Vector;

  typedef typename CGAL::Default::Get<
      Triangulation_,
      CGAL::Regular_triangulation<
          CGAL::Epick_d<DimensionTag>,
          CGAL::Triangulation_data_structure<
              typename CGAL::Epick_d<DimensionTag>::Dimension,
              CGAL::Triangulation_vertex<CGAL::Regular_triangulation_traits_adapter<CGAL::Epick_d<DimensionTag> >,
                                         Vertex_data>,
              CGAL::Triangulation_full_cell<
                  CGAL::Regular_triangulation_traits_adapter<CGAL::Epick_d<DimensionTag> > > > > >::type Triangulation;
  typedef typename Triangulation::Geom_traits Tr_traits;
  typedef typename Triangulation::Weighted_point Tr_point;
  typedef typename Tr_traits::Base::Point_d Tr_bare_point;
  typedef typename Triangulation::Vertex_handle Tr_vertex_handle;
  typedef typename Triangulation::Full_cell_handle Tr_full_cell_handle;
  typedef typename Tr_traits::Vector_d Tr_vector;

#if defined(GUDHI_USE_TBB)
  typedef tbb::mutex Mutex_for_perturb;
  typedef Vector Translation_for_perturb;
  typedef std::vector<Atomic_wrapper<FT> > Weights;
#else
  typedef Vector Translation_for_perturb;
  typedef std::vector<FT> Weights;
#endif
  typedef std::vector<Translation_for_perturb> Translations_for_perturb;

  // Store a local triangulation and a handle to its center vertex

  struct Tr_and_VH {
   public:
    Tr_and_VH() : m_tr(NULL) {}

    Tr_and_VH(int dim) : m_tr(new Triangulation(dim)) {}

    ~Tr_and_VH() { destroy_triangulation(); }

    Triangulation &construct_triangulation(int dim) {
      delete m_tr;
      m_tr = new Triangulation(dim);
      return tr();
    }

    void destroy_triangulation() {
      delete m_tr;
      m_tr = NULL;
    }

    Triangulation &tr() { return *m_tr; }

    Triangulation const &tr() const { return *m_tr; }

    Tr_vertex_handle const &center_vertex() const { return m_center_vertex; }

    Tr_vertex_handle &center_vertex() { return m_center_vertex; }

   private:
    Triangulation *m_tr;
    Tr_vertex_handle m_center_vertex;
  };

 public:
  typedef Basis<K> Tangent_space_basis;
  typedef Basis<K> Orthogonal_space_basis;
  typedef std::vector<Tangent_space_basis> TS_container;
  typedef std::vector<Orthogonal_space_basis> OS_container;

  typedef std::vector<Point> Points;

  typedef boost::container::flat_set<std::size_t> Simplex;
  typedef std::set<Simplex> Simplex_set;

 private:
  typedef sps::Kd_tree_search<K, Points> Points_ds;
  typedef typename Points_ds::KNS_range KNS_range;
  typedef typename Points_ds::INS_range INS_range;

  typedef std::vector<Tr_and_VH> Tr_container;
  typedef std::vector<Vector> Vectors;

  // An Incident_simplex is the list of the vertex indices
  // except the center vertex
  typedef boost::container::flat_set<std::size_t> Incident_simplex;
  typedef std::vector<Incident_simplex> Star;
  typedef std::vector<Star> Stars_container;

  // For transform_iterator

  static const Tr_point &vertex_handle_to_point(Tr_vertex_handle vh) { return vh->point(); }

  template <typename P, typename VH>
  static const P &vertex_handle_to_point(VH vh) {
    return vh->point();
  }

 public:
  typedef internal::Simplicial_complex Simplicial_complex;

  /** \brief Constructor from a range of points.
   * Points are copied into the instance, and a search data structure is initialized.
   * Note the complex is not computed: `compute_tangential_complex` must be called after the creation
   * of the object.
   *
   * @param[in] points Range of points (`Point_range::value_type` must be the same as `Kernel_::Point_d`).
   * @param[in] intrinsic_dimension Intrinsic dimension of the manifold.
   * @param[in] k Kernel instance.
   */
  template <typename Point_range>
  Tangential_complex(Point_range points, int intrinsic_dimension,
#ifdef GUDHI_TC_USE_ANOTHER_POINT_SET_FOR_TANGENT_SPACE_ESTIM
                     InputIterator first_for_tse, InputIterator last_for_tse,
#endif
                     const K &k = K())
      : m_k(k),
        m_intrinsic_dim(intrinsic_dimension),
        m_ambient_dim(points.empty() ? 0 : k.point_dimension_d_object()(*points.begin())),
        m_points(points.begin(), points.end()),
        m_weights(m_points.size(), FT(0))
#if defined(GUDHI_USE_TBB) && defined(GUDHI_TC_PERTURB_POSITION)
        ,
        m_p_perturb_mutexes(NULL)
#endif
        ,
        m_points_ds(m_points),
        m_last_max_perturb(0.),
        m_are_tangent_spaces_computed(m_points.size(), false),
        m_tangent_spaces(m_points.size(), Tangent_space_basis())
#ifdef GUDHI_TC_EXPORT_NORMALS
        ,
        m_orth_spaces(m_points.size(), Orthogonal_space_basis())
#endif
#ifdef GUDHI_TC_USE_ANOTHER_POINT_SET_FOR_TANGENT_SPACE_ESTIM
        ,
        m_points_for_tse(first_for_tse, last_for_tse),
        m_points_ds_for_tse(m_points_for_tse)
#endif
  {
  }

  /// Destructor
  ~Tangential_complex() {
#if defined(GUDHI_USE_TBB) && defined(GUDHI_TC_PERTURB_POSITION)
    delete[] m_p_perturb_mutexes;
#endif
  }

  /// Returns the intrinsic dimension of the manifold.
  int intrinsic_dimension() const { return m_intrinsic_dim; }

  /// Returns the ambient dimension.
  int ambient_dimension() const { return m_ambient_dim; }

  Points const &points() const { return m_points; }

  /** \brief Returns the point corresponding to the vertex given as parameter.
   *
   * @param[in] vertex Vertex handle of the point to retrieve.
   * @return The point found.
   */
  Point get_point(std::size_t vertex) const { return m_points[vertex]; }

  /** \brief Returns the perturbed position of the point corresponding to the vertex given as parameter.
   *
   * @param[in] vertex Vertex handle of the point to retrieve.
   * @return The perturbed position of the point found.
   */
  Point get_perturbed_point(std::size_t vertex) const { return compute_perturbed_point(vertex); }

  /// Returns the number of vertices.

  std::size_t number_of_vertices() const { return m_points.size(); }

  void set_weights(const Weights &weights) { m_weights = weights; }

  void set_tangent_planes(const TS_container &tangent_spaces
#ifdef GUDHI_TC_EXPORT_NORMALS
                          ,
                          const OS_container &orthogonal_spaces
#endif
  ) {
#ifdef GUDHI_TC_EXPORT_NORMALS
    GUDHI_CHECK(m_points.size() == tangent_spaces.size() && m_points.size() == orthogonal_spaces.size(),
                std::logic_error("Wrong sizes"));
#else
    GUDHI_CHECK(m_points.size() == tangent_spaces.size(), std::logic_error("Wrong sizes"));
#endif
    m_tangent_spaces = tangent_spaces;
#ifdef GUDHI_TC_EXPORT_NORMALS
    m_orth_spaces = orthogonal_spaces;
#endif
    for (std::size_t i = 0; i < m_points.size(); ++i) m_are_tangent_spaces_computed[i] = true;
  }

  /** \brief Computes the tangential complex.
   *  \exception std::invalid_argument In debug mode, if the computed star dimension is too low. Try to set a bigger
   *  maximal edge length value with `Tangential_complex::set_max_squared_edge_length` if
   *  this happens.
   */
  void compute_tangential_complex() {
#ifdef GUDHI_TC_PERFORM_EXTRA_CHECKS
    std::cerr << red << "WARNING: GUDHI_TC_PERFORM_EXTRA_CHECKS is defined. "
              << "Computation might be slower than usual.\n"
              << white;
#endif

#if defined(GUDHI_TC_PROFILING) && defined(GUDHI_USE_TBB)
    Gudhi::Clock t;
#endif

    // We need to do that because we don't want the container to copy the
    // already-computed triangulations (while resizing) since it would
    // invalidate the vertex handles stored beside the triangulations
    m_triangulations.resize(m_points.size());
    m_stars.resize(m_points.size());
    m_squared_star_spheres_radii_incl_margin.resize(m_points.size(), FT(-1));
#ifdef GUDHI_TC_PERTURB_POSITION
    if (m_points.empty())
      m_translations.clear();
    else
      m_translations.resize(m_points.size(), m_k.construct_vector_d_object()(m_ambient_dim));
#if defined(GUDHI_USE_TBB)
    delete[] m_p_perturb_mutexes;
    m_p_perturb_mutexes = new Mutex_for_perturb[m_points.size()];
#endif
#endif

#ifdef GUDHI_USE_TBB
    // Parallel
    if (boost::is_convertible<Concurrency_tag, CGAL::Parallel_tag>::value) {
      tbb::parallel_for(tbb::blocked_range<size_t>(0, m_points.size()), Compute_tangent_triangulation(*this));
    } else {
#endif  // GUDHI_USE_TBB
        // Sequential
      for (std::size_t i = 0; i < m_points.size(); ++i) compute_tangent_triangulation(i);
#ifdef GUDHI_USE_TBB
    }
#endif  // GUDHI_USE_TBB

#if defined(GUDHI_TC_PROFILING) && defined(GUDHI_USE_TBB)
    t.end();
    std::cerr << "Tangential complex computed in " << t.num_seconds() << " seconds.\n";
#endif
  }

  /// \brief Type returned by `Tangential_complex::fix_inconsistencies_using_perturbation`.
  struct Fix_inconsistencies_info {
    /// `true` if all inconsistencies could be removed, `false` if the time limit has been reached before
    bool success = false;
    /// number of steps performed
    unsigned int num_steps = 0;
    /// initial number of inconsistent stars
    std::size_t initial_num_inconsistent_stars = 0;
    /// best number of inconsistent stars during the process
    std::size_t best_num_inconsistent_stars = 0;
    /// final number of inconsistent stars
    std::size_t final_num_inconsistent_stars = 0;
  };

  /** \brief Attempts to fix inconsistencies by perturbing the point positions.
   *
   * @param[in] max_perturb Maximum length of the translations used by the perturbation.
   * @param[in] time_limit Time limit in seconds. If -1, no time limit is set.
   */
  Fix_inconsistencies_info fix_inconsistencies_using_perturbation(double max_perturb, double time_limit = -1.) {
    Fix_inconsistencies_info info;

    if (time_limit == 0.) return info;

    Gudhi::Clock t;

#ifdef GUDHI_TC_SHOW_DETAILED_STATS_FOR_INCONSISTENCIES
    std::tuple<std::size_t, std::size_t, std::size_t> stats_before = number_of_inconsistent_simplices(false);

    if (std::get<1>(stats_before) == 0) {
#ifdef DEBUG_TRACES
      std::cerr << "Nothing to fix.\n";
#endif
      info.success = false;
      return info;
    }
#endif  // GUDHI_TC_SHOW_DETAILED_STATS_FOR_INCONSISTENCIES

    m_last_max_perturb = max_perturb;

    bool done = false;
    info.best_num_inconsistent_stars = m_triangulations.size();
    info.num_steps = 0;
    while (!done) {
#ifdef GUDHI_TC_SHOW_DETAILED_STATS_FOR_INCONSISTENCIES
      std::cerr << "\nBefore fix step:\n"
                << "  * Total number of simplices in stars (incl. duplicates): " << std::get<0>(stats_before) << "\n"
                << "  * Num inconsistent simplices in stars (incl. duplicates): " << red << std::get<1>(stats_before)
                << white << " (" << 100. * std::get<1>(stats_before) / std::get<0>(stats_before) << "%)\n"
                << "  * Number of stars containing inconsistent simplices: " << red << std::get<2>(stats_before)
                << white << " (" << 100. * std::get<2>(stats_before) / m_points.size() << "%)\n";
#endif

#if defined(DEBUG_TRACES) || defined(GUDHI_TC_PROFILING)
      std::cerr << yellow << "\nAttempt to fix inconsistencies using perturbations - step #" << info.num_steps + 1
                << "... " << white;
#endif

      std::size_t num_inconsistent_stars = 0;
      std::vector<std::size_t> updated_points;

#ifdef GUDHI_TC_PROFILING
      Gudhi::Clock t_fix_step;
#endif

      // Parallel
#if defined(GUDHI_USE_TBB)
      if (boost::is_convertible<Concurrency_tag, CGAL::Parallel_tag>::value) {
        tbb::combinable<std::size_t> num_inconsistencies;
        tbb::combinable<std::vector<std::size_t> > tls_updated_points;
        tbb::parallel_for(tbb::blocked_range<size_t>(0, m_triangulations.size()),
                          Try_to_solve_inconsistencies_in_a_local_triangulation(*this, max_perturb, num_inconsistencies,
                                                                                tls_updated_points));
        num_inconsistent_stars = num_inconsistencies.combine(std::plus<std::size_t>());
        updated_points =
            tls_updated_points.combine([](std::vector<std::size_t> const &x, std::vector<std::size_t> const &y) {
              std::vector<std::size_t> res;
              res.reserve(x.size() + y.size());
              res.insert(res.end(), x.begin(), x.end());
              res.insert(res.end(), y.begin(), y.end());
              return res;
            });
      } else {
#endif  // GUDHI_USE_TBB
        // Sequential
        for (std::size_t i = 0; i < m_triangulations.size(); ++i) {
          num_inconsistent_stars +=
              try_to_solve_inconsistencies_in_a_local_triangulation(i, max_perturb, std::back_inserter(updated_points));
        }
#if defined(GUDHI_USE_TBB)
      }
#endif  // GUDHI_USE_TBB

#ifdef GUDHI_TC_PROFILING
      t_fix_step.end();
#endif

#if defined(GUDHI_TC_SHOW_DETAILED_STATS_FOR_INCONSISTENCIES) || defined(DEBUG_TRACES)
      std::cerr << "\nEncountered during fix:\n"
                << "  * Num stars containing inconsistent simplices: " << red << num_inconsistent_stars << white << " ("
                << 100. * num_inconsistent_stars / m_points.size() << "%)\n";
#endif

#ifdef GUDHI_TC_PROFILING
      std::cerr << yellow << "done in " << t_fix_step.num_seconds() << " seconds.\n" << white;
#elif defined(DEBUG_TRACES)
      std::cerr << yellow << "done.\n" << white;
#endif

      if (num_inconsistent_stars > 0) refresh_tangential_complex(updated_points);

#ifdef GUDHI_TC_PERFORM_EXTRA_CHECKS
      // Confirm that all stars were actually refreshed
      std::size_t num_inc_1 = std::get<1>(number_of_inconsistent_simplices(false));
      refresh_tangential_complex();
      std::size_t num_inc_2 = std::get<1>(number_of_inconsistent_simplices(false));
      if (num_inc_1 != num_inc_2)
        std::cerr << red << "REFRESHMENT CHECK: FAILED. (" << num_inc_1 << " vs " << num_inc_2 << ")\n" << white;
      else
        std::cerr << green << "REFRESHMENT CHECK: PASSED.\n" << white;
#endif

#ifdef GUDHI_TC_SHOW_DETAILED_STATS_FOR_INCONSISTENCIES
      std::tuple<std::size_t, std::size_t, std::size_t> stats_after = number_of_inconsistent_simplices(false);

      std::cerr << "\nAfter fix:\n"
                << "  * Total number of simplices in stars (incl. duplicates): " << std::get<0>(stats_after) << "\n"
                << "  * Num inconsistent simplices in stars (incl. duplicates): " << red << std::get<1>(stats_after)
                << white << " (" << 100. * std::get<1>(stats_after) / std::get<0>(stats_after) << "%)\n"
                << "  * Number of stars containing inconsistent simplices: " << red << std::get<2>(stats_after) << white
                << " (" << 100. * std::get<2>(stats_after) / m_points.size() << "%)\n";

      stats_before = stats_after;
#endif

      if (info.num_steps == 0) info.initial_num_inconsistent_stars = num_inconsistent_stars;

      if (num_inconsistent_stars < info.best_num_inconsistent_stars)
        info.best_num_inconsistent_stars = num_inconsistent_stars;

      info.final_num_inconsistent_stars = num_inconsistent_stars;

      done = (num_inconsistent_stars == 0);
      if (!done) {
        ++info.num_steps;
        if (time_limit > 0. && t.num_seconds() > time_limit) {
#ifdef DEBUG_TRACES
          std::cerr << red << "Time limit reached.\n" << white;
#endif
          info.success = false;
          return info;
        }
      }
    }

#ifdef DEBUG_TRACES
    std::cerr << green << "Fixed!\n" << white;
#endif
    info.success = true;
    return info;
  }

  /// \brief Type returned by `Tangential_complex::number_of_inconsistent_simplices`.
  struct Num_inconsistencies {
    /// Total number of simplices in stars (including duplicates that appear in several stars)
    std::size_t num_simplices = 0;
    /// Number of inconsistent simplices
    std::size_t num_inconsistent_simplices = 0;
    /// Number of stars containing at least one inconsistent simplex
    std::size_t num_inconsistent_stars = 0;
  };

  /// Returns the number of inconsistencies
  /// @param[in] verbose If true, outputs a message into `std::cerr`.

  Num_inconsistencies number_of_inconsistent_simplices(
#ifdef DEBUG_TRACES
      bool verbose = true
#else
      bool verbose = false
#endif
      ) const {
    Num_inconsistencies stats;

    // For each triangulation
    for (std::size_t idx = 0; idx < m_points.size(); ++idx) {
      bool is_star_inconsistent = false;

      // For each cell
      Star::const_iterator it_inc_simplex = m_stars[idx].begin();
      Star::const_iterator it_inc_simplex_end = m_stars[idx].end();
      for (; it_inc_simplex != it_inc_simplex_end; ++it_inc_simplex) {
        // Don't check infinite cells
        if (is_infinite(*it_inc_simplex)) continue;

        Simplex c = *it_inc_simplex;
        c.insert(idx);  // Add the missing index

        if (!is_simplex_consistent(c)) {
          ++stats.num_inconsistent_simplices;
          is_star_inconsistent = true;
        }

        ++stats.num_simplices;
      }
      stats.num_inconsistent_stars += is_star_inconsistent;
    }

    if (verbose) {
      std::cerr << "\n==========================================================\n"
                << "Inconsistencies:\n"
                << "  * Total number of simplices in stars (incl. duplicates): " << stats.num_simplices << "\n"
                << "  * Number of inconsistent simplices in stars (incl. duplicates): "
                << stats.num_inconsistent_simplices << " ("
                << 100. * stats.num_inconsistent_simplices / stats.num_simplices << "%)\n"
                << "  * Number of stars containing inconsistent simplices: " << stats.num_inconsistent_stars << " ("
                << 100. * stats.num_inconsistent_stars / m_points.size() << "%)\n"
                << "==========================================================\n";
    }

    return stats;
  }

  /** \brief Exports the complex into a Simplex_tree.
   *
   * \tparam Simplex_tree_ must be a `Simplex_tree`.
   *
   * @param[out] tree The result, where each `Vertex_handle` is the index of the
   *   corresponding point in the range provided to the constructor (it can also be
   *   retrieved through the `Tangential_complex::get_point` function.
   * @param[in] export_inconsistent_simplices Also export inconsistent simplices or not?
   * @return The maximal dimension of the simplices.
   */
  template <typename Simplex_tree_>
  int create_complex(Simplex_tree_ &tree,
                     bool export_inconsistent_simplices = true
                     /// \cond ADVANCED_PARAMETERS
                     ,
                     bool export_infinite_simplices = false, Simplex_set *p_inconsistent_simplices = NULL
                     /// \endcond
                     ) const {
#if defined(DEBUG_TRACES) || defined(GUDHI_TC_PROFILING)
    std::cerr << yellow << "\nExporting the TC as a Simplex_tree... " << white;
#endif
#ifdef GUDHI_TC_PROFILING
    Gudhi::Clock t;
#endif

    int max_dim = -1;

    // For each triangulation
    for (std::size_t idx = 0; idx < m_points.size(); ++idx) {
      // For each cell of the star
      Star::const_iterator it_inc_simplex = m_stars[idx].begin();
      Star::const_iterator it_inc_simplex_end = m_stars[idx].end();
      for (; it_inc_simplex != it_inc_simplex_end; ++it_inc_simplex) {
        Simplex c = *it_inc_simplex;

        // Don't export infinite cells
        if (!export_infinite_simplices && is_infinite(c)) continue;

        if (static_cast<int>(c.size()) > max_dim) max_dim = static_cast<int>(c.size());
        // Add the missing center vertex
        c.insert(idx);

        if (!export_inconsistent_simplices && !is_simplex_consistent(c)) continue;

        // Try to insert the simplex
        bool inserted = tree.insert_simplex_and_subfaces(c).second;

        // Inconsistent?
        if (p_inconsistent_simplices && inserted && !is_simplex_consistent(c)) {
          p_inconsistent_simplices->insert(c);
        }
      }
    }

#ifdef GUDHI_TC_PROFILING
    t.end();
    std::cerr << yellow << "done in " << t.num_seconds() << " seconds.\n" << white;
#elif defined(DEBUG_TRACES)
    std::cerr << yellow << "done.\n" << white;
#endif

    return max_dim;
  }

  // First clears the complex then exports the TC into it
  // Returns the max dimension of the simplices
  // check_lower_and_higher_dim_simplices : 0 (false), 1 (true), 2 (auto)
  //   If the check is enabled, the function:
  //   - won't insert the simplex if it is already in a higher dim simplex
  //   - will erase any lower-dim simplices that are faces of the new simplex
  //   "auto" (= 2) will enable the check as a soon as it encounters a
  //   simplex whose dimension is different from the previous ones.
  //   N.B.: The check is quite expensive.

  int create_complex(Simplicial_complex &complex, bool export_inconsistent_simplices = true,
                     bool export_infinite_simplices = false, int check_lower_and_higher_dim_simplices = 2,
                     Simplex_set *p_inconsistent_simplices = NULL) const {
#if defined(DEBUG_TRACES) || defined(GUDHI_TC_PROFILING)
    std::cerr << yellow << "\nExporting the TC as a Simplicial_complex... " << white;
#endif
#ifdef GUDHI_TC_PROFILING
    Gudhi::Clock t;
#endif

    int max_dim = -1;
    complex.clear();

    // For each triangulation
    for (std::size_t idx = 0; idx < m_points.size(); ++idx) {
      // For each cell of the star
      Star::const_iterator it_inc_simplex = m_stars[idx].begin();
      Star::const_iterator it_inc_simplex_end = m_stars[idx].end();
      for (; it_inc_simplex != it_inc_simplex_end; ++it_inc_simplex) {
        Simplex c = *it_inc_simplex;

        // Don't export infinite cells
        if (!export_infinite_simplices && is_infinite(c)) continue;

        if (static_cast<int>(c.size()) > max_dim) max_dim = static_cast<int>(c.size());
        // Add the missing center vertex
        c.insert(idx);

        if (!export_inconsistent_simplices && !is_simplex_consistent(c)) continue;

        // Unusual simplex dim?
        if (check_lower_and_higher_dim_simplices == 2 && max_dim != -1 && static_cast<int>(c.size()) != max_dim) {
          // Let's activate the check
          std::cerr << red
                    << "Info: check_lower_and_higher_dim_simplices ACTIVATED. "
                       "Export might be take some time...\n"
                    << white;
          check_lower_and_higher_dim_simplices = 1;
        }

        // Try to insert the simplex
        bool added = complex.add_simplex(c, check_lower_and_higher_dim_simplices == 1);

        // Inconsistent?
        if (p_inconsistent_simplices && added && !is_simplex_consistent(c)) {
          p_inconsistent_simplices->insert(c);
        }
      }
    }

#ifdef GUDHI_TC_PROFILING
    t.end();
    std::cerr << yellow << "done in " << t.num_seconds() << " seconds.\n" << white;
#elif defined(DEBUG_TRACES)
    std::cerr << yellow << "done.\n" << white;
#endif

    return max_dim;
  }

  template <typename ProjectionFunctor = CGAL::Identity<Point> >
  std::ostream &export_to_off(const Simplicial_complex &complex, std::ostream &os,
                              Simplex_set const *p_simpl_to_color_in_red = NULL,
                              Simplex_set const *p_simpl_to_color_in_green = NULL,
                              Simplex_set const *p_simpl_to_color_in_blue = NULL,
                              ProjectionFunctor const &point_projection = ProjectionFunctor()) const {
    return export_to_off(os, false, p_simpl_to_color_in_red, p_simpl_to_color_in_green, p_simpl_to_color_in_blue,
                         &complex, point_projection);
  }

  template <typename ProjectionFunctor = CGAL::Identity<Point> >
  std::ostream &export_to_off(std::ostream &os, bool color_inconsistencies = false,
                              Simplex_set const *p_simpl_to_color_in_red = NULL,
                              Simplex_set const *p_simpl_to_color_in_green = NULL,
                              Simplex_set const *p_simpl_to_color_in_blue = NULL,
                              const Simplicial_complex *p_complex = NULL,
                              ProjectionFunctor const &point_projection = ProjectionFunctor()) const {
    if (m_points.empty()) return os;

    if (m_ambient_dim < 2) {
      std::cerr << "Error: export_to_off => ambient dimension should be >= 2.\n";
      os << "Error: export_to_off => ambient dimension should be >= 2.\n";
      return os;
    }
    if (m_ambient_dim > 3) {
      std::cerr << "Warning: export_to_off => ambient dimension should be "
                   "<= 3. Only the first 3 coordinates will be exported.\n";
    }

    if (m_intrinsic_dim < 1 || m_intrinsic_dim > 3) {
      std::cerr << "Error: export_to_off => intrinsic dimension should be "
                   "between 1 and 3.\n";
      os << "Error: export_to_off => intrinsic dimension should be "
            "between 1 and 3.\n";
      return os;
    }

    std::stringstream output;
    std::size_t num_simplices, num_vertices;
    export_vertices_to_off(output, num_vertices, false, point_projection);
    if (p_complex) {
      export_simplices_to_off(*p_complex, output, num_simplices, p_simpl_to_color_in_red, p_simpl_to_color_in_green,
                              p_simpl_to_color_in_blue);
    } else {
      export_simplices_to_off(output, num_simplices, color_inconsistencies, p_simpl_to_color_in_red,
                              p_simpl_to_color_in_green, p_simpl_to_color_in_blue);
    }

#ifdef GUDHI_TC_EXPORT_NORMALS
    os << "N";
#endif

    os << "OFF \n"
       << num_vertices << " " << num_simplices << " "
       << "0 \n"
       << output.str();

    return os;
  }

 private:
  void refresh_tangential_complex() {
#if defined(DEBUG_TRACES) || defined(GUDHI_TC_PROFILING)
    std::cerr << yellow << "\nRefreshing TC... " << white;
#endif

#ifdef GUDHI_TC_PROFILING
    Gudhi::Clock t;
#endif
#ifdef GUDHI_USE_TBB
    // Parallel
    if (boost::is_convertible<Concurrency_tag, CGAL::Parallel_tag>::value) {
      tbb::parallel_for(tbb::blocked_range<size_t>(0, m_points.size()), Compute_tangent_triangulation(*this));
    } else {
#endif  // GUDHI_USE_TBB
        // Sequential
      for (std::size_t i = 0; i < m_points.size(); ++i) compute_tangent_triangulation(i);
#ifdef GUDHI_USE_TBB
    }
#endif  // GUDHI_USE_TBB

#ifdef GUDHI_TC_PROFILING
    t.end();
    std::cerr << yellow << "done in " << t.num_seconds() << " seconds.\n" << white;
#elif defined(DEBUG_TRACES)
    std::cerr << yellow << "done.\n" << white;
#endif
  }

  // If the list of perturbed points is provided, it is much faster
  template <typename Point_indices_range>
  void refresh_tangential_complex(Point_indices_range const &perturbed_points_indices) {
#if defined(DEBUG_TRACES) || defined(GUDHI_TC_PROFILING)
    std::cerr << yellow << "\nRefreshing TC... " << white;
#endif

#ifdef GUDHI_TC_PROFILING
    Gudhi::Clock t;
#endif

    // ANN tree containing only the perturbed points
    Points_ds updated_pts_ds(m_points, perturbed_points_indices);

#ifdef GUDHI_USE_TBB
    // Parallel
    if (boost::is_convertible<Concurrency_tag, CGAL::Parallel_tag>::value) {
      tbb::parallel_for(tbb::blocked_range<size_t>(0, m_points.size()),
                        Refresh_tangent_triangulation(*this, updated_pts_ds));
    } else {
#endif  // GUDHI_USE_TBB
        // Sequential
      for (std::size_t i = 0; i < m_points.size(); ++i) refresh_tangent_triangulation(i, updated_pts_ds);
#ifdef GUDHI_USE_TBB
    }
#endif  // GUDHI_USE_TBB

#ifdef GUDHI_TC_PROFILING
    t.end();
    std::cerr << yellow << "done in " << t.num_seconds() << " seconds.\n" << white;
#elif defined(DEBUG_TRACES)
    std::cerr << yellow << "done.\n" << white;
#endif
  }

  void export_inconsistent_stars_to_OFF_files(std::string const &filename_base) const {
    // For each triangulation
    for (std::size_t idx = 0; idx < m_points.size(); ++idx) {
      // We build a SC along the way in case it's inconsistent
      Simplicial_complex sc;
      // For each cell
      bool is_inconsistent = false;
      Star::const_iterator it_inc_simplex = m_stars[idx].begin();
      Star::const_iterator it_inc_simplex_end = m_stars[idx].end();
      for (; it_inc_simplex != it_inc_simplex_end; ++it_inc_simplex) {
        // Skip infinite cells
        if (is_infinite(*it_inc_simplex)) continue;

        Simplex c = *it_inc_simplex;
        c.insert(idx);  // Add the missing index

        sc.add_simplex(c);

        // If we do not already know this star is inconsistent, test it
        if (!is_inconsistent && !is_simplex_consistent(c)) is_inconsistent = true;
      }

      if (is_inconsistent) {
        // Export star to OFF file
        std::stringstream output_filename;
        output_filename << filename_base << "_" << idx << ".off";
        std::ofstream off_stream(output_filename.str().c_str());
        export_to_off(sc, off_stream);
      }
    }
  }

  class Compare_distance_to_ref_point {
   public:
    Compare_distance_to_ref_point(Point const &ref, K const &k) : m_ref(ref), m_k(k) {}

    bool operator()(Point const &p1, Point const &p2) {
      typename K::Squared_distance_d sqdist = m_k.squared_distance_d_object();
      return sqdist(p1, m_ref) < sqdist(p2, m_ref);
    }

   private:
    Point const &m_ref;
    K const &m_k;
  };

#ifdef GUDHI_USE_TBB
  // Functor for compute_tangential_complex function
  class Compute_tangent_triangulation {
    Tangential_complex &m_tc;

   public:
    // Constructor
    Compute_tangent_triangulation(Tangential_complex &tc) : m_tc(tc) {}

    // Constructor
    Compute_tangent_triangulation(const Compute_tangent_triangulation &ctt) : m_tc(ctt.m_tc) {}

    // operator()
    void operator()(const tbb::blocked_range<size_t> &r) const {
      for (size_t i = r.begin(); i != r.end(); ++i) m_tc.compute_tangent_triangulation(i);
    }
  };

  // Functor for refresh_tangential_complex function
  class Refresh_tangent_triangulation {
    Tangential_complex &m_tc;
    Points_ds const &m_updated_pts_ds;

   public:
    // Constructor
    Refresh_tangent_triangulation(Tangential_complex &tc, Points_ds const &updated_pts_ds)
        : m_tc(tc), m_updated_pts_ds(updated_pts_ds) {}

    // Constructor
    Refresh_tangent_triangulation(const Refresh_tangent_triangulation &ctt)
        : m_tc(ctt.m_tc), m_updated_pts_ds(ctt.m_updated_pts_ds) {}

    // operator()
    void operator()(const tbb::blocked_range<size_t> &r) const {
      for (size_t i = r.begin(); i != r.end(); ++i) m_tc.refresh_tangent_triangulation(i, m_updated_pts_ds);
    }
  };
#endif  // GUDHI_USE_TBB

  bool is_infinite(Simplex const &s) const { return *s.rbegin() == (std::numeric_limits<std::size_t>::max)(); }

  // Output: "triangulation" is a Regular Triangulation containing at least the
  // star of "center_pt"
  // Returns the handle of the center vertex
  Tr_vertex_handle compute_star(std::size_t i, const Point &center_pt, const Tangent_space_basis &tsb,
                                Triangulation &triangulation, bool verbose = false) {
    int tangent_space_dim = tsb.dimension();
    const Tr_traits &local_tr_traits = triangulation.geom_traits();

    // Kernel functor & objects
    typename K::Squared_distance_d k_sqdist = m_k.squared_distance_d_object();

    // Triangulation's traits functor & objects
    typename Tr_traits::Compute_weight_d point_weight = local_tr_traits.compute_weight_d_object();
    typename Tr_traits::Power_center_d power_center = local_tr_traits.power_center_d_object();

    //***************************************************
    // Build a minimal triangulation in the tangent space
    // (we only need the star of p)
    //***************************************************

    // Insert p
    Tr_point proj_wp;
    if (i == tsb.origin()) {
      // Insert {(0, 0, 0...), m_weights[i]}
      proj_wp = local_tr_traits.construct_weighted_point_d_object()(
          local_tr_traits.construct_point_d_object()(tangent_space_dim, CGAL::ORIGIN), m_weights[i]);
    } else {
      const Weighted_point &wp = compute_perturbed_weighted_point(i);
      proj_wp = project_point_and_compute_weight(wp, tsb, local_tr_traits);
    }

    Tr_vertex_handle center_vertex = triangulation.insert(proj_wp);
    center_vertex->data() = i;
    if (verbose) std::cerr << "* Inserted point #" << i << "\n";

#ifdef GUDHI_TC_VERY_VERBOSE
    std::size_t num_attempts_to_insert_points = 1;
    std::size_t num_inserted_points = 1;
#endif
    // const int NUM_NEIGHBORS = 150;
    // KNS_range ins_range = m_points_ds.k_nearest_neighbors(center_pt, NUM_NEIGHBORS);
    INS_range ins_range = m_points_ds.incremental_nearest_neighbors(center_pt);

    // While building the local triangulation, we keep the radius
    // of the sphere "star sphere" centered at "center_vertex"
    // and which contains all the
    // circumspheres of the star of "center_vertex"
    // If th the m_max_squared_edge_length is set the maximal radius of the "star sphere"
    // is at most square root of m_max_squared_edge_length
    boost::optional<FT> squared_star_sphere_radius_plus_margin = m_max_squared_edge_length;

    // Insert points until we find a point which is outside "star sphere"
    for (auto nn_it = ins_range.begin(); nn_it != ins_range.end(); ++nn_it) {
      std::size_t neighbor_point_idx = nn_it->first;

      // ith point = p, which is already inserted
      if (neighbor_point_idx != i) {
        // No need to lock the Mutex_for_perturb here since this will not be
        // called while other threads are perturbing the positions
        Point neighbor_pt;
        FT neighbor_weight;
        compute_perturbed_weighted_point(neighbor_point_idx, neighbor_pt, neighbor_weight);
        GUDHI_CHECK(!m_max_squared_edge_length ||
                    squared_star_sphere_radius_plus_margin.value() <= m_max_squared_edge_length.value(),
                    std::invalid_argument("Tangential_complex::compute_star - set a bigger value with set_max_squared_edge_length."));
        if (squared_star_sphere_radius_plus_margin &&
            k_sqdist(center_pt, neighbor_pt) > squared_star_sphere_radius_plus_margin.value()) {
          GUDHI_CHECK(triangulation.current_dimension() >= tangent_space_dim,
                      std::invalid_argument("Tangential_complex::compute_star - Dimension of the star is only " + \
                                            std::to_string(triangulation.current_dimension())));
          break;
        }

        Tr_point proj_pt = project_point_and_compute_weight(neighbor_pt, neighbor_weight, tsb, local_tr_traits);

#ifdef GUDHI_TC_VERY_VERBOSE
        ++num_attempts_to_insert_points;
#endif

        Tr_vertex_handle vh = triangulation.insert_if_in_star(proj_pt, center_vertex);
        // Tr_vertex_handle vh = triangulation.insert(proj_pt);
        if (vh != Tr_vertex_handle() && vh->data() == (std::numeric_limits<std::size_t>::max)()) {
#ifdef GUDHI_TC_VERY_VERBOSE
          ++num_inserted_points;
#endif
          if (verbose) std::cerr << "* Inserted point #" << neighbor_point_idx << "\n";

          vh->data() = neighbor_point_idx;

          // Let's recompute squared_star_sphere_radius_plus_margin
          if (triangulation.current_dimension() >= tangent_space_dim) {
            squared_star_sphere_radius_plus_margin = boost::none;
            // Get the incident cells and look for the biggest circumsphere
            std::vector<Tr_full_cell_handle> incident_cells;
            triangulation.incident_full_cells(center_vertex, std::back_inserter(incident_cells));
            for (typename std::vector<Tr_full_cell_handle>::iterator cit = incident_cells.begin();
                 cit != incident_cells.end(); ++cit) {
              Tr_full_cell_handle cell = *cit;
              if (triangulation.is_infinite(cell)) {
                squared_star_sphere_radius_plus_margin = boost::none;
                break;
              } else {
                // Note that this uses the perturbed point since it uses
                // the points of the local triangulation
                Tr_point c =
                    power_center(boost::make_transform_iterator(cell->vertices_begin(),
                                                                vertex_handle_to_point<Tr_point, Tr_vertex_handle>),
                                 boost::make_transform_iterator(cell->vertices_end(),
                                                                vertex_handle_to_point<Tr_point, Tr_vertex_handle>));

                FT sq_power_sphere_diam = 4 * point_weight(c);

                if (!squared_star_sphere_radius_plus_margin ||
                    sq_power_sphere_diam > squared_star_sphere_radius_plus_margin.value()) {
                  squared_star_sphere_radius_plus_margin = sq_power_sphere_diam;
                }
              }
            }

            // Let's add the margin, now
            // The value depends on whether we perturb weight or position
            if (squared_star_sphere_radius_plus_margin) {
              // "2*m_last_max_perturb" because both points can be perturbed
              squared_star_sphere_radius_plus_margin =
                  CGAL::square(std::sqrt(squared_star_sphere_radius_plus_margin.value()) + 2 * m_last_max_perturb);

              // Reduce the square radius to  m_max_squared_edge_length if necessary
              if (m_max_squared_edge_length && squared_star_sphere_radius_plus_margin.value() > m_max_squared_edge_length.value()) {
                squared_star_sphere_radius_plus_margin = m_max_squared_edge_length.value();
              }

              // Save it in `m_squared_star_spheres_radii_incl_margin`
              m_squared_star_spheres_radii_incl_margin[i] = squared_star_sphere_radius_plus_margin.value();
            } else {
              if (m_max_squared_edge_length) {
                squared_star_sphere_radius_plus_margin = m_max_squared_edge_length.value();
                m_squared_star_spheres_radii_incl_margin[i] = m_max_squared_edge_length.value();
              } else {
                m_squared_star_spheres_radii_incl_margin[i] = FT(-1);
              }
            }
          }
        }
      }
    }

    return center_vertex;
  }

  void refresh_tangent_triangulation(std::size_t i, Points_ds const &updated_pts_ds, bool verbose = false) {
    if (verbose) std::cerr << "** Refreshing tangent tri #" << i << " **\n";

    if (m_squared_star_spheres_radii_incl_margin[i] == FT(-1)) return compute_tangent_triangulation(i, verbose);

    Point center_point = compute_perturbed_point(i);
    // Among updated point, what is the closer from our center point?
    std::size_t closest_pt_index = updated_pts_ds.k_nearest_neighbors(center_point, 1, false).begin()->first;

    typename K::Construct_weighted_point_d k_constr_wp = m_k.construct_weighted_point_d_object();
    typename K::Power_distance_d k_power_dist = m_k.power_distance_d_object();

    // Construct a weighted point equivalent to the star sphere
    Weighted_point star_sphere = k_constr_wp(compute_perturbed_point(i), m_squared_star_spheres_radii_incl_margin[i]);
    Weighted_point closest_updated_point = compute_perturbed_weighted_point(closest_pt_index);

    // Is the "closest point" inside our star sphere?
    if (k_power_dist(star_sphere, closest_updated_point) <= FT(0)) compute_tangent_triangulation(i, verbose);
  }

  void compute_tangent_triangulation(std::size_t i, bool verbose = false) {
    if (verbose) std::cerr << "** Computing tangent tri #" << i << " **\n";
    // std::cerr << "***********************************************\n";

    // No need to lock the mutex here since this will not be called while
    // other threads are perturbing the positions
    const Point center_pt = compute_perturbed_point(i);
    Tangent_space_basis &tsb = m_tangent_spaces[i];

    // Estimate the tangent space
    if (!m_are_tangent_spaces_computed[i]) {
#ifdef GUDHI_TC_EXPORT_NORMALS
      tsb = compute_tangent_space(center_pt, i, true /*normalize*/, &m_orth_spaces[i]);
#else
      tsb = compute_tangent_space(center_pt, i);
#endif
    }

#if defined(GUDHI_TC_PROFILING) && defined(GUDHI_TC_VERY_VERBOSE)
    Gudhi::Clock t;
#endif
    int tangent_space_dim = tangent_basis_dim(i);
    Triangulation &local_tr = m_triangulations[i].construct_triangulation(tangent_space_dim);

    m_triangulations[i].center_vertex() = compute_star(i, center_pt, tsb, local_tr, verbose);

#if defined(GUDHI_TC_PROFILING) && defined(GUDHI_TC_VERY_VERBOSE)
    t.end();
    std::cerr << "  - triangulation construction: " << t.num_seconds() << " s.\n";
    t.reset();
#endif

#ifdef GUDHI_TC_VERY_VERBOSE
    std::cerr << "Inserted " << num_inserted_points << " points / " << num_attempts_to_insert_points
              << " attemps to compute the star\n";
#endif

    update_star(i);

#if defined(GUDHI_TC_PROFILING) && defined(GUDHI_TC_VERY_VERBOSE)
    t.end();
    std::cerr << "  - update_star: " << t.num_seconds() << " s.\n";
#endif
  }

  // Updates m_stars[i] directly from m_triangulations[i]

  void update_star(std::size_t i) {
    Star &star = m_stars[i];
    star.clear();
    Triangulation &local_tr = m_triangulations[i].tr();
    Tr_vertex_handle center_vertex = m_triangulations[i].center_vertex();
    int cur_dim_plus_1 = local_tr.current_dimension() + 1;

    std::vector<Tr_full_cell_handle> incident_cells;
    local_tr.incident_full_cells(center_vertex, std::back_inserter(incident_cells));

    typename std::vector<Tr_full_cell_handle>::const_iterator it_c = incident_cells.begin();
    typename std::vector<Tr_full_cell_handle>::const_iterator it_c_end = incident_cells.end();
    // For each cell
    for (; it_c != it_c_end; ++it_c) {
      // Will contain all indices except center_vertex
      Incident_simplex incident_simplex;
      for (int j = 0; j < cur_dim_plus_1; ++j) {
        std::size_t index = (*it_c)->vertex(j)->data();
        if (index != i) incident_simplex.insert(index);
      }
      GUDHI_CHECK(incident_simplex.size() == cur_dim_plus_1 - 1,
                  std::logic_error("update_star: wrong size of incident simplex"));
      star.push_back(incident_simplex);
    }
  }

  // Estimates tangent subspaces using PCA

  Tangent_space_basis compute_tangent_space(const Point &p, const std::size_t i, bool normalize_basis = true,
                                            Orthogonal_space_basis *p_orth_space_basis = NULL) {
    unsigned int num_pts_for_pca =
        (std::min)(static_cast<unsigned int>(std::pow(GUDHI_TC_BASE_VALUE_FOR_PCA, m_intrinsic_dim)),
                   static_cast<unsigned int>(m_points.size()));

    // Kernel functors
    typename K::Construct_vector_d constr_vec = m_k.construct_vector_d_object();
    typename K::Compute_coordinate_d coord = m_k.compute_coordinate_d_object();

#ifdef GUDHI_TC_USE_ANOTHER_POINT_SET_FOR_TANGENT_SPACE_ESTIM
    KNS_range kns_range = m_points_ds_for_tse.k_nearest_neighbors(p, num_pts_for_pca, false);
    const Points &points_for_pca = m_points_for_tse;
#else
    KNS_range kns_range = m_points_ds.k_nearest_neighbors(p, num_pts_for_pca, false);
    const Points &points_for_pca = m_points;
#endif

    // One row = one point
    Eigen::MatrixXd mat_points(num_pts_for_pca, m_ambient_dim);
    auto nn_it = kns_range.begin();
    for (unsigned int j = 0; j < num_pts_for_pca && nn_it != kns_range.end(); ++j, ++nn_it) {
      for (int i = 0; i < m_ambient_dim; ++i) {
        mat_points(j, i) = CGAL::to_double(coord(points_for_pca[nn_it->first], i));
      }
    }
    Eigen::MatrixXd centered = mat_points.rowwise() - mat_points.colwise().mean();
    Eigen::MatrixXd cov = centered.adjoint() * centered;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);

    Tangent_space_basis tsb(i);  // p = compute_perturbed_point(i) here

    // The eigenvectors are sorted in increasing order of their corresponding
    // eigenvalues
    for (int j = m_ambient_dim - 1; j >= m_ambient_dim - m_intrinsic_dim; --j) {
      if (normalize_basis) {
        Vector v = constr_vec(m_ambient_dim, eig.eigenvectors().col(j).data(),
                              eig.eigenvectors().col(j).data() + m_ambient_dim);
        tsb.push_back(normalize_vector(v, m_k));
      } else {
        tsb.push_back(constr_vec(m_ambient_dim, eig.eigenvectors().col(j).data(),
                                 eig.eigenvectors().col(j).data() + m_ambient_dim));
      }
    }

    if (p_orth_space_basis) {
      p_orth_space_basis->set_origin(i);
      for (int j = m_ambient_dim - m_intrinsic_dim - 1; j >= 0; --j) {
        if (normalize_basis) {
          Vector v = constr_vec(m_ambient_dim, eig.eigenvectors().col(j).data(),
                                eig.eigenvectors().col(j).data() + m_ambient_dim);
          p_orth_space_basis->push_back(normalize_vector(v, m_k));
        } else {
          p_orth_space_basis->push_back(constr_vec(m_ambient_dim, eig.eigenvectors().col(j).data(),
                                                   eig.eigenvectors().col(j).data() + m_ambient_dim));
        }
      }
    }

    m_are_tangent_spaces_computed[i] = true;

    return tsb;
  }

  // Compute the space tangent to a simplex (p1, p2, ... pn)
  // TODO(CJ): Improve this?
  // Basically, it takes all the neighbor points to p1, p2... pn and runs PCA
  // on it. Note that most points are duplicated.

  Tangent_space_basis compute_tangent_space(const Simplex &s, bool normalize_basis = true) {
    unsigned int num_pts_for_pca =
        (std::min)(static_cast<unsigned int>(std::pow(GUDHI_TC_BASE_VALUE_FOR_PCA, m_intrinsic_dim)),
                   static_cast<unsigned int>(m_points.size()));

    // Kernel functors
    typename K::Construct_vector_d constr_vec = m_k.construct_vector_d_object();
    typename K::Compute_coordinate_d coord = m_k.compute_coordinate_d_object();
    typename K::Squared_length_d sqlen = m_k.squared_length_d_object();
    typename K::Scaled_vector_d scaled_vec = m_k.scaled_vector_d_object();
    typename K::Scalar_product_d scalar_pdct = m_k.scalar_product_d_object();
    typename K::Difference_of_vectors_d diff_vec = m_k.difference_of_vectors_d_object();

    // One row = one point
    Eigen::MatrixXd mat_points(s.size() * num_pts_for_pca, m_ambient_dim);
    unsigned int current_row = 0;

    for (Simplex::const_iterator it_index = s.begin(); it_index != s.end(); ++it_index) {
      const Point &p = m_points[*it_index];

#ifdef GUDHI_TC_USE_ANOTHER_POINT_SET_FOR_TANGENT_SPACE_ESTIM
      KNS_range kns_range = m_points_ds_for_tse.k_nearest_neighbors(p, num_pts_for_pca, false);
      const Points &points_for_pca = m_points_for_tse;
#else
      KNS_range kns_range = m_points_ds.k_nearest_neighbors(p, num_pts_for_pca, false);
      const Points &points_for_pca = m_points;
#endif

      auto nn_it = kns_range.begin();
      for (; current_row < num_pts_for_pca && nn_it != kns_range.end(); ++current_row, ++nn_it) {
        for (int i = 0; i < m_ambient_dim; ++i) {
          mat_points(current_row, i) = CGAL::to_double(coord(points_for_pca[nn_it->first], i));
        }
      }
    }
    Eigen::MatrixXd centered = mat_points.rowwise() - mat_points.colwise().mean();
    Eigen::MatrixXd cov = centered.adjoint() * centered;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);

    Tangent_space_basis tsb;

    // The eigenvectors are sorted in increasing order of their corresponding
    // eigenvalues
    for (int j = m_ambient_dim - 1; j >= m_ambient_dim - m_intrinsic_dim; --j) {
      if (normalize_basis) {
        Vector v = constr_vec(m_ambient_dim, eig.eigenvectors().col(j).data(),
                              eig.eigenvectors().col(j).data() + m_ambient_dim);
        tsb.push_back(normalize_vector(v, m_k));
      } else {
        tsb.push_back(constr_vec(m_ambient_dim, eig.eigenvectors().col(j).data(),
                                 eig.eigenvectors().col(j).data() + m_ambient_dim));
      }
    }

    return tsb;
  }

  // Returns the dimension of the ith local triangulation

  int tangent_basis_dim(std::size_t i) const { return m_tangent_spaces[i].dimension(); }

  Point compute_perturbed_point(std::size_t pt_idx) const {
#ifdef GUDHI_TC_PERTURB_POSITION
    return m_k.translated_point_d_object()(m_points[pt_idx], m_translations[pt_idx]);
#else
    return m_points[pt_idx];
#endif
  }

  void compute_perturbed_weighted_point(std::size_t pt_idx, Point &p, FT &w) const {
#ifdef GUDHI_TC_PERTURB_POSITION
    p = m_k.translated_point_d_object()(m_points[pt_idx], m_translations[pt_idx]);
#else
    p = m_points[pt_idx];
#endif
    w = m_weights[pt_idx];
  }

  Weighted_point compute_perturbed_weighted_point(std::size_t pt_idx) const {
    typename K::Construct_weighted_point_d k_constr_wp = m_k.construct_weighted_point_d_object();

    Weighted_point wp = k_constr_wp(
#ifdef GUDHI_TC_PERTURB_POSITION
        m_k.translated_point_d_object()(m_points[pt_idx], m_translations[pt_idx]),
#else
        m_points[pt_idx],
#endif
        m_weights[pt_idx]);

    return wp;
  }

  Point unproject_point(const Tr_point &p, const Tangent_space_basis &tsb, const Tr_traits &tr_traits) const {
    typename K::Translated_point_d k_transl = m_k.translated_point_d_object();
    typename K::Scaled_vector_d k_scaled_vec = m_k.scaled_vector_d_object();
    typename Tr_traits::Compute_coordinate_d coord = tr_traits.compute_coordinate_d_object();

    Point global_point = compute_perturbed_point(tsb.origin());
    for (int i = 0; i < m_intrinsic_dim; ++i) global_point = k_transl(global_point, k_scaled_vec(tsb[i], coord(p, i)));

    return global_point;
  }

  // Project the point in the tangent space
  // Resulting point coords are expressed in tsb's space
  Tr_bare_point project_point(const Point &p, const Tangent_space_basis &tsb, const Tr_traits &tr_traits) const {
    typename K::Scalar_product_d scalar_pdct = m_k.scalar_product_d_object();
    typename K::Difference_of_points_d diff_points = m_k.difference_of_points_d_object();

    Vector v = diff_points(p, compute_perturbed_point(tsb.origin()));

    std::vector<FT> coords;
    // Ambiant-space coords of the projected point
    coords.reserve(tsb.dimension());
    for (std::size_t i = 0; i < m_intrinsic_dim; ++i) {
      // Local coords are given by the scalar product with the vectors of tsb
      FT coord = scalar_pdct(v, tsb[i]);
      coords.push_back(coord);
    }

    return tr_traits.construct_point_d_object()(static_cast<int>(coords.size()), coords.begin(), coords.end());
  }

  // Project the point in the tangent space
  // The weight will be the squared distance between p and the projection of p
  // Resulting point coords are expressed in tsb's space

  Tr_point project_point_and_compute_weight(const Weighted_point &wp, const Tangent_space_basis &tsb,
                                            const Tr_traits &tr_traits) const {
    typename K::Point_drop_weight_d k_drop_w = m_k.point_drop_weight_d_object();
    typename K::Compute_weight_d k_point_weight = m_k.compute_weight_d_object();
    return project_point_and_compute_weight(k_drop_w(wp), k_point_weight(wp), tsb, tr_traits);
  }

  // Same as above, with slightly different parameters
  Tr_point project_point_and_compute_weight(const Point &p, const FT w, const Tangent_space_basis &tsb,
                                            const Tr_traits &tr_traits) const {
    const int point_dim = m_k.point_dimension_d_object()(p);

    typename K::Construct_point_d constr_pt = m_k.construct_point_d_object();
    typename K::Scalar_product_d scalar_pdct = m_k.scalar_product_d_object();
    typename K::Difference_of_points_d diff_points = m_k.difference_of_points_d_object();
    typename K::Compute_coordinate_d coord = m_k.compute_coordinate_d_object();
    typename K::Construct_cartesian_const_iterator_d ccci = m_k.construct_cartesian_const_iterator_d_object();

    Point origin = compute_perturbed_point(tsb.origin());
    Vector v = diff_points(p, origin);

    // Same dimension? Then the weight is 0
    bool same_dim = (point_dim == tsb.dimension());

    std::vector<FT> coords;
    // Ambiant-space coords of the projected point
    std::vector<FT> p_proj(ccci(origin), ccci(origin, 0));
    coords.reserve(tsb.dimension());
    for (int i = 0; i < tsb.dimension(); ++i) {
      // Local coords are given by the scalar product with the vectors of tsb
      FT c = scalar_pdct(v, tsb[i]);
      coords.push_back(c);

      // p_proj += c * tsb[i]
      if (!same_dim) {
        for (int j = 0; j < point_dim; ++j) p_proj[j] += c * coord(tsb[i], j);
      }
    }

    // Same dimension? Then the weight is 0
    FT sq_dist_to_proj_pt = 0;
    if (!same_dim) {
      Point projected_pt = constr_pt(point_dim, p_proj.begin(), p_proj.end());
      sq_dist_to_proj_pt = m_k.squared_distance_d_object()(p, projected_pt);
    }

    return tr_traits.construct_weighted_point_d_object()(
        tr_traits.construct_point_d_object()(static_cast<int>(coords.size()), coords.begin(), coords.end()),
        w - sq_dist_to_proj_pt);
  }

  // Project all the points in the tangent space

  template <typename Indexed_point_range>
  std::vector<Tr_point> project_points_and_compute_weights(const Indexed_point_range &point_indices,
                                                           const Tangent_space_basis &tsb,
                                                           const Tr_traits &tr_traits) const {
    std::vector<Tr_point> ret;
    for (typename Indexed_point_range::const_iterator it = point_indices.begin(), it_end = point_indices.end();
         it != it_end; ++it) {
      ret.push_back(project_point_and_compute_weight(compute_perturbed_weighted_point(*it), tsb, tr_traits));
    }
    return ret;
  }

  // A simplex here is a local tri's full cell handle

  bool is_simplex_consistent(Tr_full_cell_handle fch, int cur_dim) const {
    Simplex c;
    for (int i = 0; i < cur_dim + 1; ++i) {
      std::size_t data = fch->vertex(i)->data();
      c.insert(data);
    }
    return is_simplex_consistent(c);
  }

  // A simplex here is a list of point indices
  // TODO(CJ): improve it like the other "is_simplex_consistent" below

  bool is_simplex_consistent(Simplex const &simplex) const {
    // Check if the simplex is in the stars of all its vertices
    Simplex::const_iterator it_point_idx = simplex.begin();
    // For each point p of the simplex, we parse the incidents cells of p
    // and we check if "simplex" is among them
    for (; it_point_idx != simplex.end(); ++it_point_idx) {
      std::size_t point_idx = *it_point_idx;
      // Don't check infinite simplices
      if (point_idx == (std::numeric_limits<std::size_t>::max)()) continue;

      Star const &star = m_stars[point_idx];

      // What we're looking for is "simplex" \ point_idx
      Incident_simplex is_to_find = simplex;
      is_to_find.erase(point_idx);

      // For each cell
      if (std::find(star.begin(), star.end(), is_to_find) == star.end()) return false;
    }

    return true;
  }

  // A simplex here is a list of point indices
  // "s" contains all the points of the simplex except "center_point"
  // This function returns the points whose star doesn't contain the simplex
  // N.B.: the function assumes that the simplex is contained in
  //       star(center_point)

  template <typename OutputIterator>  // value_type = std::size_t
  bool is_simplex_consistent(std::size_t center_point,
                             Incident_simplex const &s,  // without "center_point"
                             OutputIterator points_whose_star_does_not_contain_s,
                             bool check_also_in_non_maximal_faces = false) const {
    Simplex full_simplex = s;
    full_simplex.insert(center_point);

    // Check if the simplex is in the stars of all its vertices
    Incident_simplex::const_iterator it_point_idx = s.begin();
    // For each point p of the simplex, we parse the incidents cells of p
    // and we check if "simplex" is among them
    for (; it_point_idx != s.end(); ++it_point_idx) {
      std::size_t point_idx = *it_point_idx;
      // Don't check infinite simplices
      if (point_idx == (std::numeric_limits<std::size_t>::max)()) continue;

      Star const &star = m_stars[point_idx];

      // What we're looking for is full_simplex \ point_idx
      Incident_simplex is_to_find = full_simplex;
      is_to_find.erase(point_idx);

      if (check_also_in_non_maximal_faces) {
        // For each simplex "is" of the star, check if ic_to_simplex is
        // included in "is"
        bool found = false;
        for (Star::const_iterator is = star.begin(), is_end = star.end(); !found && is != is_end; ++is) {
          if (std::includes(is->begin(), is->end(), is_to_find.begin(), is_to_find.end())) found = true;
        }

        if (!found) *points_whose_star_does_not_contain_s++ = point_idx;
      } else {
        // Does the star contain is_to_find?
        if (std::find(star.begin(), star.end(), is_to_find) == star.end())
          *points_whose_star_does_not_contain_s++ = point_idx;
      }
    }

    return true;
  }

  // A simplex here is a list of point indices
  // It looks for s in star(p).
  // "s" contains all the points of the simplex except p.
  bool is_simplex_in_star(std::size_t p, Incident_simplex const &s, bool check_also_in_non_maximal_faces = true) const {
    Star const &star = m_stars[p];

    if (check_also_in_non_maximal_faces) {
      // For each simplex "is" of the star, check if ic_to_simplex is
      // included in "is"
      bool found = false;
      for (Star::const_iterator is = star.begin(), is_end = star.end(); !found && is != is_end; ++is) {
        if (std::includes(is->begin(), is->end(), s.begin(), s.end())) found = true;
      }

      return found;
    } else {
      return !(std::find(star.begin(), star.end(), s) == star.end());
    }
  }

#ifdef GUDHI_USE_TBB
  // Functor for try_to_solve_inconsistencies_in_a_local_triangulation function
  class Try_to_solve_inconsistencies_in_a_local_triangulation {
    Tangential_complex &m_tc;
    double m_max_perturb;
    tbb::combinable<std::size_t> &m_num_inconsistencies;
    tbb::combinable<std::vector<std::size_t> > &m_updated_points;

   public:
    // Constructor
    Try_to_solve_inconsistencies_in_a_local_triangulation(Tangential_complex &tc, double max_perturb,
                                                          tbb::combinable<std::size_t> &num_inconsistencies,
                                                          tbb::combinable<std::vector<std::size_t> > &updated_points)
        : m_tc(tc),
          m_max_perturb(max_perturb),
          m_num_inconsistencies(num_inconsistencies),
          m_updated_points(updated_points) {}

    // Constructor
    Try_to_solve_inconsistencies_in_a_local_triangulation(
        const Try_to_solve_inconsistencies_in_a_local_triangulation &tsilt)
        : m_tc(tsilt.m_tc),
          m_max_perturb(tsilt.m_max_perturb),
          m_num_inconsistencies(tsilt.m_num_inconsistencies),
          m_updated_points(tsilt.m_updated_points) {}

    // operator()
    void operator()(const tbb::blocked_range<size_t> &r) const {
      for (size_t i = r.begin(); i != r.end(); ++i) {
        m_num_inconsistencies.local() += m_tc.try_to_solve_inconsistencies_in_a_local_triangulation(
            i, m_max_perturb, std::back_inserter(m_updated_points.local()));
      }
    }
  };
#endif  // GUDHI_USE_TBB

  void perturb(std::size_t point_idx, double max_perturb) {
    const Tr_traits &local_tr_traits = m_triangulations[point_idx].tr().geom_traits();
    typename Tr_traits::Compute_coordinate_d coord = local_tr_traits.compute_coordinate_d_object();
    typename K::Translated_point_d k_transl = m_k.translated_point_d_object();
    typename K::Construct_vector_d k_constr_vec = m_k.construct_vector_d_object();
    typename K::Scaled_vector_d k_scaled_vec = m_k.scaled_vector_d_object();

    CGAL::Random_points_in_ball_d<Tr_bare_point> tr_point_in_ball_generator(
        m_intrinsic_dim, m_random_generator.get_double(0., max_perturb));

    Tr_point local_random_transl =
        local_tr_traits.construct_weighted_point_d_object()(*tr_point_in_ball_generator++, 0);
    Translation_for_perturb global_transl = k_constr_vec(m_ambient_dim);
    const Tangent_space_basis &tsb = m_tangent_spaces[point_idx];
    for (int i = 0; i < m_intrinsic_dim; ++i) {
      global_transl = k_transl(global_transl, k_scaled_vec(tsb[i], coord(local_random_transl, i)));
    }
    // Parallel
#if defined(GUDHI_USE_TBB)
    m_p_perturb_mutexes[point_idx].lock();
    m_translations[point_idx] = global_transl;
    m_p_perturb_mutexes[point_idx].unlock();
    // Sequential
#else
    m_translations[point_idx] = global_transl;
#endif
  }

  // Return true if inconsistencies were found
  template <typename OutputIt>
  bool try_to_solve_inconsistencies_in_a_local_triangulation(
      std::size_t tr_index, double max_perturb, OutputIt perturbed_pts_indices = CGAL::Emptyset_iterator()) {
    bool is_inconsistent = false;

    Star const &star = m_stars[tr_index];

    // For each incident simplex
    Star::const_iterator it_inc_simplex = star.begin();
    Star::const_iterator it_inc_simplex_end = star.end();
    for (; it_inc_simplex != it_inc_simplex_end; ++it_inc_simplex) {
      const Incident_simplex &incident_simplex = *it_inc_simplex;

      // Don't check infinite cells
      if (is_infinite(incident_simplex)) continue;

      Simplex c = incident_simplex;
      c.insert(tr_index);  // Add the missing index

      // Perturb the center point
      if (!is_simplex_consistent(c)) {
        is_inconsistent = true;

        std::size_t idx = tr_index;

        perturb(tr_index, max_perturb);
        *perturbed_pts_indices++ = idx;

        // We will try the other cells next time
        break;
      }
    }

    return is_inconsistent;
  }

  // 1st line: number of points
  // Then one point per line
  std::ostream &export_point_set(std::ostream &os, bool use_perturbed_points = false,
                                 const char *coord_separator = " ") const {
    if (use_perturbed_points) {
      std::vector<Point> perturbed_points;
      perturbed_points.reserve(m_points.size());
      for (std::size_t i = 0; i < m_points.size(); ++i) perturbed_points.push_back(compute_perturbed_point(i));

      return export_point_set(m_k, perturbed_points, os, coord_separator);
    } else {
      return export_point_set(m_k, m_points, os, coord_separator);
    }
  }

  template <typename ProjectionFunctor = CGAL::Identity<Point> >
  std::ostream &export_vertices_to_off(std::ostream &os, std::size_t &num_vertices, bool use_perturbed_points = false,
                                       ProjectionFunctor const &point_projection = ProjectionFunctor()) const {
    if (m_points.empty()) {
      num_vertices = 0;
      return os;
    }

    // If m_intrinsic_dim = 1, we output each point two times
    // to be able to export each segment as a flat triangle with 3 different
    // indices (otherwise, Meshlab detects degenerated simplices)
    const int N = (m_intrinsic_dim == 1 ? 2 : 1);

    // Kernel functors
    typename K::Compute_coordinate_d coord = m_k.compute_coordinate_d_object();

#ifdef GUDHI_TC_EXPORT_ALL_COORDS_IN_OFF
    int num_coords = m_ambient_dim;
#else
    int num_coords = (std::min)(m_ambient_dim, 3);
#endif

#ifdef GUDHI_TC_EXPORT_NORMALS
    OS_container::const_iterator it_os = m_orth_spaces.begin();
#endif
    typename Points::const_iterator it_p = m_points.begin();
    typename Points::const_iterator it_p_end = m_points.end();
    // For each point p
    for (std::size_t i = 0; it_p != it_p_end; ++it_p, ++i) {
      Point p = point_projection(use_perturbed_points ? compute_perturbed_point(i) : *it_p);
      for (int ii = 0; ii < N; ++ii) {
        int j = 0;
        for (; j < num_coords; ++j) os << CGAL::to_double(coord(p, j)) << " ";
        if (j == 2) os << "0";

#ifdef GUDHI_TC_EXPORT_NORMALS
        for (j = 0; j < num_coords; ++j) os << " " << CGAL::to_double(coord(*it_os->begin(), j));
#endif
        os << "\n";
      }
#ifdef GUDHI_TC_EXPORT_NORMALS
      ++it_os;
#endif
    }

    num_vertices = N * m_points.size();
    return os;
  }

  std::ostream &export_simplices_to_off(std::ostream &os, std::size_t &num_OFF_simplices,
                                        bool color_inconsistencies = false,
                                        Simplex_set const *p_simpl_to_color_in_red = NULL,
                                        Simplex_set const *p_simpl_to_color_in_green = NULL,
                                        Simplex_set const *p_simpl_to_color_in_blue = NULL) const {
    // If m_intrinsic_dim = 1, each point is output two times
    // (see export_vertices_to_off)
    num_OFF_simplices = 0;
    std::size_t num_maximal_simplices = 0;
    std::size_t num_inconsistent_maximal_simplices = 0;
    std::size_t num_inconsistent_stars = 0;
    typename Tr_container::const_iterator it_tr = m_triangulations.begin();
    typename Tr_container::const_iterator it_tr_end = m_triangulations.end();
    // For each triangulation
    for (std::size_t idx = 0; it_tr != it_tr_end; ++it_tr, ++idx) {
      bool is_star_inconsistent = false;

      Triangulation const &tr = it_tr->tr();

      if (tr.current_dimension() < m_intrinsic_dim) continue;

      // Color for this star
      std::stringstream color;
      // color << rand()%256 << " " << 100+rand()%156 << " " << 100+rand()%156;
      color << 128 << " " << 128 << " " << 128;

      // Gather the triangles here, with an int telling its color
      typedef std::vector<std::pair<Simplex, int> > Star_using_triangles;
      Star_using_triangles star_using_triangles;

      // For each cell of the star
      Star::const_iterator it_inc_simplex = m_stars[idx].begin();
      Star::const_iterator it_inc_simplex_end = m_stars[idx].end();
      for (; it_inc_simplex != it_inc_simplex_end; ++it_inc_simplex) {
        Simplex c = *it_inc_simplex;
        c.insert(idx);
        std::size_t num_vertices = c.size();
        ++num_maximal_simplices;

        int color_simplex = -1;  // -1=no color, 0=yellow, 1=red, 2=green, 3=blue
        if (color_inconsistencies && !is_simplex_consistent(c)) {
          ++num_inconsistent_maximal_simplices;
          color_simplex = 0;
          is_star_inconsistent = true;
        } else {
          if (p_simpl_to_color_in_red && std::find(p_simpl_to_color_in_red->begin(), p_simpl_to_color_in_red->end(),
                                                   c) != p_simpl_to_color_in_red->end()) {
            color_simplex = 1;
          } else if (p_simpl_to_color_in_green &&
                     std::find(p_simpl_to_color_in_green->begin(), p_simpl_to_color_in_green->end(), c) !=
                         p_simpl_to_color_in_green->end()) {
            color_simplex = 2;
          } else if (p_simpl_to_color_in_blue &&
                     std::find(p_simpl_to_color_in_blue->begin(), p_simpl_to_color_in_blue->end(), c) !=
                         p_simpl_to_color_in_blue->end()) {
            color_simplex = 3;
          }
        }

        // If m_intrinsic_dim = 1, each point is output two times,
        // so we need to multiply each index by 2
        // And if only 2 vertices, add a third one (each vertex is duplicated in
        // the file when m_intrinsic dim = 2)
        if (m_intrinsic_dim == 1) {
          Simplex tmp_c;
          Simplex::iterator it = c.begin();
          for (; it != c.end(); ++it) tmp_c.insert(*it * 2);
          if (num_vertices == 2) tmp_c.insert(*tmp_c.rbegin() + 1);

          c = tmp_c;
        }

        if (num_vertices <= 3) {
          star_using_triangles.push_back(std::make_pair(c, color_simplex));
        } else {
          // num_vertices >= 4: decompose the simplex in triangles
          std::vector<bool> booleans(num_vertices, false);
          std::fill(booleans.begin() + num_vertices - 3, booleans.end(), true);
          do {
            Simplex triangle;
            Simplex::iterator it = c.begin();
            for (int i = 0; it != c.end(); ++i, ++it) {
              if (booleans[i]) triangle.insert(*it);
            }
            star_using_triangles.push_back(std::make_pair(triangle, color_simplex));
          } while (std::next_permutation(booleans.begin(), booleans.end()));
        }
      }

      // For each cell
      Star_using_triangles::const_iterator it_simplex = star_using_triangles.begin();
      Star_using_triangles::const_iterator it_simplex_end = star_using_triangles.end();
      for (; it_simplex != it_simplex_end; ++it_simplex) {
        const Simplex &c = it_simplex->first;

        // Don't export infinite cells
        if (is_infinite(c)) continue;

        int color_simplex = it_simplex->second;

        std::stringstream sstr_c;

        Simplex::const_iterator it_point_idx = c.begin();
        for (; it_point_idx != c.end(); ++it_point_idx) {
          sstr_c << *it_point_idx << " ";
        }

        os << 3 << " " << sstr_c.str();
        if (color_inconsistencies || p_simpl_to_color_in_red || p_simpl_to_color_in_green || p_simpl_to_color_in_blue) {
          switch (color_simplex) {
            case 0:
              os << " 255 255 0";
              break;
            case 1:
              os << " 255 0 0";
              break;
            case 2:
              os << " 0 255 0";
              break;
            case 3:
              os << " 0 0 255";
              break;
            default:
              os << " " << color.str();
              break;
          }
        }
        ++num_OFF_simplices;
        os << "\n";
      }
      if (is_star_inconsistent) ++num_inconsistent_stars;
    }

#ifdef DEBUG_TRACES
    std::cerr << "\n==========================================================\n"
              << "Export from list of stars to OFF:\n"
              << "  * Number of vertices: " << m_points.size() << "\n"
              << "  * Total number of maximal simplices: " << num_maximal_simplices << "\n";
    if (color_inconsistencies) {
      std::cerr << "  * Number of inconsistent stars: " << num_inconsistent_stars << " ("
                << (m_points.size() > 0 ? 100. * num_inconsistent_stars / m_points.size() : 0.) << "%)\n"
                << "  * Number of inconsistent maximal simplices: " << num_inconsistent_maximal_simplices << " ("
                << (num_maximal_simplices > 0 ? 100. * num_inconsistent_maximal_simplices / num_maximal_simplices : 0.)
                << "%)\n";
    }
    std::cerr << "==========================================================\n";
#endif

    return os;
  }

 public:
  std::ostream &export_simplices_to_off(const Simplicial_complex &complex, std::ostream &os,
                                        std::size_t &num_OFF_simplices,
                                        Simplex_set const *p_simpl_to_color_in_red = NULL,
                                        Simplex_set const *p_simpl_to_color_in_green = NULL,
                                        Simplex_set const *p_simpl_to_color_in_blue = NULL) const {
    typedef Simplicial_complex::Simplex Simplex;
    typedef Simplicial_complex::Simplex_set Simplex_set;

    // If m_intrinsic_dim = 1, each point is output two times
    // (see export_vertices_to_off)
    num_OFF_simplices = 0;
    std::size_t num_maximal_simplices = 0;

    typename Simplex_set::const_iterator it_s = complex.simplex_range().begin();
    typename Simplex_set::const_iterator it_s_end = complex.simplex_range().end();
    // For each simplex
    for (; it_s != it_s_end; ++it_s) {
      Simplex c = *it_s;
      ++num_maximal_simplices;

      int color_simplex = -1;  // -1=no color, 0=yellow, 1=red, 2=green, 3=blue
      if (p_simpl_to_color_in_red && std::find(p_simpl_to_color_in_red->begin(), p_simpl_to_color_in_red->end(), c) !=
                                         p_simpl_to_color_in_red->end()) {
        color_simplex = 1;
      } else if (p_simpl_to_color_in_green &&
                 std::find(p_simpl_to_color_in_green->begin(), p_simpl_to_color_in_green->end(), c) !=
                     p_simpl_to_color_in_green->end()) {
        color_simplex = 2;
      } else if (p_simpl_to_color_in_blue &&
                 std::find(p_simpl_to_color_in_blue->begin(), p_simpl_to_color_in_blue->end(), c) !=
                     p_simpl_to_color_in_blue->end()) {
        color_simplex = 3;
      }

      // Gather the triangles here
      typedef std::vector<Simplex> Triangles;
      Triangles triangles;

      int num_vertices = static_cast<int>(c.size());
      // Do not export smaller dimension simplices
      if (num_vertices < m_intrinsic_dim + 1) continue;

      // If m_intrinsic_dim = 1, each point is output two times,
      // so we need to multiply each index by 2
      // And if only 2 vertices, add a third one (each vertex is duplicated in
      // the file when m_intrinsic dim = 2)
      if (m_intrinsic_dim == 1) {
        Simplex tmp_c;
        Simplex::iterator it = c.begin();
        for (; it != c.end(); ++it) tmp_c.insert(*it * 2);
        if (num_vertices == 2) tmp_c.insert(*tmp_c.rbegin() + 1);

        c = tmp_c;
      }

      if (num_vertices <= 3) {
        triangles.push_back(c);
      } else {
        // num_vertices >= 4: decompose the simplex in triangles
        std::vector<bool> booleans(num_vertices, false);
        std::fill(booleans.begin() + num_vertices - 3, booleans.end(), true);
        do {
          Simplex triangle;
          Simplex::iterator it = c.begin();
          for (int i = 0; it != c.end(); ++i, ++it) {
            if (booleans[i]) triangle.insert(*it);
          }
          triangles.push_back(triangle);
        } while (std::next_permutation(booleans.begin(), booleans.end()));
      }

      // For each cell
      Triangles::const_iterator it_tri = triangles.begin();
      Triangles::const_iterator it_tri_end = triangles.end();
      for (; it_tri != it_tri_end; ++it_tri) {
        // Don't export infinite cells
        if (is_infinite(*it_tri)) continue;

        os << 3 << " ";
        Simplex::const_iterator it_point_idx = it_tri->begin();
        for (; it_point_idx != it_tri->end(); ++it_point_idx) {
          os << *it_point_idx << " ";
        }

        if (p_simpl_to_color_in_red || p_simpl_to_color_in_green || p_simpl_to_color_in_blue) {
          switch (color_simplex) {
            case 0:
              os << " 255 255 0";
              break;
            case 1:
              os << " 255 0 0";
              break;
            case 2:
              os << " 0 255 0";
              break;
            case 3:
              os << " 0 0 255";
              break;
            default:
              os << " 128 128 128";
              break;
          }
        }

        ++num_OFF_simplices;
        os << "\n";
      }
    }

#ifdef DEBUG_TRACES
    std::cerr << "\n==========================================================\n"
              << "Export from complex to OFF:\n"
              << "  * Number of vertices: " << m_points.size() << "\n"
              << "  * Total number of maximal simplices: " << num_maximal_simplices << "\n"
              << "==========================================================\n";
#endif

    return os;
  }

  /** \brief Sets the maximal possible squared edge length for the edges in the triangulations.
   *
   * @param[in] max_squared_edge_length Maximal possible squared edge length.
   *
   * If the maximal edge length value is too low `Tangential_complex::compute_tangential_complex` will throw an
   * exception in debug mode.
   */
  void set_max_squared_edge_length(FT max_squared_edge_length) { m_max_squared_edge_length = max_squared_edge_length; }

 private:
  const K m_k;
  const int m_intrinsic_dim;
  const int m_ambient_dim;

  Points m_points;
  Weights m_weights;
#ifdef GUDHI_TC_PERTURB_POSITION
  Translations_for_perturb m_translations;
#if defined(GUDHI_USE_TBB)
  Mutex_for_perturb *m_p_perturb_mutexes;
#endif
#endif

  Points_ds m_points_ds;
  double m_last_max_perturb;
  std::vector<bool> m_are_tangent_spaces_computed;
  TS_container m_tangent_spaces;
#ifdef GUDHI_TC_EXPORT_NORMALS
  OS_container m_orth_spaces;
#endif
  Tr_container m_triangulations;  // Contains the triangulations
  // and their center vertex
  Stars_container m_stars;
  std::vector<FT> m_squared_star_spheres_radii_incl_margin;
  boost::optional<FT> m_max_squared_edge_length;

#ifdef GUDHI_TC_USE_ANOTHER_POINT_SET_FOR_TANGENT_SPACE_ESTIM
  Points m_points_for_tse;
  Points_ds m_points_ds_for_tse;
#endif

  mutable CGAL::Random m_random_generator;
};  // /class Tangential_complex

}  // end namespace tangential_complex
}  // end namespace Gudhi

#endif  // TANGENTIAL_COMPLEX_H_

/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2015 Inria
 *
 *    Modification(s):
 *      - 2019/08 Vincent Rouvreau: Fix issue #10 for CGAL and Eigen3
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef ALPHA_COMPLEX_H_
#define ALPHA_COMPLEX_H_

#include <gudhi/Alpha_complex/Alpha_kernel_d.h>
#include <gudhi/Debug_utils.h>
// to construct Alpha_complex from a OFF file of points
#include <gudhi/Points_off_io.h>

#include <cmath>  // isnan, fmax
#include <memory>  // for std::unique_ptr
#include <cstddef>  // for std::size_t

#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Regular_triangulation.h>  // aka. Weighted Delaunay triangulation
#include <CGAL/Epeck_d.h>  // For EXACT or SAFE version
#include <CGAL/Epick_d.h>  // For FAST version
#include <CGAL/Spatial_sort_traits_adapter_d.h>
#include <CGAL/property_map.h>  // for CGAL::Identity_property_map
#include <CGAL/version.h>  // for CGAL_VERSION_NR
#include <CGAL/NT_converter.h>

#include <Eigen/src/Core/util/Macros.h>  // for EIGEN_VERSION_AT_LEAST

#include <boost/range/size.hpp>
#include <boost/range/combine.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include <iostream>
#include <vector>
#include <string>
#include <limits>  // NaN
#include <map>
#include <utility>  // std::pair
#include <stdexcept>
#include <numeric>  // for std::iota
#include <algorithm>  // for std::sort

// Make compilation fail - required for external projects - https://github.com/GUDHI/gudhi-devel/issues/10
#if CGAL_VERSION_NR < 1041101000
# error Alpha_complex is only available for CGAL >= 4.11
#endif

#if !EIGEN_VERSION_AT_LEAST(3,1,0)
# error Alpha_complex is only available for Eigen3 >= 3.1.0 installed with CGAL
#endif

namespace Gudhi {

namespace alpha_complex {

template<typename D> struct Is_Epeck_D { static const bool value = false; };
template<typename D> struct Is_Epeck_D<CGAL::Epeck_d<D>> { static const bool value = true; };

/**
 * \class Alpha_complex Alpha_complex.h gudhi/Alpha_complex.h
 * \brief Alpha complex data structure.
 * 
 * \ingroup alpha_complex
 * 
 * \details
 * The data structure is constructing a CGAL Delaunay triangulation (for more information on CGAL Delaunay 
 * triangulation, please refer to the corresponding chapter in page http://doc.cgal.org/latest/Triangulation/) from a
 * range of points or from an OFF file (cf. Points_off_reader).
 * 
 * Please refer to \ref alpha_complex for examples.
 *
 * The complex is a template class requiring an <a target="_blank"
 * href="https://doc.cgal.org/latest/Kernel_d/structCGAL_1_1Epeck__d.html">CGAL::Epeck_d</a>,
 * or an <a target="_blank"
 * href="https://doc.cgal.org/latest/Kernel_d/structCGAL_1_1Epick__d.html">CGAL::Epick_d</a> <a target="_blank"
 * href="http://doc.cgal.org/latest/Kernel_d/index.html#Chapter_dD_Geometry_Kernel">dD Geometry Kernel</a>
 * \cite cgal:s-gkd-19b from CGAL as template, default value is <a target="_blank"
 * href="https://doc.cgal.org/latest/Kernel_d/structCGAL_1_1Epeck__d.html">CGAL::Epeck_d</a>
 * < <a target="_blank" href="http://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Dynamic__dimension__tag.html">
 * CGAL::Dynamic_dimension_tag </a> >
 * 
 * \remark
 * - When Alpha_complex is constructed with an infinite value of alpha, the complex is a Delaunay complex.
 * - Using the default `CGAL::Epeck_d` makes the construction safe. If you pass exact=true to create_complex, the
 * filtration values are the exact ones converted to the filtration value type of the simplicial complex. This can be
 * very slow. If you pass exact=false (the default), the filtration values are only guaranteed to have a small
 * multiplicative error compared to the exact value, see <code><a class="el" target="_blank"
 * href="https://doc.cgal.org/latest/Number_types/classCGAL_1_1Lazy__exact__nt.html">
 * CGAL::Lazy_exact_nt<NT>::set_relative_precision_of_to_double</a></code> for details. A drawback, when computing
 * persistence, is that an empty exact interval [10^12,10^12] may become a non-empty approximate interval
 * [10^12,10^12+10^6]. Using `CGAL::Epick_d` makes the computations slightly faster, and the combinatorics are still
 * exact, but the computation of filtration values can exceptionally be arbitrarily bad. In all cases, we still
 * guarantee that the output is a valid filtration (faces have a filtration value no larger than their cofaces).
 * - For performances reasons, it is advised to use `Alpha_complex` with \ref cgal &ge; 5.0.0.
 */
template<class Kernel = CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>, bool Weighted = false>
class Alpha_complex {
 private:
  // Vertex_handle internal type (required by triangulation_ and vertices_).
  using Internal_vertex_handle = std::ptrdiff_t;

 public:
  /** \brief Geometric traits class that provides the geometric types and predicates needed by the triangulations.*/
  using Geom_traits = std::conditional_t<Weighted, CGAL::Regular_triangulation_traits_adapter<Kernel>, Kernel>;

  // Add an int in TDS to save point index in the structure
  using TDS = CGAL::Triangulation_data_structure<typename Geom_traits::Dimension,
                                                 CGAL::Triangulation_vertex<Geom_traits, Internal_vertex_handle>,
                                                 CGAL::Triangulation_full_cell<Geom_traits> >;

  /** \brief A (Weighted or not) Delaunay triangulation of a set of points in \f$ \mathbb{R}^D\f$.*/
  using Triangulation = std::conditional_t<Weighted, CGAL::Regular_triangulation<Kernel, TDS>,
                                                     CGAL::Delaunay_triangulation<Kernel, TDS>>;

  /** \brief CGAL kernel container for computations in function of the weighted or not characteristics.*/
  using A_kernel_d = Alpha_kernel_d<Kernel, Weighted>;

  // Numeric type of coordinates in the kernel
  using FT = typename A_kernel_d::FT;

  /** \brief Sphere is a std::pair<Kernel::Point_d, Kernel::FT> (aka. circurmcenter and squared radius).
   * If Weighted, Sphere is a Kernel::Weighted_point_d (aka. circurmcenter and the weight value is the squared radius).
  */
  using Sphere = typename A_kernel_d::Sphere;

  /** \brief A point, or a weighted point in Euclidean space.*/
  using Point_d = typename Geom_traits::Point_d;

 private:
  // Vertex_iterator type from CGAL.
  using CGAL_vertex_iterator = typename Triangulation::Vertex_iterator;

  // Structure to switch from simplex tree vertex handle to CGAL vertex iterator.
  using Vector_vertex_iterator = std::vector< CGAL_vertex_iterator >;

 private:
  /** \brief Vertex iterator vector to switch from simplex tree vertex handle to CGAL vertex iterator.
   * Vertex handles are inserted sequentially, starting at 0.*/
  Vector_vertex_iterator vertex_handle_to_iterator_;
  /** \brief Pointer on the CGAL Delaunay triangulation.*/
  std::unique_ptr<Triangulation> triangulation_;
  /** \brief Kernel for triangulation_ functions access.*/
  A_kernel_d kernel_;
  /** \brief Vertices to be inserted first by the create_complex method to avoid quadratic complexity.
   * It isn't just [0, n) if some points have multiplicity (only one copy appears in the complex).
  */
  std::vector<Internal_vertex_handle> vertices_;

  /** \brief Cache for geometric constructions: circumcenter and squared radius of a simplex.*/
  std::vector<Sphere> cache_, old_cache_;

 public:
  /** \brief Alpha_complex constructor from an OFF file name.
   * 
   * Uses the Points_off_reader to construct the Delaunay triangulation required to initialize 
   * the Alpha_complex.
   * 
   * Duplicate points are inserted once in the Alpha_complex. This is the reason why the vertices may be not contiguous.
   *
   * @param[in] off_file_name OFF file [path and] name.
   */
  Alpha_complex(const std::string& off_file_name) {
    Gudhi::Points_off_reader<Point_d> off_reader(off_file_name);
    if (!off_reader.is_valid()) {
      std::cerr << "Alpha_complex - Unable to read file " << off_file_name << "\n";
      exit(-1);  // ----- >>
    }

    init_from_range(off_reader.get_point_cloud());
  }

  /** \brief Alpha_complex constructor from a list of points.
   *
   * The vertices may be not contiguous as some points may be discarded in the triangulation (duplicate points,
   * weighted hidden point, ...).
   * 
   * @param[in] points Range of points to triangulate. Points must be in Kernel::Point_d or Kernel::Weighted_point_d.
   * 
   * The type InputPointRange must be a range for which std::begin and std::end return input iterators on a
   * Kernel::Point_d or Kernel::Weighted_point_d.
   */
  template<typename InputPointRange >
  Alpha_complex(const InputPointRange& points) {
    init_from_range(points);
  }

  /** \brief Alpha_complex constructor from a list of points and weights.
   *
   * The vertices may be not contiguous as some points may be discarded in the triangulation (duplicate points,
   * weighted hidden point, ...).
   * 
   * @param[in] points Range of points to triangulate. Points must be in Kernel::Point_d.
   * 
   * @param[in] weights Range of points weights. Weights must be in Kernel::FT.
   * 
   * The type InputPointRange must be a range for which std::begin and std::end return input iterators on a
   * Kernel::Point_d.
   */
  template <typename InputPointRange, typename WeightRange>
  Alpha_complex(const InputPointRange& points, WeightRange weights) {
    static_assert(Weighted, "This constructor is not available for non-weighted versions of Alpha_complex");
    // FIXME: this test is only valid if we have a forward range
    GUDHI_CHECK(boost::size(weights) == boost::size(points),
                std::invalid_argument("Points number in range different from weights range number"));
    auto weighted_points = boost::range::combine(points, weights)
      | boost::adaptors::transformed([](auto const&t){return Point_d(boost::get<0>(t), boost::get<1>(t));});
    init_from_range(weighted_points);
  }

  // Forbid copy/move constructor/assignment operator
  Alpha_complex(const Alpha_complex& other) = delete;
  Alpha_complex& operator= (const Alpha_complex& other) = delete;
  Alpha_complex (Alpha_complex&& other) = delete;
  Alpha_complex& operator= (Alpha_complex&& other) = delete;

  /** \brief Returns the number of finite vertices in the triangulation.
   */
  std::size_t num_vertices() const {
    if (triangulation_ == nullptr)
      return 0;
    else
      return triangulation_->number_of_vertices();
  }

  /** \brief get_point returns the point corresponding to the vertex given as parameter.
   *
   * @param[in] vertex Vertex handle of the point to retrieve.
   * @return The point found.
   * @exception std::out_of_range In case vertex is not found (cf. std::vector::at).
   */
  const Point_d& get_point(std::size_t vertex) const {
    auto it = vertex_handle_to_iterator_.at(vertex);
    if (it == nullptr) throw std::out_of_range("This vertex is missing, maybe hidden by a duplicate or another heavier point.");
    return it->point();
  }

 private:
  template<typename InputPointRange >
  void init_from_range(const InputPointRange& points) {
    #if CGAL_VERSION_NR < 1050000000
    if (Is_Epeck_D<Kernel>::value)
      std::cerr << "It is strongly advised to use a CGAL version >= 5.0 with Epeck_d Kernel for performance reasons."
                << std::endl;
    #endif

#if CGAL_VERSION_NR < 1050101000
    // Make compilation fail if weighted and CGAL < 5.1
    static_assert(!Weighted, "Weighted Alpha_complex is only available for CGAL >= 5.1");
#endif

    auto first = std::begin(points);
    auto last = std::end(points);

    if (first != last) {
      // Delaunay triangulation init with point dimension.
      triangulation_ = std::make_unique<Triangulation>(kernel_.get_dimension(*first));

      std::vector<Point_d> point_cloud(first, last);

      // Creates a vector {0, 1, ..., N-1}
      std::vector<Internal_vertex_handle> indices(boost::counting_iterator<Internal_vertex_handle>(0),
                                                  boost::counting_iterator<Internal_vertex_handle>(point_cloud.size()));

      using Point_property_map = boost::iterator_property_map<typename std::vector<Point_d>::iterator,
                                                              CGAL::Identity_property_map<Internal_vertex_handle>>;
      using Search_traits_d = CGAL::Spatial_sort_traits_adapter_d<Geom_traits, Point_property_map>;

      CGAL::spatial_sort(indices.begin(), indices.end(), Search_traits_d(std::begin(point_cloud)));

      typename Triangulation::Full_cell_handle hint;
      for (auto index : indices) {
        typename Triangulation::Vertex_handle pos = triangulation_->insert(point_cloud[index], hint);
        if (pos != nullptr) {
          // Save index value as data to retrieve it after insertion
          pos->data() = index;
          hint = pos->full_cell();
        }
      }
      // --------------------------------------------------------------------------------------------
      // structure to retrieve CGAL points from vertex handle - one vertex handle per point.
      // Needs to be constructed before as vertex handles arrives in no particular order.
      vertex_handle_to_iterator_.resize(point_cloud.size());
      // List of sorted unique vertices in the triangulation. We take advantage of the existing loop to construct it
      // Vertices list avoids quadratic complexity with the Simplex_tree. We should not fill it up with Toplex_map e.g.
      vertices_.reserve(triangulation_->number_of_vertices());
      // Loop on triangulation vertices list
      for (CGAL_vertex_iterator vit = triangulation_->vertices_begin(); vit != triangulation_->vertices_end(); ++vit) {
        if (!triangulation_->is_infinite(*vit)) {
#ifdef DEBUG_TRACES
          std::clog << "Vertex insertion - " << vit->data() << " -> " << vit->point() << std::endl;
#endif  // DEBUG_TRACES
          vertex_handle_to_iterator_[vit->data()] = vit;
          vertices_.push_back(vit->data());
        }
      }
      std::sort(vertices_.begin(), vertices_.end());
      // --------------------------------------------------------------------------------------------
    }
  }

  /** \brief get_point_ returns the point corresponding to the vertex given as parameter.
   * Only for internal use for faster access.
   *
   * @param[in] vertex Vertex handle of the point to retrieve.
   * @return The point found.
   */
  const Point_d& get_point_(std::size_t vertex) const {
    return vertex_handle_to_iterator_[vertex]->point();
  }

  /// Return a reference to the circumcenter and circumradius, writing them in the cache if necessary.
  template<class SimplicialComplexForAlpha>
  auto& get_cache(SimplicialComplexForAlpha& cplx, typename SimplicialComplexForAlpha::Simplex_handle s) {
    auto k = cplx.key(s);
    if(k==cplx.null_key()){
      k = cache_.size();
      cplx.assign_key(s, k);
      // Using a transform_range is slower, currently.
      thread_local std::vector<Point_d> v;
      v.clear();
      for (auto vertex : cplx.simplex_vertex_range(s))
        v.push_back(get_point_(vertex));
      cache_.emplace_back(kernel_.get_sphere(v.cbegin(), v.cend()));
    }
    return cache_[k];
  }

  /// Return the circumradius, either from the old cache or computed, without writing to the cache.
  template<class SimplicialComplexForAlpha>
  auto radius(SimplicialComplexForAlpha& cplx, typename SimplicialComplexForAlpha::Simplex_handle s) {
    auto k = cplx.key(s);
    if(k!=cplx.null_key())
      return kernel_.get_squared_radius(old_cache_[k]);
    // Using a transform_range is slower, currently.
    thread_local std::vector<Point_d> v;
    v.clear();
    for (auto vertex : cplx.simplex_vertex_range(s))
      v.push_back(get_point_(vertex));
    return kernel_.get_squared_radius(v.cbegin(), v.cend());
  }

 public:
  /** \brief Inserts all Delaunay triangulation into the simplicial complex.
   * It also computes the filtration values accordingly to the \ref createcomplexalgorithm if default_filtration_value
   * is not set.
   *
   * \tparam SimplicialComplexForAlpha must meet `SimplicialComplexForAlpha` concept.
   * 
   * @param[in] complex SimplicialComplexForAlpha to be created.
   * @param[in] max_alpha_square maximum for alpha square value. Default value is +\f$\infty\f$, and there is very
   * little point using anything else since it does not save time. Useless if `default_filtration_value` is set to
   * `true`.
   * @param[in] exact Exact filtration values computation. Not exact if `Kernel` is not <a target="_blank"
   * href="https://doc.cgal.org/latest/Kernel_d/structCGAL_1_1Epeck__d.html">CGAL::Epeck_d</a>.
   * @param[in] default_filtration_value Set this value to `true` if filtration values are not needed to be computed
   * (will be set to `NaN`).
   * Default value is `false` (which means compute the filtration values).
   *
   * @return true if creation succeeds, false otherwise.
   * 
   * @pre Delaunay triangulation must be already constructed with dimension strictly greater than 0.
   * @pre The simplicial complex must be empty (no vertices)
   * 
   * Initialization can be launched once.
   */
  template <typename SimplicialComplexForAlpha,
            typename Filtration_value = typename SimplicialComplexForAlpha::Filtration_value>
  bool create_complex(SimplicialComplexForAlpha& complex,
                      Filtration_value max_alpha_square = std::numeric_limits<Filtration_value>::infinity(),
                      bool exact = false,
                      bool default_filtration_value = false) {
    // Filtration_value must be capable to represent the special value "Not-A-Number"
    static_assert(std::numeric_limits<Filtration_value>::has_quiet_NaN);
    // To support more general types for Filtration_value
    using std::isnan;

    // From SimplicialComplexForAlpha type required to insert into a simplicial complex (with or without subfaces).
    using Vertex_handle = typename SimplicialComplexForAlpha::Vertex_handle;
    using Simplex_handle = typename SimplicialComplexForAlpha::Simplex_handle;
    using Vector_vertex = std::vector<Vertex_handle>;

    if (triangulation_ == nullptr) {
      std::cerr << "Alpha_complex cannot create_complex from a NULL triangulation\n";
      return false;  // ----- >>
    }
    if (triangulation_->maximal_dimension() < 1) {
      std::cerr << "Alpha_complex cannot create_complex from a zero-dimension triangulation\n";
      return false;  // ----- >>
    }
    if (complex.num_vertices() > 0) {
      std::cerr << "Alpha_complex create_complex - complex is not empty\n";
      return false;  // ----- >>
    }

    // --------------------------------------------------------------------------------------------
    // Simplex_tree construction from loop on triangulation finite full cells list
    if (num_vertices() > 0) {
      std::vector<Vertex_handle> one_vertex(1);
      for (auto vertex : vertices_) {
#ifdef DEBUG_TRACES
        std::clog << "SimplicialComplex insertion " << vertex << std::endl;
#endif  // DEBUG_TRACES
        one_vertex[0] = vertex;
        complex.insert_simplex_and_subfaces(one_vertex, std::numeric_limits<Filtration_value>::quiet_NaN());
      }

      for (auto cit = triangulation_->finite_full_cells_begin();
           cit != triangulation_->finite_full_cells_end();
           ++cit) {
        Vector_vertex vertexVector;
#ifdef DEBUG_TRACES
        std::clog << "SimplicialComplex insertion ";
#endif  // DEBUG_TRACES
        for (auto vit = cit->vertices_begin(); vit != cit->vertices_end(); ++vit) {
          if (*vit != nullptr) {
#ifdef DEBUG_TRACES
            std::clog << " " << (*vit)->data();
#endif  // DEBUG_TRACES
            // Vector of vertex construction for simplex_tree structure
            vertexVector.push_back((*vit)->data());
          }
        }
#ifdef DEBUG_TRACES
        std::clog << std::endl;
#endif  // DEBUG_TRACES
        // Insert each simplex and its subfaces in the simplex tree - filtration is NaN
        complex.insert_simplex_and_subfaces(vertexVector, std::numeric_limits<Filtration_value>::quiet_NaN());
      }
    }
    // --------------------------------------------------------------------------------------------

    if (!default_filtration_value) {
      CGAL::NT_converter<FT, Filtration_value> cgal_converter;
      // --------------------------------------------------------------------------------------------
      // ### For i : d -> 0
      for (int decr_dim = triangulation_->maximal_dimension(); decr_dim >= 0; decr_dim--) {
        // ### Foreach Sigma of dim i
        for (Simplex_handle f_simplex : complex.skeleton_simplex_range(decr_dim)) {
          int f_simplex_dim = complex.dimension(f_simplex);
          if (decr_dim == f_simplex_dim) {
            // ### If filt(Sigma) is NaN : filt(Sigma) = alpha(Sigma)
            if (isnan(complex.filtration(f_simplex))) {
              Filtration_value alpha_complex_filtration = 0.0;
              // No need to compute squared_radius on a non-weighted single point - alpha is 0.0
              if (Weighted || f_simplex_dim > 0) {
                auto const& sqrad = radius(complex, f_simplex);
#if CGAL_VERSION_NR >= 1050000000
                if(exact) CGAL::exact(sqrad);
#endif
                alpha_complex_filtration = cgal_converter(sqrad);
              }
              complex.assign_filtration(f_simplex, alpha_complex_filtration);
#ifdef DEBUG_TRACES
              std::clog << "filt(Sigma) is NaN : filt(Sigma) =" << complex.filtration(f_simplex) << std::endl;
#endif  // DEBUG_TRACES
            }
            // No need to propagate further, unweighted points all have value 0
            if (decr_dim > !Weighted)
              propagate_alpha_filtration(complex, f_simplex);
          }
        }
        old_cache_ = std::move(cache_);
        cache_.clear();
      }
      // --------------------------------------------------------------------------------------------
  
      // --------------------------------------------------------------------------------------------
      if (!exact)
        // As Alpha value is an approximation, we have to make filtration non decreasing while increasing the dimension
        // Only in not exact version, cf. https://github.com/GUDHI/gudhi-devel/issues/57
        complex.make_filtration_non_decreasing();
      // Remove all simplices that have a filtration value greater than max_alpha_square
      complex.prune_above_filtration(max_alpha_square);
      // --------------------------------------------------------------------------------------------
    }
    return true;
  }

 private:
  template <typename SimplicialComplexForAlpha, typename Simplex_handle>
  void propagate_alpha_filtration(SimplicialComplexForAlpha& complex, Simplex_handle f_simplex) {
    // From SimplicialComplexForAlpha type required to assign filtration values.
    using Filtration_value = typename SimplicialComplexForAlpha::Filtration_value;
    // To support more general types for Filtration_value
    using std::isnan;

    // ### Foreach Tau face of Sigma
    for (auto face_opposite_vertex : complex.boundary_opposite_vertex_simplex_range(f_simplex)) {
      auto f_boundary = face_opposite_vertex.first;
#ifdef DEBUG_TRACES
      std::clog << " | --------------------------------------------------\n";
      std::clog << " | Tau ";
      for (auto vertex : complex.simplex_vertex_range(f_boundary)) {
        std::clog << vertex << " ";
      }
      std::clog << "is a face of Sigma\n";
      std::clog << " | isnan(complex.filtration(Tau)=" << isnan(complex.filtration(f_boundary)) << std::endl;
#endif  // DEBUG_TRACES
      // ### If filt(Tau) is not NaN
      if (!isnan(complex.filtration(f_boundary))) {
        // ### filt(Tau) = fmin(filt(Tau), filt(Sigma))
        Filtration_value alpha_complex_filtration = fmin(complex.filtration(f_boundary),
                                                                             complex.filtration(f_simplex));
        complex.assign_filtration(f_boundary, alpha_complex_filtration);
#ifdef DEBUG_TRACES
        std::clog << " | filt(Tau) = fmin(filt(Tau), filt(Sigma)) = " << complex.filtration(f_boundary) << std::endl;
#endif  // DEBUG_TRACES
        // ### Else
      } else {
        auto const& cache=get_cache(complex, f_boundary);
        bool is_gab = kernel_.is_gabriel(cache, get_point_(face_opposite_vertex.second));
#ifdef DEBUG_TRACES
        std::clog << " | Tau is_gabriel(Sigma)=" << is_gab << " - vertexForGabriel=" << face_opposite_vertex.second << std::endl;
#endif  // DEBUG_TRACES
        // ### If Tau is not Gabriel of Sigma
        if (false == is_gab) {
          // ### filt(Tau) = filt(Sigma)
          Filtration_value alpha_complex_filtration = complex.filtration(f_simplex);
          complex.assign_filtration(f_boundary, alpha_complex_filtration);
#ifdef DEBUG_TRACES
          std::clog << " | filt(Tau) = filt(Sigma) = " << complex.filtration(f_boundary) << std::endl;
#endif  // DEBUG_TRACES
        }
      }
    }
  }
};

}  // namespace alpha_complex

namespace alphacomplex = alpha_complex;

}  // namespace Gudhi

#endif  // ALPHA_COMPLEX_H_

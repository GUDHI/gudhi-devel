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

#include <gudhi/Debug_utils.h>
// to construct Alpha_complex from a OFF file of points
#include <gudhi/Points_off_io.h>

#include <stdlib.h>
#include <math.h>  // isnan, fmax

#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Epeck_d.h>  // For EXACT or SAFE version
#include <CGAL/Epick_d.h>  // For FAST version
#include <CGAL/Spatial_sort_traits_adapter_d.h>
#include <CGAL/property_map.h>  // for CGAL::Identity_property_map
#include <CGAL/version.h>  // for CGAL_VERSION_NR
#include <CGAL/NT_converter.h>

#include <Eigen/src/Core/util/Macros.h>  // for EIGEN_VERSION_AT_LEAST

#include <iostream>
#include <vector>
#include <string>
#include <limits>  // NaN
#include <map>
#include <utility>  // std::pair
#include <stdexcept>
#include <numeric>  // for std::iota

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
 * The data structure is constructing a CGAL Delaunay triangulation (for more informations on CGAL Delaunay 
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
template<class Kernel = CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>>
class Alpha_complex {
 public:
  // Add an int in TDS to save point index in the structure
  typedef CGAL::Triangulation_data_structure<typename Kernel::Dimension,
                              CGAL::Triangulation_vertex<Kernel, std::ptrdiff_t>,
                              CGAL::Triangulation_full_cell<Kernel> > TDS;
  /** \brief A Delaunay triangulation of a set of points in \f$ \mathbb{R}^D\f$.*/
  typedef CGAL::Delaunay_triangulation<Kernel, TDS> Delaunay_triangulation;

  /** \brief A point in Euclidean space.*/
  typedef typename Kernel::Point_d Point_d;
  /** \brief Geometric traits class that provides the geometric types and predicates needed by Delaunay
   * triangulations.*/
  typedef Kernel Geom_traits;

 private:
  typedef typename Kernel::Compute_squared_radius_d Squared_Radius;
  typedef typename Kernel::Side_of_bounded_sphere_d Is_Gabriel;
  typedef typename Kernel::Point_dimension_d        Point_Dimension;

  // Type required to compute squared radius, or side of bounded sphere on a vector of points.
  typedef typename std::vector<Point_d> Vector_of_CGAL_points;

  // Vertex_iterator type from CGAL.
  typedef typename Delaunay_triangulation::Vertex_iterator CGAL_vertex_iterator;

  // size_type type from CGAL.
  typedef typename Delaunay_triangulation::size_type size_type;

  // Structure to switch from simplex tree vertex handle to CGAL vertex iterator.
  typedef typename std::vector< CGAL_vertex_iterator > Vector_vertex_iterator;

 private:
  /** \brief Vertex iterator vector to switch from simplex tree vertex handle to CGAL vertex iterator.
   * Vertex handles are inserted sequentially, starting at 0.*/
  Vector_vertex_iterator vertex_handle_to_iterator_;
  /** \brief Pointer on the CGAL Delaunay triangulation.*/
  Delaunay_triangulation* triangulation_;
  /** \brief Kernel for triangulation_ functions access.*/
  Kernel kernel_;

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
  Alpha_complex(const std::string& off_file_name)
      : triangulation_(nullptr) {
    Gudhi::Points_off_reader<Point_d> off_reader(off_file_name);
    if (!off_reader.is_valid()) {
      std::cerr << "Alpha_complex - Unable to read file " << off_file_name << "\n";
      exit(-1);  // ----- >>
    }

    init_from_range(off_reader.get_point_cloud());
  }

  /** \brief Alpha_complex constructor from a list of points.
   *
   * Duplicate points are inserted once in the Alpha_complex. This is the reason why the vertices may be not contiguous.
   * 
   * @param[in] points Range of points to triangulate. Points must be in Kernel::Point_d
   * 
   * The type InputPointRange must be a range for which std::begin and
   * std::end return input iterators on a Kernel::Point_d.
   */
  template<typename InputPointRange >
  Alpha_complex(const InputPointRange& points)
      : triangulation_(nullptr) {
    init_from_range(points);
  }

  /** \brief Alpha_complex destructor deletes the Delaunay triangulation.
   */
  ~Alpha_complex() {
    delete triangulation_;
  }

  // Forbid copy/move constructor/assignment operator
  Alpha_complex(const Alpha_complex& other) = delete;
  Alpha_complex& operator= (const Alpha_complex& other) = delete;
  Alpha_complex (Alpha_complex&& other) = delete;
  Alpha_complex& operator= (Alpha_complex&& other) = delete;

  /** \brief get_point returns the point corresponding to the vertex given as parameter.
   *
   * @param[in] vertex Vertex handle of the point to retrieve.
   * @return The point found.
   * @exception std::out_of_range In case vertex is not found (cf. std::vector::at).
   */
  const Point_d& get_point(std::size_t vertex) const {
    return vertex_handle_to_iterator_.at(vertex)->point();
  }

 private:
  template<typename InputPointRange >
  void init_from_range(const InputPointRange& points) {
    #if CGAL_VERSION_NR < 1050000000
    if (Is_Epeck_D<Kernel>::value)
      std::cerr << "It is strongly advised to use a CGAL version >= 5.0 with Epeck_d Kernel for performance reasons."
                << std::endl;
    #endif

    auto first = std::begin(points);
    auto last = std::end(points);

    if (first != last) {
      // point_dimension function initialization
      Point_Dimension point_dimension = kernel_.point_dimension_d_object();

      // Delaunay triangulation is point dimension.
      triangulation_ = new Delaunay_triangulation(point_dimension(*first));

      std::vector<Point_d> point_cloud(first, last);

      // Creates a vector {0, 1, ..., N-1}
      std::vector<std::ptrdiff_t> indices(boost::counting_iterator<std::ptrdiff_t>(0),
                                          boost::counting_iterator<std::ptrdiff_t>(point_cloud.size()));

      typedef boost::iterator_property_map<typename std::vector<Point_d>::iterator,
                                           CGAL::Identity_property_map<std::ptrdiff_t>> Point_property_map;
      typedef CGAL::Spatial_sort_traits_adapter_d<Kernel, Point_property_map> Search_traits_d;

      CGAL::spatial_sort(indices.begin(), indices.end(), Search_traits_d(std::begin(point_cloud)));

      typename Delaunay_triangulation::Full_cell_handle hint;
      for (auto index : indices) {
        typename Delaunay_triangulation::Vertex_handle pos = triangulation_->insert(point_cloud[index], hint);
        // Save index value as data to retrieve it after insertion
        pos->data() = index;
        hint = pos->full_cell();
      }
      // --------------------------------------------------------------------------------------------
      // structure to retrieve CGAL points from vertex handle - one vertex handle per point.
      // Needs to be constructed before as vertex handles arrives in no particular order.
      vertex_handle_to_iterator_.resize(point_cloud.size());
      // Loop on triangulation vertices list
      for (CGAL_vertex_iterator vit = triangulation_->vertices_begin(); vit != triangulation_->vertices_end(); ++vit) {
        if (!triangulation_->is_infinite(*vit)) {
#ifdef DEBUG_TRACES
          std::clog << "Vertex insertion - " << vit->data() << " -> " << vit->point() << std::endl;
#endif  // DEBUG_TRACES
          vertex_handle_to_iterator_[vit->data()] = vit;
        }
      }
      // --------------------------------------------------------------------------------------------
    }
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
   * @param[in] default_filtration_value Set this value to `true` if filtration values are not needed to be computed.
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
    // From SimplicialComplexForAlpha type required to insert into a simplicial complex (with or without subfaces).
    typedef typename SimplicialComplexForAlpha::Vertex_handle Vertex_handle;
    typedef typename SimplicialComplexForAlpha::Simplex_handle Simplex_handle;
    typedef std::vector<Vertex_handle> Vector_vertex;

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
    if (triangulation_->number_of_vertices() > 0) {
      for (auto cit = triangulation_->finite_full_cells_begin();
           cit != triangulation_->finite_full_cells_end();
           ++cit) {
        Vector_vertex vertexVector;
#ifdef DEBUG_TRACES
        std::clog << "Simplex_tree insertion ";
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
        complex.insert_simplex_and_subfaces(vertexVector, std::numeric_limits<double>::quiet_NaN());
      }
    }
    // --------------------------------------------------------------------------------------------

    if (!default_filtration_value) {
      // --------------------------------------------------------------------------------------------
      // Will be re-used many times
      Vector_of_CGAL_points pointVector;
      // ### For i : d -> 0
      for (int decr_dim = triangulation_->maximal_dimension(); decr_dim >= 0; decr_dim--) {
        // ### Foreach Sigma of dim i
        for (Simplex_handle f_simplex : complex.skeleton_simplex_range(decr_dim)) {
          int f_simplex_dim = complex.dimension(f_simplex);
          if (decr_dim == f_simplex_dim) {
            pointVector.clear();
  #ifdef DEBUG_TRACES
            std::clog << "Sigma of dim " << decr_dim << " is";
  #endif  // DEBUG_TRACES
            for (auto vertex : complex.simplex_vertex_range(f_simplex)) {
              pointVector.push_back(get_point(vertex));
  #ifdef DEBUG_TRACES
              std::clog << " " << vertex;
  #endif  // DEBUG_TRACES
            }
  #ifdef DEBUG_TRACES
            std::clog << std::endl;
  #endif  // DEBUG_TRACES
            // ### If filt(Sigma) is NaN : filt(Sigma) = alpha(Sigma)
            if (std::isnan(complex.filtration(f_simplex))) {
              Filtration_value alpha_complex_filtration = 0.0;
              // No need to compute squared_radius on a single point - alpha is 0.0
              if (f_simplex_dim > 0) {
                // squared_radius function initialization
                Squared_Radius squared_radius = kernel_.compute_squared_radius_d_object();
  
                CGAL::NT_converter<typename Geom_traits::FT, Filtration_value> cv;
                auto sqrad = squared_radius(pointVector.begin(), pointVector.end());
  #if CGAL_VERSION_NR >= 1050000000
                if(exact) CGAL::exact(sqrad);
  #endif
                alpha_complex_filtration = cv(sqrad);
              }
              complex.assign_filtration(f_simplex, alpha_complex_filtration);
  #ifdef DEBUG_TRACES
              std::clog << "filt(Sigma) is NaN : filt(Sigma) =" << complex.filtration(f_simplex) << std::endl;
  #endif  // DEBUG_TRACES
            }
            // No need to propagate further, unweighted points all have value 0
            if (decr_dim > 1)
              propagate_alpha_filtration(complex, f_simplex);
          }
        }
      }
      // --------------------------------------------------------------------------------------------
  
      // --------------------------------------------------------------------------------------------
      // As Alpha value is an approximation, we have to make filtration non decreasing while increasing the dimension
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
    typedef typename SimplicialComplexForAlpha::Filtration_value Filtration_value;
#ifdef DEBUG_TRACES
    typedef typename SimplicialComplexForAlpha::Vertex_handle Vertex_handle;
#endif  // DEBUG_TRACES

    // ### Foreach Tau face of Sigma
    for (auto f_boundary : complex.boundary_simplex_range(f_simplex)) {
#ifdef DEBUG_TRACES
      std::clog << " | --------------------------------------------------\n";
      std::clog << " | Tau ";
      for (auto vertex : complex.simplex_vertex_range(f_boundary)) {
        std::clog << vertex << " ";
      }
      std::clog << "is a face of Sigma\n";
      std::clog << " | isnan(complex.filtration(Tau)=" << std::isnan(complex.filtration(f_boundary)) << std::endl;
#endif  // DEBUG_TRACES
      // ### If filt(Tau) is not NaN
      if (!std::isnan(complex.filtration(f_boundary))) {
        // ### filt(Tau) = fmin(filt(Tau), filt(Sigma))
        Filtration_value alpha_complex_filtration = fmin(complex.filtration(f_boundary),
                                                                             complex.filtration(f_simplex));
        complex.assign_filtration(f_boundary, alpha_complex_filtration);
#ifdef DEBUG_TRACES
        std::clog << " | filt(Tau) = fmin(filt(Tau), filt(Sigma)) = " << complex.filtration(f_boundary) << std::endl;
#endif  // DEBUG_TRACES
        // ### Else
      } else {
        // insert the Tau points in a vector for is_gabriel function
        Vector_of_CGAL_points pointVector;
#ifdef DEBUG_TRACES
        Vertex_handle vertexForGabriel = Vertex_handle();
#endif  // DEBUG_TRACES
        for (auto vertex : complex.simplex_vertex_range(f_boundary)) {
          pointVector.push_back(get_point(vertex));
        }
        // Retrieve the Sigma point that is not part of Tau - parameter for is_gabriel function
        Point_d point_for_gabriel;
        for (auto vertex : complex.simplex_vertex_range(f_simplex)) {
          point_for_gabriel = get_point(vertex);
          if (std::find(pointVector.begin(), pointVector.end(), point_for_gabriel) == pointVector.end()) {
#ifdef DEBUG_TRACES
            // vertex is not found in Tau
            vertexForGabriel = vertex;
#endif  // DEBUG_TRACES
            // No need to continue loop
            break;
          }
        }
        // is_gabriel function initialization
        Is_Gabriel is_gabriel = kernel_.side_of_bounded_sphere_d_object();
        bool is_gab = is_gabriel(pointVector.begin(), pointVector.end(), point_for_gabriel)
          != CGAL::ON_BOUNDED_SIDE;
#ifdef DEBUG_TRACES
        std::clog << " | Tau is_gabriel(Sigma)=" << is_gab << " - vertexForGabriel=" << vertexForGabriel << std::endl;
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

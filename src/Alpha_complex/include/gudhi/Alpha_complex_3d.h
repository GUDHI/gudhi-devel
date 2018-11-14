/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2018 Inria
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ALPHA_COMPLEX_3D_H_
#define ALPHA_COMPLEX_3D_H_

#include <boost/version.hpp>
#include <boost/variant.hpp>

#if BOOST_VERSION >= 105400
#include <boost/container/static_vector.hpp>
#endif

#include <gudhi/Debug_utils.h>
#include <gudhi/Alpha_complex_options.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_traits_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>
#include <CGAL/Periodic_3_regular_triangulation_traits_3.h>
#include <CGAL/Periodic_3_regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#include <CGAL/Alpha_shape_vertex_base_3.h>

#include <CGAL/Object.h>
#include <CGAL/tuple.h>
#include <CGAL/iterator.h>
#include <CGAL/version.h>

#include <iostream>
#include <vector>
#include <unordered_map>
#include <stdexcept>
#include <cstddef>
#include <memory>       // for std::unique_ptr
#include <type_traits>  // for std::conditional and std::enable_if

#if CGAL_VERSION_NR < 1041101000
// Make compilation fail - required for external projects - https://gitlab.inria.fr/GUDHI/gudhi-devel/issues/10
static_assert(false, "Alpha_complex_3d is only available for CGAL >= 4.11");
#endif

namespace Gudhi {

namespace alpha_complex {

#ifdef GUDHI_CAN_USE_CXX11_THREAD_LOCAL
thread_local
#endif  // GUDHI_CAN_USE_CXX11_THREAD_LOCAL
    double RELATIVE_PRECISION_OF_TO_DOUBLE = 0.00001;

// Value_from_iterator returns the filtration value from an iterator on alpha shapes values
//
//                        FAST                         SAFE                         EXACT
// not weighted and       *iterator                    Specific case due to CGAL    CGAL::to_double(iterator->exact())
// not periodic                                        issue # 3153
//
// otherwise              *iterator                    CGAL::to_double(*iterator)   CGAL::to_double(iterator->exact())

template <complexity Complexity, bool Weighted_or_periodic>
struct Value_from_iterator {
  template <typename Iterator>
  static double perform(Iterator it) {
    // Default behaviour is to return the value pointed by the given iterator
    return *it;
  }
};

template <>
struct Value_from_iterator<complexity::SAFE, true> {
  template <typename Iterator>
  static double perform(Iterator it) {
    // In SAFE mode, we are with Epick with EXACT value set to CGAL::Tag_true.
    return CGAL::to_double(*it);
  }
};

template <>
struct Value_from_iterator<complexity::SAFE, false> {
  template <typename Iterator>
  static double perform(Iterator it) {
    // In SAFE mode, we are with Epeck with EXACT value set to CGAL::Tag_true.
    // Specific case due to CGAL issue https://github.com/CGAL/cgal/issues/3153
    auto approx = it->approx();
    double r;
    if (CGAL::fit_in_double(approx, r)) return r;

    // If it's precise enough, then OK.
    if (CGAL::has_smaller_relative_precision(approx, RELATIVE_PRECISION_OF_TO_DOUBLE)) return CGAL::to_double(approx);

    it->exact();
    return CGAL::to_double(it->approx());
  }
};

template <>
struct Value_from_iterator<complexity::EXACT, true> {
  template <typename Iterator>
  static double perform(Iterator it) {
    // In EXACT mode, we are with Epeck or Epick with EXACT value set to CGAL::Tag_true.
    return CGAL::to_double(it->exact());
  }
};

template <>
struct Value_from_iterator<complexity::EXACT, false> {
  template <typename Iterator>
  static double perform(Iterator it) {
    // In EXACT mode, we are with Epeck or Epick with EXACT value set to CGAL::Tag_true.
    return CGAL::to_double(it->exact());
  }
};

/**
 * \class Alpha_complex_3d
 * \brief Alpha complex data structure for 3d specific case.
 *
 * \ingroup alpha_complex
 *
 * \details
 * The data structure is constructing a <a href="https://doc.cgal.org/latest/Alpha_shapes_3/index.html">CGAL 3D Alpha
 * Shapes</a> from a range of points (can be read from an OFF file, cf. Points_off_reader).
 * Duplicate points are inserted once in the Alpha_complex. This is the reason why the vertices may be not contiguous.
 *
 * \tparam Complexity shall be `Gudhi::alpha_complex::complexity`. Default value is
 * `Gudhi::alpha_complex::complexity::FAST`.
 *
 * \tparam Weighted Boolean used to set/unset the weighted version of Alpha_complex_3d. Default value is false.
 *
 * \tparam Periodic Boolean used to set/unset the periodic version of Alpha_complex_3d. Default value is false.
 *
 * For the weighted version, weights values are explained on CGAL
 * <a href="https://doc.cgal.org/latest/Alpha_shapes_3/index.html#title0">Alpha shapes 3d</a> and
 * <a href="https://doc.cgal.org/latest/Triangulation_3/index.html#Triangulation3secclassRegulartriangulation">Regular
 * triangulation</a> documentation.
 *
 * For the periodic version, refer to the
 * <a href="https://doc.cgal.org/latest/Periodic_3_triangulation_3/index.html">CGAL’s 3D Periodic Triangulations User
 * Manual </a> for more details.
 * The periodicity is defined by an iso-oriented cuboid with diagonal opposite vertices (x_min, y_min, z_min) and
 * (x_max, y_max, z_max).
 *
 * Please refer to \ref alpha_complex for examples.
 *
 * \remark When Alpha_complex_3d is constructed with an infinite value of alpha (default value), the complex is a
 * 3d Delaunay complex.
 *
 */
template <complexity Complexity = complexity::FAST, bool Weighted = false, bool Periodic = false>
class Alpha_complex_3d {
  // Epick = Exact_predicates_inexact_constructions_kernel
  // Epeck = Exact_predicates_exact_constructions_kernel
  // ExactAlphaComparisonTag = exact version of CGAL Alpha_shape_3 and of its objects (Alpha_shape_vertex_base_3 and
  //                           Alpha_shape_cell_base_3). Not available if weighted or periodic.
  //                           Can be CGAL::Tag_false or CGAL::Tag_true
  //                           cf. https://doc.cgal.org/latest/Alpha_shapes_3/classCGAL_1_1Alpha__shape__3.html
  //
  //
  //                        FAST                         SAFE                         EXACT
  // not weighted and       Epick + CGAL::Tag_false      Epick + CGAL::Tag_true       Epick + CGAL::Tag_true
  // not periodic
  //
  // otherwise              Epick + CGAL::Tag_false      Epeck                        Epeck
  using Predicates = typename std::conditional<((!Weighted && !Periodic) || (Complexity == complexity::FAST)),
                                               CGAL::Exact_predicates_inexact_constructions_kernel,
                                               CGAL::Exact_predicates_exact_constructions_kernel>::type;

  // The other way to do a conditional type. Here there are 3 possibilities
  template <typename Predicates, bool Weighted_version, bool Periodic_version>
  struct Kernel_3 {};

  template <typename Predicates>
  struct Kernel_3<Predicates, false, false> {
    using Kernel = Predicates;
  };
  template <typename Predicates>
  struct Kernel_3<Predicates, true, false> {
    using Kernel = Predicates;
  };
  template <typename Predicates>
  struct Kernel_3<Predicates, false, true> {
    using Kernel = CGAL::Periodic_3_Delaunay_triangulation_traits_3<Predicates>;
  };
  template <typename Predicates>
  struct Kernel_3<Predicates, true, true> {
    using Kernel = CGAL::Periodic_3_regular_triangulation_traits_3<Predicates>;
  };

  using Kernel = typename Kernel_3<Predicates, Weighted, Periodic>::Kernel;

  using Exact_tag = typename std::conditional<(Complexity == complexity::FAST), CGAL::Tag_false, CGAL::Tag_true>::type;

  using TdsVb = typename std::conditional<Periodic, CGAL::Periodic_3_triangulation_ds_vertex_base_3<>,
                                          CGAL::Triangulation_ds_vertex_base_3<>>::type;

  using Tvb = typename std::conditional<Weighted, CGAL::Regular_triangulation_vertex_base_3<Kernel, TdsVb>,
                                        CGAL::Triangulation_vertex_base_3<Kernel, TdsVb>>::type;

  using Vb = CGAL::Alpha_shape_vertex_base_3<Kernel, Tvb, Exact_tag>;

  using TdsCb = typename std::conditional<Periodic, CGAL::Periodic_3_triangulation_ds_cell_base_3<>,
                                          CGAL::Triangulation_ds_cell_base_3<>>::type;

  using Tcb = typename std::conditional<Weighted, CGAL::Regular_triangulation_cell_base_3<Kernel, TdsCb>,
                                        CGAL::Triangulation_cell_base_3<Kernel, TdsCb>>::type;

  using Cb = CGAL::Alpha_shape_cell_base_3<Kernel, Tcb, Exact_tag>;
  using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb>;

  // The other way to do a conditional type. Here there 4 possibilities, cannot use std::conditional
  template <typename Kernel, typename Tds, bool Weighted_version, bool Periodic_version>
  struct Triangulation {};

  template <typename Kernel, typename Tds>
  struct Triangulation<Kernel, Tds, false, false> {
    using Triangulation_3 = CGAL::Delaunay_triangulation_3<Kernel, Tds>;
  };
  template <typename Kernel, typename Tds>
  struct Triangulation<Kernel, Tds, true, false> {
    using Triangulation_3 = CGAL::Regular_triangulation_3<Kernel, Tds>;
  };
  template <typename Kernel, typename Tds>
  struct Triangulation<Kernel, Tds, false, true> {
    using Triangulation_3 = CGAL::Periodic_3_Delaunay_triangulation_3<Kernel, Tds>;
  };
  template <typename Kernel, typename Tds>
  struct Triangulation<Kernel, Tds, true, true> {
    using Triangulation_3 = CGAL::Periodic_3_regular_triangulation_3<Kernel, Tds>;
  };

 public:
  using Triangulation_3 = typename Triangulation<Kernel, Tds, Weighted, Periodic>::Triangulation_3;

  using Alpha_shape_3 = CGAL::Alpha_shape_3<Triangulation_3, Exact_tag>;

  using Point_3 = typename Kernel::Point_3;

 private:
  using Alpha_value_type = typename Alpha_shape_3::FT;
  using Dispatch =
      CGAL::Dispatch_output_iterator<CGAL::cpp11::tuple<CGAL::Object, Alpha_value_type>,
                                     CGAL::cpp11::tuple<std::back_insert_iterator<std::vector<CGAL::Object>>,
                                                        std::back_insert_iterator<std::vector<Alpha_value_type>>>>;

  using Cell_handle = typename Alpha_shape_3::Cell_handle;
  using Facet = typename Alpha_shape_3::Facet;
  using Edge = typename Alpha_shape_3::Edge;
  using Alpha_vertex_handle = typename Alpha_shape_3::Vertex_handle;
#if BOOST_VERSION >= 105400
  using Vertex_list = boost::container::static_vector<Alpha_vertex_handle, 4>;
#else
  using Vertex_list = std::vector<Alpha_vertex_handle>;
#endif

 public:
  /** \brief Alpha_complex constructor from a list of points.
   *
   * @param[in] points Range of points to triangulate. Points must be in `Alpha_complex_3d::Point_3` or
   * `Alpha_complex_3d::Triangulation_3::Weighted_point`.
   *
   * @pre Available if Alpha_complex_3d is not Periodic.
   *
   * The type InputPointRange must be a range for which std::begin and std::end return input iterators on a
   * `Alpha_complex_3d::Point_3` or a `Alpha_complex_3d::Triangulation_3::Weighted_point`.
   */
  template <typename InputPointRange>
  Alpha_complex_3d(const InputPointRange& points) {
    static_assert(!Periodic, "This constructor is not available for periodic versions of Alpha_complex_3d");

    alpha_shape_3_ptr_ = std::unique_ptr<Alpha_shape_3>(
        new Alpha_shape_3(std::begin(points), std::end(points), 0, Alpha_shape_3::GENERAL));
  }

  /** \brief Alpha_complex constructor from a list of points and associated weights.
   *
   * @exception std::invalid_argument In debug mode, if points and weights do not have the same size.
   *
   * @param[in] points Range of points to triangulate. Points must be in `Alpha_complex_3d::Point_3`
   * @param[in] weights Range of weights on points. Weights shall be in `Alpha_complex_3d::Alpha_shape_3::FT`
   *
   * @pre Available if Alpha_complex_3d is Weighted and not Periodic.
   *
   * The type InputPointRange must be a range for which std::begin and
   * std::end return input iterators on a `Alpha_complex_3d::Point_3`.
   * The type WeightRange must be a range for which std::begin and
   * std::end return an input iterator on a `Alpha_complex_3d::Alpha_shape_3::FT`.
   */
  template <typename InputPointRange, typename WeightRange>
  Alpha_complex_3d(const InputPointRange& points, WeightRange weights) {
    static_assert(Weighted, "This constructor is not available for non-weighted versions of Alpha_complex_3d");
    static_assert(!Periodic, "This constructor is not available for periodic versions of Alpha_complex_3d");
    GUDHI_CHECK((weights.size() == points.size()),
                std::invalid_argument("Points number in range different from weights range number"));

    using Weighted_point_3 = typename Triangulation_3::Weighted_point;
    std::vector<Weighted_point_3> weighted_points_3;

    std::size_t index = 0;
    weighted_points_3.reserve(points.size());
    while ((index < weights.size()) && (index < points.size())) {
      weighted_points_3.push_back(Weighted_point_3(points[index], weights[index]));
      index++;
    }

    alpha_shape_3_ptr_ = std::unique_ptr<Alpha_shape_3>(
        new Alpha_shape_3(std::begin(weighted_points_3), std::end(weighted_points_3), 0, Alpha_shape_3::GENERAL));
  }

  /** \brief Alpha_complex constructor from a list of points and an iso-cuboid coordinates.
   *
   * @exception std::invalid_argument In debug mode, if the size of the cuboid in every directions is not the same.
   *
   * @param[in] points Range of points to triangulate. Points must be in `Alpha_complex_3d::Point_3` or
   * `Alpha_complex_3d::Triangulation_3::Weighted_point`.
   * @param[in] x_min Iso-oriented cuboid x_min.
   * @param[in] y_min Iso-oriented cuboid y_min.
   * @param[in] z_min Iso-oriented cuboid z_min.
   * @param[in] x_max Iso-oriented cuboid x_max.
   * @param[in] y_max Iso-oriented cuboid y_max.
   * @param[in] z_max Iso-oriented cuboid z_max.
   *
   * @pre Available if Alpha_complex_3d is Periodic.
   *
   * The type InputPointRange must be a range for which std::begin and std::end return input iterators on a
   * `Alpha_complex_3d::Point_3` or a `Alpha_complex_3d::Triangulation_3::Weighted_point`.
   *
   * @note In weighted version, please check weights are greater than zero, and lower than 1/64*cuboid length
   * squared.
   */
  template <typename InputPointRange>
  Alpha_complex_3d(const InputPointRange& points, Alpha_value_type x_min, Alpha_value_type y_min,
                   Alpha_value_type z_min, Alpha_value_type x_max, Alpha_value_type y_max, Alpha_value_type z_max) {
    static_assert(Periodic, "This constructor is not available for non-periodic versions of Alpha_complex_3d");
    // Checking if the cuboid is the same in x,y and z direction. If not, CGAL will not process it.
    GUDHI_CHECK(
        (x_max - x_min == y_max - y_min) && (x_max - x_min == z_max - z_min) && (z_max - z_min == y_max - y_min),
        std::invalid_argument("The size of the cuboid in every directions is not the same."));

    // Define the periodic cube
    Triangulation_3 pdt(typename Kernel::Iso_cuboid_3(x_min, y_min, z_min, x_max, y_max, z_max));
    // Heuristic for inserting large point sets (if pts is reasonably large)
    pdt.insert(std::begin(points), std::end(points), true);
    // As pdt won't be modified anymore switch to 1-sheeted cover if possible
    if (!pdt.is_triangulation_in_1_sheet()) {
      throw std::invalid_argument("Unable to construct a triangulation within a single periodic domain.");
    }
    pdt.convert_to_1_sheeted_covering();

    // alpha shape construction from points. CGAL has a strange behavior in REGULARIZED mode. This is the default mode
    // Maybe need to set it to GENERAL mode
    alpha_shape_3_ptr_ = std::unique_ptr<Alpha_shape_3>(new Alpha_shape_3(pdt, 0, Alpha_shape_3::GENERAL));
  }

  /** \brief Alpha_complex constructor from a list of points, associated weights and an iso-cuboid coordinates.
   *
   * @exception std::invalid_argument In debug mode, if points and weights do not have the same size.
   * @exception std::invalid_argument In debug mode, if the size of the cuboid in every directions is not the same.
   * @exception std::invalid_argument In debug mode, if a weight is negative, zero, or greater than 1/64*cuboid length
   * squared.
   *
   * @param[in] points Range of points to triangulate. Points must be in `Alpha_complex_3d::Point_3`
   * @param[in] weights Range of weights on points. Weights shall be in `Alpha_complex_3d::Alpha_shape_3::FT`
   * @param[in] x_min Iso-oriented cuboid x_min.
   * @param[in] y_min Iso-oriented cuboid y_min.
   * @param[in] z_min Iso-oriented cuboid z_min.
   * @param[in] x_max Iso-oriented cuboid x_max.
   * @param[in] y_max Iso-oriented cuboid y_max.
   * @param[in] z_max Iso-oriented cuboid z_max.
   *
   * @pre Available if Alpha_complex_3d is Weighted and Periodic.
   *
   * The type InputPointRange must be a range for which std::begin and
   * std::end return input iterators on a `Alpha_complex_3d::Point_3`.
   * The type WeightRange must be a range for which std::begin and
   * std::end return an input iterator on a `Alpha_complex_3d::Alpha_shape_3::FT`.
   * The type of x_min, y_min, z_min, x_max, y_max and z_max is `Alpha_complex_3d::Alpha_shape_3::FT`.
   */
  template <typename InputPointRange, typename WeightRange>
  Alpha_complex_3d(const InputPointRange& points, WeightRange weights, Alpha_value_type x_min, Alpha_value_type y_min,
                   Alpha_value_type z_min, Alpha_value_type x_max, Alpha_value_type y_max, Alpha_value_type z_max) {
    static_assert(Weighted, "This constructor is not available for non-weighted versions of Alpha_complex_3d");
    static_assert(Periodic, "This constructor is not available for non-periodic versions of Alpha_complex_3d");
    GUDHI_CHECK((weights.size() == points.size()),
                std::invalid_argument("Points number in range different from weights range number"));
    // Checking if the cuboid is the same in x,y and z direction. If not, CGAL will not process it.
    GUDHI_CHECK(
        (x_max - x_min == y_max - y_min) && (x_max - x_min == z_max - z_min) && (z_max - z_min == y_max - y_min),
        std::invalid_argument("The size of the cuboid in every directions is not the same."));

    using Weighted_point_3 = typename Triangulation_3::Weighted_point;
    std::vector<Weighted_point_3> weighted_points_3;

    std::size_t index = 0;
    weighted_points_3.reserve(points.size());

#ifdef GUDHI_DEBUG
    // Defined in GUDHI_DEBUG to avoid unused variable warning for GUDHI_CHECK
    Alpha_value_type maximal_possible_weight = 0.015625 * (x_max - x_min) * (x_max - x_min);
#endif

    while ((index < weights.size()) && (index < points.size())) {
      GUDHI_CHECK((weights[index] < maximal_possible_weight) && (weights[index] >= 0),
                  std::invalid_argument("Invalid weight at index " + std::to_string(index + 1) +
                                        ". Must be positive and less than maximal possible weight = 1/64*cuboid length "
                                        "squared, which is not an acceptable input."));
      weighted_points_3.push_back(Weighted_point_3(points[index], weights[index]));
      index++;
    }

    // Define the periodic cube
    Triangulation_3 pdt(typename Kernel::Iso_cuboid_3(x_min, y_min, z_min, x_max, y_max, z_max));
    // Heuristic for inserting large point sets (if pts is reasonably large)
    pdt.insert(std::begin(weighted_points_3), std::end(weighted_points_3), true);
    // As pdt won't be modified anymore switch to 1-sheeted cover if possible
    if (!pdt.is_triangulation_in_1_sheet()) {
      throw std::invalid_argument("Unable to construct a triangulation within a single periodic domain.");
    }
    pdt.convert_to_1_sheeted_covering();

    // alpha shape construction from points. CGAL has a strange behavior in REGULARIZED mode. This is the default mode
    // Maybe need to set it to GENERAL mode
    alpha_shape_3_ptr_ = std::unique_ptr<Alpha_shape_3>(new Alpha_shape_3(pdt, 0, Alpha_shape_3::GENERAL));
  }

  /** \brief Inserts all Delaunay triangulation into the simplicial complex.
   * It also computes the filtration values accordingly to the \ref createcomplexalgorithm
   *
   * \tparam SimplicialComplexForAlpha3d must meet `SimplicialComplexForAlpha3d` concept.
   *
   * @param[in] complex SimplicialComplexForAlpha3d to be created.
   * @param[in] max_alpha_square maximum for alpha square value. Default value is +\f$\infty\f$.
   *
   * @return true if creation succeeds, false otherwise.
   *
   * @pre The simplicial complex must be empty (no vertices)
   *
   * Initialization can be launched once.
   *
   */
  template <typename SimplicialComplexForAlpha3d,
            typename Filtration_value = typename SimplicialComplexForAlpha3d::Filtration_value>
  bool create_complex(SimplicialComplexForAlpha3d& complex,
                      Filtration_value max_alpha_square = std::numeric_limits<Filtration_value>::infinity()) {
    if (complex.num_vertices() > 0) {
      std::cerr << "Alpha_complex_3d create_complex - complex is not empty\n";
      return false;  // ----- >>
    }

    // using Filtration_value = typename SimplicialComplexForAlpha3d::Filtration_value;
    using Complex_vertex_handle = typename SimplicialComplexForAlpha3d::Vertex_handle;
    using Alpha_shape_simplex_tree_map = std::unordered_map<Alpha_vertex_handle, Complex_vertex_handle>;
    using Simplex_tree_vector_vertex = std::vector<Complex_vertex_handle>;

#ifdef DEBUG_TRACES
    std::size_t count_vertices = 0;
    std::size_t count_edges = 0;
    std::size_t count_facets = 0;
    std::size_t count_cells = 0;
#endif  // DEBUG_TRACES
    std::vector<CGAL::Object> objects;
    std::vector<Alpha_value_type> alpha_values;

    Dispatch dispatcher = CGAL::dispatch_output<CGAL::Object, Alpha_value_type>(std::back_inserter(objects),
                                                                                std::back_inserter(alpha_values));

    alpha_shape_3_ptr_->filtration_with_alpha_values(dispatcher);
#ifdef DEBUG_TRACES
    std::cout << "filtration_with_alpha_values returns : " << objects.size() << " objects" << std::endl;
#endif  // DEBUG_TRACES

    Alpha_shape_simplex_tree_map map_cgal_simplex_tree;
    using Alpha_value_iterator = typename std::vector<Alpha_value_type>::const_iterator;
    Alpha_value_iterator alpha_value_iterator = alpha_values.begin();
    for (auto object_iterator : objects) {
      Vertex_list vertex_list;

      // Retrieve Alpha shape vertex list from object
      if (const Cell_handle* cell = CGAL::object_cast<Cell_handle>(&object_iterator)) {
        for (auto i = 0; i < 4; i++) {
#ifdef DEBUG_TRACES
          std::cout << "from cell[" << i << "]=" << (*cell)->vertex(i)->point() << std::endl;
#endif  // DEBUG_TRACES
          vertex_list.push_back((*cell)->vertex(i));
        }
#ifdef DEBUG_TRACES
        count_cells++;
#endif  // DEBUG_TRACES
      } else if (const Facet* facet = CGAL::object_cast<Facet>(&object_iterator)) {
        for (auto i = 0; i < 4; i++) {
          if ((*facet).second != i) {
#ifdef DEBUG_TRACES
            std::cout << "from facet=[" << i << "]" << (*facet).first->vertex(i)->point() << std::endl;
#endif  // DEBUG_TRACES
            vertex_list.push_back((*facet).first->vertex(i));
          }
        }
#ifdef DEBUG_TRACES
        count_facets++;
#endif  // DEBUG_TRACES
      } else if (const Edge* edge = CGAL::object_cast<Edge>(&object_iterator)) {
        for (auto i : {(*edge).second, (*edge).third}) {
#ifdef DEBUG_TRACES
          std::cout << "from edge[" << i << "]=" << (*edge).first->vertex(i)->point() << std::endl;
#endif  // DEBUG_TRACES
          vertex_list.push_back((*edge).first->vertex(i));
        }
#ifdef DEBUG_TRACES
        count_edges++;
#endif  // DEBUG_TRACES
      } else if (const Alpha_vertex_handle* vertex = CGAL::object_cast<Alpha_vertex_handle>(&object_iterator)) {
#ifdef DEBUG_TRACES
        count_vertices++;
        std::cout << "from vertex=" << (*vertex)->point() << std::endl;
#endif  // DEBUG_TRACES
        vertex_list.push_back((*vertex));
      }
      // Construction of the vector of simplex_tree vertex from list of alpha_shapes vertex
      Simplex_tree_vector_vertex the_simplex;
      for (auto the_alpha_shape_vertex : vertex_list) {
        auto the_map_iterator = map_cgal_simplex_tree.find(the_alpha_shape_vertex);
        if (the_map_iterator == map_cgal_simplex_tree.end()) {
          // alpha shape not found
          Complex_vertex_handle vertex = map_cgal_simplex_tree.size();
#ifdef DEBUG_TRACES
          std::cout << "vertex [" << the_alpha_shape_vertex->point() << "] not found - insert " << vertex << std::endl;
#endif  // DEBUG_TRACES
          the_simplex.push_back(vertex);
          map_cgal_simplex_tree.emplace(the_alpha_shape_vertex, vertex);
        } else {
          // alpha shape found
          Complex_vertex_handle vertex = the_map_iterator->second;
#ifdef DEBUG_TRACES
          std::cout << "vertex [" << the_alpha_shape_vertex->point() << "] found in " << vertex << std::endl;
#endif  // DEBUG_TRACES
          the_simplex.push_back(vertex);
        }
      }
      // Construction of the simplex_tree
      Filtration_value filtr = Value_from_iterator<Complexity, (Weighted || Periodic)>::perform(alpha_value_iterator);

#ifdef DEBUG_TRACES
      std::cout << "filtration = " << filtr << std::endl;
#endif  // DEBUG_TRACES
      complex.insert_simplex(the_simplex, static_cast<Filtration_value>(filtr));
      GUDHI_CHECK(alpha_value_iterator != alpha_values.end(), "CGAL provided more simplices than values");
      ++alpha_value_iterator;
    }

#ifdef DEBUG_TRACES
    std::cout << "vertices \t" << count_vertices << std::endl;
    std::cout << "edges \t\t" << count_edges << std::endl;
    std::cout << "facets \t\t" << count_facets << std::endl;
    std::cout << "cells \t\t" << count_cells << std::endl;
#endif  // DEBUG_TRACES
    // --------------------------------------------------------------------------------------------
    // As Alpha value is an approximation, we have to make filtration non decreasing while increasing the dimension
    complex.make_filtration_non_decreasing();
    // Remove all simplices that have a filtration value greater than max_alpha_square
    complex.prune_above_filtration(max_alpha_square);
    // --------------------------------------------------------------------------------------------
    return true;
  }

 private:
  // use of a unique_ptr on cgal Alpha_shape_3, as copy and default constructor is not available - no need to be freed
  std::unique_ptr<Alpha_shape_3> alpha_shape_3_ptr_;
};

}  // namespace alpha_complex

}  // namespace Gudhi

#endif  // ALPHA_COMPLEX_3D_H_

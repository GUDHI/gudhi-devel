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
#include <gudhi/Alpha_complex_3d_options.h>

#include <CGAL/Object.h>
#include <CGAL/tuple.h>
#include <CGAL/iterator.h>


namespace Gudhi {

namespace alpha_complex {

template<typename AlphaComplex3dOptions>
class Alpha_complex_3d {
private:
  using Alpha_shape_3 = typename AlphaComplex3dOptions::Alpha_shape_3;
  // filtration with alpha values needed type definition
  using Alpha_value_type = typename Alpha_shape_3::FT;
  using Object = CGAL::Object;
  using Dispatch = CGAL::Dispatch_output_iterator<CGAL::cpp11::tuple<Object, Alpha_value_type>,
      CGAL::cpp11::tuple<std::back_insert_iterator<std::vector<Object> >,
          std::back_insert_iterator<std::vector<Alpha_value_type> > > >;
  using Cell_handle = typename Alpha_shape_3::Cell_handle;
  using Facet = typename Alpha_shape_3::Facet;
  using Edge_3 = typename Alpha_shape_3::Edge;
  using Alpha_vertex_handle = typename Alpha_shape_3::Vertex_handle;

#if BOOST_VERSION >= 105400
  using Vertex_list = boost::container::static_vector<Alpha_vertex_handle, 4>;
#else
  using Vertex_list = std::vector<Alpha_vertex_handle>;
#endif

public:
  /** \brief Alpha_complex constructor from a list of points.
  *
  * Duplicate points are inserted once in the Alpha_complex. This is the reason why the vertices may be not contiguous.
  *
  * @param[in] points Range of points to triangulate. Points must be in AlphaComplex3dOptions::Point_3
  *
  * The type InputPointRange must be a range for which std::begin and
  * std::end return input iterators on a AlphaComplex3dOptions::Point_3.
  */
  template<typename InputPointRange >
  Alpha_complex_3d(const InputPointRange& points) {
    static_assert(!AlphaComplex3dOptions::weighted, "This constructor is not available for weighted versions of Alpha_complex_3d");
    static_assert(!AlphaComplex3dOptions::periodic, "This constructor is not available for periodic versions of Alpha_complex_3d");
    std::cout << points[0] << std::endl;
    Alpha_shape_3 alpha_shape_3(std::begin(points), std::end(points), 0, Alpha_shape_3::GENERAL);

    Dispatch disp = CGAL::dispatch_output<Object, Alpha_value_type>(std::back_inserter(the_objects),
                                                                    std::back_inserter(the_alpha_values));

    alpha_shape_3.filtration_with_alpha_values(disp);
#ifdef DEBUG_TRACES
    std::cout << "filtration_with_alpha_values returns : " << the_objects.size() << " objects" << std::endl;
#endif  // DEBUG_TRACES

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
  template<typename InputPointRange , typename WeightRange>
  Alpha_complex_3d(const InputPointRange& points, WeightRange weights) {
    static_assert(AlphaComplex3dOptions::weighted,
                  "This constructor is not available for non-weighted versions of Alpha_complex_3d");
    static_assert(!AlphaComplex3dOptions::periodic,
                  "This constructor is not available for periodic versions of Alpha_complex_3d");
    GUDHI_CHECK((weights.size() == points.size()),
                std::invalid_argument("Alpha_complex_3d constructor with weights requires points number to be equal with points number"));

    using Weighted_point_3 = typename AlphaComplex3dOptions::Weighted_point_3;
    std::vector<Weighted_point_3> weighted_points_3;

    std::size_t index = 0;
    weighted_points_3.reserve(points.size());
    while ((index < weights.size()) && (index < points.size())) {
      weighted_points_3.push_back(Weighted_point_3(points[index], weights[index]));
      index++;
    }

    Alpha_shape_3 alpha_shape_3(std::begin(weighted_points_3), std::end(weighted_points_3), 0, Alpha_shape_3::GENERAL);

    Dispatch disp = CGAL::dispatch_output<Object, Alpha_value_type>(std::back_inserter(the_objects),
                                                                    std::back_inserter(the_alpha_values));

    alpha_shape_3.filtration_with_alpha_values(disp);
#ifdef DEBUG_TRACES
    std::cout << "filtration_with_alpha_values returns : " << the_objects.size() << " objects" << std::endl;
#endif  // DEBUG_TRACES
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
  template<typename InputPointRange>
  Alpha_complex_3d(const InputPointRange& points,
                   Alpha_value_type x_min, Alpha_value_type y_min, Alpha_value_type z_min,
                   Alpha_value_type x_max, Alpha_value_type y_max, Alpha_value_type z_max) {
    static_assert(!AlphaComplex3dOptions::weighted,
                  "This constructor is not available for weighted versions of Alpha_complex_3d");
    static_assert(AlphaComplex3dOptions::periodic,
                  "This constructor is not available for non-periodic versions of Alpha_complex_3d");
    // Checking if the cuboid is the same in x,y and z direction. If not, CGAL will not process it.
    GUDHI_CHECK((x_max - x_min != y_max - y_min) || (x_max - x_min != z_max - z_min) || (z_max - z_min != y_max - y_min),
                std::invalid_argument("The size of the cuboid in every directions is not the same."));

    using Periodic_delaunay_triangulation_3 = typename AlphaComplex3dOptions::Periodic_delaunay_triangulation_3;
    using Iso_cuboid_3 = typename AlphaComplex3dOptions::Iso_cuboid_3;
    // Define the periodic cube
    Periodic_delaunay_triangulation_3 pdt(Iso_cuboid_3(x_min, y_min, z_min, x_max, y_max, z_max));
    // Heuristic for inserting large point sets (if pts is reasonably large)
    pdt.insert(std::begin(points), std::end(points), true);
    // As pdt won't be modified anymore switch to 1-sheeted cover if possible
    GUDHI_CHECK(pdt.is_triangulation_in_1_sheet(),
                std::invalid_argument("Uable to construct a triangulation within a single periodic domain."));

    // alpha shape construction from points. CGAL has a strange behavior in REGULARIZED mode. This is the default mode
    // Maybe need to set it to GENERAL mode
    Alpha_shape_3 as(pdt, 0, Alpha_shape_3::GENERAL);

    // filtration with alpha values from alpha shape
    std::vector<Object> the_objects;
    std::vector<Alpha_value_type> the_alpha_values;

    Dispatch disp = CGAL::dispatch_output<Object, Alpha_value_type>(std::back_inserter(the_objects),
                                                                    std::back_inserter(the_alpha_values));

    as.filtration_with_alpha_values(disp);
#ifdef DEBUG_TRACES
    std::cout << "filtration_with_alpha_values returns : " << the_objects.size() << " objects" << std::endl;
#endif  // DEBUG_TRACES

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
  template<typename InputPointRange , typename WeightRange>
  Alpha_complex_3d(const InputPointRange& points, WeightRange weights,
                   Alpha_value_type x_min, Alpha_value_type y_min, Alpha_value_type z_min,
                   Alpha_value_type x_max, Alpha_value_type y_max, Alpha_value_type z_max) {
    static_assert(AlphaComplex3dOptions::weighted,
                  "This constructor is not available for non-weighted versions of Alpha_complex_3d");
    static_assert(AlphaComplex3dOptions::periodic,
                  "This constructor is not available for non-periodic versions of Alpha_complex_3d");
    GUDHI_CHECK((weights.size() == points.size()),
                std::invalid_argument("Alpha_complex_3d constructor with weights requires points number to be equal with points number"));
    // Checking if the cuboid is the same in x,y and z direction. If not, CGAL will not process it.
    GUDHI_CHECK((x_max - x_min != y_max - y_min) || (x_max - x_min != z_max - z_min) || (z_max - z_min != y_max - y_min),
                std::invalid_argument("The size of the cuboid in every directions is not the same."));

    double maximal_possible_weight = 0.015625 * (x_max - x_min) * (x_max - x_min);

    using Weighted_point_3 = typename AlphaComplex3dOptions::Weighted_point_3;
    std::vector<Weighted_point_3> weighted_points_3;

    std::size_t index = 0;
    weighted_points_3.reserve(points.size());
    while ((index < weights.size()) && (index < points.size())) {
      GUDHI_CHECK((weights[index] >= maximal_possible_weight) || (weights[index] < 0),
                  std::invalid_argument("Invalid weight at line" + std::to_string(index + 1) +
                                        ". Must be positive and less than maximal possible weight = 1/64*cuboid length "
                                        "squared, which is not an acceptable input."));
      weighted_points_3.push_back(Weighted_point_3(points[index], weights[index]));
      index++;
    }

    using Periodic_delaunay_triangulation_3 = typename AlphaComplex3dOptions::Periodic_delaunay_triangulation_3;
    using Iso_cuboid_3 = typename AlphaComplex3dOptions::Iso_cuboid_3;
    // Define the periodic cube
    Periodic_delaunay_triangulation_3 pdt(Iso_cuboid_3(x_min, y_min, z_min, x_max, y_max, z_max));
    // Heuristic for inserting large point sets (if pts is reasonably large)
    pdt.insert(std::begin(weighted_points_3), std::end(weighted_points_3), true);
    // As pdt won't be modified anymore switch to 1-sheeted cover if possible
    GUDHI_CHECK(pdt.is_triangulation_in_1_sheet(),
                std::invalid_argument("Uable to construct a triangulation within a single periodic domain."));

    // alpha shape construction from points. CGAL has a strange behavior in REGULARIZED mode. This is the default mode
    // Maybe need to set it to GENERAL mode
    Alpha_shape_3 as(pdt, 0, Alpha_shape_3::GENERAL);

    // filtration with alpha values from alpha shape
    std::vector<Object> the_objects;
    std::vector<Alpha_value_type> the_alpha_values;

    Dispatch disp = CGAL::dispatch_output<Object, Alpha_value_type>(std::back_inserter(the_objects),
                                                                    std::back_inserter(the_alpha_values));

    as.filtration_with_alpha_values(disp);
#ifdef DEBUG_TRACES
    std::cout << "filtration_with_alpha_values returns : " << the_objects.size() << " objects" << std::endl;
#endif  // DEBUG_TRACES
  }

  template <typename SimplicialComplexForAlpha3d>
  void create_complex(SimplicialComplexForAlpha3d& complex) {
    using Filtration_value = typename SimplicialComplexForAlpha3d::Filtration_value;
    create_complex(complex, std::numeric_limits<Filtration_value>::infinity());
  }

  /** \brief Inserts all Delaunay triangulation into the simplicial complex.
   * It also computes the filtration values accordingly to the \ref createcomplexalgorithm
   *
   * \tparam SimplicialComplexForAlpha must meet `SimplicialComplexForAlpha` concept.
   *
   * @param[in] complex SimplicialComplexForAlpha to be created.
   * @param[in] max_alpha_square maximum for alpha square value. Default value is +\f$\infty\f$.
   *
   * @return true if creation succeeds, false otherwise.
   *
   * @pre Delaunay triangulation must be already constructed with dimension strictly greater than 0.
   * @pre The simplicial complex must be empty (no vertices)
   *
   * Initialization can be launched once.
   */
  template <typename SimplicialComplexForAlpha3d>
  void create_complex(SimplicialComplexForAlpha3d& complex,
                      typename SimplicialComplexForAlpha3d::Filtration_value max_alpha_square) {
    using Filtration_value = typename SimplicialComplexForAlpha3d::Filtration_value;
    using Vertex_handle = typename SimplicialComplexForAlpha3d::Vertex_handle;
#ifdef DEBUG_TRACES
    Alpha_shape_3::size_type count_vertices = 0;
    Alpha_shape_3::size_type count_edges = 0;
    Alpha_shape_3::size_type count_facets = 0;
    Alpha_shape_3::size_type count_cells = 0;
#endif  // DEBUG_TRACES

    // Loop on objects vector
    Vertex_list vertex_list;
    SimplicialComplexForAlpha3d simplex_tree;
    std::map<Alpha_vertex_handle, Vertex_handle> map_cgal_simplex_tree;
    typename std::vector<Alpha_value_type>::iterator the_alpha_value_iterator = the_alpha_values.begin();
    for (auto object_iterator : the_objects) {
      // Retrieve Alpha shape vertex list from object
      if (const Cell_handle *cell = CGAL::object_cast<Cell_handle>(&object_iterator)) {
        vertex_list = from_cell<Vertex_list, Cell_handle>(*cell);
#ifdef DEBUG_TRACES
        count_cells++;
#endif  // DEBUG_TRACES
      } else if (const Facet *facet = CGAL::object_cast<Facet>(&object_iterator)) {
        vertex_list = from_facet<Vertex_list, Facet>(*facet);
#ifdef DEBUG_TRACES
        count_facets++;
#endif  // DEBUG_TRACES
      } else if (const Edge_3 *edge = CGAL::object_cast<Edge_3>(&object_iterator)) {
        vertex_list = from_edge<Vertex_list, Edge_3>(*edge);
#ifdef DEBUG_TRACES
        count_edges++;
#endif  // DEBUG_TRACES
      } else if (const Vertex_handle *vertex = CGAL::object_cast<Vertex_handle>(&object_iterator)) {
        vertex_list = from_vertex<Vertex_list, Vertex_handle>(*vertex);
#ifdef DEBUG_TRACES
        count_vertices++;
#endif  // DEBUG_TRACES
      }
      // Construction of the vector of simplex_tree vertex from list of alpha_shapes vertex
      std::vector<Vertex_handle> the_simplex;
      for (auto the_alpha_shape_vertex : vertex_list) {
        typename std::map<Alpha_vertex_handle, Vertex_handle>::iterator the_map_iterator =
            map_cgal_simplex_tree.find(the_alpha_shape_vertex);
        if (the_map_iterator == map_cgal_simplex_tree.end()) {
          // alpha shape not found
          Vertex_handle vertex = map_cgal_simplex_tree.size();
#ifdef DEBUG_TRACES
          std::cout << "vertex [" << the_alpha_shape_vertex->point() << "] not found - insert " << vertex << std::endl;
#endif  // DEBUG_TRACES
          the_simplex.push_back(vertex);
          map_cgal_simplex_tree.emplace(the_alpha_shape_vertex, vertex);
        } else {
          // alpha shape found
          Vertex_handle vertex = the_map_iterator->second;
#ifdef DEBUG_TRACES
          std::cout << "vertex [" << the_alpha_shape_vertex->point() << "] found in " << vertex << std::endl;
#endif  // DEBUG_TRACES
          the_simplex.push_back(vertex);
        }
      }
      // Construction of the simplex_tree
      Filtration_value filtr = /*std::sqrt*/ (*the_alpha_value_iterator);
#ifdef DEBUG_TRACES
      std::cout << "filtration = " << filtr << std::endl;
#endif  // DEBUG_TRACES
      simplex_tree.insert_simplex(the_simplex, filtr);
      GUDHI_CHECK(the_alpha_value_iterator != the_alpha_values.end(), "CGAL provided more simplices than values");
      ++the_alpha_value_iterator;
    }

#ifdef DEBUG_TRACES
    std::cout << "vertices \t\t" << count_vertices << std::endl;
    std::cout << "edges \t\t" << count_edges << std::endl;
    std::cout << "facets \t\t" << count_facets << std::endl;
    std::cout << "cells \t\t" << count_cells << std::endl;

    std::cout << "Information of the Simplex Tree: " << std::endl;
    std::cout << "  Number of vertices = " << simplex_tree.num_vertices() << " ";
    std::cout << "  Number of simplices = " << simplex_tree.num_simplices() << std::endl << std::endl;
    std::cout << "  Dimension = " << simplex_tree.dimension() << " ";
#endif  // DEBUG_TRACES

#ifdef DEBUG_TRACES
    std::cout << "Iterator on vertices: " << std::endl;
    for (auto vertex : simplex_tree.complex_vertex_range()) {
      std::cout << vertex << " ";
    }
#endif  // DEBUG_TRACES

  }

private:
  template <class Vertex_list, class Cell_handle>
  Vertex_list from_cell(const Cell_handle& ch) {
    Vertex_list the_list;
    for (auto i = 0; i < 4; i++) {
#ifdef DEBUG_TRACES
      std::cout << "from cell[" << i << "]=" << ch->vertex(i)->point() << std::endl;
#endif  // DEBUG_TRACES
      the_list.push_back(ch->vertex(i));
    }
    return the_list;
  }

  template <class Vertex_list, class Facet>
  Vertex_list from_facet(const Facet& fct) {
    Vertex_list the_list;
    for (auto i = 0; i < 4; i++) {
      if (fct.second != i) {
#ifdef DEBUG_TRACES
        std::cout << "from facet=[" << i << "]" << fct.first->vertex(i)->point() << std::endl;
#endif  // DEBUG_TRACES
        the_list.push_back(fct.first->vertex(i));
      }
    }
    return the_list;
  }

  template <class Vertex_list, class Edge_3>
  Vertex_list from_edge(const Edge_3& edg) {
    Vertex_list the_list;
    for (auto i : {edg.second, edg.third}) {
#ifdef DEBUG_TRACES
      std::cout << "from edge[" << i << "]=" << edg.first->vertex(i)->point() << std::endl;
#endif  // DEBUG_TRACES
      the_list.push_back(edg.first->vertex(i));
    }
    return the_list;
  }

  template <class Vertex_list, class Vertex_handle>
  Vertex_list from_vertex(const Vertex_handle& vh) {
    Vertex_list the_list;
#ifdef DEBUG_TRACES
    std::cout << "from vertex=" << vh->point() << std::endl;
#endif  // DEBUG_TRACES
    the_list.push_back(vh);
    return the_list;
  }

private:
  std::vector<Object> the_objects;
  std::vector<Alpha_value_type> the_alpha_values;

};

}  // namespace alpha_complex

}  // namespace Gudhi

#endif  // ALPHA_COMPLEX_3D_H_

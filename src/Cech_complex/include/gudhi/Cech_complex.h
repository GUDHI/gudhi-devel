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

#ifndef CECH_COMPLEX_H_
#define CECH_COMPLEX_H_

#include <gudhi/distance_functions.h>        // for Gudhi::Minimal_enclosing_ball_radius
#include <gudhi/graph_simplicial_complex.h>  // for Gudhi::Proximity_graph
#include <gudhi/Debug_utils.h>               // for GUDHI_CHECK
#include <gudhi/Cech_complex_blocker.h>      // for Gudhi::cech_complex::Cech_blocker

#include <iostream>
#include <stdexcept>  // for exception management

namespace Gudhi {

namespace cech_complex {

/**
 * \class Cech_complex
 * \brief Cech complex data structure.
 * 
 * \ingroup cech_complex
 * 
 * \details
 * The data structure is a proximity graph, containing edges when the edge length is less or equal
 * to a given max_radius. Edge length is computed from `Gudhi::Minimal_enclosing_ball_radius` distance function.
 *
 * \tparam SimplicialComplexForProximityGraph furnishes `Vertex_handle` and `Filtration_value` type definition required
 * by `Gudhi::Proximity_graph`.
 *
 * \tparam ForwardPointRange must be a range for which `std::begin()` and `std::end()` methods return input
 * iterators on a point. `std::begin()` and `std::end()` methods are also required for a point.
 */
template<typename SimplicialComplexForProximityGraph, typename ForwardPointRange>
class Cech_complex {
 private:
  using Vertex_handle = typename SimplicialComplexForProximityGraph::Vertex_handle;
  using Filtration_value = typename SimplicialComplexForProximityGraph::Filtration_value;
  using Proximity_graph = Gudhi::Proximity_graph<SimplicialComplexForProximityGraph>;
  using Point_iterator = typename boost::range_const_iterator<ForwardPointRange>::type;
  using Point= typename std::iterator_traits<Point_iterator>::value_type;
  using Point_cloud = std::vector<Point>;

 public:
  /** \brief Cech_complex constructor from a list of points.
   *
   * @param[in] points Range of points.
   * @param[in] max_radius Maximal radius value.
   *
   * \tparam ForwardPointRange must be a range for which `std::begin()` and `std::end()` methods return input
   * iterators on a point. `std::begin()` and `std::end()` methods are also required for a point.
   *
   */
  Cech_complex(const ForwardPointRange& points, Filtration_value max_radius)
    : max_radius_(max_radius),
      point_cloud_(std::begin(points), std::end(points)) {
    cech_skeleton_graph_ =
        Gudhi::compute_proximity_graph<SimplicialComplexForProximityGraph>(point_cloud_,
                                                                           max_radius_,
                                                                           Gudhi::Minimal_enclosing_ball_radius());
  }

  /** \brief Initializes the simplicial complex from the proximity graph and expands it until a given maximal
   * dimension, using the Cech blocker oracle.
   *
   * @param[in] complex SimplicialComplexForCech to be created.
   * @param[in] dim_max graph expansion until this given maximal dimension.
   * @exception std::invalid_argument In debug mode, if `complex.num_vertices()` does not return 0.
   *
   */
  template<typename SimplicialComplexForCechComplex >
  void create_complex(SimplicialComplexForCechComplex& complex, int dim_max) {
    GUDHI_CHECK(complex.num_vertices() == 0,
                std::invalid_argument("Cech_complex::create_complex - simplicial complex is not empty"));

    // insert the proximity graph in the simplicial complex
    complex.insert_graph(cech_skeleton_graph_);
    // expand the graph until dimension dim_max
    complex.expansion_with_blockers(dim_max,
                                    Cech_blocker<SimplicialComplexForCechComplex, ForwardPointRange>(complex, this));
  }

  /** @return max_radius value given at construction. */
  Filtration_value max_radius() const {
    return max_radius_;
  }

  /** @param[in] vertex Point position in the range.
   * @return A const iterator on the point.
   * @exception std::out_of_range In debug mode, if point position in the range is out.
   */
  typename ForwardPointRange::const_iterator point_iterator(std::size_t vertex) const {
    GUDHI_CHECK((std::begin(point_cloud_) + vertex) < std::end(point_cloud_),
                std::out_of_range("Cech_complex::point - simplicial complex is not empty"));
    return (std::begin(point_cloud_) + vertex);
  }

 private:
  Proximity_graph cech_skeleton_graph_;
  Filtration_value max_radius_;
  Point_cloud point_cloud_;
};

}  // namespace cech_complex

}  // namespace Gudhi

#endif  // CECH_COMPLEX_H_

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
 * \tparam InputPointRange must be a range for which `std::begin()` and `std::end()` methods return input
 * iterators on a point. `std::begin()` and `std::end()` methods are also required for a point.
 */
template<typename SimplicialComplexForProximityGraph, typename InputPointRange>
class Cech_complex {
 private:
  // Required by compute_proximity_graph
  using Vertex_handle = typename SimplicialComplexForProximityGraph::Vertex_handle;
  using Filtration_value = typename SimplicialComplexForProximityGraph::Filtration_value;
  using Proximity_graph = Gudhi::Proximity_graph<SimplicialComplexForProximityGraph>;

  // Retrieve Coordinate type from InputPointRange
  using Point_from_range_iterator = typename boost::range_const_iterator<InputPointRange>::type;
  using Point_from_range = typename std::iterator_traits<Point_from_range_iterator>::value_type;
  using Coordinate_iterator = typename boost::range_const_iterator<Point_from_range>::type;
  using Coordinate= typename std::iterator_traits<Coordinate_iterator>::value_type;

 public:
  // Point and Point_cloud type definition
  using Point = std::vector<Coordinate>;
  using Point_cloud = std::vector<Point>;

 public:
  /** \brief Cech_complex constructor from a list of points.
   *
   * @param[in] points Range of points.
   * @param[in] max_radius Maximal radius value.
   *
   * \tparam InputPointRange must be a range of Point. Point must be a range of <b>copyable</b> Cartesian coordinates.
   *
   */
  Cech_complex(const InputPointRange& points, Filtration_value max_radius)
    : max_radius_(max_radius) {
    // Point cloud deep copy
    auto points_begin_itr = std::begin(points);
    auto points_end_itr = std::end(points);

    point_cloud_.reserve(points_end_itr - points_begin_itr);
    for (auto point_itr = points_begin_itr; point_itr < points_end_itr; point_itr++) {
      point_cloud_.push_back(Point(std::begin(*point_itr), std::end(*point_itr)));
    }

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
                                    Cech_blocker<SimplicialComplexForCechComplex, Cech_complex>(&complex, this));
  }

  /** @return max_radius value given at construction. */
  Filtration_value max_radius() const {
    return max_radius_;
  }

  /** @param[in] vertex Point position in the range.
   * @return The point.
   */
  const Point& get_point(Vertex_handle vertex) const{
    return point_cloud_[vertex];
  }

 private:
  Proximity_graph cech_skeleton_graph_;
  Filtration_value max_radius_;
  Point_cloud point_cloud_;
};

}  // namespace cech_complex

}  // namespace Gudhi

#endif  // CECH_COMPLEX_H_

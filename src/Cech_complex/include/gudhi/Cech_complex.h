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

#include <gudhi/graph_simplicial_complex.h>  // for Gudhi::Proximity_graph
#include <gudhi/Debug_utils.h>               // for GUDHI_CHECK
#include <gudhi/Cech_complex_blocker.h>      // for Gudhi::cech_complex::Cech_blocker

#include <iostream>
#include <vector>
#include <cstddef>
#include <iterator>  // for std::size

namespace Gudhi {

namespace cech_complex {

/**
 * \class Cech_complex
 * \brief Cech complex data structure.
 * 
 * \ingroup Cech_complex
 * 
 * \details
 * The data structure is a one skeleton graph, or Rips graph, containing edges when the edge length is less or equal
 * to a given threshold. Edge length is computed from a user given point cloud with a given distance function, or a
 * distance matrix.
 * 
 * \tparam Filtration_value is the type used to store the filtration values of the simplicial complex.
 */
template<typename SimplicialComplexForCechComplex, typename ForwardPointRange>
class Cech_complex {
 private:
  using Vertex_handle = typename SimplicialComplexForCechComplex::Vertex_handle;
  using Filtration_value = typename SimplicialComplexForCechComplex::Filtration_value;
  using Proximity_graph = Gudhi::Proximity_graph<SimplicialComplexForCechComplex>;

 public:
  /** \brief Cech_complex constructor from a list of points.
   *
   * @param[in] points Range of points.
   * @param[in] threshold Rips value.
   * @param[in] distance distance function that returns a `Filtration_value` from 2 given points.
   * 
   * \tparam ForwardPointRange must be a range for which `std::begin` and `std::end` return input iterators on a
   * point.
   *
   * \tparam Distance furnishes `operator()(const Point& p1, const Point& p2)`, where
   * `Point` is a point from the `ForwardPointRange`, and that returns a `Filtration_value`.
   */
  template<typename Distance >
  Cech_complex(const ForwardPointRange& points, Filtration_value threshold, Distance distance)
    : threshold_(threshold),
      point_cloud_(points) {
    GUDHI_CHECK(std::size(points) > 0,
                std::invalid_argument("Cech_complex::create_complex - point cloud is empty"));
    dimension_ = points[0].size();
    cech_skeleton_graph_ = Gudhi::compute_proximity_graph<SimplicialComplexForCechComplex>(point_cloud_, threshold_, distance);
  }

  /** \brief Initializes the simplicial complex from the Rips graph and expands it until a given maximal
   * dimension.
   *
   * @param[in] complex SimplicialComplexForCech to be created.
   * @param[in] dim_max graph expansion for Rips until this given maximal dimension.
   * @exception std::invalid_argument In debug mode, if `complex.num_vertices()` does not return 0.
   *
   */
  void create_complex(SimplicialComplexForCechComplex& complex, int dim_max) {
    GUDHI_CHECK(complex.num_vertices() == 0,
                std::invalid_argument("Cech_complex::create_complex - simplicial complex is not empty"));

    // insert the proximity graph in the simplicial complex
    complex.insert_graph(cech_skeleton_graph_);
    // expand the graph until dimension dim_max
    complex.expansion_with_blockers(dim_max,
                                    Cech_blocker<SimplicialComplexForCechComplex, ForwardPointRange>(complex, this));
  }

  Filtration_value threshold() const {
    return threshold_;
  }

  std::size_t dimension() const {
    return dimension_;
  }

  auto point(std::size_t vertex) const -> decltype(point_cloud_.begin() + vertex) {
    GUDHI_CHECK((point_cloud_.begin() + vertex) < point_cloud_.end(),
                std::invalid_argument("Cech_complex::point - simplicial complex is not empty"));
    return (point_cloud_.begin() + vertex);
  }

 private:
  Proximity_graph cech_skeleton_graph_;
  Filtration_value threshold_;
  ForwardPointRange point_cloud_;
  std::size_t dimension_;
};

}  // namespace cech_complex

}  // namespace Gudhi

#endif  // CECH_COMPLEX_H_

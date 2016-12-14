/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clément Maria, Pawel Dlotko, Vincent Rouvreau
 *
 *    Copyright (C) 2016  INRIA
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

#ifndef RIPS_COMPLEX_H_
#define RIPS_COMPLEX_H_

#include <gudhi/Debug_utils.h>
#include <gudhi/graph_simplicial_complex.h>

#include <boost/graph/adjacency_list.hpp>

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <limits>  // for numeric_limits
#include <utility>  // for pair<>


namespace Gudhi {

namespace rips_complex {

/**
 * \class Rips_complex
 * \brief Rips complex data structure.
 * 
 * \ingroup rips_complex
 * 
 * \details
 * The data structure is a one skeleton graph, or Rips graph, constructed from a point cloud, containing edges when
 * the edge length is less or equal to a given threshold. Edge length is computed from a user given function.
 * 
 * The complex is a template class requiring a Filtration_value type.
 * 
 * \tparam Filtration_value must meet `SimplicialComplexForRips` concept.
 */
template<typename Filtration_value>
class Rips_complex {
 private:
  typedef typename boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS
  , boost::property < vertex_filtration_t, Filtration_value >
  , boost::property < edge_filtration_t, Filtration_value >> Graph_t;

  typedef int Vertex_handle;

 public:
  /** \brief Rips_complex constructor from a list of points.
   *
   * @param[in] points Range of points.
   * @param[in] threshold rips value.
   * @param[in] distance distance function that returns a Filtration_value from 2 given points.
   * 
   * The type InputPointRange must be a range for which std::begin and std::end return input iterators on a point.
   */
  template<typename InputPointRange, typename Distance >
  Rips_complex(const InputPointRange& points, Filtration_value threshold, Distance distance) {
    compute_proximity_graph<InputPointRange, Distance >(points, threshold, distance);
  }

  /** \brief Rips_complex constructor from a distance matrix.
   *
   * @param[in] distance_matrix Range of distances.
   * @param[in] threshold rips value.
   * 
   * The type InputDistanceRange must have a \code size() \endcode method and on which distance_matrix[i][j] returns
   * the distance between points \f$i\f$ and \f$j\f$ as long as \f$ 0 \subseteq i \subseteq j \subseteq
   * distance_matrix.size().\f$
   */
  template<typename InputDistanceRange>
  Rips_complex(const InputDistanceRange& distance_matrix, Filtration_value threshold) {
    compute_proximity_graph(boost::irange((size_t)0, distance_matrix.size()), threshold,
                            [&](size_t i, size_t j){return distance_matrix[j][i];});
  }

  /** \brief Initializes the simplicial complex from the Rips graph and expands it until a given maximal
   * dimension.
   *
   * \tparam SimplicialComplexForRips must meet `SimplicialComplexForRips` concept.
   * 
   * @param[in] complex SimplicialComplexForRips to be created.
   * @param[in] dim_max graph expansion for rips until this given maximal dimension.
   * @exception std::invalid_argument In debug mode, if \code complex.num_vertices() \endcode does not return 0.
   *
   */
  template <typename SimplicialComplexForRips>
  void create_complex(SimplicialComplexForRips& complex, int dim_max) {
    GUDHI_CHECK(complex.num_vertices() == 0,
                std::invalid_argument("Rips_complex::create_complex - simplicial complex is not empty"));

    // insert the proximity graph in the simplicial complex
    complex.insert_graph(rips_skeleton_graph_);
    // expand the graph until dimension dim_max
    complex.expansion(dim_max);
  }

 private:
  /** \brief Computes the proximity graph of the points.
   *
   * If points contains n elements, the proximity graph is the graph with n vertices, and an edge [u,v] iff the
   * distance function between points u and v is smaller than threshold.
   *
   * \tparam The type InputPointRange furnishes .begin() and .end() methods, that return iterators with
   * value_type Point.
   */
  template< typename InputPointRange, typename Distance >
  void compute_proximity_graph(const InputPointRange& points, Filtration_value threshold,
               Distance distance) {
    std::vector< std::pair< Vertex_handle, Vertex_handle > > edges;
    std::vector< Filtration_value > edges_fil;

    // Compute the proximity graph of the points.
    // If points contains n elements, the proximity graph is the graph with n vertices, and an edge [u,v] iff the
    // distance function between points u and v is smaller than threshold.
    // --------------------------------------------------------------------------------------------
    // Creates the vector of edges and its filtration values (returned by distance function)
    Vertex_handle idx_u = 0;
    for (auto it_u = std::begin(points); it_u != std::end(points); ++it_u) {
      Vertex_handle idx_v = idx_u + 1;
      for (auto it_v = it_u + 1; it_v != std::end(points); ++it_v, ++idx_v) {
        Filtration_value fil = distance(*it_u, *it_v);
        if (fil <= threshold) {
          edges.emplace_back(idx_u, idx_v);
          edges_fil.push_back(fil);
        }
      }
      ++idx_u;
    }

    // --------------------------------------------------------------------------------------------
    // Creates the proximity graph from edges and sets the property with the filtration value.
    // Number of points is labeled from 0 to idx_u-1
    // --------------------------------------------------------------------------------------------
    // Do not use : rips_skeleton_graph_ = Graph_t(...) -> deep copy of the graph (boost graph is not move-enabled)
    rips_skeleton_graph_.~Graph_t();
    new(&rips_skeleton_graph_)Graph_t(edges.begin(), edges.end(), edges_fil.begin(), idx_u);

    auto vertex_prop = boost::get(vertex_filtration_t(), rips_skeleton_graph_);

    using vertex_iterator = typename boost::graph_traits<Graph_t>::vertex_iterator;
    vertex_iterator vi, vi_end;
    for (std::tie(vi, vi_end) = boost::vertices(rips_skeleton_graph_);
         vi != vi_end; ++vi) {
      boost::put(vertex_prop, *vi, 0.);
    }
  }

 private:
  Graph_t rips_skeleton_graph_;
};

}  // namespace rips_complex

}  // namespace Gudhi

#endif  // RIPS_COMPLEX_H_

/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria, Pawel Dlotko, Vincent Rouvreau
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
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
 * The data structure is a one skeleton graph, or Rips graph, containing edges when the edge length is less or equal
 * to a given threshold. Edge length is computed from a user given point cloud with a given distance function, or a
 * distance matrix.
 * 
 * \tparam Filtration_value is the type used to store the filtration values of the simplicial complex.
 */
template<typename Filtration_value>
class Rips_complex {
 public:
  /**
   * \brief Type of the one skeleton graph stored inside the Rips complex structure.
   */
  typedef typename boost::adjacency_list < boost::vecS, boost::vecS, boost::directedS
  , boost::property < vertex_filtration_t, Filtration_value >
  , boost::property < edge_filtration_t, Filtration_value >> OneSkeletonGraph;

 private:
  typedef int Vertex_handle;

 public:
  /** \brief Rips_complex constructor from a list of points.
   *
   * @param[in] points Range of points.
   * @param[in] threshold Maximal edge length. All edges of the graph strictly greater than `threshold` are not
   * inserted in the graph.
   * 
   * @param[in] distance distance function that returns a `Filtration_value` from 2 given points.
   * 
   * \tparam ForwardPointRange must be a range for which `std::begin` and `std::end` return input iterators on a
   * point.
   *
   * \tparam Distance furnishes `operator()(const Point& p1, const Point& p2)`, where
   * `Point` is a point from the `ForwardPointRange`, and that returns a `Filtration_value`.
   */
  template<typename ForwardPointRange, typename Distance >
  Rips_complex(const ForwardPointRange& points, Filtration_value threshold, Distance distance) {
    compute_proximity_graph(points, threshold, distance);
  }

  /** \brief Rips_complex constructor from a distance matrix.
   *
   * @param[in] distance_matrix Range of distances.
   * @param[in] threshold Maximal edge length. All edges of the graph strictly greater than `threshold` are not
   * inserted in the graph.
   * 
   * \tparam DistanceMatrix must have a `size()` method and on which `distance_matrix[i][j]` returns
   * the distance between points \f$i\f$ and \f$j\f$ as long as \f$ 0 \leqslant j < i \leqslant
   * distance\_matrix.size().\f$
   */
  template<typename DistanceMatrix>
  Rips_complex(const DistanceMatrix& distance_matrix, Filtration_value threshold) {
    compute_proximity_graph(boost::irange((size_t)0, distance_matrix.size()), threshold,
                            [&](size_t i, size_t j){return distance_matrix[j][i];});
  }

  /** \brief Initializes the simplicial complex from the Rips graph and expands it until a given maximal
   * dimension.
   *
   * \tparam SimplicialComplexForRips must meet `SimplicialComplexForRips` concept.
   * 
   * @param[in] complex SimplicialComplexForRips to be created.
   * @param[in] dim_max graph expansion for Rips until this given maximal dimension.
   * @exception std::invalid_argument In debug mode, if `complex.num_vertices()` does not return 0.
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
   * \tparam ForwardPointRange furnishes `.begin()` and `.end()`
   * methods.
   *
   * \tparam Distance furnishes `operator()(const Point& p1, const Point& p2)`, where
   * `Point` is a point from the `ForwardPointRange`, and that returns a `Filtration_value`.
   */
  template< typename ForwardPointRange, typename Distance >
  void compute_proximity_graph(const ForwardPointRange& points, Filtration_value threshold,
               Distance distance) {
    std::vector< std::pair< Vertex_handle, Vertex_handle > > edges;
    std::vector< Filtration_value > edges_fil;

    // Compute the proximity graph of the points.
    // If points contains n elements, the proximity graph is the graph with n vertices, and an edge [u,v] iff the
    // distance function between points u and v is smaller than threshold.
    // --------------------------------------------------------------------------------------------
    // Creates the vector of edges and its filtration values (returned by distance function)
    Vertex_handle idx_u = 0;
    for (auto it_u = std::begin(points); it_u != std::end(points); ++it_u, ++idx_u) {
      Vertex_handle idx_v = idx_u + 1;
      for (auto it_v = it_u + 1; it_v != std::end(points); ++it_v, ++idx_v) {
        Filtration_value fil = distance(*it_u, *it_v);
        if (fil <= threshold) {
          edges.emplace_back(idx_u, idx_v);
          edges_fil.push_back(fil);
        }
      }
    }

    // --------------------------------------------------------------------------------------------
    // Creates the proximity graph from edges and sets the property with the filtration value.
    // Number of points is labeled from 0 to idx_u-1
    // --------------------------------------------------------------------------------------------
    // Do not use : rips_skeleton_graph_ = OneSkeletonGraph(...) -> deep copy of the graph (boost graph is not
    // move-enabled)
    rips_skeleton_graph_.~OneSkeletonGraph();
    new(&rips_skeleton_graph_)OneSkeletonGraph(edges.begin(), edges.end(), edges_fil.begin(), idx_u);

    auto vertex_prop = boost::get(vertex_filtration_t(), rips_skeleton_graph_);

    using vertex_iterator = typename boost::graph_traits<OneSkeletonGraph>::vertex_iterator;
    vertex_iterator vi, vi_end;
    for (std::tie(vi, vi_end) = boost::vertices(rips_skeleton_graph_);
         vi != vi_end; ++vi) {
      boost::put(vertex_prop, *vi, 0.);
    }
  }

 private:
  OneSkeletonGraph rips_skeleton_graph_;
};

}  // namespace rips_complex

}  // namespace Gudhi

#endif  // RIPS_COMPLEX_H_

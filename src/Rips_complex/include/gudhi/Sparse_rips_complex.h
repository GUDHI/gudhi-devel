/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s): Marc Glisse
 *
 *    Copyright (C) 2018 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef SPARSE_RIPS_COMPLEX_H_
#define SPARSE_RIPS_COMPLEX_H_

#include <gudhi/Debug_utils.h>
#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/choose_n_farthest_points.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/range/metafunctions.hpp>

#include <vector>

namespace Gudhi {

namespace rips_complex {

// The whole interface is copied on Rips_complex. A redesign should be discussed with all complex creation classes in
// mind.

/**
 * \class Sparse_rips_complex
 * \brief Sparse Rips complex data structure.
 *
 * \ingroup rips_complex
 *
 * \details
 * This class is used to construct a sparse \f$(1+O(\epsilon))\f$-approximation of `Rips_complex`, i.e. a filtered
 * simplicial complex that is multiplicatively
 * \f$(1+O(\epsilon))\f$-interleaved with the Rips filtration. More precisely,
 * this is a \f$(1,\frac{1}{1-\epsilon}\f$-interleaving.
 *
 * \tparam Filtration_value is the type used to store the filtration values of the simplicial complex.
 */
template <typename Filtration_value>
class Sparse_rips_complex {
 private:
  // TODO(MG): use a different graph where we know we can safely insert in parallel.
  typedef typename boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
                                         boost::property<vertex_filtration_t, Filtration_value>,
                                         boost::property<edge_filtration_t, Filtration_value>>
      Graph;

  typedef int Vertex_handle;

 public:
  /** \brief Sparse_rips_complex constructor from a list of points.
   *
   * @param[in] points Range of points.
   * @param[in] distance Distance function that returns a `Filtration_value` from 2 given points.
   * @param[in] epsilon Approximation parameter. epsilon must be positive.
   * @param[in] mini Minimal filtration value. Ignore anything below this scale. This is a less efficient version of `Gudhi::subsampling::sparsify_point_set()`.
   * @param[in] maxi Maximal filtration value. Ignore anything above this scale.
   *
   */
  template <typename RandomAccessPointRange, typename Distance>
  Sparse_rips_complex(const RandomAccessPointRange& points, Distance distance, double epsilon, Filtration_value mini=-std::numeric_limits<Filtration_value>::infinity(), Filtration_value maxi=std::numeric_limits<Filtration_value>::infinity())
      : epsilon_(epsilon) {
    GUDHI_CHECK(epsilon > 0, "epsilon must be positive");
    auto dist_fun = [&](Vertex_handle i, Vertex_handle j) { return distance(points[i], points[j]); };
    Ker<decltype(dist_fun)> kernel(dist_fun);
    subsampling::choose_n_farthest_points(kernel, boost::irange<Vertex_handle>(0, boost::size(points)), -1, -1,
                                          std::back_inserter(sorted_points), std::back_inserter(params));
    compute_sparse_graph(dist_fun, epsilon, mini, maxi);
  }

  /** \brief Sparse_rips_complex constructor from a distance matrix.
   *
   * @param[in] distance_matrix Range of range of distances.
   * `distance_matrix[i][j]` returns the distance between points \f$i\f$ and
   * \f$j\f$ as long as \f$ 0 \leqslant j < i \leqslant
   * distance\_matrix.size().\f$
   * @param[in] epsilon Approximation parameter. epsilon must be positive.
   * @param[in] mini Minimal filtration value. Ignore anything below this scale. This is a less efficient version of `Gudhi::subsampling::sparsify_point_set()`.
   * @param[in] maxi Maximal filtration value. Ignore anything above this scale.
   */
  template <typename DistanceMatrix>
  Sparse_rips_complex(const DistanceMatrix& distance_matrix, double epsilon, Filtration_value mini=-std::numeric_limits<Filtration_value>::infinity(), Filtration_value maxi=std::numeric_limits<Filtration_value>::infinity())
      : Sparse_rips_complex(boost::irange<Vertex_handle>(0, boost::size(distance_matrix)),
                            [&](Vertex_handle i, Vertex_handle j) { return (i==j) ? 0 : (i<j) ? distance_matrix[j][i] : distance_matrix[i][j]; },
                            epsilon, mini, maxi) {}

  /** \brief Fills the simplicial complex with the sparse Rips graph and
   * expands it with all the cliques, stopping at a given maximal dimension.
   *
   * \tparam SimplicialComplexForRips must meet `SimplicialComplexForRips` concept.
   *
   * @param[in] complex the complex to fill
   * @param[in] dim_max maximal dimension of the simplicial complex.
   * @exception std::invalid_argument In debug mode, if `complex.num_vertices()` does not return 0.
   *
   */
  template <typename SimplicialComplexForRips>
  void create_complex(SimplicialComplexForRips& complex, int dim_max) {
    GUDHI_CHECK(complex.num_vertices() == 0,
                std::invalid_argument("Sparse_rips_complex::create_complex - simplicial complex is not empty"));

    complex.insert_graph(graph_);
    if(epsilon_ >= 1) {
      complex.expansion(dim_max);
      return;
    }
    const int n = boost::size(params);
    std::vector<Filtration_value> lambda(n);
    // lambda[original_order]=params[sorted_order]
    for(int i=0;i<n;++i)
      lambda[sorted_points[i]] = params[i];
    double cst = epsilon_ * (1 - epsilon_) / 2;
    auto block = [cst,&complex,&lambda](typename SimplicialComplexForRips::Simplex_handle sh){
      auto filt = complex.filtration(sh);
      auto mini = filt * cst;
      for(auto v : complex.simplex_vertex_range(sh)){
        if(lambda[v] < mini)
          return true; // v died before this simplex could be born
      }
      return false;
    };
    complex.expansion_with_blockers(dim_max, block);
  }

 private:
  // choose_n_farthest_points wants the distance function in this form...
  template <class Distance>
  struct Ker {
    typedef std::size_t Point_d;  // index into point range
    Ker(Distance& d) : dist(d) {}
    // Despite the name, this is not squared...
    typedef Distance Squared_distance_d;
    Squared_distance_d& squared_distance_d_object() const { return dist; }
    Distance& dist;
  };

  // PointRange must be random access.
  template <typename Distance>
  void compute_sparse_graph(Distance& dist, double epsilon, Filtration_value mini, Filtration_value maxi) {
    const auto& points = sorted_points; // convenience alias
    const int n = boost::size(points);
    double cst = epsilon * (1 - epsilon) / 2;
    graph_.~Graph();
    new (&graph_) Graph(n);
    // for(auto v : vertices(g)) // doesn't work :-(
    typename boost::graph_traits<Graph>::vertex_iterator v_i, v_e;
    for (std::tie(v_i, v_e) = vertices(graph_); v_i != v_e; ++v_i) {
      auto v = *v_i;
      // This whole loop might not be necessary, leave it until someone investigates if it is safe to remove.
      put(vertex_filtration_t(), graph_, v, 0);
    }

    // TODO(MG):
    // - make it parallel
    // - only test near-enough neighbors
    for (int i = 0; i < n; ++i) {
      auto&& pi = points[i];
      auto li = params[i];
      if (li < mini) break;
      for (int j = i + 1; j < n; ++j) {
        auto&& pj = points[j];
        auto d = dist(pi, pj);
        auto lj = params[j];
        if (lj < mini) break;
        GUDHI_CHECK(lj <= li, "Bad furthest point sorting");
        Filtration_value alpha;

        // The paper has d/2 and d-lj/e to match the Cech, but we use doubles to match the Rips
        if (d * epsilon <= 2 * lj)
          alpha = d;
        else if (d * epsilon > li + lj)
          continue;
        else {
          alpha = (d - lj / epsilon) * 2;
          // Keep the test exactly the same as in block to avoid inconsistencies
          if (epsilon < 1 && alpha * cst > lj)
            continue;
        }

        if (alpha <= maxi)
          add_edge(pi, pj, alpha, graph_);
      }
    }
  }

  Graph graph_;
  double epsilon_;
  // Because of the arbitrary split between constructor and create_complex
  // sorted_points[sorted_order]=original_order
  std::vector<Vertex_handle> sorted_points;
  // params[sorted_order]=distance to previous points
  std::vector<Filtration_value> params;
};

}  // namespace rips_complex

}  // namespace Gudhi

#endif  // SPARSE_RIPS_COMPLEX_H_

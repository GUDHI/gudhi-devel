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

#include <boost/graph/graph_traits.hpp>
#include <boost/range/metafunctions.hpp>
#include <boost/iterator/counting_iterator.hpp>

#include <vector>

namespace Gudhi {
namespace rips_complex {
// A custom graph class, because boost::adjacency_list does not conveniently allow to choose vertex descriptors
template <class Vertex_handle, class Filtration_value>
struct Graph {
  typedef std::vector<Vertex_handle> VList;
  typedef std::vector<std::tuple<Vertex_handle, Vertex_handle, Filtration_value>> EList;
  typedef typename VList::const_iterator vertex_iterator;
  typedef boost::counting_iterator<std::size_t> edge_iterator;
  VList vlist;
  EList elist;
};
template <class Vertex_handle, class Filtration_value>
void add_vertex(Vertex_handle v, Graph<Vertex_handle, Filtration_value>&g) { g.vlist.push_back(v); }
template <class Vertex_handle, class Filtration_value>
void add_edge(Vertex_handle u, Vertex_handle v, Filtration_value f, Graph<Vertex_handle, Filtration_value>&g) { g.elist.emplace_back(u, v, f); }
template <class Vertex_handle, class Filtration_value>
std::size_t num_vertices(Graph<Vertex_handle, Filtration_value> const&g) { return g.vlist.size(); }
template <class Vertex_handle, class Filtration_value>
std::size_t num_edges(Graph<Vertex_handle, Filtration_value> const&g) { return g.elist.size(); }
template <class Vertex_handle, class Filtration_value, class Iter = typename Graph<Vertex_handle, Filtration_value>::vertex_iterator>
std::pair<Iter, Iter>
vertices(Graph<Vertex_handle, Filtration_value> const&g) {
  return { g.vlist.begin(), g.vlist.end() };
}
template <class Vertex_handle, class Filtration_value>
std::pair<boost::counting_iterator<std::size_t>, boost::counting_iterator<std::size_t>>
edges(Graph<Vertex_handle, Filtration_value> const&g) {
  typedef boost::counting_iterator<std::size_t> I;
  return { I(0), I(g.elist.size()) };
}
template <class Vertex_handle, class Filtration_value>
Vertex_handle source(std::size_t e, Graph<Vertex_handle, Filtration_value> const&g) { return std::get<0>(g.elist[e]); }
template <class Vertex_handle, class Filtration_value>
Vertex_handle target(std::size_t e, Graph<Vertex_handle, Filtration_value> const&g) { return std::get<1>(g.elist[e]); }
template <class Vertex_handle, class Filtration_value>
Filtration_value get(vertex_filtration_t, Graph<Vertex_handle, Filtration_value> const&, Vertex_handle) { return 0; }
template <class Vertex_handle, class Filtration_value>
Filtration_value get(edge_filtration_t, Graph<Vertex_handle, Filtration_value> const&g, std::size_t e) { return std::get<2>(g.elist[e]); }
}  // namespace rips_complex
}  // namespace Gudhi
namespace boost {
template <class Vertex_handle, class Filtration_value>
struct graph_traits<Gudhi::rips_complex::Graph<Vertex_handle, Filtration_value>> {
  typedef Gudhi::rips_complex::Graph<Vertex_handle, Filtration_value> G;
  struct traversal_category : vertex_list_graph_tag, edge_list_graph_tag {};
  typedef Vertex_handle vertex_descriptor;
  typedef typename G::vertex_iterator vertex_iterator;
  typedef std::size_t vertices_size_type;
  typedef std::size_t edge_descriptor;
  typedef typename G::edge_iterator edge_iterator;
  typedef std::size_t edges_size_type;
  typedef directed_tag directed_category;
  typedef disallow_parallel_edge_tag edge_parallel_category;
};
// Etc, since we don't expose this graph to the world, we know we are not going to query property_traits for instance.
}

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
  typedef int Vertex_handle;
  typedef rips_complex::Graph<Vertex_handle, Filtration_value> Graph;

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
  Sparse_rips_complex(const RandomAccessPointRange& points, Distance distance, double const epsilon, Filtration_value const mini=-std::numeric_limits<Filtration_value>::infinity(), Filtration_value const maxi=std::numeric_limits<Filtration_value>::infinity())
      : epsilon_(epsilon) {
    GUDHI_CHECK(epsilon > 0, "epsilon must be positive");
    auto dist_fun = [&](Vertex_handle i, Vertex_handle j) { return distance(points[i], points[j]); };
    // TODO: stop choose_n_farthest_points once it reaches mini or 0?
    subsampling::choose_n_farthest_points_metric(dist_fun, boost::irange<Vertex_handle>(0, boost::size(points)), -1, -1,
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
  Sparse_rips_complex(const DistanceMatrix& distance_matrix, double const epsilon, Filtration_value const mini=-std::numeric_limits<Filtration_value>::infinity(), Filtration_value const maxi=std::numeric_limits<Filtration_value>::infinity())
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
  void create_complex(SimplicialComplexForRips& complex, int const dim_max) {
    GUDHI_CHECK(complex.num_vertices() == 0,
                std::invalid_argument("Sparse_rips_complex::create_complex - simplicial complex is not empty"));

    complex.insert_graph(graph_);
    if(epsilon_ >= 1) {
      complex.expansion(dim_max);
      return;
    }
    const std::size_t n = num_vertices(graph_);
    std::vector<Filtration_value> lambda(max_v + 1);
    // lambda[original_order]=params[sorted_order]
    for(std::size_t i=0;i<n;++i)
      lambda[sorted_points[i]] = params[i];
    double cst = epsilon_ * (1 - epsilon_) / 2;
    auto block = [cst,&complex,&lambda](typename SimplicialComplexForRips::Simplex_handle sh){
      auto filt = complex.filtration(sh);
      auto min_f = filt * cst;
      for(auto v : complex.simplex_vertex_range(sh)){
        if(lambda[v] < min_f)
          return true; // v died before this simplex could be born
      }
      return false;
    };
    complex.expansion_with_blockers(dim_max, block);
  }

 private:
  // PointRange must be random access.
  template <typename Distance>
  void compute_sparse_graph(Distance& dist, double const epsilon, Filtration_value const mini, Filtration_value const maxi) {
    const auto& points = sorted_points; // convenience alias
    std::size_t n = boost::size(points);
    double cst = epsilon * (1 - epsilon) / 2;
    max_v = -1; // Useful for the size of the map lambda.
    for (std::size_t i = 0; i < n; ++i) {
      if ((params[i] < mini || params[i] <= 0) && i != 0) break;
      // The parameter of the first point is not very meaningful, it is supposed to be infinite,
      // but if the type does not support it...
      // It would be better to do this reduction of the number of points earlier, around choose_n_farthest_points.
      add_vertex(points[i], graph_);
      max_v = std::max(max_v, points[i]);
    }
    n = num_vertices(graph_);

    // TODO(MG):
    // - make it parallel
    // - only test near-enough neighbors
    for (std::size_t i = 0; i < n; ++i) {
      auto&& pi = points[i];
      auto li = params[i];
      // If we inserted all the points, points with multiplicity would get connected to their first representative,
      // no need to handle the redundant ones in the outer loop.
      // if (li <= 0 && i != 0) break;
      for (std::size_t j = i + 1; j < n; ++j) {
        auto&& pj = points[j];
        auto d = dist(pi, pj);
        auto lj = params[j];
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
  Vertex_handle max_v;
  // Because of the arbitrary split between constructor and create_complex
  // sorted_points[sorted_order]=original_order
  std::vector<Vertex_handle> sorted_points;
  // params[sorted_order]=distance to previous points
  std::vector<Filtration_value> params;
};

}  // namespace rips_complex

}  // namespace Gudhi

#endif  // SPARSE_RIPS_COMPLEX_H_

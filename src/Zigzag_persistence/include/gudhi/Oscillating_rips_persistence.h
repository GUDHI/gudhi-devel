/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file Oscillating_rips_persistence.h
 * @author Hannah Schreiber
 * @brief Contains the implementation of the @ref Gudhi::zigzag_persistence::compute_oscillating_rips_persistence
 * method.
 */

#ifndef OSCILLATING_RIPS_PERSISTENCE_H_
#define OSCILLATING_RIPS_PERSISTENCE_H_

#include <boost/container/map.hpp>

#include <gudhi/Zigzag_persistence/oscillating_rips_iterators.h>
#include <gudhi/filtered_zigzag_persistence.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/distance_functions.h>

namespace Gudhi {
namespace zigzag_persistence {

/**
 * @ingroup zigzag_persistence
 *
 * @brief Determines the type of the oscillating Rips edge range.
 *
 * The two options are usually quite equivalent, as the base algorithm behind them is the same.
 * They mostly differ by their memory management, so it can be worth it to try both out.
 */
enum Edge_range_type {
  VECTOR,     /**< The edges are computed all at once and stored in a vector. */
  BOOST_RANGE /**< The edges are computed one by one at each increment of the
                   range iterator and therefore not stored. */
};

/**
 * @ingroup zigzag_persistence
 *
 * @brief Model of @ref SimplexTreeOptions, as expected from @ref Oscillating_rips_simplex_range.
 * The values of `stable_simplex_handles`, `store_key` and `store_filtration` are mandatory.
 * `Simplex_key` has to be a signed integer type if @ref Zigzag_persistence is used
 * for the range. Otherwise, the options can be readapted.
 *
 * Maximum number of simplices to compute persistence is <CODE>std::numeric_limits<std::int64_t>::max()</CODE>
 * (way enough simplices). If less are needed, just inherit from this structure and redefine
 * the type `Simplex_key`.
 */
struct Simplex_tree_options_oscillating_rips {
  typedef linear_indexing_tag Indexing_tag;
  typedef int Vertex_handle;
  typedef double Filtration_value;
  typedef std::int64_t Simplex_key;
  static const bool store_key = true;
  static const bool store_filtration = true;
  static const bool contiguous_vertices = true;
  static const bool link_nodes_by_label = true;
  static const bool stable_simplex_handles = true;
};

template <class StableFilteredComplex>
struct Default_oscillating_rips_zigzag_options : Default_filtered_zigzag_options {
  using Face_key = typename StableFilteredComplex::Simplex_handle;
  using Filtration_value = typename StableFilteredComplex::Filtration_value;
  using Dimension = int;  // it is `int` in the simplex tree
  struct Hash {
    std::size_t operator()(const Face_key& sh) const {
      return sh->second.key();
    }
  };
  struct KeyEqual {
    bool operator()(const Face_key& sh1, const Face_key& sh2) const {
      return sh1->second.key() == sh2->second.key();
    }
  };
};

/**
 * @brief Computes the oscillating Rips filtration based on the given parameters and computes the
 * corresponding zigzag barcode.
 *
 * @ingroup zigzag_persistence
 *
 * @details It is thought to be a helper function to easily compute oscillating Rips persistence with the
 * possibility to switch out a few common options. But additional options exists for the iterators
 * (e.g. replacing the Euclidean distance by another one, or using your own epsilon values; see
 * the documentation of @ref Oscillating_rips_edge_range and @ref Oscillating_rips_simplex_range for
 * more information). One can easily create their own method based on this one.
 *
 * @tparam PointRange Range containing the point cloud.
 * @tparam edge_range_type Either Edge_range_type::BOOST_RANGE or Edge_range_type::VECTOR.
 * Default value: Edge_range_type::BOOST_RANGE.
 * @tparam OscillatingRipsSimplexRange Type of the oscillating Rips simplex range.
 * Default value: Oscillating_rips_simplex_range with templates depending on @p edge_range_type.
 * @tparam ZigzagPersistenceOptions Options for the matrix used by the zigzag persistence algorithm.
 * Default value: @ref Gudhi::persistence_matrix::Zigzag_options<>.
 * @param points Point cloud.
 * @param nu Lower multiplier.
 * @param mu Upper multiplier.
 * @param maxDim Maximum dimension to which to expand the Rips complex. If set to -1, there is no limit.
 * @param p Order policy for the points.
 * Can be either @ref Oscillating_rips_edge_range::Order_policy::FARTHEST_POINT_ORDERING,
 * @ref Oscillating_rips_edge_range::Order_policy::ALREADY_ORDERED or
 * @ref Oscillating_rips_edge_range::Order_policy::RANDOM_POINT_ORDERING.
 * Default value: @ref Oscillating_rips_edge_range::Order_policy::FARTHEST_POINT_ORDERING.
 * @return The persistence diagram of the oscillating Rips filtration.
 */
template <typename PointRange,
          typename F,
          Edge_range_type edge_range_type = Edge_range_type::BOOST_RANGE,
          class StableFilteredComplex = Gudhi::Simplex_tree<Simplex_tree_options_oscillating_rips> >
void compute_oscillating_rips_persistence(
    const PointRange& points,
    double nu,
    double mu,
    int maxDim,
    F&& outStream,
    Oscillating_rips_edge_order_policy p = Oscillating_rips_edge_order_policy::FARTHEST_POINT_ORDERING)
{
  using ZPOptions = Default_oscillating_rips_zigzag_options<StableFilteredComplex>;
  using ZP = Filtered_zigzag_persistence<ZPOptions>;
  using Filtration_value = typename ZPOptions::Filtration_value;
  using Dimension = typename ZPOptions::Dimension;  // always `int` for now
  using Bar = Gudhi::persistence_matrix::Persistence_interval<Dimension, Filtration_value>;
  using EdgeRange = Oscillating_rips_edge_range<Filtration_value>;
  using VectorEdgeRange = std::vector<Zigzag_edge<Filtration_value> >;
  using EdgeRangeIterator = typename std::conditional<edge_range_type == Edge_range_type::BOOST_RANGE,
                                                      typename EdgeRange::Oscillating_rips_edge_iterator,
                                                      typename VectorEdgeRange::const_iterator>::type;
  using OscillatingRipsSimplexRange = Oscillating_rips_simplex_range<StableFilteredComplex, EdgeRangeIterator>;

  StableFilteredComplex st;
  ZP zp(outStream);

  EdgeRangeIterator start, end;
  VectorEdgeRange vec;

  if constexpr (edge_range_type == Edge_range_type::BOOST_RANGE) {
    start = EdgeRange::begin(nu, mu, points, Gudhi::Euclidean_distance(), p);
    end = EdgeRange::end();
  } else {
    vec = EdgeRange::compute_vector_range(nu, mu, points, Gudhi::Euclidean_distance(), p);
    start = vec.begin();
    end = vec.end();
  }

  for (const auto& t : OscillatingRipsSimplexRange::get_iterator_range(start, end, st, maxDim)) {
    if (std::get<2>(t))
      zp.insert_face(std::get<0>(t),
                     st.boundary_simplex_range(std::get<0>(t)),
                     st.dimension(std::get<0>(t)),
                     std::get<1>(t));
    else
      zp.remove_face(std::get<0>(t), std::get<1>(t));
  }

  zp.get_current_infinite_intervals([&](Dimension dim, Filtration_value birth) { outStream(dim, birth, Bar::inf); });
}

/**
 * @brief Computes the oscillating Rips filtration based on the given parameters and computes the
 * corresponding zigzag barcode.
 *
 * @ingroup zigzag_persistence
 *
 * @details It is thought to be a helper function to easily compute oscillating Rips persistence with the
 * possibility to switch out a few common options. But additional options exists for the iterators
 * (e.g. replacing the Euclidean distance by another one, or using your own epsilon values; see
 * the documentation of @ref Oscillating_rips_edge_range and @ref Oscillating_rips_simplex_range for
 * more information). One can easily create their own method based on this one.
 *
 * @tparam PointRange Range containing the point cloud.
 * @tparam edge_range_type Either Edge_range_type::BOOST_RANGE or Edge_range_type::VECTOR.
 * Default value: Edge_range_type::BOOST_RANGE.
 * @tparam OscillatingRipsSimplexRange Type of the oscillating Rips simplex range.
 * Default value: Oscillating_rips_simplex_range with templates depending on @p edge_range_type.
 * @tparam ZigzagPersistenceOptions Options for the matrix used by the zigzag persistence algorithm.
 * Default value: @ref Gudhi::persistence_matrix::Zigzag_options<>.
 * @param points Point cloud.
 * @param nu Lower multiplier.
 * @param mu Upper multiplier.
 * @param maxDim Maximum dimension to which to expand the Rips complex. If set to -1, there is no limit.
 * @param p Order policy for the points.
 * Can be either @ref Oscillating_rips_edge_range::Order_policy::FARTHEST_POINT_ORDERING,
 * @ref Oscillating_rips_edge_range::Order_policy::ALREADY_ORDERED or
 * @ref Oscillating_rips_edge_range::Order_policy::RANDOM_POINT_ORDERING.
 * Default value: @ref Oscillating_rips_edge_range::Order_policy::FARTHEST_POINT_ORDERING.
 * @return The persistence diagram of the oscillating Rips filtration.
 */
template <typename PointRange,
          Edge_range_type edge_range_type = Edge_range_type::BOOST_RANGE,
          class StableFilteredComplex = Gudhi::Simplex_tree<Simplex_tree_options_oscillating_rips> >
std::vector<Gudhi::persistence_matrix::Persistence_interval<int, typename StableFilteredComplex::Filtration_value> >
compute_oscillating_rips_persistence(
    const PointRange& points,
    double nu,
    double mu,
    int maxDim,
    Oscillating_rips_edge_order_policy p = Oscillating_rips_edge_order_policy::FARTHEST_POINT_ORDERING)
{
  using ZPOptions = Default_oscillating_rips_zigzag_options<StableFilteredComplex>;
  using ZP = Filtered_zigzag_persistence_with_storage<ZPOptions>;
  using Filtration_value = typename ZPOptions::Filtration_value;
  using EdgeRange = Oscillating_rips_edge_range<Filtration_value>;
  using VectorEdgeRange = std::vector<Zigzag_edge<Filtration_value> >;
  using EdgeRangeIterator = typename std::conditional<edge_range_type == Edge_range_type::BOOST_RANGE,
                                                      typename EdgeRange::Oscillating_rips_edge_iterator,
                                                      typename VectorEdgeRange::const_iterator>::type;
  using OscillatingRipsSimplexRange = Oscillating_rips_simplex_range<StableFilteredComplex, EdgeRangeIterator>;

  StableFilteredComplex st;
  ZP zp;

  EdgeRangeIterator start, end;
  VectorEdgeRange vec;

  if constexpr (edge_range_type == Edge_range_type::BOOST_RANGE) {
    start = EdgeRange::begin(nu, mu, points, Gudhi::Euclidean_distance(), p);
    end = EdgeRange::end();
  } else {
    vec = EdgeRange::compute_vector_range(nu, mu, points, Gudhi::Euclidean_distance(), p);
    start = vec.begin();
    end = vec.end();
  }

  for (const auto& t : OscillatingRipsSimplexRange::get_iterator_range(start, end, st, maxDim)) {
    if (std::get<2>(t))
      zp.insert_face(std::get<0>(t),
                     st.boundary_simplex_range(std::get<0>(t)),
                     st.dimension(std::get<0>(t)),
                     std::get<1>(t));
    else
      zp.remove_face(std::get<0>(t), std::get<1>(t));
  }

  return zp.get_persistence_diagram();
}

}  // namespace zigzag_persistence

}  // namespace Gudhi

#endif  // OSCILLATING_RIPS_PERSISTENCE_H_

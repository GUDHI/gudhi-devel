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

#include <gudhi/Zigzag_persistence/oscillating_rips_iterators.h>
#include <gudhi/Zigzag_persistence.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/options.h>
#include <gudhi/distance_functions.h>

namespace Gudhi {
namespace zigzag_persistence {

/**
 * @brief Determines the type of the oscillanting Rips edge range.
 *
 * The two options are usually quite equivalent, as the base algorithm behind them is the same.
 * They mostly differ by their memory managment, so it can be worth it to try both out.
 */
enum Edge_range_type {
  VECTOR,     /**< The edges are computed all at once and stored in a vector. */
  BOOST_RANGE /**< The edges are computed one by one at each incrementation of
               * the range iterator and therefore not stored. */
};

/**
 * @brief Computes the oscillating Rips filtration based on the given parameters and computes the 
 * corresponding zigzag barcode.
 *
 * @ingroup zigzag_persistence
 *
 * @details It is thought to be a helper function to easily compute oscillating Rips persistence with the 
 * possibility to switch out a few common options. But additional options exists for the iterators 
 * (e.g. replacing the Euclidian distance by another one, or using your own epsilon values; see
 * the documention of @ref Oscillating_rips_edge_range and @ref Oscillating_rips_simplex_range for
 * more information). One can easily create their own method based on this one.
 * 
 * @tparam PointRange Range containing the point cloud.
 * @tparam edge_range_type Either Edge_range_type::BOOST_RANGE or Edge_range_type::VECTOR. 
 * Default value: Edge_range_type::BOOST_RANGE.
 * @tparam OscillatingRipsSimplexRange Type of the oscilliating Rips simplex range. 
 * Default value: Oscillating_rips_simplex_range with templates depending on @p edge_range_type.
 * @tparam ZigzagPersistenceOptions Options for the matrix used by the zigzag persistence algorithm.
 * Default value: @ref Gudhi::persistence_matrix::Zigzag_options<>.
 * @param points Point cloud.
 * @param nu Lower multiplicator.
 * @param mu Upper multiplicator.
 * @param maxDim Maximum dimension to which to expand the Rips complex. If set to -1, there is no limit.
 * @param p Order policy for the points. 
 * Can be either @ref Oscillating_rips_edge_range::Order_policy::FARTHEST_POINT_ORDERING, 
 * @ref Oscillating_rips_edge_range::Order_policy::ALREADY_ORDERED or 
 * @ref Oscillating_rips_edge_range::Order_policy::RANDOM_POINT_ORDERING. 
 * Default value: @ref Oscillating_rips_edge_range::Order_policy::FARTHEST_POINT_ORDERING.
 * @return The persistence diagram of the oscillating Rips filtration.
 */
template <typename PointRange, Edge_range_type edge_range_type = Edge_range_type::BOOST_RANGE,
          class OscillatingRipsSimplexRange = Oscillating_rips_simplex_range<
              Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_oscillating_rips>,
              typename std::conditional<
                  edge_range_type == Edge_range_type::BOOST_RANGE,
                  Oscillating_rips_edge_range<Gudhi::Simplex_tree<
                      Gudhi::Simplex_tree_options_oscillating_rips>::Filtration_value>::Oscillating_rips_edge_iterator,
                  std::vector<Zigzag_edge<Gudhi::Simplex_tree<
                      Gudhi::Simplex_tree_options_oscillating_rips>::Filtration_value> >::iterator>::type>,
          typename ZigzagPersistenceOptions = Gudhi::persistence_matrix::Zigzag_options<> >
std::vector<typename Zigzag_persistence<Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_oscillating_rips>,
                                        ZigzagPersistenceOptions>::filtration_value_interval>
compute_oscillating_rips_persistence(
    const PointRange& points, double nu, double mu, int maxDim,
    typename Oscillating_rips_edge_range<
        Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_oscillating_rips>::Filtration_value>::Order_policy p =
        Oscillating_rips_edge_range<Gudhi::Simplex_tree<
            Gudhi::Simplex_tree_options_oscillating_rips>::Filtration_value>::Order_policy::FARTHEST_POINT_ORDERING) 
{
  using ST = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_oscillating_rips>;
  using Filtration_value = ST::Filtration_value;
  using EdgeRange = Oscillating_rips_edge_range<Filtration_value>;

  Zigzag_persistence<ST> zp;
  ST& st = zp.get_complex();

  if constexpr (edge_range_type == Edge_range_type::BOOST_RANGE) {
    auto start = EdgeRange::begin(nu, mu, points, Gudhi::Euclidean_distance(), p);
    auto end = EdgeRange::end();
    for (const auto& t : OscillatingRipsSimplexRange::get_iterator_range(start, end, st, maxDim)) {
      if (std::get<2>(t))
        zp.insert_simplex_handle(std::get<0>(t));
      else
        zp.remove_simplex_handle(std::get<0>(t), std::get<1>(t));
    }
  } else {
    auto vec = EdgeRange::compute_vector_range(nu, mu, points, Gudhi::Euclidean_distance(), p);
    auto start = vec.begin();
    auto end = vec.end();
    for (const auto& t : OscillatingRipsSimplexRange::get_iterator_range(start, end, st, maxDim)) {
      if (std::get<2>(t))
        zp.insert_simplex_handle(std::get<0>(t));
      else
        zp.remove_simplex_handle(std::get<0>(t), std::get<1>(t));
    }
  }

  return zp.get_persistence_diagram();
}

}  // namespace zigzag_persistence

}  // namespace Gudhi

#endif  // OSCILLATING_RIPS_PERSISTENCE_H_

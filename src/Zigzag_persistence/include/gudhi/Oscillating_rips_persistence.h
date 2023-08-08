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

#include <gudhi/Debug_utils.h>
#include <gudhi/Zigzag_persistence/oscillating_rips_iterators.h>
#include <gudhi/Zigzag_persistence.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/options.h>
#include <gudhi/distance_functions.h>
#include <type_traits>

namespace Gudhi {
namespace zigzag_persistence {

enum Edge_range_type { VECTOR, BOOST_RANGE };

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
    //   std::cout << (std::get<2>(t) ? "in" : "out") << ": ";
    //   for (auto v : st.simplex_vertex_range(std::get<0>(t))) std::cout << v << " ";
    //   std::cout << " - ";
    //   for (auto sh : st.boundary_simplex_range(std::get<0>(t))) std::cout << st.key(sh) << " ";
    //   std::cout << "\n";
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

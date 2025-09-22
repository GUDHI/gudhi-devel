/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - 2025/06 Hannah Schreiber: Various small bug fixes (missing `inline`s, `DEBUG_TRACES`s etc.)
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PERSISTENCE_INTERVALS_WITH_DISTANCES_H_
#define PERSISTENCE_INTERVALS_WITH_DISTANCES_H_

#ifdef DEBUG_TRACES
#include <iostream>  // std::cerr
#endif
#include <limits>     // std::numeric_limits
#include <stdexcept>  // std::logic_error

#include <gudhi/Persistence_intervals.h>
#include <gudhi/Bottleneck.h>
#include <gudhi/Debug_utils.h>

namespace Gudhi {
namespace Persistence_representations {

// TODO: it would have been better to have this file in a subfolder "Persistence_representations"
// to avoid including it with "<gudhi/Persistence_intervals_with_distances.h>" which makes it sound universal within
// gudhi even though it is only used in this format within this module.
// How critical would it be for retro-compatibility to change that? It does not seem to appear in the documentation.

/**
 * @class Persistence_intervals_with_distances Persistence_intervals_with_distances.h \
 * gudhi/Persistence_intervals_with_distances.h
 * @brief This class implements the following concepts: Vectorized_topological_data, Topological_data_with_distances,
 * Real_valued_topological_data
 *
 * @ingroup Persistence_representations
 **/
class Persistence_intervals_with_distances : public Persistence_intervals
{
 public:
  using Persistence_intervals::Persistence_intervals;
  using Base = Persistence_intervals;

  /**
   * @brief Computes the distance to the persistence diagram given as a parameter.
   *
   * @param second Diagram to which to compute the distance.
   * @param power It is here in case we would like to compute \f$ power^{th} \f$ Wasserstein distance. At the
   * moment, this method only implement Bottleneck distance, which is infinity Wasserstein distance.
   * Therefore any power which is not the default `std::numeric_limits<double>::max()` will be ignored and an
   * exception will be thrown.
   * @param tolerance It is the additive error of the approximation, set by default to zero.
   */
  double distance(const Persistence_intervals_with_distances& second,
                  double power = std::numeric_limits<double>::max(),
                  double tolerance = (std::numeric_limits<double>::min)()) const
  {
    if (power >= std::numeric_limits<double>::max()) {
      return Gudhi::persistence_diagram::bottleneck_distance(Base::intervals_, second.intervals_, tolerance);
    } else {
#ifdef DEBUG_TRACES
      std::cerr << "At the moment Gudhi do not support Wasserstein distances. We only support Bottleneck distance."
                << std::endl;
#endif
      throw std::logic_error(
          "At the moment Gudhi do not support Wasserstein distances. We only support Bottleneck distance.");
    }
  }
};

}  // namespace Persistence_representations
}  // namespace Gudhi

#endif  // PERSISTENCE_INTERVALS_WITH_DISTANCES_H_

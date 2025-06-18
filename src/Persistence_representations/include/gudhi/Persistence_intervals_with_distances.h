/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - 2025/06 Hannah Schreiber: Divers small bug fixes (missing `inline`s, `GUDHI_DEBUG`s etc.)
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PERSISTENCE_INTERVALS_WITH_DISTANCES_H_
#define PERSISTENCE_INTERVALS_WITH_DISTANCES_H_

#ifdef GUDHI_DEBUG
#include <iostream>  // std::cerr
#endif
#include <limits>     // std::numeric_limits
#include <stdexcept>  // std::logic_error

#include <gudhi/Persistence_intervals.h>
#include <gudhi/Bottleneck.h>

namespace Gudhi {
namespace Persistence_representations {

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
#ifdef GUDHI_DEBUG
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

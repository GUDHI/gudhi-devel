/*    This file is part of the Gudhi hiLibrary. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PERSISTENCE_INTERVALS_WITH_DISTANCES_H_
#define PERSISTENCE_INTERVALS_WITH_DISTANCES_H_

#include <gudhi/Persistence_intervals.h>
#include <gudhi/Bottleneck.h>

#include <limits>

namespace Gudhi {
namespace Persistence_representations {

class Persistence_intervals_with_distances : public Persistence_intervals {
 public:
  using Persistence_intervals::Persistence_intervals;

  /**
   *Computations of distance from the current persistnce diagram to the persistence diagram given as a parameter of this
   *function.
   *The last but one parameter, power, is here in case we would like to compute p=th Wasserstein distance. At the
   *moment, this method only implement Bottleneck distance,
   * which is infinity Wasserstein distance. Therefore any power which is not the default std::numeric_limits< double
   *>::max() will be ignored and an
   * exception will be thrown.
   * The last parameter, tolerance, it is an additiv error of the approimation, set by default to zero.
  **/
  double distance(const Persistence_intervals_with_distances& second, double power = std::numeric_limits<double>::max(),
                  double tolerance = (std::numeric_limits<double>::min)()) const {
    if (power >= std::numeric_limits<double>::max()) {
      return Gudhi::persistence_diagram::bottleneck_distance(this->intervals, second.intervals, tolerance);
    } else {
      std::cerr << "At the moment Gudhi do not support Wasserstein distances. We only support Bottleneck distance."
                << std::endl;
      throw "At the moment Gudhi do not support Wasserstein distances. We only support Bottleneck distance.";
    }
  }
};

}  // namespace Persistence_representations
}  // namespace Gudhi

#endif  // PERSISTENCE_INTERVALS_WITH_DISTANCES_H_

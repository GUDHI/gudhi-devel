/*    This file is part of the Gudhi Library. The Gudhi library
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

#ifndef CONCEPT_TOPOLOGICAL_DATA_WITH_DISTANCES_H_
#define CONCEPT_TOPOLOGICAL_DATA_WITH_DISTANCES_H_

namespace Gudhi {

namespace Persistence_representations {

/** \brief The concept Topological_data_with_distances describes the requirements
  * for a type to implement a container that allows computations of distance to another contained of that type.
  * \details
  * The second parameter of the distance function allow to declare power of a distance. The exact meaning of that
  * number will be different for different distances. A few examples are given below:
  * In case of p-Wasserstein distance, the power is equal to p. power = std::limit<double>::max() for bottleneck
  * distance.
  *
  * In case of L^p landscape distance, the power is equal to p. s
  */
class Topological_data_with_distances {
 public:
  double distance(const Topological_data_with_distances& second, double power = 1);
};

}  // namespace Persistence_representations

}  // namespace Gudhi

#endif  // CONCEPT_TOPOLOGICAL_DATA_WITH_DISTANCES_H_

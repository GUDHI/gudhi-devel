/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2017  INRIA (France)
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
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

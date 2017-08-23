/*    This file is part of the Gudhi hiLibrary. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2017 Swansea University
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

#ifndef PERSISTENCE_INTERVALS_WITH_DISTANCES_INTERFACE_H_
#define PERSISTENCE_INTERVALS_WITH_DISTANCES_INTERFACE_H_

#include <gudhi/Persistence_intervals_with_distances.h>

namespace Gudhi {
namespace Persistence_representations {

class Persistence_intervals_with_distances_interface : public Persistence_intervals_with_distances {
 public:  
  double distance_interface(const Persistence_intervals_with_distances& second, double power = std::numeric_limits<double>::max(),
                  double tolerance = 0) const 
  {
	  return this->distance( second, power, tolerance );
  }  
};

}  // namespace Persistence_representations
}  // namespace Gudhi

#endif  // PERSISTENCE_INTERVALS_WITH_DISTANCES_INTERFACE_H_

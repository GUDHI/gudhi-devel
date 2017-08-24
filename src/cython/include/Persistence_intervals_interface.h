/*    This file is part of the Gudhi Library. The Gudhi library
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

#ifndef INCLUDE_PERSISTENCE_REPRESENTATIONS_INTERVALS_
#define INCLUDE_PERSISTENCE_REPRESENTATIONS_INTERVALS_

#include <gudhi/Persistence_intervals.h>

#include <iostream>
#include <vector>
#include <string>

namespace Gudhi {
	
//But if we want to have the same names of classes in C++ and cyton side we ned this interface, because othervise we will have a name conflict. And we want to have the same names on the 
//C++ and python side for various reasonc (clarity, documentantions etc.).
//If the C++ class we inherid from are template class, we are inherid from concretization, for instance Persistence_intervals<double>.
//Also in this class, we create an interface functions that will be used in the python side. That will allow to have the same name of the functions in the C++ and python side. 

namespace Persistence_representations {

class Persistence_intervals_interface : public Persistence_intervals 
{
 public: 
	Persistence_intervals_interface(const char* filename, unsigned dimension = std::numeric_limits<unsigned>::max())
	: Persistence_intervals(filename, dimension) {
  }

	Persistence_intervals_interface(const std::vector<std::pair<double, double> >& intervals)
	: Persistence_intervals(intervals) {
	}

};

}  // namespace Persistence_representations

}  // namespace Gudhi

#endif  // INCLUDE_PERSISTENCE_REPRESENTATIONS_DIAGRAMS_


/*    This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author:       François Godi
 *
 *    Copyright (C) 2015  INRIA
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

#ifndef CONCEPT_BOTTLENECK_DISTANCE_PERSISTENCE_DIAGRAM_H_
#define CONCEPT_BOTTLENECK_DISTANCE_PERSISTENCE_DIAGRAM_H_

namespace Gudhi {

namespace bottleneck_distance {

/** \brief Concept of Diagram_point. std::get<0>(point) must return the birth of the corresponding component and std::get<1>(point) its death.
 * A valid implementation of this concept is std::pair<double,double>.
 * Death should be larger than birth, death can be std::numeric_limits<double>::infinity() for components which stay alive.
 *
 * \ingroup bottleneck_distance
 */
typename Diagram_point;

/** \brief Concept of persistence diagram. It's a range of Diagram_point.
 * std::begin(diagram) and std::end(diagram) must return corresponding iterators.
 *
 * \ingroup bottleneck_distance
 */
typename Persistence_Diagram;

}  // namespace bottleneck_distance

}  // namespace Gudhi

#endif  // CONCEPT_BOTTLENECK_DISTANCE_PERSISTENCE_DIAGRAM_H_

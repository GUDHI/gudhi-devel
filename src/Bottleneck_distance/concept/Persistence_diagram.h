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

/** \brief Concept of persistence diagram point. get<0>() must return the birth of the component and get<1>() its death.
 *
 * \ingroup bottleneck_distance
 */
struct Diagram_point{
    double get<int>();
};

/** \brief Concept of persistence diagram. It's a range of Diagram_point.
 *
 * \ingroup bottleneck_distance
 */
struct Persistence_Diagram
{
    iterator<Diagram_point> begin();
    iterator<Diagram_point> end();
};

}  // namespace bottleneck_distance

}  // namespace Gudhi

#endif  // CONCEPT_BOTTLENECK_DISTANCE_PERSISTENCE_DIAGRAM_H_

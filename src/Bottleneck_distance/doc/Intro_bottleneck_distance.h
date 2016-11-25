/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author:       François Godi
 *
 *    Copyright (C) 2015  INRIA (France)
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

#ifndef DOC_BOTTLENECK_DISTANCE_H_
#define DOC_BOTTLENECK_DISTANCE_H_

// needs namespace for Doxygen to link on classes
namespace Gudhi {
// needs namespace for Doxygen to link on classes
namespace bottleneck_distance {

/**  \defgroup bottleneck_distance Bottleneck distance
 * 
 * \author    Fran&ccedil;ois Godi
 * @{
 * 
 * \section bottleneckdefinition Definition
 * 
 * Bottleneck distance measures the similarity between two persistence diagrams. 
 * It's the shortest distance b for which there exists a perfect matching between 
 * the points of the two diagrams (and also all the points of the diagonal) such that
 * any couple of matched points are at distance at most b.
 *
 * \image html perturb_pd.png Bottleneck distance is the length of the longest edge on this picture.
 *
 */
/** @} */  // end defgroup bottleneck_distance

}  // namespace bottleneck_distance

}  // namespace Gudhi


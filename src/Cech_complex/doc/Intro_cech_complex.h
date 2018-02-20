/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2018  Inria
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

#ifndef DOC_CECH_COMPLEX_INTRO_CECH_COMPLEX_H_
#define DOC_CECH_COMPLEX_INTRO_CECH_COMPLEX_H_

namespace Gudhi {

namespace cech_complex {

/**  \defgroup cech_complex Cech complex
 * 
 * \author    Clément Maria, Pawel Dlotko, Vincent Rouvreau
 * 
 * @{
 * 
 * \section cechdefinition Cech complex definition
 * 
 * Cech_complex
 * <a target="_blank" href="https://en.wikipedia.org/wiki/%C4%8Cech_cohomology">(Wikipedia)</a> is a
 * proximity graph that allows to construct a
 * <a target="_blank" href="https://en.wikipedia.org/wiki/Simplicial_complex">simplicial complex</a>
 * from it.
 * The input can be a point cloud with a given distance function.
 * 
 * The filtration value of each edge is computed from a user-given distance function.
 * 
 * All edges that have a filtration value strictly greater than a given threshold value are not inserted into
 * the complex.
 * 
 * When creating a simplicial complex from this proximity graph, Cech inserts the proximity graph into the data
 * structure, and then expands the simplicial complex when required.
 *
 * Vertex name correspond to the index of the point in the given range (aka. the point cloud).
 * 
 * \image html "cech_complex_representation.png" "Cech complex proximity graph representation"
 * 
 * On this example, as edges (4,5), (4,6) and (5,6) are in the complex, simplex (4,5,6) is added with the filtration
 * value set with \f$max(filtration(4,5), filtration(4,6), filtration(5,6))\f$.
 * And so on for simplex (0,1,2,3).
 * 
 * If the Cech_complex interfaces are not detailed enough for your need, please refer to
 * <a href="_persistent_cohomology_2cech_persistence_step_by_step_8cpp-example.html">
 * cech_persistence_step_by_step.cpp</a> example, where the graph construction over the Simplex_tree is more detailed.
 *
 * \section cechpointsdistance Point cloud and distance function
 * 
 * \subsection cechpointscloudexample Example from a point cloud and a distance function
 * 
 * This example builds the proximity graph from the given points, threshold value, and distance function.
 * Then it creates a `Simplex_tree` with it.
 * 
 * Then, it is asked to display information about the simplicial complex.
 * 
 * \include Cech_complex/cech_complex_example_from_points.cpp
 * 
 * When launching (Cech maximal distance between 2 points is 7.1, is expanded until dimension 2):
 * 
 * \code $> ./Cech_complex_example_from_points
 * \endcode
 *
 * the program output is:
 * 
 * \include Cech_complex/cech_complex_example_from_points_for_doc.txt
 * 
 */
/** @} */  // end defgroup cech_complex

}  // namespace cech_complex

}  // namespace Gudhi

#endif  // DOC_CECH_COMPLEX_INTRO_CECH_COMPLEX_H_

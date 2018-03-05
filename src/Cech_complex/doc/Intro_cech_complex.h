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
 * \author    Vincent Rouvreau
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
 * The input shall be a point cloud in an Euclidean space.
 * 
 * The filtration value of each edge of the `Gudhi::Proximity_graph` is computed from `Gudhi::Radius_distance` function.
 * 
 * All edges that have a filtration value strictly greater than a user given maximal radius value, \f$max\_radius\f$,
 * are not inserted into the complex.
 * 
 * Vertex name correspond to the index of the point in the given range (aka. the point cloud).
 * 
 * \image html "cech_one_skeleton.png" "Cech complex proximity graph representation"
 * 
 * When creating a simplicial complex from this proximity graph, Cech inserts the proximity graph into the simplicial
 * complex data structure, and then expands the simplicial complex when required.
 *
 * On this example, as edges \f$(x,y)\f$, \f$(y,z)\f$ and \f$(z,y)\f$ are in the complex, the minimal ball radius
 * containing the points \f$(x,y,z)\f$ is computed.
 *
 * \f$(x,y,z)\f$ is inserted to the simplicial complex with the filtration value set with
 * \f$mini\_ball\_radius(x,y,z))\f$ iff \f$mini\_ball\_radius(x,y,z)) \leq max\_radius\f$.
 *
 * And so on for higher dimensions.
 *
 * \image html "cech_complex_representation.png" "Cech complex expansion"
 *
 * The minimal ball radius computation is insured by
 * <a target="_blank" href="https://people.inf.ethz.ch/gaertner/subdir/software/miniball.html">
 * the miniball software (V3.0)</a> - Smallest Enclosing Balls of Points - and distributed with GUDHI.
 *
 * Please refer to
 * <a target="_blank" href="https://people.inf.ethz.ch/gaertner/subdir/texts/own_work/esa99_final.pdf">
 * the miniball software design description</a> for more information about this computation.
 *
 * If the Cech_complex interfaces are not detailed enough for your need, please refer to
 * <a href="_cech_complex_2cech_complex_step_by_step_8cpp-example.html">
 * cech_complex_step_by_step.cpp</a> example, where the graph construction over the Simplex_tree is more detailed.
 *
 * \section cechpointsdistance Point cloud
 * 
 * \subsection cechpointscloudexample Example from a point cloud
 * 
 * This example builds the proximity graph from the given points, and maximal radius values.
 * Then it creates a `Simplex_tree` with it.
 * 
 * Then, it is asked to display information about the simplicial complex.
 * 
 * \include Cech_complex/cech_complex_example_from_points.cpp
 * 
 * When launching (Cech maximal distance between 2 points is 1., is expanded until dimension 2):
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

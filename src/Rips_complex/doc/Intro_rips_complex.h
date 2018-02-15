/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clément Maria, Pawel Dlotko, Vincent Rouvreau
 *
 *    Copyright (C) 2016  INRIA
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

#ifndef DOC_RIPS_COMPLEX_INTRO_RIPS_COMPLEX_H_
#define DOC_RIPS_COMPLEX_INTRO_RIPS_COMPLEX_H_

namespace Gudhi {

namespace rips_complex {

/**  \defgroup rips_complex Rips complex
 * 
 * \author    Clément Maria, Pawel Dlotko, Vincent Rouvreau
 * 
 * @{
 * 
 * \section ripsdefinition Rips complex definition
 * 
 * Rips_complex
 * <a target="_blank" href="https://en.wikipedia.org/wiki/Vietoris%E2%80%93Rips_complex">(Wikipedia)</a> is a
 * one skeleton graph that allows to construct a
 * <a target="_blank" href="https://en.wikipedia.org/wiki/Simplicial_complex">simplicial complex</a>
 * from it.
 * The input can be a point cloud with a given distance function, or a distance matrix.
 * 
 * The filtration value of each edge is computed from a user-given distance function, or directly from the distance
 * matrix.
 * 
 * All edges that have a filtration value strictly greater than a given threshold value are not inserted into
 * the complex.
 * 
 * When creating a simplicial complex from this one skeleton graph, Rips inserts the one skeleton graph into the data
 * structure, and then expands the simplicial complex when required.
 *
 * Vertex name correspond to the index of the point in the given range (aka. the point cloud).
 * 
 * \image html "rips_complex_representation.png" "Rips-complex one skeleton graph representation"
 * 
 * On this example, as edges (4,5), (4,6) and (5,6) are in the complex, simplex (4,5,6) is added with the filtration
 * value set with \f$max(filtration(4,5), filtration(4,6), filtration(5,6))\f$.
 * And so on for simplex (0,1,2,3).
 * 
 * If the Rips_complex interfaces are not detailed enough for your need, please refer to
 * <a href="_persistent_cohomology_2rips_persistence_step_by_step_8cpp-example.html">
 * rips_persistence_step_by_step.cpp</a> example, where the graph construction over the Simplex_tree is more detailed.
 *
 * \section ripspointsdistance Point cloud and distance function
 * 
 * \subsection ripspointscloudexample Example from a point cloud and a distance function
 * 
 * This example builds the one skeleton graph from the given points, threshold value, and distance function.
 * Then it creates a `Simplex_tree` with it.
 * 
 * Then, it is asked to display information about the simplicial complex.
 * 
 * \include Rips_complex/example_one_skeleton_rips_from_points.cpp
 * 
 * When launching (Rips maximal distance between 2 points is 12.0, is expanded until dimension 1 - one skeleton graph
 * in other words):
 * 
 * \code $> ./Rips_complex_example_one_skeleton_from_points
 * \endcode
 *
 * the program output is:
 * 
 * \include Rips_complex/one_skeleton_rips_for_doc.txt
 * 
 * \subsection ripsoffexample Example from OFF file
 * 
 * This example builds the Rips_complex from the given points in an OFF file, threshold value, and distance
 * function.
 * Then it creates a `Simplex_tree` with it.
 * 
 * 
 * Then, it is asked to display information about the Rips complex.
 * 
 * \include Rips_complex/example_rips_complex_from_off_file.cpp
 * 
 * When launching:
 * 
 * \code $> ./Rips_complex_example_from_off ../../data/points/alphacomplexdoc.off 12.0 3
 * \endcode
 *
 * the program output is:
 * 
 * \include Rips_complex/full_skeleton_rips_for_doc.txt
 * 
 * 
 * 
 * \section ripsdistancematrix Distance matrix
 * 
 * \subsection ripsdistancematrixexample Example from a distance matrix
 * 
 * This example builds the one skeleton graph from the given distance matrix and threshold value.
 * Then it creates a `Simplex_tree` with it.
 * 
 * Then, it is asked to display information about the simplicial complex.
 * 
 * \include Rips_complex/example_one_skeleton_rips_from_distance_matrix.cpp
 * 
 * When launching (Rips maximal distance between 2 points is 1.0, is expanded until dimension 1 - one skeleton graph
 * with other words):
 * 
 * \code $> ./Rips_complex_example_one_skeleton_from_distance_matrix
 * \endcode
 *
 * the program output is:
 * 
 * \include Rips_complex/one_skeleton_rips_for_doc.txt
 * 
 * \subsection ripscsvdistanceexample Example from a distance matrix read in a csv file
 * 
 * This example builds the one skeleton graph from the given distance matrix read in a csv file and threshold value.
 * Then it creates a `Simplex_tree` with it.
 * 
 * 
 * Then, it is asked to display information about the Rips complex.
 * 
 * \include Rips_complex/example_rips_complex_from_csv_distance_matrix_file.cpp
 * 
 * When launching:
 * 
 * \code $> ./Rips_complex_example_from_csv_distance_matrix ../../data/distance_matrix/full_square_distance_matrix.csv 1.0 3
 * \endcode
 *
 * the program output is:
 * 
 * \include Rips_complex/full_skeleton_rips_for_doc.txt
 * 
 */
/** @} */  // end defgroup rips_complex

}  // namespace rips_complex

}  // namespace Gudhi

#endif  // DOC_RIPS_COMPLEX_INTRO_RIPS_COMPLEX_H_

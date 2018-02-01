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
 * \author    Clément Maria, Pawel Dlotko, Vincent Rouvreau, Marc Glisse
 * 
 * @{
 * 
 * \section ripsdefinition Rips complex definition
 * 
 * The Vietoris-Rips complex
 * <a target="_blank" href="https://en.wikipedia.org/wiki/Vietoris%E2%80%93Rips_complex">(Wikipedia)</a> is an abstract simplicial complex
 * defined on a finite metric space, where each simplex corresponds to a subset
 * of point whose diameter is smaller that some given threshold.
 * Varying the threshold, we can also see the Rips complex as a filtration of
 * the \f$(n-1)-\f$dimensional simplex, where the filtration value of each
 * simplex is the diameter of the corresponding subset of points.
 *
 * This filtered complex is most often used as an approximation of the
 * Čech complex. After rescaling (Rips using the length of the edges and Čech
 * the half-length), they share the same 1-skeleton and are multiplicatively
 * 2-interleaved or better. While it is slightly bigger, it is also much
 * easier to compute.
 *
 * The number of simplices in the full Rips complex is exponential in the
 * number of vertices, it is thus usually restricted, by excluding all the
 * simplices with filtration value larger than some threshold, and keeping only
 * the dim_max-skeleton.
 *
 * In order to build this complex, the algorithm first builds the graph.
 * The filtration value of each edge is computed from a user-given distance
 * function, or directly read from the distance matrix.
 * In a second step, this graph is inserted in a simplicial complex, which then
 * gets expanded to a flag complex.
 * 
 * The input can be given as a range of points and a distance function, or as a
 * distance matrix.
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
 * rips_persistence_step_by_step.cpp</a> example, where the constructions of the graph and the Simplex_tree are more detailed.
 *
 * \section sparserips Sparse Rips complex
 *
 * Even truncated in filtration value and dimension, the Rips complex remains
 * quite large. However, it is possible to approximate it by a much smaller
 * filtered simplicial complex (linear size, with constants that depend on
 * &epsilon; and the doubling dimension of the space) that is
 * \f$(1+O(\epsilon))-\f$interleaved with it (in particular, their persistence
 * diagrams are at log-bottleneck distance at most &epsilon;).
 *
 * The sparse Rips filtration was introduced by Don Sheehy
 * \cite sheehy13linear. We are using the version from \cite buchet16efficient
 * (except that we multiply all filtration values by 2, to match the usual
 * Rips complex).
 * A more intuitive presentation of the idea is available in
 * \cite cavanna15geometric, and in a video \cite cavanna15visualizing.
 *
 * The interface of `Sparse_rips_complex` is similar to the one for the usual
 * `Rips_complex`, except that one has to specify the approximation factor, and
 * there is no option to limit the maximum filtration value (the way the
 * approximation is done means that larger filtration values are much cheaper
 * to handle than low filtration values, so the gain would be too small).
 *
 * Theoretical guarantees are only available for \f$\epsilon<1\f$. The
 * construction accepts larger values of &epsilon;, and the size of the complex
 * keeps decreasing, but there is no guarantee on the quality of the result.
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
 * \subsection sparseripspointscloudexample Example of a sparse Rips from a point cloud
 * 
 * This example builds the full sparse Rips of a set of 2D Euclidean points, then prints some minimal information about the complex.
 * 
 * \include Rips_complex/example_sparse_rips.cpp
 * 
 * When launching:
 * 
 * \code $> ./Rips_complex_example_sparse
 * \endcode
 *
 * the program output may be (the exact output varies from one run to the next):
 *
 * \code Sparse Rips complex is of dimension 2 - 19 simplices - 7 vertices.
 * \endcode
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

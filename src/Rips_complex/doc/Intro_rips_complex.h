/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clément Maria & Vincent Rouvreau
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

#ifndef DOC_RIPS_COMPLEX_INTRO_RIPS_COMPLEX_H_
#define DOC_RIPS_COMPLEX_INTRO_RIPS_COMPLEX_H_

namespace Gudhi {

namespace rips_complex {

/**  \defgroup rips_complex Rips complex
 * 
 * \author    Clément Maria & Vincent Rouvreau
 * 
 * @{
 * 
 * \section definition Definition
 * 
 * Rips_complex
 * <a target="_blank" href="https://en.wikipedia.org/wiki/Vietoris%E2%80%93Rips_complex">(Wikipedia)</a> is a
 * <a target="_blank" href="https://en.wikipedia.org/wiki/Simplicial_complex">simplicial complex</a>
 * constructed from a one skeleton graph.
 * 
 * The filtration value of each edge is computed from a user-given distance function.
 * 
 * All edges that have a filtration value strictly greater than a given threshold value are not inserted into
 * the complex.
 * 
 * When creating a simplicial complex from this one skeleton graph, rips inserts the one skeleton graph into the data
 * structure, and then expands the simplicial when required.
 * 
 * \image html "rips_complex_representation.png" "Rips-complex one skeleton graph representation"
 * 
 * On this example, as edges (4,5), (4,6) and (5,6) are in the complex, simplex (4,5,6) is added with the filtration
 * value set with \f$max(filtration(4,5), filtration(4,6), filtration(5,6))\f$.
 * And so on for simplex (0,1,2,3).
 * 
 * \section ripspointsexample Example from points
 * 
 * This example builds the one skeleton graph from the given points, threshold value, and distance function.
 * Then it creates a `Simplex_tree` with it.
 * 
 * Then, it is asked to display information about the simplicial complex.
 * 
 * \include Rips_complex/example_rips_complex_from_points.cpp
 * 
 * When launching:
 * 
 * \code $> ./ripspoints
 * \endcode
 *
 * the program output is:
 * 
 * \include Rips_complex/rips_points_for_doc_12_2.txt
 * 
 * \section offexample Example from OFF file
 * 
 * This example builds the one skeleton graph from the given points in an OFF file, threshold value, and distance
 * function.
 * Then it creates a `Simplex_tree` with it.
 * 
 * 
 * Then, it is asked to display information about the rips complex.
 * 
 * \include Rips_complex/example_rips_complex_from_off_file.cpp
 * 
 * When launching:
 * 
 * \code $> ./ripsoffreader ../../data/points/alphacomplexdoc.off 12.0 3
 * \endcode
 *
 * the program output is:
 * 
 * \include Rips_complex/rips_points_for_doc_12_3.txt
 * 
 * \copyright GNU General Public License v3.                         
 * \verbatim  Contact: gudhi-users@lists.gforge.inria.fr \endverbatim
 */
/** @} */  // end defgroup rips_complex

}  // namespace rips_complex

}  // namespace Gudhi

#endif  // DOC_RIPS_COMPLEX_INTRO_RIPS_COMPLEX_H_

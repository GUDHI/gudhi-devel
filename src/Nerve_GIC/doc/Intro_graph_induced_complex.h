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

#ifndef DOC_GRAPH_INDUCED_COMPLEX_INTRO_GRAPH_INDUCED_COMPLEX_H_
#define DOC_GRAPH_INDUCED_COMPLEX_INTRO_GRAPH_INDUCED_COMPLEX_H_

namespace Gudhi {

namespace graph_induced_complex {

/**  \defgroup graph_induced_complex Graph induced complex
 * 
 * \author    Mathieu Carrière
 * 
 * @{
 * 
 * \section complexes Graph induced complexes (GIC) and Nerves
 * 
 * GIC and Nerves are simplicial complexes built on top of a point cloud P.
 *
 * \subsection nervedefinition Nerve definition
 *
 * Assume you are given a cover C of your point cloud P, that is a set of subsets of P
 * whose union is P itself. Then, the Nerve of this cover
 * is the simplicial complex that has one k-simplex per k-fold intersection of cover elements.
 * See also <a target="_blank" href="https://en.wikipedia.org/wiki/Nerve_of_a_covering"> Wikipedia </a>.
 *
 * \subsection nerveexample Example
 *
 * This example builds the Nerve of a point cloud sampled on a 3D human shape.
 * The cover C comes from the preimages of intervals covering the height function.
 * All intervals have the resolution (either the length or the number of the intervals)
 * and gain (overlap percentage).
 *
 * \include
 *
 * When launching:
 *
 * \code $>
 * \endcode
 *
 * the program output is:
 *
 * \include
 *
 * \section gicdefinition GIC definition
 *
 * Again, assume you are given a cover C of your point cloud P. Moreover, assume
 * you are also given a graph G built on top of P. Then, for any clique in G
 * whose nodes all belong to different elements of C, the GIC includes a corresponding
 * simplex, whose dimension is the number of nodes in the clique minus one.
 *
 * \subsection gicexample Example
 *
 * This example builds the GIC of a point cloud sampled on a 3D human shape.
 * The cover C comes from the preimages of intervals covering the height function,
 * and the graph G comes from a Rips complex built with a threshold parameter.
 * Note that if the gain is too big, the number of cliques increases a lot,
 * which make the computation time much larger.
 *
 * \include
 *
 * When launching:
 *
 * \code $>
 * \endcode
 *
 * the program output is:
 *
 * \include
 *
 * \subsection mapperdeltadefinition Mapper Delta
 *
 * If one restricts to the cliques in G whose nodes all belong to preimages of consecutive
 * intervals (assuming the cover of the height function is minimal, i.e. no more than
 * two intervals can intersect at a time), the GIC is of dimension one, i.e. a graph.
 * We call this graph the Mapper Delta, since it is related to the usual Mapper (see
 * <a target="_blank" href="https://arxiv.org/abs/1511.05823"> this article </a>).
 *
 * \subsection mapperdeltaexample Example
 *
 * Mapper Delta comes with optimal selection for the Rips threshold,
 * the resolution and the gain of the function cover. In this example,
 * we compute the Mapper Delta of a point cloud sampled on a 3D human shape,
 * where the graph G comes from a Rips complex with optimal threshold,
 * and the cover C comes from the preimages of intervals covering the height function,
 * with optimal resolution and gain. Note that optimal threshold, resolution and gain
 * also exist for the Nerve of this cover.
 *
 * \include
 *
 * When launching:
 *
 * \code $>
 * \endcode
 *
 * the program output is:
 *
 * \include
 *
 * 
 * \copyright GNU General Public License v3.                         
 * \verbatim  Contact: gudhi-users@lists.gforge.inria.fr \endverbatim
 */
/** @} */  // end defgroup graph_induced_complex

}  // namespace graph_induced_complex

}  // namespace Gudhi

#endif  // DOC_GRAPH_INDUCED_COMPLEX_INTRO_GRAPH_INDUCED_COMPLEX_H_

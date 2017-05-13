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
 * \section covers Covers
 *
 * Nerves and Graph Induced Complexes require a cover C of the input point cloud P,
 * that is a set of subsets of P whose union is P itself.
 * Very often, this cover is obtained from the preimage of a family of intervals covering
 * the image of some scalar-valued function f defined on P. This family is parameterized
 * by its resolution, which can be either the number or the length of the intervals,
 * and its gain, which is the overlap percentage between consecutive intervals (ordered by their first values).
 *
 * \section nerves Nerves
 *
 * \subsection nervedefinition Nerve definition
 *
 * Assume you are given a cover C of your point cloud P. Then, the Nerve of this cover
 * is the simplicial complex that has one k-simplex per k-fold intersection of cover elements.
 * See also <a target="_blank" href="https://en.wikipedia.org/wiki/Nerve_of_a_covering"> Wikipedia </a>.
 *
 * \image html "nerve.png" "Nerve of a double torus"
 *
 * \subsection nerveexample Example
 *
 * This example builds the Nerve of a point cloud sampled on a 3D human shape (human.off).
 * The cover C comes from the preimages of intervals (10 intervals with gain 0.3)
 * covering the height function (coordinate 2),
 * which are then refined into their connected components using the triangulation of the .OFF file.
 *
 * \include Nerve_GIC/Nerve.cpp
 *
 * When launching:
 *
 * \code $> ./Nerve ../../../data/points/human.off 2 10 0.3
 * \endcode
 *
 * the program output is:
 *
 * \include Nerve_GIC/Nerve.txt
 *
 * The first three lines are requirements for visualization with Kepler-Mapper.
 * The fourth line contains the number of vertices nv and edges ne of the Nerve.
 * The next nv lines represent the vertices. Each line contains the vertex ID,
 * the number of data points it contains, and their average color function value.
 * Finally, the next ne lines represent the edges, characterized by the ID of their vertices.
 *
 *
 * \section gic Graph Induced Complexes (GIC)
 *
 * \subsection gicdefinition GIC definition
 *
 * Again, assume you are given a cover C of your point cloud P. Moreover, assume
 * you are also given a graph G built on top of P. Then, for any clique in G
 * whose nodes all belong to different elements of C, the GIC includes a corresponding
 * simplex, whose dimension is the number of nodes in the clique minus one.
 * See <a target="_blank" href="https://arxiv.org/abs/1304.0662"> this article </a>
 * for more details.
 *
 * \image html "gic_complex.png" "GIC of a point cloud."
 *
 * \subsection gicexample Example
 *
 * This example builds the GIC of a point cloud sampled on a 3D human shape (human.off).
 * The cover C comes from the preimages of intervals (with length 0.075 and gain 0)
 * covering the height function (coordinate 2),
 * and the graph G comes from a Rips complex built with threshold 0.075.
 * Note that if the gain is too big, the number of cliques increases a lot,
 * which make the computation time much larger.
 *
 * \include Nerve_GIC/GIC.cpp
 *
 * When launching:
 *
 * \code $> ./GIC ../../../data/points/human.off 0.075 2 0.075 0
 * \endcode
 *
 * the program output is:
 *
 * \include Nerve_GIC/GIC.txt
 *
 * \subsection mapperdeltadefinition Mapper Delta
 *
 * If one restricts to the cliques in G whose nodes all belong to preimages of consecutive
 * intervals (assuming the cover of the height function is minimal, i.e. no more than
 * two intervals can intersect at a time), the GIC is of dimension one, i.e. a graph.
 * We call this graph the Mapper Delta, since it is related to the usual Mapper. See
 * <a target="_blank" href="https://arxiv.org/abs/1511.05823"> this article </a> for more details.
 *
 * \subsection mapperdeltaexample Example
 *
 * Mapper Delta comes with optimal selection for the Rips threshold,
 * the resolution and the gain of the function cover. In this example,
 * we compute the Mapper Delta of a point cloud sampled on a 3D human shape (human.off),
 * where the graph G comes from a Rips complex with optimal threshold,
 * and the cover C comes from the preimages of intervals covering the height function (coordinate 2),
 * with optimal resolution and gain. Note that optimal threshold, resolution and gain
 * also exist for the Nerve of this cover.
 *
 * \include Nerve_GIC/MapperDeltaCoord.cpp
 *
 * When launching:
 *
 * \code $> ./MapperDeltaCoord ../../../data/points/human.off 2
 * \endcode
 *
 * the program output is:
 *
 * \include MapperDeltaCoord.txt
 *
 * We also provide an example on a set of 72 pictures taken around the same object (lucky_cat.off).
 * The function is now the first eigenfunction given by PCA, whose values
 * are written in a file (lucky_cat_PCA1). Threshold, resolution and gain are automatically selected as before.
 *
 * \include Nerve_GIC/MapperDeltaFunc.cpp
 *
 * When launching:
 *
 * \code $> ./MapperDeltaFunc ../../../data/points/COIL_database/lucky_cat.off ../../../data/points/COIL_database/lucky_cat_PCA1
 * \endcode
 *
 * the program output is:
 *
 * \include MapperDeltaFunc.txt
 *
 * If you have python and firefox, all the previous .txt files can then be displayed in a browser.
 * We provide a .py file called visu.py that comes from the
 * <a target="_blank" href="https://github.com/MLWave/kepler-mapper"> Kepler-Mapper </a> library.
 * One can visualize data by launching:
 *
 * \code python visu.py && firefox SC.html
 * \endcode
 *
 * \copyright GNU General Public License v3.                         
 * \verbatim  Contact: gudhi-users@lists.gforge.inria.fr \endverbatim
 */
/** @} */  // end defgroup graph_induced_complex

}  // namespace graph_induced_complex

}  // namespace Gudhi

#endif  // DOC_GRAPH_INDUCED_COMPLEX_INTRO_GRAPH_INDUCED_COMPLEX_H_

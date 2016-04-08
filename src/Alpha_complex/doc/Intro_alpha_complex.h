/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2015  INRIA Saclay (France)
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

#ifndef DOC_ALPHA_COMPLEX_INTRO_ALPHA_COMPLEX_H_
#define DOC_ALPHA_COMPLEX_INTRO_ALPHA_COMPLEX_H_

// needs namespace for Doxygen to link on classes
namespace Gudhi {
// needs namespace for Doxygen to link on classes
namespace alphacomplex {

/**  \defgroup alpha_complex Alpha complex
 * 
 * \author    Vincent Rouvreau
 * 
 * @{
 * 
 * \section definition Definition
 * 
 * Alpha_complex is a <a target="_blank" href="https://en.wikipedia.org/wiki/Simplicial_complex">simplicial complex</a>
 * constructed from the finite cells of a Delaunay Triangulation.
 * 
 * The filtration value of each simplex is computed as the square of the circumradius of the simplex if the circumsphere is empty (the simplex is then said to be Gabriel), and as the minimum of the filtration
 * values of the codimension 1 cofaces that make it not Gabriel otherwise.
 * 
 * All simplices that have a filtration value strictly greater than a given alpha squared value are not inserted into
 * the complex.
 * 
 * \image html "alpha_complex_representation.png" "Alpha-complex representation"
 * 
 * Alpha_complex is constructing a `Simplex_tree` using <a target="_blank"
 * href="http://doc.cgal.org/latest/Triangulation/index.html#Chapter_Triangulations">Delaunay Triangulation</a>
 * \cite cgal:hdj-t-15b from <a target="_blank" href="http://www.cgal.org/">CGAL</a> (the Computational Geometry
 * Algorithms Library \cite cgal:eb-15b).
 * 
 * The complex is a template class requiring an Epick_d <a target="_blank"
 * href="http://doc.cgal.org/latest/Kernel_d/index.html#Chapter_dD_Geometry_Kernel">dD Geometry Kernel</a>
 * \cite cgal:s-gkd-15b from CGAL as template parameter.
 * 
 * \remark When Alpha_complex is constructed with an infinite value of alpha, the complex is a Delaunay complex.
 * 
 * \section pointsexample Example from points
 * 
 * This example builds the Delaunay triangulation from the given points in a 2D static kernel, and initializes the
 * alpha complex with it.
 * 
 * Then, it is asked to display information about the alpha complex.
 * 
 * \include Alpha_complex/Alpha_complex_from_points.cpp
 * 
 * When launching:
 * 
 * \code $> ./alphapoints
 * \endcode
 *
 * the program output is:
 * 
 * \include Alpha_complex/alphaoffreader_for_doc_60.txt
 * 
 * \section algorithm Algorithm
 * 
 * \subsection datastructure Data structure
 * 
 * In order to build the alpha complex, first, a Simplex tree is built from the cells of a Delaunay Triangulation.
 * (The filtration value is set to NaN, which stands for unknown value):
 * \image html "alpha_complex_doc.png" "Simplex tree structure construction example"
 *
 * \subsection filtrationcomputation Filtration value computation algorithm
 * 
 * \f{algorithm}{ 
 * \caption{Filtration value computation algorithm}\label{alpha}
 * \begin{algorithmic}
 * \For{i : dimension $\rightarrow$ 0}
 *   \ForAll{$\sigma$ of dimension i}
 *     \If {filtration($\sigma$) is NaN}
 *       \State filtration($\sigma$) = $\alpha^2(\sigma)$
 *     \EndIf
 *     \ForAll{$\tau$ face of $\sigma$} \Comment{propagate alpha filtration value}
 *       \If {filtration($\tau$) is not NaN} 
 *         \State filtration($\tau$) = min (filtration($\tau$), filtration($\sigma$))
 *       \Else
 *         \If {$\tau$ is not Gabriel for $\sigma$} 
 *           \State filtration($\tau$) = filtration($\sigma$)
 *         \EndIf
 *       \EndIf
 *     \EndFor
 *   \EndFor
 * \EndFor
 * \State make\_filtration\_non\_decreasing()
 * \State prune\_above\_filtration()
 * \end{algorithmic}
 * \f}
 * 
 * \subsubsection dimension2 Dimension 2
 * 
 * From the example above, it means the algorithm looks into each triangle ([0,1,2], [0,2,4], [1,2,3], ...),
 * computes the filtration value of the triangle, and then propagates the filtration value as described
 * here :
 * \image html "alpha_complex_doc_420.png" "Filtration value propagation example"
 * 
 * \subsubsection dimension1 Dimension 1
 * 
 * Then, the algorithm looks into each edge ([0,1], [0,2], [1,2], ...),
 * computes the filtration value of the edge (in this case, propagation will have no effect).
 * 
 * \subsubsection dimension0 Dimension 0
 * 
 * Finally, the algorithm looks into each vertex ([0], [1], [2], [3], [4], [5] and [6]) and
 * sets the filtration value (0 in case of a vertex - propagation will have no effect).
 * 
 * \subsubsection nondecreasing Non decreasing filtration values
 * 
 * As the squared radii computed by CGAL are an approximation, it might happen that these alpha squared values do not quite define a proper filtration (i.e. non-decreasing with respect to inclusion).
 * We fix that up by calling `Simplex_tree::make_filtration_non_decreasing()`.
 * 
 * \subsubsection pruneabove Prune above given filtration value
 * 
 * The simplex tree is pruned from the given maximum alpha squared value (cf. `Simplex_tree::prune_above_filtration()`).
 * In the following example, the value is given by the user as argument of the program.
 * 
 * 
 * \section offexample Example from OFF file
 * 
 * This example builds the Delaunay triangulation in a dynamic kernel, and initializes the alpha complex with it.
 * 
 * 
 * Then, it is asked to display information about the alpha complex.
 * 
 * \include Alpha_complex/Alpha_complex_from_off.cpp
 * 
 * When launching:
 * 
 * \code $> ./alphaoffreader ../../data/points/alphacomplexdoc.off 32.0
 * \endcode
 *
 * the program output is:
 * 
 * \include Alpha_complex/alphaoffreader_for_doc_32.txt
 * 
 * \copyright GNU General Public License v3.                         
 * \verbatim  Contact: gudhi-users@lists.gforge.inria.fr \endverbatim
 */
/** @} */  // end defgroup alpha_complex

}  // namespace alphacomplex

}  // namespace Gudhi

#endif  // DOC_ALPHA_COMPLEX_INTRO_ALPHA_COMPLEX_H_

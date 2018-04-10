/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
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

#ifndef DOC_ALPHA_COMPLEX_INTRO_ALPHA_COMPLEX_H_
#define DOC_ALPHA_COMPLEX_INTRO_ALPHA_COMPLEX_H_

// needs namespace for Doxygen to link on classes
namespace Gudhi {
// needs namespace for Doxygen to link on classes
namespace alpha_complex {

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
 * The filtration value of each simplex is computed as the square of the circumradius of the simplex if the
 * circumsphere is empty (the simplex is then said to be Gabriel), and as the minimum of the filtration
 * values of the codimension 1 cofaces that make it not Gabriel otherwise.
 * 
 * All simplices that have a filtration value strictly greater than a given alpha squared value are not inserted into
 * the complex.
 * 
 * \image html "alpha_complex_representation.png" "Alpha-complex representation"
 * 
 * Alpha_complex is constructing a <a target="_blank"
 * href="http://doc.cgal.org/latest/Triangulation/index.html#Chapter_Triangulations">Delaunay Triangulation</a>
 * \cite cgal:hdj-t-15b from <a target="_blank" href="http://www.cgal.org/">CGAL</a> (the Computational Geometry
 * Algorithms Library \cite cgal:eb-15b) and is able to create a `SimplicialComplexForAlpha`.
 * 
 * The complex is a template class requiring an Epick_d <a target="_blank"
 * href="http://doc.cgal.org/latest/Kernel_d/index.html#Chapter_dD_Geometry_Kernel">dD Geometry Kernel</a>
 * \cite cgal:s-gkd-15b from CGAL as template parameter.
 * 
 * \remark
 * - When the simplicial complex is constructed with an infinite value of alpha, the complex is a Delaunay
 * complex.
 * - For people only interested in the topology of the \ref alpha_complex (for instance persistence),
 * \ref alpha_complex is equivalent to the \ref cech_complex and much smaller if you do not bound the radii.
 * \ref cech_complex can still make sense in higher dimension precisely because you can bound the radii.
 *
 * \section pointsexample Example from points
 * 
 * This example builds the Delaunay triangulation from the given points in a 2D static kernel, and creates a
 * `Simplex_tree` with it.
 * 
 * Then, it is asked to display information about the simplicial complex.
 * 
 * \include Alpha_complex/Alpha_complex_from_points.cpp
 * 
 * When launching:
 * 
 * \code $> ./Alpha_complex_example_from_points
 * \endcode
 *
 * the program output is:
 * 
 * \include Alpha_complex/alphaoffreader_for_doc_60.txt
 * 
 * \section createcomplexalgorithm Create complex algorithm
 * 
 * \subsection datastructure Data structure
 * 
 * In order to create the simplicial complex, first, it is built from the cells of the Delaunay Triangulation.
 * The filtration values are set to NaN, which stands for unknown value.
 * 
 * In example, :
 * \image html "alpha_complex_doc.png" "Simplicial complex structure construction example"
 *
 * \subsection filtrationcomputation Filtration value computation algorithm
 *
 *
 *
 * <ul>
 *   <li style="list-style-type: none;">\f$ \textbf{for } i : dimension \rightarrow 0 \textbf{ do} \f$
 *     <ul>
 *       <li style="list-style-type: none;">\f$\textbf{for all } \sigma of dimension i \f$
 *         <ul>
 *           <li style="list-style-type: none;">\f$\textbf{if } filtration( \sigma ) is NaN \textbf{ then} \f$
 *             <ul>
 *               <li style="list-style-type: none;">\f$ filtration( \sigma ) = \alpha^2( \sigma ) \f$
 *               </li>
 *             </ul>
 *           </li>
 *           <li style="list-style-type: none;">\f$\textbf{end if}\f$
 *           </li>
 *           <li style="list-style-type: none;">\f$\textbf{for all } \tau face of \sigma \textbf{ do} \f$
 *             &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;// propagate alpha filtration value
 *             <ul>
 *               <li style="list-style-type: none;">\f$\textbf{if } filtration( \tau ) is not NaN \textbf{ then} \f$
 *                 <ul>
 *                   <li style="list-style-type: none;">\f$ filtration( \tau ) = min ( filtration( \tau ), filtration( \sigma ) ) \f$
 *                   </li>
 *                 </ul>
 *               </li>
 *               <li style="list-style-type: none;">\f$\textbf{else}\f$
 *                 <ul>
 *                   <li style="list-style-type: none;">\f$\textbf{if } \tau is not Gabriel for \sigma \textbf{ then} \f$
 *                     <ul>
 *                       <li style="list-style-type: none;">\f$ filtration( \tau ) = filtration( \sigma ) \f$
 *                       </li>
 *                     </ul>
 *                   </li>
 *                   <li style="list-style-type: none;">\f$\textbf{end if}\f$
 *                   </li>
 *                 </ul>
 *               </li>
 *               <li style="list-style-type: none;">\f$\textbf{end if}\f$
 *               </li>
 *             </ul>
 *           </li>
 *           <li style="list-style-type: none;">\f$\textbf{end for}\f$
 *           </li>
 *         </ul>
 *       </li>
 *       <li style="list-style-type: none;">\f$\textbf{end for}\f$
 *       </li>
 *     </ul>
 *   </li>
 *   <li style="list-style-type: none;">\f$\textbf{end for}\f$
 *   </li>
 *   <li style="list-style-type: none;">\f$make\_filtration\_non\_decreasing()\f$
 *   </li>
 *   <li style="list-style-type: none;">\f$prune\_above\_filtration()\f$
 *   </li>
 * </ul>
 *
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
 * As the squared radii computed by CGAL are an approximation, it might happen that these alpha squared values do not
 * quite define a proper filtration (i.e. non-decreasing with respect to inclusion).
 * We fix that up by calling `SimplicialComplexForAlpha::make_filtration_non_decreasing()`.
 * 
 * \subsubsection pruneabove Prune above given filtration value
 * 
 * The simplex tree is pruned from the given maximum alpha squared value (cf.
 * `SimplicialComplexForAlpha::prune_above_filtration()`).
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
 * \code $> ./Alpha_complex_example_from_off ../../data/points/alphacomplexdoc.off 32.0
 * \endcode
 *
 * the program output is:
 * 
 * \include Alpha_complex/alphaoffreader_for_doc_32.txt
 * 
 */
/** @} */  // end defgroup alpha_complex

}  // namespace alpha_complex

namespace alphacomplex = alpha_complex;

}  // namespace Gudhi

#endif  // DOC_ALPHA_COMPLEX_INTRO_ALPHA_COMPLEX_H_

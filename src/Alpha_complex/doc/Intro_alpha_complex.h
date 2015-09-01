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
 * Alpha_complex is a Simplex_tree constructed from each finite cell of a Delaunay Triangulation.
 * 
 * The filtration value of each simplex is computed from the alpha square value of the simplex if it is Gabriel or
 * from the alpha value of the simplex coface that makes the simplex not Gabriel.
 * 
 * Please refer to \cite AlphaShapesDefinition for a more complete alpha complex definition.
 * 
 * Alpha complex are interesting because it looks like an \ref alpha-shape "Alpha shape" as described in
 * \cite AlphaShapesIntroduction (an alpha complex concept vulgarization).
 * 
 * \section example Example
 * 
 * This example loads points from an OFF file, builds the Delaunay triangulation from the points, and finally
 * initialize the alpha complex with it.
 * 
 * Then, it is asked to display information about the alpha complex.
 * 
 * \include Alpha_complex_from_off.cpp
 * 
 * When launching:
 * 
 * \code $> ./alphaoffreader ../../data/points/alphacomplexdoc.off 60.0
 * \endcode
 *
 * the program output is:
 * 
 * \include alphaoffreader_for_doc.txt
 * 
 * \section algorithm Algorithm
 * 
 * <b>Data structure</b>
 * 
 * In order to build the alpha complex, first, a Simplex tree is build from the cells of a Delaunay Triangulation.
 * (The filtration value is set to NaN, which stands for unknown value):
 * \image html "alpha_complex_doc.png" "Simplex tree structure construction example"
 *
 * <b>Filtration value computation algorithm</b>
 *
 * \f{algorithm}{ 
 * \caption{Filtration value computation algorithm}\label{alpha}
 * \begin{algorithmic}
 * \For{i : dimension $\rightarrow$ 1}
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
 * \end{algorithmic}
 * \f}
 * 
 * From the example above, it means the algorithm will look into each triangulation ([1,2,3], [2,3,4], [1,3,5], ...),
 * will compute the filtration value of the triangulation, and then will propagate the filtration value as described
 * here :
 * \image html "alpha_complex_doc_135.png" "Filtration value propagation example"
 * Then, the algorithm will look into each edge ([1,2], [2,3], [1,3], ...),
 * will compute the filtration value of the edge (in this case, propagation will have no effect).
 * 
 * Finally, the algorithm will look into each vertex ([1], [2], [3], [4], [5], [6] and [7]),
 * will set the filtration value (0 in case of a vertex - propagation will have no effect).
 * 
 * \section alpha-shape Alpha shape
 * 
 * In the example above, the alpha shape of \f$\alpha^2_{74} < \alpha^2 < \alpha^2_{73}\f$ is the alpha complex where the 
 * \f$\alpha^2_{74} <\f$ filtration value \f$< \alpha^2_{73}\f$ as described in \cite AlphaShapesIntroduction
 * 
 * \image html "alpha_complex_doc_alpha_shape.png" "Alpha shape example"
 * \copyright GNU General Public License v3.                         
 * \verbatim  Contact: gudhi-users@lists.gforge.inria.fr \endverbatim
 */
/** @} */  // end defgroup alpha_complex

} // namespace alphacomplex

} // namespace Gudhi

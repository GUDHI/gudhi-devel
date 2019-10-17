/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author:       François Godi
 *
 *    Copyright (C) 2017  INRIA
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

#ifndef DOC_TOPLEX_MAP_H_
#define DOC_TOPLEX_MAP_H_

// needs namespace for Doxygen to link on classes
namespace Gudhi {

/**  \defgroup toplex_map Toplex Map
 * 
 * \author    Fran&ccedil;ois Godi
 * @{
 * 
 * \section toplexmapdefinition Definition
 * 
 * Let's consider a simplicial complex, denote by \f$d\f$ its dimension
 * and by \f$k\f$ its number of maximal simplices.
 * Furthermore, denote by \f$\gamma_0\f$ the maximal number of toplices, i.e. maximal simplices,
 * that contain a same vertex.
 *
 * The goal of the Toplex Map is both to represent the complex in optimal
 * O(kd) space and to provide fast standard operations such as : insertion, removal
 * and membership of a simplex, contraction of an edge, collapses. The time needed
 * for these operation is linear or quadratic in \f$\gamma_0\f$ and \f$d\f$.
 *
 * Toplex map is composed firstly of a raw storage of toplices and secondly of a
 * map which associate any vertex to a set of pointers toward all toplices
 * containing this vertex.
 *
 * \image html map.png
 *
 * The performances are a lot better than in simplex tree as soon you use maximal simplices and not simplices,
 * here the construction of a strong witness complex of a point set with growing parameter :
 *
 * \image html graph.png
 *
 */
/** @} */  // end defgroup toplex_map

}  // namespace Gudhi

#endif  // DOC_TOPLEX_MAP_H_

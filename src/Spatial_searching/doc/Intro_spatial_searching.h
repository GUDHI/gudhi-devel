/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clement Jamin
 *
 *    Copyright (C) 2016 INRIA
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

#ifndef DOC_SPATIAL_SEARCHING_INTRO_SPATIAL_SEARCHING_H_
#define DOC_SPATIAL_SEARCHING_INTRO_SPATIAL_SEARCHING_H_

// needs namespaces for Doxygen to link on classes
namespace Gudhi {
namespace spatial_searching {

/**  \defgroup spatial_searching Spatial_searching
 * 
 * \author Cl&eacute;ment Jamin
 * 
 * @{
 * 
 * \section introduction Introduction
 * 
 * This Gudhi component is a wrapper around 
 * <a target="_blank" href="http://doc.cgal.org/latest/Spatial_searching/index.html">CGAL dD spatial searching algorithms</a>.
 * It provides a simplified API to perform (approximate) neighbor searches. Contrary to CGAL default behavior, the tree
 * does not store the points themselves, but stores indices.
 *
 * For more details about the data structure or the algorithms, or for more advanced usages, reading 
 * <a target="_blank" href="http://doc.cgal.org/latest/Spatial_searching/index.html">CGAL documentation</a>
 * is highly recommended.
 * 
 * \section spatial_searching_examples Example
 * 
 * This example generates 500 random points, then performs queries for nearest and farthest points using different methods.
 * 
 * \include Spatial_searching/example_spatial_searching.cpp
 * 
 * \copyright GNU General Public License v3.                         
 * \verbatim  Contact: gudhi-users@lists.gforge.inria.fr \endverbatim
 */
/** @} */  // end defgroup spatial_searching

}  // namespace spatial_searching

}  // namespace Gudhi

#endif  // DOC_SPATIAL_SEARCHING_INTRO_SPATIAL_SEARCHING_H_

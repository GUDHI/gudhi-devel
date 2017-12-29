/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Mathieu Carriere
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

#ifndef DOC_KERNEL_INTRO_KERNEL_H_
#define DOC_KERNEL_INTRO_KERNEL_H_

namespace Gudhi {

namespace kernel {

/**  \defgroup kernel Kernels
 * 
 * \author    Mathieu Carrière
 * 
 * @{
 * 
 * Kernels are generalized scalar products. They take the form of functions whose evaluations on pairs of persistence diagrams are equal
 * to the scalar products of the images of the diagrams under some feature map into a (generally unknown and infinite dimensional)
 * Hilbert space. Kernels are
 * very useful to handle any type of data for algorithms that require at least a Hilbert structure, such as Principal Component Analysis
 * or Support Vector Machines. In this package, we implement three kernels for persistence diagrams: the Persistence Scale Space kernel,
 * the Persistence Weighted Gaussian kernel and the Sliced Wasserstein kernel.
 *
 *
 * When launching:
 *
 * \code $> ./BasicEx
 * \endcode
 *
 * the program output is:
 *
 *
 * \copyright GNU General Public License v3.                         
 * \verbatim  Contact: gudhi-users@lists.gforge.inria.fr \endverbatim
 */
/** @} */  // end defgroup kernel

}  // namespace kernel

}  // namespace Gudhi

#endif  // DOC_KERNEL_INTRO_KERNEL_H_

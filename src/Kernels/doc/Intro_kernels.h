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
 * or Support Vector Machines. In this package, we implement three kernels for persistence diagrams:
 * the Persistence Scale Space Kernel (PSSK)---see \cite Reininghaus_Huber_ALL_PSSK,
 * the Persistence Weighted Gaussian Kernel (PWGK)---see \cite Kusano_Fukumizu_Hiraoka_PWGK,
 * and the Sliced Wasserstein Kernel (SWK)---see \cite pmlr-v70-carriere17a.
 *
 * \section  pwg Persistence Weighted Gaussian Kernel and Persistence Scale Space Kernel
 *
 * The PWGK is built with Gaussian Kernel Mean Embedding, meaning that each persistence diagram is first
 * sent to the Hilbert space of a Gaussian kernel with bandwidth parameter \f$\sigma >0\f$ using a weighted mean embedding \f$\Phi\f$:
 *
 * \f$ \Phi\,:\,D\,\rightarrow\,\sum_{p\in D}\,w(p)\,{\rm exp}\left(-\frac{\|p-\cdot\|_2^2}{2\sigma^2}\right) \f$,
 *
 * Usually, the weight function is chosen to be an arctan function of the distance of the point to the diagonal:
 * \f$w(p) = {\rm arctan}(C\,|y-x|^\alpha)\f$, for some parameters \f$C,\alpha >0\f$.
 * Then, either their scalar product in this space is
 * computed (Linear Persistence Weighted Gaussian Kernel):
 *
 * \f$ LPWGK(D_1,D_2)=\langle\Phi(D_1),\Phi(D_2)\rangle
 * \,=\,\sum_{p\in D_1}\,\sum_{q\in D_2}\,w(p)\,w(q)\,{\rm exp}\left(-\frac{\|p-q\|_2^2}{2\sigma^2}\right)\f$,
 *
 * or a second Gaussian kernel with bandwidth parameter \f$\tau >0\f$ is applied to their distance in this space
 * (Gaussian Persistence Weighted Gaussian Kernel):
 *
 * \f$ GPWGK(D_1,D_2)={\rm exp}\left(-\frac{\|\Phi(D_1)-\Phi(D_2)\|^2}{2\tau^2} \right)\f$,
 * where \f$\|\Phi(D_1)-\Phi(D_2)\|^2 = \langle\Phi(D_1)-\Phi(D_2),\Phi(D_1)-\Phi(D_2)\rangle\f$.
 *
 * It follows that the computation time is \f$O(n^2)\f$ where \f$n\f$ is the number of points
 * in the diagrams. This time can be improved by computing approximations of the kernel
 * with \f$m\f$ Fourier features \cite Rahimi07randomfeatures. In that case, the computation time becomes \f$O(mn)\f$.
 *
 * The PSSK is a Linear Persistence Weighted Gaussian Kernel between modified diagrams:
 * the symmetric of each point with respect to the diagonal is first added in each diagram, and then the weight function
 * is set to be +1 if the point is above the diagonal and -1 otherwise.
 *
 * \section sw Sliced Wasserstein Kernel
 *
 * The Sliced Wasserstein Kernel is defined as a Gaussian-like Kernel between persistence diagrams, where the distance used for
 * comparison is the Sliced Wasserstein distance \f$SW\f$ between persistence diagrams, defined as the integral of the 1-norm
 * between the sorted projections of the diagrams onto all lines passing through the origin:
 *
 * \f$ SW(D_1,D_2)=\int_{\theta\in\mathbb{S}}\,\|\pi_\theta(D_1\cup\pi_\Delta(D_2))-\pi_\theta(D_2\cup\pi_\Delta(D_1))\|_1{\rm d}\theta\f$,
 *
 * where \f$\pi_\theta\f$ is the projection onto the line defined with angle \f$\theta\f$ in the unit circle \f$\mathbb{S}\f$,
 * and \f$\pi_\Delta\f$ is the projection onto the diagonal.
 * The integral can be either computed exactly in \f$O(n^2{\rm log}(n))\f$ time, where \f$n\f$ is the number of points
 * in the diagrams, or approximated by sampling \f$m\f$ lines in the circle in \f$O(mn{\rm log}(n))\f$ time. The SWK is then computed as:
 *
 * \f$ SWK(D_1,D_2) = {\rm exp}\left(-\frac{SW(D_1,D_2)}{2\sigma^2}\right).\f$
 *
 * When launching:
 *
 * \code $> ./BasicEx ../../../../data/persistence_diagram/PD1 ../../../../data/persistence_diagram/PD2
 * \endcode
 *
 * the program output is:
 *
 * \include Kernels/kernel.txt
 *
 *
 * \copyright GNU General Public License v3.                         
 * \verbatim  Contact: gudhi-users@lists.gforge.inria.fr \endverbatim
 */
/** @} */  // end defgroup kernel

}  // namespace kernel

}  // namespace Gudhi

#endif  // DOC_KERNEL_INTRO_KERNEL_H_

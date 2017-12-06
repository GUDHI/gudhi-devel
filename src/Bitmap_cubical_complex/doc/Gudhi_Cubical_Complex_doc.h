/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2015  INRIA Sophia-Saclay (France)
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


#ifndef DOC_GUDHI_CUBICAL_COMPLEX_COMPLEX_H_
#define DOC_GUDHI_CUBICAL_COMPLEX_COMPLEX_H_

namespace Gudhi {

namespace cubical_complex {

/**  \defgroup cubical_complex Cubical complex
 *
 * \author   Pawel Dlotko
 *
 * @{
 *

 * Bitmap_cubical_complex is an example of a structured complex useful in computational mathematics (specially rigorous
 * numerics) and image analysis. The presented implementation of cubical complexes is based on the following
 * definition.
 *
 * An <em>elementary interval</em> is an interval of a form \f$ [n,n+1] \f$, or \f$[n,n]\f$, for \f$ n \in \mathcal{Z}
 * \f$. The first one is called <em>non-degenerate</em>, while the second one is \a degenerate interval. A
 * <em>boundary of a elementary interval</em> is a chain  \f$\partial [n,n+1] = [n+1,n+1]-[n,n] \f$ in case of
 * non-degenerated elementary interval and \f$\partial [n,n] = 0 \f$ in case of degenerate elementary interval. An
 * <em>elementary cube</em> \f$ C \f$ is a product of elementary intervals, \f$C=I_1 \times \ldots \times I_n\f$.
 * <em>Embedding dimension</em> of a cube is n, the number of elementary intervals (degenerate or not) in the product.
 * A <em>dimension of a cube</em> \f$C=I_1 \times ... \times I_n\f$ is the number of non degenerate elementary
 * intervals in the product. A <em>boundary of a cube</em> \f$C=I_1 \times \ldots \times I_n\f$ is a chain obtained
 * in the following way:
 * \f[\partial C = (\partial I_1 \times \ldots \times I_n) + (I_1 \times \partial I_2 \times \ldots \times I_n) +
 * \ldots + (I_1 \times I_2 \times \ldots \times \partial I_n).\f]
 * A <em>cubical complex</em> \f$\mathcal{K}\f$ is a collection of cubes closed under operation of taking boundary
 * (i.e. boundary of every cube from the collection is in the collection). A cube \f$C\f$ in cubical complex
 * \f$\mathcal{K}\f$ is <em>maximal</em> if it is not in a boundary of any other cube in \f$\mathcal{K}\f$. A \a
 * support of a cube \f$C\f$ is the set in \f$\mathbb{R}^n\f$ occupied by \f$C\f$ (\f$n\f$ is the embedding dimension
 * of \f$C\f$).
 *
 * Cubes may be equipped with a filtration values in which case we have filtered cubical complex. All the cubical
 * complexes considered in this implementation are filtered cubical complexes (although, the range of a filtration may
 * be a set of two elements).
 *
 * For further details and theory of cubical complexes, please consult \cite kaczynski2004computational as well as the
 * following paper \cite peikert2012topological .
 *
 * \section cubicalcomplexdatastructure Data structure
 *
 * The implementation of Cubical complex provides a representation of complexes that occupy a rectangular region in
 * \f$\mathbb{R}^n\f$. This extra assumption allows for a memory efficient way of storing cubical complexes in a form
 * of so called bitmaps. Let \f$R = [b_1,e_1] \times \ldots \times [b_n,e_n]\f$, for \f$b_1,...b_n,e_1,...,e_n \in
 * \mathbb{Z}\f$, \f$b_i \leq d_i\f$ be the considered rectangular region and let \f$\mathcal{K}\f$ be a filtered
 * cubical complex having the rectangle \f$R\f$ as its support. Note that the structure of the coordinate system gives
 * a way a lexicographical ordering of cells of \f$\mathcal{K}\f$. This ordering is a base of the presented
 * bitmap-based implementation. In this implementation, the whole cubical complex is stored as a vector of the values
 * of filtration. This, together with dimension of \f$\mathcal{K}\f$ and the sizes of \f$\mathcal{K}\f$ in all
 * directions, allows to determine, dimension, neighborhood, boundary and coboundary of every cube \f$C \in
 * \mathcal{K}\f$.
 *
 * \image html "Cubical_complex_representation.png" Cubical complex.
 *
 * Note that the cubical complex in the figure above is, in a natural way, a product of one dimensional cubical
 * complexes in \f$\mathbb{R}\f$. The number of all cubes in each direction is equal \f$2n+1\f$, where \f$n\f$ is the
 * number of maximal cubes in the considered direction. Let us consider a cube at the position \f$k\f$ in the bitmap.
 * Knowing the sizes of the bitmap, by a series of modulo operation, we can determine which elementary intervals are
 * present in the product that gives the cube \f$C\f$. In a similar way, we can compute boundary and the coboundary of
 * each cube. Further details can be found in the literature.
 *
 * \section inputformat Input Format
 *
 * In the current implantation, filtration is given at the maximal cubes, and it is then extended by the lower star
 * filtration to all cubes. There are a number of constructors that can be used to construct cubical complex by users
 * who want to use the code directly. They can be found in the \a Bitmap_cubical_complex class.
 * Currently one input from a text file is used. It uses a format used already in Perseus software
 * (http://www.sas.upenn.edu/~vnanda/perseus/) by Vidit Nanda. The file format is described here: \ref FileFormatsPerseus.
 *
 * \section PeriodicBoundaryConditions Periodic boundary conditions
 * Often one would like to impose periodic boundary conditions to the cubical complex. Let \f$ I_1\times ... \times
 * I_n \f$ be a box that is decomposed with a cubical complex \f$ \mathcal{K} \f$. Imposing periodic boundary
 * conditions in the direction i, means that the left and the right side of a complex \f$ \mathcal{K} \f$ are
 * considered the same. In particular, if for a bitmap \f$ \mathcal{K} \f$ periodic boundary conditions are imposed
 * in all directions, then complex \f$ \mathcal{K} \f$ became n-dimensional torus. One can use various constructors
 * from the file Bitmap_cubical_complex_periodic_boundary_conditions_base.h to construct cubical complex with periodic
 * boundary conditions. One can also use Perseus style input files (see \ref FileFormatsPerseus).
 *
 * \section BitmapExamples Examples
 * End user programs are available in example/Bitmap_cubical_complex and utilities/Bitmap_cubical_complex folders.
 * 
 * \copyright GNU General Public License v3.
 */
/** @} */  // end defgroup cubical_complex

}  // namespace cubical_complex

namespace Cubical_complex = cubical_complex;

}  // namespace Gudhi

#endif  // DOC_GUDHI_CUBICAL_COMPLEX_COMPLEX_H_

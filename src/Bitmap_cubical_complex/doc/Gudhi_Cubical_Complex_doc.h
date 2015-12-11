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


#pragma once

namespace Gudhi
{

namespace Cubical_complex
{

/**  \defgroup cubical_complex Cubical complex
*
* \author   Pawel Dlotko
*
* @{
*

*Cubical complex is an example of a structured complex useful in computational mathematics (specially rigorous numerics) and image analysis. The presented implementation of cubical complexes is based on the following definition.
*
* An <em>elementary interval</em> is an interval of a form \f$ [n,n+1] \f$, or \f$[n,n]\f$, for \f$ n \in \mathcal{Z} \f$. The first one is called <em>non-degenerated</em>, while the second one is \a degenerated interval. A <em>boundary of a elementary
*interval</em> is a chain  \f$\partial [n,n+1] = [n+1,n+1]-[n,n] \f$ in case of non-degenerated elementary interval and \f$\partial [n,n] = 0 \f$ in case of degenerated elementary interval. An <em>elementary cube</em> \f$ C \f$ is a

*product of elementary intervals, \f$C=I_1 \times \ldots \times I_n\f$. <em>Embedding dimension</em> of a cube is n, the number of elementary intervals (degenerated or not) in the product. A <em>dimension of a cube</em> \f$C=I_1 \times ... \times I_n\f$ is the
*number of non degenerated elementary intervals in the product. A <em>boundary of a cube</em> \f$C=I_1 \times \ldots \times I_n\f$ is a chain obtained in the following way:
*\f[\partial C = (\partial I_1 \times \ldots \times I_n) + (I_1 \times \partial I_2 \times \ldots \times I_n) + \ldots + (I_1 \times I_2 \times \ldots \times \partial I_n).\f]
*A <em>cubical complex</em> \f$\mathcal{K}\f$ is a collection of cubes closed under operation of taking boundary (i.e. boundary of every cube from the collection is in the collection). A cube \f$C\f$ in cubical complex \f$\mathcal{K}\f$ is <em>maximal</em> if it is not in
*a boundary of any other cube in \f$\mathcal{K}\f$. A \a support of a cube \f$C\f$ is the set in \f$\mathbb{R}^n\f$ occupied by \f$C\f$ (\f$n\f$ is the embedding dimension of \f$C\f$).
*
*Cubes may be equipped with a filtration values in which case we have filtered cubical complex. All the cubical complexes considered in this implementation are filtered cubical complexes (although, the range of a filtration may be a set of two elements).
*
*For further details and theory of cubical complexes, please consult \cite kaczynski2004computational .
*
*as well as the following paper \cite peikert2012topological .
*
*\section datastructure Data structure.
*
*The implementation of Cubical complex provides a representation of complexes that occupy a rectangular region in \f$\mathbb{R}^n\f$. This extra
*assumption allows for a memory efficient way of storing cubical complexes in a form of so called bitmaps. Let \f$R = [b_1,e_1] \times \ldots \times [b_n,e_n]\f$, for \f$b_1,...b_n,e_1,...,e_n \in \mathbb{Z}\f$
*, \f$b_i \leq d_i\f$ be the considered rectangular region and let \f$\mathcal{K}\f$ be a filtered cubical complex having the rectangle \f$R\f$ as its support. Note that the structure of the coordinate system gives a way
*a lexicographical ordering of cells of \f$\mathcal{K}\f$. This ordering is a base of the presented bitmap-based implementation. In this implementation, the whole cubical complex is stored as a vector
*of the values of filtration. This, together with dimension of \f$\mathcal{K}\f$ and the sizes of \f$\mathcal{K}\f$ in all directions, allows to determine, dimension, neighborhood, boundary and coboundary of every cube \f$C \in \mathcal{K}\f$.
*
*\image html "bitmapAllCubes.png" "Cubical complex.
*
*Note that the cubical complex in the figure above is, in a natural way, a product of one dimensional cubical complexes in \f$\mathbb{R}\f$. The number of all cubes in each direction is
*equal \f$2n+1\f$, where \f$n\f$ is the number of maximal cubes in the considered direction. Let us consider a cube at the position \f$k\f$ in the bitmap. Knowing the sizes of the bitmap,
*by a series of modulo operation, we can determine which elementary intervals are present in the product that gives the cube \f$C\f$. In a similar way, we can compute boundary
*and the coboundary of each cube. Further details can be found in the literature.
*
*\section inputformat Input Format.
*
*In the current implantation, filtration is given at the maximal cubes, and it is then extended by the lower star filtration to all cubes. There are a number of constructors
*that can be used to construct cubical complex by users who want to use the code directly. They can be found in the \a Bitmap_cubical_complex class.
*Currently one input from a text file is used. It uses a format used already in Perseus software (http://www.sas.upenn.edu/~vnanda/perseus/) by Vidit Nanda.
*Below we are providing a description of the format.
*
*
*\image html "exampleBitmap.png" "Example of a input data."
*
*The input file for the following complex is:
*\verbatim
2
3
3
1
2
3
8
20
4
7
6
5
\endverbatim

*/
/** @} */  // end defgroup cubical_complex

*@}//end of the group
}
}

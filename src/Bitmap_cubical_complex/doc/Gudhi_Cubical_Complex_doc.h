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
* An \emph{elementary interval} is an interval of a form $[n,n+1]$, or $[n,n]$, for $n \in \mathcal{Z}$. The first one is called \emph{non-degenerated}, while the second one is \emph{degenerated} interval. A \emph{boundary of a elementary
*interval} is a chain  $\partial [n,n+1] = [n+1,n+1]-[n,n]$ in case of non-degenerated elementary interval and $\partial [n,n] = 0$ in case of degenerated elementary interval. An \emph{elementary cube} $C$ is a
*product of elementary intervals, $C=I_1 \times \ldots \times I_n$. \emph{Embedding dimension} of a cube is n, the number of elementary intervals (degenerated or not) in the product. A \emph{dimension of a cube} $C=I_1 \times ... \times I_n$ is the
*number of non degenerated elementary intervals in the product. A \emph{boundary of a cube} $C=I_1 \times \ldots \times I_n$ is a chain obtained in the following way:
*\[\partial C = (\partial I_1 \times \ldots \times I_n) + (I_1 \times \partial I_2 \times \ldots \times I_n) + \ldots + (I_1 \times I_2 \times \ldots \times \partial I_n).\]
*A \emph{cubical complex} $\mathcal{K}$ is a collection of cubes closed under operation of taking boundary (i.e. boundary of every cube from the collection is in the collection). A cube $C$ in cubical complex $\mathcal{K}$ is \emph{maximal} if it is not in
*a boundary of any other cube in $\mathcal{K}$. A \emph{support} of a cube $C$ is the set in $\mathbb{R}^n$ occupied by $C$ ($n$ is the embedding dimension of $C$).
*
*Cubes may be equipped with a filtration values in which case we have filtered cubical complex. All the cubical complexes considered in this implementation are filtered cubical complexes (although, the range of a filtration may be a set of two elements).
*
*For further details and theory of cubical complexes, please consult a book:\\
*Computational homology, by Tomasz Kaczynski, Konstantin Mischaikow, and Marion Mrozek, Appl. Math. Sci., vol. 157, Springer-Verlag, New York, 2004
*
*as well as the paper:
*Efficient computation of persistent homology for cubical data by Hubert Wagner, Chao Chen, Erald Vuçini (published in the proceedings of Workshop on Topology-based Methods in Data
*Analysis and Visualization)
*
*\section{Data structure.}
*
*The implementation of Cubical complex provides a representation of complexes that occupy a rectangular region in $\mathbb{R}^n$. This extra
*assumption allows for a memory efficient way of storing cubical complexes in a form of so called bitmaps. Let $R = [b_1,e_1] \times \ldots \times [b_n,e_n]$, for $b_1,...b_n,e_1,...,e_n \in \mathbb{Z}$
*, $b_i \leq d_i$ be the considered rectangular region and let $\mathcal{K}$ be a filtered cubical complex having the rectangle $R$ as its support. Note that the structure of the coordinate system gives a way
*a lexicographical ordering of cells of $\mathcal{K}$. This ordering is a base of the presented bitmap-based implementation. In this implementation, the whole cubical complex is stored as a vector
*of the values of filtration. This, together with dimension of $\mathcal{K}$ and the sizes of $\mathcal{K}$ in all directions, allows to determine, dimension, neighborhood, boundary and coboundary of every cube $C \in \mathcal{K}$.
*
*\image html "bitmapAllCubes.pdf" "Cubical complex in $\mathbb{R}^2".
*
*Note that the cubical complex in the figure above is, in a natural way, a product of one dimensional cubical complexes in $\mathbb{R}$. The number of all cubes in each direction is
*equal $2n+1$, where $n$ is the number of maximal cubes in the considered direction. Let us consider a cube at the position $k$ in the bitmap. Knowing the sizes of the bitmap,
*by a series of modulo operation, we can determine which elementary intervals are present in the product that gives the cube $C$. In a similar way, we can compute boundary
*and the coboundary of each cube. Further details can be found in the literature.
*
*\section{Input Format.}
*
*In the current implantation, filtration is given at the maximal cubes, and it is then extended by the lower star filtration to all cubes. There are a number of constructors
*that can be used to construct cubical complex by users who want to use the code directly. They can be found in the \emph{Bitmap\_cubical\_complex} class.
*Currently one input from a text file is used. It uses a format used already in Perseus software $(http://www.sas.upenn.edu/~vnanda/perseus/)$ by Vidit Nanda.
*Below we are providing a description of the format.
*
*\begin{enumerate}
*\item The first line of the file is $d$, the embedding dimension of a complex.
*\item The next $d$ lines consist of positive numbers being the numbers of top dimensional cubes in the given direction. Let us call those numbers $n_1,\ldots,n_d$.
*\item Later there is a sequence of $n_1 \dot \ldots \dot n_d$ numbers in a lexicographical ordering. Those numbers are filtrations of top dimensional cubes.
*\end{enumerate}
*
*\image html "exampleBitmap.pdf" "Example of a input data."
*
*The input file for the following complex is:
*\begin{verbatim}
*2
*3
*3
*1
*2
*3
*8
*20
*4
*7
*6
*5
*\end{verbatim}
*
*
}
}

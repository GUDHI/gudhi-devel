/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014  INRIA Sophia Antipolis-Méditerranée (France)
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

#ifndef DOC_SIMPLEX_TREE_INTRO_SIMPLEX_TREE_H_
#define DOC_SIMPLEX_TREE_INTRO_SIMPLEX_TREE_H_

// needs namespace for Doxygen to link on classes
namespace Gudhi {


/** \defgroup simplex_tree Filtered Complexes
 * \author    Cl&eacute;ment Maria
 *
 * A simplicial complex \f$\mathbf{K}\f$ on a set of vertices \f$V = \{1, \cdots ,|V|\}\f$ is a collection of
 * implices \f$\{\sigma\}\f$, \f$\sigma \subseteq V\f$ such that
 * \f$\tau \subseteq \sigma \in \mathbf{K} \rightarrow \tau \in \mathbf{K}\f$. The dimension \f$n=|\sigma|-1\f$ of
 * \f$\sigma\f$ is its number of elements minus \f$1\f$.
 *
 * A filtration of a simplicial complex is a function \f$f:\mathbf{K} \rightarrow \mathbb{R}\f$ satisfying
 * \f$f(\tau)\leq f(\sigma)\f$ whenever \f$\tau \subseteq \sigma\f$. Ordering the simplices by increasing filtration
 * values (breaking ties so as a simplex appears after its subsimplices of same filtration value) provides an
 * indexing scheme.
 *
 * \section filteredcomplexesimplementation Implementations
 * \subsection filteredcomplexessimplextree Simplex tree
 * There are two implementation of complexes. The first on is the Simplex_tree data structure. The simplex tree is an
 * efficient and flexible data structure for representing general (filtered) simplicial complexes. The data structure
 * is described in \cite boissonnatmariasimplextreealgorithmica
 * \image html "Simplex_tree_representation.png" "Simplex tree representation"
 * 
 * \subsubsection filteredcomplexessimplextreeexamples Examples
 * 
 * Here is a list of simplex tree examples :
 * \li <a href="_simplex_tree_2simple_simplex_tree_8cpp-example.html">
 * Simplex_tree/simple_simplex_tree.cpp</a> - Simple simplex tree construction and basic function use.
 *
 * \li <a href="_simplex_tree_2simplex_tree_from_cliques_of_graph_8cpp-example.html">
 * Simplex_tree/simplex_tree_from_cliques_of_graph.cpp</a> - Simplex tree construction from cliques of graph read in
 * a file.
 * 
 * Simplex tree construction with \f$\mathbb{Z}/3\mathbb{Z}\f$ coefficients on weighted graph Klein bottle file:
 * \code $> ./simplex_tree_from_cliques_of_graph ../../data/points/Klein_bottle_complex.txt 3 \endcode
 * \code Insert the 1-skeleton in the simplex tree in 0.000404 s. 
max_dim = 3
Expand the simplex tree in 3.8e-05 s. 
Information of the Simplex Tree: 
Number of vertices = 10   Number of simplices = 98 \endcode
 * 
 * \li <a href="_simplex_tree_2simplex_tree_from_alpha_shapes_3_8cpp-example.html">
 * Simplex_tree/simplex_tree_from_alpha_shapes_3.cpp</a> - Simplex tree is computed and displayed from a 3D alpha
 * complex (Requires CGAL, GMP and GMPXX to be installed)
 * 
 * 
 * \subsection filteredcomplexeshassecomplex Hasse complex
 * The second one is the Hasse_complex. The Hasse complex is a data structure representing explicitly all co-dimension
 * 1 incidence relations in a complex. It is consequently faster when accessing the boundary of a simplex, but is less
 * compact and harder to construct from scratch.
 * 
 * \copyright GNU General Public License v3.
 * @{
 */

}  // namespace Gudhi

#endif  // DOC_SIMPLEX_TREE_INTRO_SIMPLEX_TREE_H_

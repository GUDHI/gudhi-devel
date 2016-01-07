 /*    This file is part of the Gudhi Library. The Gudhi library
  *    (Geometric Understanding in Higher Dimensions) is a generic C++
  *    library for computational topology.
  *
  *    Author(s):       David Salinas
  *
  *    Copyright (C) 2014  INRIA Sophia Antipolis-Mediterranee (France)
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

#ifndef SKELETON_BLOCKER_H_
#define SKELETON_BLOCKER_H_

#include <gudhi/Skeleton_blocker_complex.h>
#include <gudhi/Skeleton_blocker_geometric_complex.h>
#include <gudhi/Skeleton_blocker_simplifiable_complex.h>
#include <gudhi/Skeleton_blocker/Skeleton_blocker_off_io.h>

#include <gudhi/Skeleton_blocker/Skeleton_blocker_simple_traits.h>
#include <gudhi/Skeleton_blocker/Skeleton_blocker_simple_geometric_traits.h>

#include <gudhi/Utils.h>  // xxx

namespace Gudhi {

namespace skbl {

/** \defgroup skbl Skeleton-Blocker 

\author David Salinas

\section Introduction
The Skeleton-Blocker data-structure proposes a light encoding for simplicial complexes by storing only an *implicit* representation of its
simplices 
\cite socg_blockers_2011,\cite blockers2012.
Intuitively, it just stores the 1-skeleton of a simplicial complex with a graph and the set of its "missing faces" that
is very small in practice (see next section for a formal definition).
This data-structure handles all  simplicial complexes operations such as
 as simplex enumeration or simplex removal but operations that are particularly efficient 
 are operations that do not require simplex enumeration such as edge iteration, link computation or simplex contraction.


\section Definitions

We recall briefly classical definitions of simplicial complexes  
 \cite Munkres-elementsalgtop1984.
An abstract simplex is a finite non-empty set and its dimension is its number of elements minus 1.
Whenever \f$\tau \subset \sigma\f$ and \f$\tau \neq \emptyset \f$, \f$ \tau \f$ is called a face of 
\f$ \sigma\f$  and \f$ \sigma\f$ is called a coface of \f$ \tau \f$ . Furthermore,
when \f$ \tau \neq \sigma\f$ we say that \f$ \tau\f$ is a proper-face of \f$ \sigma\f$.
An abstract simplicial complex is a set of simplices that contains all the faces of its simplices.
The 1-skeleton of a simplicial complex (or its graph) consists of its elements of dimension lower than 2.

*\image html "ds_representation.png" "Skeleton-blocker representation" width=20cm


To encode, a simplicial complex, one can encodes all its simplices. 
In case when this number gets too large,
a lighter and implicit version consists of encoding only its graph plus some elements called missing faces or blockers.
A blocker is a simplex of dimension greater than 1 
that does not belong to the complex but whose all proper faces does.


Remark that for a clique complex (i.e. a simplicial complex whose simplices are cliques of its graph), the set of blockers
is empty and the data-structure is then particularly sparse. 
One famous example of clique-complex is the Rips complex which is intensively used
in topological data-analysis.
In practice, the set of blockers of a simplicial complex 
remains also small when simplifying a Rips complex with edge contractions 
but also for most of the simplicial complexes used in topological data-analysis such as Delaunay, Cech or Witness complexes. 
For instance, the numbers of blockers is depicted for random 3-dimensional spheres embedded into \f$R^4\f$ 
in next figure. Storing the graph and blockers of such simplicial complexes is much compact in this case than storing 
their simplices.


*\image html "blockers_curve.png" "Number of blockers of random triangulations of 3-spheres" width=10cm




\section API

\subsection Overview

Two main classes of this package are Skeleton_blocker_complex and Skeleton_blocker_geometric_complex.
The first one can be used to represent an abstract simplicial complex and supports most used
operations in a simplicial complex such as :

\li vertex/edge/simplex enumeration
\li simplifications operations such as remove star, add star (e.g. general form of collapse),
edge contractions

The class Skeleton_blocker_geometric_complex supports the same methods as Skeleton_blocker_complex
and point access in addition.



\subsection Visitor

The class Skeleton_blocker_complex has a visitor that is called when usual operations such as adding an edge or remove a vertex are called.
You may want to use this visitor to compute statistics or to update another data-structure (for instance this visitor is heavily used in the \ref contr package).




\section Example

 
\subsection Iterating Iterating through vertices, edges, blockers and simplices	
 
Iteration through vertices, edges, simplices or blockers is straightforward with c++11 for range loops.
Note that simplex iteration with this implicit data-structure just takes
a few more time compared to iteration via an explicit representation 
such as the Simplex Tree. The following example computes the Euler Characteristic
of a simplicial complex.

  \code{.cpp}
	typedef Skeleton_blocker_complex<Skeleton_blocker_simple_traits> Complex;
	typedef Complex::Vertex_handle Vertex_handle;
	typedef Complex::Simplex Simplex;

  	const int n = 15;

	// build a full complex with 10 vertices and 2^n-1 simplices
	Complex complex;
	for(int i=0;i<n;i++)
		complex.add_vertex();
	for(int i=0;i<n;i++)
		for(int j=0;j<i;j++)
			complex.add_edge_without_blockers(Vertex_handle(i),Vertex_handle(j));

	// this is just to illustrate iterators, to count number of vertices
	// or edges, complex.num_vertices() and complex.num_edges() are
	// more appropriated!
	unsigned num_vertices = 0;
	for(auto v : complex.vertex_range()){
		++num_vertices;
	}

	unsigned num_edges = 0;
	for(auto e : complex.edge_range())
		++num_edges;

	unsigned euler = 0;
	unsigned num_simplices = 0;
	// we use a reference to a simplex instead of a copy
	// value here because a simplex is a set of integers
	// and copying it cost time
	for(const Simplex & s : complex.star_simplex_range()){
		++num_simplices;
		if(s.dimension()%2 == 0) 
			euler += 1;
		else 
			euler -= 1;
	}
	std::cout << "Saw "<<num_vertices<<" vertices, "<<num_edges<<" edges and "<<num_simplices<<" simplices"<<std::endl;
	std::cout << "The Euler Characteristic is "<<euler<<std::endl;
  \endcode


\verbatim
./SkeletonBlockerIteration
Saw 15 vertices, 105 edges and 32767 simplices
The Euler Characteristic is 1
 0.537302s wall, 0.530000s user + 0.000000s system = 0.530000s CPU (98.6%)
\endverbatim


\subsection s Constructing a skeleton-blockers from a list of maximal faces or from a list of faces

  \code{.cpp}
	std::vector<Simplex> simplices;

	//add 4 triangles of a tetrahedron 0123
	simplices.push_back(Simplex(Vertex_handle(0),Vertex_handle(1),Vertex_handle(2)));
	simplices.push_back(Simplex(Vertex_handle(1),Vertex_handle(2),Vertex_handle(3)));
	simplices.push_back(Simplex(Vertex_handle(3),Vertex_handle(0),Vertex_handle(2)));
	simplices.push_back(Simplex(Vertex_handle(3),Vertex_handle(0),Vertex_handle(1)));

	Complex complex;
	//get complex from top faces
	make_complex_from_top_faces(complex,simplices.begin(),simplices.end());

	std::cout << "Simplices:"<<std::endl;
	for(const Simplex & s : complex.star_simplex_range())
		std::cout << s << " ";
	std::cout << std::endl;

	//One blocker as simplex 0123 is not in the complex but all its proper faces are.
	std::cout << "Blockers: "<<complex.blockers_to_string()<<std::endl;

	//now build a complex from its full list of simplices
	simplices.clear();
	simplices.push_back(Simplex(Vertex_handle(0)));
	simplices.push_back(Simplex(Vertex_handle(1)));
	simplices.push_back(Simplex(Vertex_handle(2)));
	simplices.push_back(Simplex(Vertex_handle(0),Vertex_handle(1)));
	simplices.push_back(Simplex(Vertex_handle(1),Vertex_handle(2)));
	simplices.push_back(Simplex(Vertex_handle(2),Vertex_handle(0)));
	complex = Complex(simplices.begin(),simplices.end());

	std::cout << "Simplices:"<<std::endl;
	for(const Simplex & s : complex.star_simplex_range())
		std::cout << s << " ";
	std::cout << std::endl;

	//One blocker as simplex 012 is not in the complex but all its proper faces are.
	std::cout << "Blockers: "<<complex.blockers_to_string()<<std::endl;
 \endcode
\verbatim
./SkeletonBlockerFromSimplices
Simplices:
{0} {0,1} {0,2} {0,3} {0,1,2} {0,1,3} {0,2,3} {1} {1,2} {1,3} {1,2,3} {2} {2,3} {3} 
Blockers: {0,1,2,3}

Simplices:
{0} {0,1} {0,2} {1} {1,2} {2} 
Blockers: {0,1,2}
\endverbatim


\section Acknowledgements
The author wishes to thank Dominique Attali and André Lieutier for 
their collaboration to write the two initial papers 
\cite socg_blockers_2011,\cite blockers2012
 about this data-structure
 and also Dominique for leaving him use a prototype. 


\copyright GNU General Public License v3.                         
\verbatim  Contact: gudhi-users@lists.gforge.inria.fr \endverbatim
*/
/** @} */  // end defgroup

}  // namespace skbl

}  // namespace Gudhi

#endif  // SKELETON_BLOCKER_H_

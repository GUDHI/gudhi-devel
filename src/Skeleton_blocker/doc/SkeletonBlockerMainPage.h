
/*! \mainpage Skeleton blockers data-structure

\author David Salinas

\section Introduction
The Skeleton-blocker data-structure had been introduced in the two papers 
\cite skbl_socg2011 \cite skbl_ijcga2012. 
It proposes a light encoding for simplicial complexes by storing only an *implicit* representation of its
simplices.
Intuitively, it just stores the 1-skeleton of a simplicial complex with a graph and the set of its "missing faces" that
is very small in practice (see next section for a formal definition).
This data-structure handles every classical operations used for simplicial complexes such as
 as simplex enumeration or simplex removal but operations that are particularly efficient 
 are operations that do not require simplex enumeration such as edge iteration, link computation or simplex contraction.


\todo{image wont include}

\section Definitions

We recall briefly classical definitions of simplicial complexes  \cite Munkres.
An abstract simplex is a finite non-empty set and its dimension is its number of elements minus 1.
Whenever \f$\tau \subset \sigma\f$ and \f$\tau \neq \emptyset \f$, \f$ \tau \f$ is called a face of 
\f$ \sigma\f$  and \f$ \sigma\f$ is called a coface of \f$ \tau \f$ . Furthermore,
when \f$ \tau \neq \sigma\f$ we say that \f$ \tau\f$ is a proper-face of \f$ \sigma\f$.
An abstract simplicial complex is a set of simplices that contains all the faces of its simplices.
The 1-skeleton of a simplicial complex (or its graph) consists of its elements of dimension lower than 2.


\image latex "images/ds_representation.eps" "My application" width=10cm


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
For instance, the numbers of blockers is depicted for a random 3 dimensional sphere embedded into \f$R^4\f$ 
in figure X.


\image latex "images/blockers_curve.eps" "My application" width=10cm




\section API

\subsection Overview

Four classes define  a simplicial complex namely :

\li <Code>Skeleton_blocker_complex</Code> : a simplicial complex with basic operations such as vertex/edge/simplex enumeration and construction
\li <Code>Skeleton_blocker_link_complex</Code> : the link of a simplex in a parent complex. It is represented as a sub complex
of the parent complex
\li <Code>Skeleton_blocker_simplifiable_complex</Code> : a simplicial complex with simplification operations such as edge contraction or simplex collapse
\li <Code>Skeleton_blocker_geometric_complex</Code> : a simplicial complex who has access to geometric points in  \f$R^d\f$ 

The two last classes are derived classes from the <Code>Skeleton_blocker_complex</Code> class. The class <Code>Skeleton_blocker_link_complex</Code> inheritates from a template passed parameter
that may be either <Code>Skeleton_blocker_complex</Code> or <Code>Skeleton_blocker_geometric_complex</Code> (a link may store points coordinates or not).

\todo{include links}

\subsection Visitor

The class <Code>Skeleton_blocker_complex</Code> has a visitor that is called when usual operations such as adding an edge or remove a vertex are called.
You may want to use this visitor to compute statistics or to update another data-structure (for instance this visitor is heavily used in the 
<Code>Contraction</Code> package).



\section Example

 
\subsection s Iterating through vertices, edges, blockers and simplices	
 
Iteration through vertices, edges, simplices or blockers is straightforward with c++11 for range loops.
Note that simplex iteration with this implicit data-structure just takes
a few more time compared to iteration via an explicit representation 
such as the Simplex Tree. The following example computes the Euler Characteristic
of a simplicial complex.

  \code{.cpp}
	typedef Skeleton_blocker_complex<Skeleton_blocker_simple_traits> Complex;
	typedef Complex::Vertex_handle Vertex_handle;
	typedef Complex::Simplex_handle Simplex;

  	const int n = 15;

	// build a full complex with 10 vertices and 2^n-1 simplices
	Complex complex;
	for(int i=0;i<n;i++)
		complex.add_vertex();
	for(int i=0;i<n;i++)
		for(int j=0;j<i;j++)
			//note that add_edge adds the edge and all its cofaces
			complex.add_edge(Vertex_handle(i),Vertex_handle(j));

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
	for(const Simplex & s : complex.simplex_range()){
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



\subsection Acknowledgements
The author wishes to thank Dominique Attali and André Lieutier for 
their collaboration to write the two initial papers about this data-structure
 and also Dominique for leaving him use a prototype. 


\copyright GNU General Public License v3.                         
\verbatim  Contact: David Salinas,     david.salinas@inria.fr \endverbatim

*/

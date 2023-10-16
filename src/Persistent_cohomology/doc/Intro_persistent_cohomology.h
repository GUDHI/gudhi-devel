/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef DOC_PERSISTENT_COHOMOLOGY_INTRO_PERSISTENT_COHOMOLOGY_H_
#define DOC_PERSISTENT_COHOMOLOGY_INTRO_PERSISTENT_COHOMOLOGY_H_

// needs namespace for Doxygen to link on classes
namespace Gudhi {
// needs namespace for Doxygen to link on classes
namespace persistent_cohomology {

/** \defgroup persistent_cohomology Persistent Cohomology

  \author    Clément Maria

  Computation of persistent cohomology using the algorithm of 
  \cite DBLP:journals/dcg/SilvaMV11 and \cite DBLP:conf/compgeom/DeyFW14
  and the Compressed Annotation Matrix 
  implementation of \cite DBLP:conf/esa/BoissonnatDM13 . Some \ref specialcasepersistence "special cases"
  using different algorithms are also provided at the end of this page.
       
  The theory of homology consists in attaching to a topological space a sequence of 
  (homology) groups, 
  capturing global topological features 
  like connected components, holes, cavities, etc. Persistent homology studies the evolution 
  -- birth, life and death -- of 
  these features when the topological space is changing. Consequently, the theory is essentially 
  composed of three elements: 
  topological spaces, their homology groups and an evolution scheme.

 \section persistencetopolocalspaces Topological Spaces
 Topological spaces are represented by simplicial complexes.
 Let \f$V = \{1, \cdots ,|V|\}\f$ be a set of <EM>vertices</EM>.
 A <EM>simplex</EM> \f$\sigma\f$ is a subset of vertices
 \f$\sigma \subseteq V\f$. A <EM>simplicial complex</EM> \f$\mathbf{K}\f$
 on \f$V\f$ is a collection of simplices \f$\{\sigma\}\f$,
 \f$\sigma \subseteq V\f$, such that \f$\tau \subseteq \sigma \in \mathbf{K}
 \Rightarrow \tau \in \mathbf{K}\f$. The dimension \f$n=|\sigma|-1\f$ of \f$\sigma\f$
 is its number of elements minus 1. A <EM>filtration</EM> of a simplicial complex is
 a function \f$f:\mathbf{K} \rightarrow \mathbb{R}\f$ satisfying \f$f(\tau)\leq
 f(\sigma)\f$ whenever \f$\tau \subseteq \sigma\f$.

 We define the concept FilteredComplex which enumerates the requirements for a class
 to represent a filtered complex from which persistent homology may be computed.
 We use the vocabulary of simplicial complexes, but the concept
 is valid for any type of cell complex. The main requirements
 are the definition of:
 \li type <CODE>Indexing_tag</CODE>, which is a model of the concept
 <CODE>IndexingTag</CODE>,
 describing the nature of the indexing scheme,
 \li type Simplex_handle to manipulate simplices,
 \li method <CODE>int dimension(Simplex_handle)</CODE> returning
 the dimension of a simplex,
 \li type and method <CODE>Boundary_simplex_range
 boundary_simplex_range(Simplex_handle)</CODE> that returns
 a range giving access to the codimension 1 subsimplices of the
 input simplex, as-well-as the coefficients \f$(-1)^i\f$ in the
 definition of the operator \f$\partial\f$. The iterators have
 value type <CODE>Simplex_handle</CODE>,
 \li type and method
 <CODE>Filtration_simplex_range filtration_simplex_range ()</CODE>
 that returns a range giving
 access to all the simplices of the complex read in the order
 assigned by the indexing scheme,
 \li type and method
 <CODE>Filtration_value filtration (Simplex_handle)</CODE> that returns the value of
 the filtration on the simplex represented by the handle.

 \section persistencehomology Homology
 For a ring \f$\mathcal{R}\f$, the group of <EM>n-chains</EM>,
 denoted \f$\mathbf{C}_n(\mathbf{K},\mathcal{R})\f$, of \f$\mathbf{K}\f$ is the
 group of formal sums of
 n-simplices with \f$\mathcal{R}\f$ coefficients. The <EM>boundary operator</EM> is a
 linear operator
 \f$\partial_n: \mathbf{C}_n(\mathbf{K},\mathcal{R}) \rightarrow \mathbf{C}_{n-1}(\mathbf{K},\mathcal{R})\f$
 such that \f$\partial_n \sigma = \partial_n [v_0, \cdots , v_n] =
 \sum_{i=0}^n (-1)^{i}[v_0,\cdots ,\widehat{v_i}, \cdots,v_n]\f$,
 where \f$\widehat{v_i}\f$ means \f$v_i\f$ is omitted from the list. The chain
 groups form a sequence:

 \f[\cdots \ \ \mathbf{C}_n(\mathbf{K},\mathcal{R}) \xrightarrow{\ \partial_n\ } \mathbf{C}_{n-1}(\mathbf{K},\mathcal{R})
 \xrightarrow{\partial_{n-1}} \cdots \xrightarrow{\ \partial_2 \ }
 \mathbf{C}_1(\mathbf{K},\mathcal{R}) \xrightarrow{\ \partial_1 \ }  \mathbf{C}_0(\mathbf{K},\mathcal{R}) \f]

 of finitely many groups \f$\mathbf{C}_n(\mathbf{K},\mathcal{R})\f$ and homomorphisms
 \f$\partial_n\f$, indexed by the dimension \f$n \geq 0\f$.
 The boundary operators satisfy the property \f$\partial_n \circ \partial_{n+1}=0\f$
 for every \f$n > 0\f$
 and we define the homology groups:

 \f[\mathbf{H}_n(\mathbf{K},\mathcal{R}) = \ker \partial_n / \mathrm{im} \  \partial_{n+1}\f]

 We refer to \cite Munkres-elementsalgtop1984 for an introduction to homology
 theory and to \cite DBLP:books/daglib/0025666 for an introduction to persistent homology.

 \section persistenceindexingscheme Indexing Scheme
 "Changing" a simplicial complex consists in applying a simplicial map.
 An <EM>indexing scheme</EM> is a directed graph together with a traversal
 order, such that two
 consecutive nodes in the graph are connected by an arrow (either forward or backward).
 The nodes represent simplicial complexes and the directed edges simplicial maps.

 From the computational point of view, there are two types of indexing schemes of
 interest
 in persistent homology: <EM>linear</EM> ones
 \f$\bullet \longrightarrow \bullet \longrightarrow \cdots \longrightarrow \bullet
 \longrightarrow \bullet\f$
 in persistent homology \cite DBLP:journals/dcg/ZomorodianC05 ,
 and <EM>zigzag</EM> ones
 \f$\bullet \longrightarrow \bullet \longleftarrow \cdots
 \longrightarrow \bullet
 \longleftarrow \bullet \f$ in zigzag persistent
 homology \cite DBLP:journals/focm/CarlssonS10.
 These indexing schemes have a natural left-to-right traversal order, and we
 describe them with ranges and iterators.
 In the current release of the Gudhi library, only the linear case is implemented.

 In the following, we consider the case where the indexing scheme is induced
 by a filtration.
 Ordering the simplices
 by increasing filtration values (breaking ties so as a simplex appears after
 its subsimplices of same filtration value) provides an indexing scheme.

\section pcohexamples Examples

We provide several example files: run these examples with -h for details on their use, and read the README file.

\li \gudhi_example_link{Rips_complex,rips_persistence.cpp} computes the Rips complex of a point cloud and outputs its persistence
diagram.
\code $> ./rips_persistence ../../data/points/tore3D_1307.off -r 0.25 -m 0.5 -d 3 -p 3 \endcode
\code The complex contains 177838 simplices 
   and has dimension 3 
3  0 0 inf 
3  1 0.0983494 inf 
3  1 0.104347 inf 
3  2 0.138335 inf \endcode

More details on the <a href="../../ripscomplex/">Rips complex utilities</a> dedicated page.

\li \gudhi_example_link{Persistent_cohomology,rips_multifield_persistence.cpp} computes the Rips complex of a point cloud and outputs its
persistence diagram with a family of field coefficients.

\li \gudhi_example_link{Rips_complex,rips_distance_matrix_persistence.cpp} computes the Rips complex of a distance matrix and
outputs its persistence diagram.

The file should contain square or lower triangular distance matrix with semicolons as separators.
The code do not check if it is dealing with a distance matrix. It is the user responsibility to provide a valid input.
Please refer to data/distance_matrix/lower_triangular_distance_matrix.csv for an example of a file.

More details on the <a href="../../ripscomplex/">Rips complex utilities</a> dedicated page.

\li \gudhi_example_link{Rips_complex,rips_correlation_matrix_persistence.cpp}
computes the Rips complex of a correlation matrix and outputs its persistence diagram.

Note that no check is performed if the matrix given as the input is a correlation matrix.
It is the user responsibility to ensure that this is the case. The input is to be given either as a square or a lower
triangular matrix.
Please refer to data/correlation_matrix/lower_triangular_correlation_matrix.csv for an example of a file.

More details on the <a href="../../ripscomplex/">Rips complex utilities</a> dedicated page.

\li \gudhi_example_link{Alpha_complex,alpha_complex_3d_persistence.cpp} computes the persistent homology with
\f$\mathbb{Z}/2\mathbb{Z}\f$ coefficients of the alpha complex on points sampling from an OFF file.
\code $> ./alpha_complex_3d_persistence ../../data/points/tore3D_300.off -p 2 -m 0.45 \endcode
\code Simplex_tree dim: 3
2  0 0 inf 
2  1 0.0682162 1.0001 
2  1 0.0934117 1.00003 
2  2 0.56444 1.03938 \endcode

More details on the <a href="../../alphacomplex/">Alpha complex utilities</a> dedicated page.

CGAL can be forced to compute the exact values, it is slower, but it is necessary when points are on a grid
for instance (the fast version `--fast` would give incorrect values).
\code $> ./alpha_complex_3d_persistence ../../data/points/sphere3D_pts_on_grid.off --exact -p 2 -m 0.1 \endcode
\code Simplex_tree dim: 3
2  0 0 inf
2  2 0.0002 0.2028 \endcode

It can also compute the persistent homology with
\f$\mathbb{Z}/2\mathbb{Z}\f$ coefficients of the weighted alpha complex on points sampling from an OFF file
and a weights file.
\code $> ./alpha_complex_3d_persistence ../../data/points/tore3D_300.off
--weight-file ../../data/points/tore3D_300.weights -p 2 -m 0.45 \endcode
\code Simplex_tree dim: 3
2  0 -1 inf
2  1 -0.931784 0.000103311
2  1 -0.906588 2.60165e-05
2  2 -0.43556 0.0393798 \endcode

One can also compute the persistent homology with
\f$\mathbb{Z}/2\mathbb{Z}\f$ coefficients of the periodic alpha complex on points sampling from an OFF file.
The second parameter is a \ref FileFormatsIsoCuboid file with coordinates of the periodic cuboid.
Note that the lengths of the sides of the periodic cuboid have to be the same.
\code $> ./alpha_complex_3d_persistence ../../data/points/grid_10_10_10_in_0_1.off
--cuboid-file ../../data/points/iso_cuboid_3_in_0_1.txt -p 3 -m 1.0 \endcode
\code Simplex_tree dim: 3
3  0 0 inf 
3  1 0.0025 inf 
3  1 0.0025 inf 
3  1 0.0025 inf 
3  2 0.005 inf 
3  2 0.005 inf 
3  2 0.005 inf 
3  3 0.0075 inf \endcode

In order to compute the persistent homology with
\f$\mathbb{Z}/2\mathbb{Z}\f$ coefficients of the periodic alpha complex on weighted points from an OFF file. The
additional parameters of this program are:<br>
(a) The file with the weights of points. The file consist of a sequence of numbers (as many as points).
Note that the weight of each single point have to be bounded by 1/64 times the square of the cuboid edge length.<br>
(b) A \ref FileFormatsIsoCuboid file with coordinates of the periodic cuboid.
Note that the lengths of the sides of the periodic cuboid have to be the same.<br>
\code $> ./alpha_complex_3d_persistence ../../data/points/shifted_sphere.off 
--weight-file ../../data/points/shifted_sphere.weights 
--cuboid-file  ../../data/points/iso_cuboid_3_in_0_10.txt -p 3 -m 1.0 \endcode
\code Simplex_tree dim: 3
3  0 -0.0001 inf
3  1 16.0264 inf
3  1 16.0273 inf
3  1 16.0303 inf
3  2 36.8635 inf
3  2 36.8704 inf
3  2 36.8838 inf
3  3 58.6783 inf  \endcode

\li \gudhi_example_link{Alpha_complex,alpha_complex_persistence.cpp} computes the persistent homology with
\f$\mathbb{Z}/p\mathbb{Z}\f$ coefficients of the alpha complex on points sampling from an OFF file.
\code $> ./alpha_complex_persistence -r 32 -p 2 -m 0.45 ../../data/points/tore3D_300.off \endcode
\code Alpha complex is of dimension 3 - 9273 simplices - 300 vertices.
Simplex_tree dim: 3
2  0 0 inf 
2  1 0.0682162 1.0001 
2  1 0.0934117 1.00003 
2  2 0.56444 1.03938 \endcode

More details on the <a href="../../alphacomplex/">Alpha complex utilities</a> dedicated page.

\li \gudhi_example_link{Persistent_cohomology,plain_homology.cpp} computes the plain homology of a simple simplicial complex without
filtration values.

\section specialcasepersistence Special cases

In order to compute the persistent homology of a piecewise-linear function on \f$\mathbb{R}\f$, the standard strategy would be to
create a path complex (special case of cubical or simplicial complex), define a lower-star filtration on it, and finally
apply a general persistence algorithm. compute_persistence_of_function_on_line() computes this persistence diagram
directly on the function in linear time using the algorithm of \cite arxivpers1d .
\include persistence_on_a_line.cpp
outputs
\code
0.7 1.5
0 3.1
-1 inf
\endcode
 */

}  // namespace persistent_cohomology

}  // namespace Gudhi

#endif  // DOC_PERSISTENT_COHOMOLOGY_INTRO_PERSISTENT_COHOMOLOGY_H_

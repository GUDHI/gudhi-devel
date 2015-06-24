/**
\mainpage 

\image html "Gudhi_banner.jpg" "" width=20cm

The Gudhi library (Geometric Understanding in Higher Dimensions) is a generic C++ library for 
topological analysis of high-dimensional data whose goal is to provide robust, efficient, flexible and easy to use 
implementations of 
state-of-the-art algorithms and data structures for computational topology.
This library is part of the <a href="https://project.inria.fr/gudhi/">Gudhi project</a>.

The current release of the library allows to use several data-structures for simplicial complexes :
simplex tree, Hasse diagram or skeleton-blocker. Several  operations can then be done on top of these
representations such as persistent homology computation or simplification. 
All data-structures are generic and several of their aspects (such as stored elements, policies) 
can be parameterized via template classes.
We refer to 
\cite gudhilibrary_ICMS14
for a detailed description of the design of the library.

\section installation Gudhi installation

As Gudhi is a header only library, there is no need to install the library.

Examples of Gudhi headers inclusion can be found in \ref demos.


\section compiling Compiling

The library uses c++11 and requires <a href="http://www.boost.org/">Boost</a> with version 1.48.0 or more recent.
It is a multi-platform library and compiles on Linux, Mac OSX and Visual Studio 2013. 


\subsection gmp GMP:
The multi-field persistent homology algorithm requires GMP which is a free library for arbitrary-precision
arithmetic, operating on signed integers, rational numbers, and floating point numbers.

The following examples require the <a href="http://gmplib.org/">GNU Multiple Precision Arithmetic Library</a> (GMP) 
and will not be built if GMP is not installed:
  - Persistent_cohomology/rips_multifield_persistence
  - Simplex_tree/simplex_tree_from_alpha_shapes_3

Having GMP version 4.2 or higher installed is recommended.

\subsection cgal CGAL:
CGAL is a C++ library which provides easy access to efficient and reliable geometric algorithms.

The following example requires the <a href="http://www.cgal.org/">Computational Geometry Algorithms Library</a> (CGAL) 
and will not be built if CGAL is not installed:
  - Simplex_tree/simplex_tree_from_alpha_shapes_3
  - Alpha_complex/Alpha_complex_from_off
  - Alpha_complex/Alpha_complex_from_points

Having CGAL version 4.7 or higher installed is recommended. The procedure to install this library according to
your operating system is detailed here http://doc.cgal.org/latest/Manual/installation.html

\subsection demos Demos and examples

To build the demos and libraries, run the following commands in a terminal:

\verbatim
cd /path-to-gudhi/
mkdir build
cd build/
cmake ..
make
\endverbatim

\details 

\copyright GNU General Public License v3.                         
\verbatim  Contact: gudhi-users@lists.gforge.inria.fr \endverbatim

*/

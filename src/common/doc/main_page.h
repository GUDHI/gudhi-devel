/**
\mainpage 

\image html "Gudhi_banner.jpg" "" width=20cm

The Gudhi library (Geometric Understanding in Higher Dimensions) is a generic C++ library for 
topological analysis of high-dimensional data whose goal is to provide robust, efficient, flexible and easy to use 
implementations of 
state-of-the-art algorithms and data structures for computational topology. 

The current release of the library allows to use several data-structures for simplicial complexes :
simplex tree, Hasse diagram or skeleton-blocker. Several  operations can then be done on top of these
representations such a persistent homology computation or simplification.

All data-structures are generic and several of their aspects (such as stored elements, policies) 
can be parameterized via template classes.

We refer to 
\cite gudhilibrary_ICMS14
for a detailed description of the design of the library.


\section compiling Compiling

The library uses c++11 and requires Boost with version 1.48.0 or more recent :  http://www.boost.org/.

The library compiles in Linux and Mac OSX.

\subsection gmp GMP:
The multi-field persistent homology algorithm requires GMP which is a free library for arbitrary-precision
arithmetic, operating on signed integers, rational numbers, and floating point numbers

The following examples require The GNU Multiple Precision Arithmetic Library (GMP) http://gmplib.org/
and will not be built if GMP is not installed:
  - Persistent_cohomology/rips_multifield_persistence
  - Simplex_tree/simplex_tree_from_alpha_shapes_3

Having GMP version 4.2 or higher installed is recommended. This library can be obtained from http://gmplib.org/

\subsection cgal CGAL:
CGAL is a C++ library which provides easy access to efficient and reliable geometric algorithms.

The following example requires CGAL https://www.cgal.org/ and will not be built if CGAL is not installed:
  - Simplex_tree/simplex_tree_from_alpha_shapes_3

Having CGAL version 4.5 or higher installed is recommended. The procedure to install this library according to
your operating system is detailed here http://doc.cgal.org/latest/Manual/installation.html

\section demos Demos and Examples

To build the library, run the following in a terminal:

\verbatim
cd /path-to-gudhi/
mkdir build
cd build/
cmake -DCMAKE_BUILD_TYPE=Release ..
make
\endverbatim





\details 

\copyright GNU General Public License v3.                         
\verbatim  Contact: gudhi-users@lists.gforge.inria.fr \endverbatim

*/

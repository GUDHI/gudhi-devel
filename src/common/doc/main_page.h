/**
\mainpage 

\image html "Gudhi_banner.jpg" "" width=20cm

The Gudhi library (Geometric Understanding in Higher Dimensions) is a generic C++ library for 
topological analysis of high-dimensional data whose goal is to provide robust, efficient, flexible and easy to use 
implementations of 
state-of-the-art algorithms and data structures for computational topology. 

The current release of the library allows to use several data-structures for simplicial complexes :
simplex tree, Hasse diagram or skeleton-blocker. Several  operations can then be done on top of these
representations persistent homology computation, or simplification. 

All data-structures are generic and several of their aspects (such as stored elements, policies) 
can be parametrized via template classes.

We refer to 
\cite gudhilibrary_ICMS14
for a detailed description of the design of the library.


\section Compiling

The library uses c++11 and requires Boost with version 1.48.0 or more recent :  http://www.boost.org/.
Some packages require additional libraries :
\li the multi-field persistent homology algorithm has an optional dependency with GMP
\li Qt demos require CGAL, Qt4 and QGLViewer

The procedure to install these packages according to your operating system is 
detailled here http://doc.cgal.org/latest/Manual/installation.html

The library compiles in Linux, Mac OSX and visual studio 2013.

\section d Demos and Examples

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
\verbatim  Contact: Clément Maria,     clement.maria@inria.fr \endverbatim

*/
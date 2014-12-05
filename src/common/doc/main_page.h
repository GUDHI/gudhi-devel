/**
\mainpage 

The Gudhi library (Geometric Understanding in Higher Dimensions) is a generic C++ library for 
computational topology. Its goal is to provide robust, efficient, flexible and easy to use 
implementations of 
state-of-the-art algorithms and data structures for computational topology. We refer to 
\cite gudhilibrary_ICMS14
for a detailed description of the design of the library.

The current release of the library allows the user to construct representations of simplicial complexes -- 
simplex tree or Hasse diagram -- from a point 
cloud (Rips complex) or 
a list of simplices, and to compute their persistent homology with coefficients in a field 
\f$\mathbb{Z}/p\mathbb{Z}\f$ (for an arbitrary prime \f$p\f$), or simultaneously with coefficients 
in a family of fields (multi-field persistent homology).


To build the library, run the following in a terminal:

\verbatim
cd /path-to-gudhi/
mkdir build
cd build/
cmake -DCMAKE_BUILD_TYPE=Release ..
make
\endverbatim

The library has dependencies with Boost 1.48.0 or more recent (required):  http://www.boost.org/
and with GMP: https://gmplib.org/ The dependency with GMP is optional, and is used only for the 
multi-field persistent homology algorithm.


We provide example files: run these examples with -h for details on their use, and read the README file.

\li <CODE>rips_persistence.cpp</CODE> computes the Rips complex of a point cloud and its persistence diagram.

\li <CODE>rips_multifield_persistence.cpp</CODE> computes the Rips complex of a point cloud and its persistence diagram 
with a family of field coefficients.

\li <CODE>performance_rips_persistence.cpp</CODE> provides timings for the construction of the Rips complex on a set of 
points sampling a Klein bottle in \f$\mathbb{R}^5\f$ with a simplex tree, its conversion to a 
Hasse diagram and the computation of persistent homology and multi-field persistent homology for the 
different representations.


\details 

\copyright GNU General Public License v3.                         
\verbatim  Contact: Clément Maria,     clement.maria@inria.fr \endverbatim

*/
/**
\mainpage 

\image html "Gudhi_banner.jpg" "" width=20cm

The Gudhi library (Geometric Understanding in Higher Dimensions) is a generic C++ library for 
topological analysis of high-dimensional data whose goal is to provide robust, efficient, flexible and easy to use 
implementations of 
state-of-the-art algorithms and data structures for computational topology.
This library is part of the <a class="el" href="https://project.inria.fr/gudhi/">Gudhi project</a>.

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

The library uses c++11 and requires <a target="_blank" href="http://www.boost.org/">Boost</a> with version 1.48.0 or more recent.
It is a multi-platform library and compiles on Linux, Mac OSX and Visual Studio 2013. 


\subsection gmp GMP:
The multi-field persistent homology algorithm requires GMP which is a free library for arbitrary-precision
arithmetic, operating on signed integers, rational numbers, and floating point numbers.

The following examples require the <a target="_blank" href="http://gmplib.org/">GNU Multiple Precision Arithmetic Library</a> (GMP) 
and will not be built if GMP is not installed:
  - Persistent_cohomology/rips_multifield_persistence

Having GMP version 4.2 or higher installed is recommended.

\subsection cgal CGAL:
CGAL is a C++ library which provides easy access to efficient and reliable geometric algorithms.

The following example requires the <a target="_blank" href="http://www.cgal.org/">Computational Geometry Algorithms Library</a> (CGAL) 
and will not be built if CGAL is not installed:
  - GudhUI
  - Persistent_cohomology/alpha_shapes_persistence
  - Simplex_tree/simplex_tree_from_alpha_shapes_3

Having CGAL version 4.4 or higher installed is recommended. The procedure to install this library according to
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

\subsection testsuites Test suites

To test your build, run the following command in a terminal:

\verbatim
make test
\endverbatim

\details 

\copyright GNU General Public License v3.                         
\verbatim  Contact: gudhi-users@lists.gforge.inria.fr \endverbatim

*/

/*! \page Software Software
 * \tableofcontents
 * \section SoftwareIntroduction Introduction
 * The GUDHI library is a C++ open source library **intended to provide** the central data structures and algorithms
 * that underly applications in Geometric and Topological Data Analysis
 * (<a class="el" target="_blank" href="https://en.wikipedia.org/wiki/Topological_data_analysis">TDA</a>). The GUDHI
 * library is developed as part of the <a class="el" target="_blank" href="https://project.inria.fr/gudhi/">GUDHI
 * project</a> supported by the European Research Council. The GUDHI library can both help the development of new
 * algorithmic solutions  and to facilitate the transfer of state of the art results and new applications of TDA.
 * 
 * The current release of the GUDHI library includes:
 * 
 * \li Data structures to represent, construct and manipulate simplicial complexes.
 * \li Algorithms to compute persistent homology and multi-field persistent homology.
 * \li Simplification methods via implicit representations.
 * 
 *
 * The library is available <a class="el" target="_blank" href="https://gforge.inria.fr/frs/?group_id=3865">here</a>
 * and the documentation is available at this <a class="el" href="http://gudhi.gforge.inria.fr/doc/latest/">
 * webpage</a>.
 * 
 * The library comes with data sets, \ref demos and \ref testsuites.
 * 
 * \subsection People People
 * 
 * The development of the GUDHI library is steered by an Editorial Board, which is responsible for guiding the
 * development of the library, developers, and the user community.
 * 
 * The Editorial board is composed of:
 * 
 * \li <a class="el" target="_blank" href="http://www-sop.inria.fr/members/Jean-Daniel.Boissonnat/">
 * Jean-Daniel Boissonnat</a> | INRIA Sophia Antipolis - Méditerranée
 * \li <a class="el" target="_blank" href="http://www.loria.fr/~glisse/">Marc Glisse</a> | INRIA Saclay - Ile de France
 * \li Clément Jamin | INRIA Sophia Antipolis - Méditerranée
 * \li Vincent Rouvreau | INRIA Saclay - Ile de France
 * 
 * \section Contributions Bug reports and contributions
 * Please help us improving the quality of the GUDHI library. You may report bugs or suggestions to:
 * \verbatim  Contact: gudhi-users@lists.gforge.inria.fr \endverbatim
 * 
 * Gudhi is **open** to external contributions. If you want to join our development team, please contact us.
 * 
 * 
 * \section ReleaseHistory Release history
 *
 * \li 24-10-2015; release v.1.2.0, GudhUI (Gudhi Qt demo), Simplex tree coface function, Clang build issue fix.
 * \li 18-12-2014; release v.1.1, Skeleton-Blocker data-structure, simplification package, additional examples for topological persistence.
 * \li 08-12-2014; release v. 1.0.2, initialize simplex keys in initialize_filtration in Simplex_tree
 * \li 07-11-2014: release v. 1.0.1, bug fix in summing columns in Persistent_cohomology
 * \li 23-06-2014: release v. 1.0
 * 
 * \section Upcoming Upcoming
 *
 * The library is under active development. New packages to be released next include:
 * \li Alpha complex.
 * \li Bottleneck distance.
 * \li Zig zag persistence.
 * \li Witness complex.
 * \li Tangential complex.
 * \li Clustering.
*/

/*! \page Citation Acknowledging the GUDHI library
 * We kindly ask users to cite the GUDHI library as appropriately as possible in their papers, and to mention the use
 * of the GUDHI library on the web pages of their projects using GUDHI and provide us with links to these web pages.
 * Feel free to contact us in case you have any question or remark on this topic.
 * 
 * We provide \ref GudhiBibtex entries for the modules of the User and Reference Manual, as well as for publications
 * directly related to the GUDHI library.
 * \section GudhiBibtex GUDHI bibtex
 * \verbinclude  biblio/how_to_cite_gudhi.bib
*/


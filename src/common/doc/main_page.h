/*! \mainpage
 * \image html "Gudhi_banner.jpg" "" width=20cm
 * 
 * \section Introduction Introduction
 * The Gudhi library (Geometric Understanding in Higher Dimensions) is a generic open source C++ library for
 * Computational Topology and Topological Data Analysis
 * (<a class="el" target="_blank" href="https://en.wikipedia.org/wiki/Topological_data_analysis">TDA</a>).
 * The GUDHI library is developed as part of the
 * <a class="el" target="_blank" href="https://project.inria.fr/gudhi/">GUDHI project</a> supported by the European
 * Research Council. The GUDHI library intends  to help the development of new algorithmic solutions in TDA and their
 * transfer to applications. It provides robust, efficient, flexible and easy to use implementations of
 * state-of-the-art algorithms and data structures.
 * 
 * The current release of the GUDHI library includes:
 * 
 * \li Data structures to represent, construct and manipulate simplicial complexes.
 * \li Algorithms to compute persistent homology and multi-field persistent homology.
 * \li Simplication of simplicial complexes by edge contraction.
 * 
 * All data-structures are generic and several of their aspects can be parameterized via template classes.
 * We refer to \cite gudhilibrary_ICMS14 for a detailed description of the design of the library.
 *
 * The library is available <a class="el" target="_blank" href="https://gforge.inria.fr/frs/?group_id=3865">here</a>
 * and the documentation is available at this <a class="el" href="http://gudhi.gforge.inria.fr/doc/latest/">
 * webpage</a>.
 * 
 * The library comes with data sets, \ref demos and \ref testsuites.
 * 
 * Gudhi is also accessible though the
 * <a class="el" target="_blank" href="https://cran.r-project.org/web/packages/TDA/index.html">R package TDA</a>
 * (Statistical Tools for Topological Data Analysis).
 * 
 * The development of the GUDHI library is steered by an Editorial Board composed of:
 * 
 * \li <a class="el" target="_blank" href="http://www-sop.inria.fr/members/Jean-Daniel.Boissonnat/">
 * Jean-Daniel Boissonnat</a> | INRIA Sophia Antipolis - Méditerranée
 * \li <a class="el" target="_blank" href="http://geometrica.saclay.inria.fr/team/Marc.Glisse/">Marc Glisse</a> | INRIA Saclay - Ile de France
 * \li Clément Jamin | INRIA Sophia Antipolis - Méditerranée
 * \li Vincent Rouvreau | INRIA Saclay - Ile de France
 * 
*/

/*! \page installation Gudhi installation
 * As Gudhi is a header only library, there is no need to install the library.
 * 
 * Examples of Gudhi headers inclusion can be found in \ref demos.
 * 
 * \section compiling Compiling
 * The library uses c++11 and requires <a target="_blank" href="http://www.boost.org/">Boost</a> with version 1.48.0 or
 * more recent. It is a multi-platform library and compiles on Linux, Mac OSX and Visual Studio 2013.
 * 
 * \subsection gmp GMP:
 * The multi-field persistent homology algorithm requires GMP which is a free library for arbitrary-precision
 * arithmetic, operating on signed integers, rational numbers, and floating point numbers.
 * 
 * The following example requires the <a target="_blank" href="http://gmplib.org/">GNU Multiple Precision Arithmetic
 * Library</a> (GMP) and will not be built if GMP is not installed:
 * \li Persistent_cohomology/rips_multifield_persistence
 *
 * Having GMP version 4.2 or higher installed is recommended.
 * 
 * \subsection cgal CGAL:
 * CGAL is a C++ library which provides easy access to efficient and reliable geometric algorithms.
 * 
 * The following examples require the <a target="_blank" href="http://www.cgal.org/">Computational Geometry Algorithms
 * Library</a> (CGAL) and will not be built if CGAL is not installed:
 * \li GudhUI
 * \li Persistent_cohomology/alpha_shapes_persistence
 * \li Simplex_tree/simplex_tree_from_alpha_shapes_3
 * \li Alpha_complex/Alpha_complex_from_off
 * \li Alpha_complex/Alpha_complex_from_points
 * 
 * Having CGAL version 4.4 or higher installed is recommended. The procedure to install this library according to
 * your operating system is detailed here http://doc.cgal.org/latest/Manual/installation.html
 * 
 * \subsection demos Demos and examples
 * To build the demos and libraries, run the following commands in a terminal:
\verbatim  cd /path-to-gudhi/
mkdir build
cd build/
cmake ..
make \endverbatim
 * 
 * \subsection testsuites Test suites
 * To test your build, run the following command in a terminal:
 * \verbatim  make test \endverbatim
 * 
 * \section Contributions Bug reports and contributions
 * Please help us improving the quality of the GUDHI library. You may report bugs or suggestions to:
 * \verbatim  Contact: gudhi-users@lists.gforge.inria.fr \endverbatim
 * 
 * Gudhi is open to external contributions. If you want to join our development team, please contact us.
 * 
*/

/*! \page Upcoming Upcoming
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


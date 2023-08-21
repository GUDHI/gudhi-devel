/*! \page installation GUDHI installation
 *  \tableofcontents
 * As GUDHI is a header only library, there is no need to install the library.
 * 
 * Examples of GUDHI headers inclusion can be found in \ref utilities.
 * 
 * \section compiling Compiling
 * The library uses c++17 and requires <a target="_blank" href="https://www.boost.org/">Boost</a>  &ge; 1.71.0
 * and <a target="_blank" href="https://cmake.org/">CMake</a> &ge; 3.8.
 * It is a multi-platform library and compiles on Linux, Mac OSX and Visual Studio 2017.
 * 
 * \subsection utilities Utilities and examples
 * To build the utilities, run the following commands in a terminal:
\verbatim  cd /path-to-gudhi/
mkdir build
cd build/
cmake -DCMAKE_BUILD_TYPE=Release ..
make \endverbatim
 * By default, examples are disabled. You can activate their compilation with
 * <a href="https://cmake.org/cmake/help/latest/manual/ccmake.1.html">ccmake</a> (on Linux and Mac OSX),
 * <a href="https://cmake.org/cmake/help/latest/manual/cmake-gui.1.html">cmake-gui</a> (on Windows) or by modifying the
 * cmake command as follows :
\verbatim  cmake -DCMAKE_BUILD_TYPE=Release -DWITH_GUDHI_EXAMPLE=ON ..
make \endverbatim
 * A list of utilities and examples is available <a href="examples.html">here</a>.
 * 
 * \subsection libraryinstallation Installation
 * To install the library (headers and activated utilities), run the following command in a terminal:
 * \verbatim  make install \endverbatim
 * This action may require to be in the sudoer or administrator of the machine in function of the operating system and
 * of <a href="https://cmake.org/cmake/help/latest/variable/CMAKE_INSTALL_PREFIX.html">CMAKE_INSTALL_PREFIX</a>.
 *
 * \note Python module will be compiled by the `make` command, but `make install` will not install it. Please refer to
 * the <a href="https://gudhi.inria.fr/python/latest/installation.html#gudhi-python-module-installation">Python
 * module installation documentation</a>.
 * 
 * \subsection testsuites Test suites
 * To test your build, run the following command in a terminal:
 * \verbatim make test \endverbatim
 * `make test` is using <a href="https://cmake.org/cmake/help/latest/manual/ctest.1.html">Ctest</a> (CMake test driver
 * program). If some of the tests are failing, please send us the result of the following command:
 * \verbatim ctest --output-on-failure \endverbatim
 * Testing fetching datasets feature requires the use of the internet and is disabled by default. If you want to include this test, set WITH_GUDHI_REMOTE_TEST to ON when building in the previous step (note that this test is included in the python module):
 * \verbatim cmake -DCMAKE_BUILD_TYPE=Release -DWITH_GUDHI_TEST=ON -DWITH_GUDHI_REMOTE_TEST=ON --DWITH_GUDHI_PYTHON=ON .. \endverbatim
 * 
 * \subsection documentationgeneration C++ documentation
 * To generate the C++ documentation, the <a target="_blank" href="http://www.doxygen.org/">doxygen</a> program
 * is required (version &ge; 1.9.5 is advised). Run the following command in a terminal:
 * \verbatim make doxygen \endverbatim
 * Documentation will be generated in a folder named <code>html</code>.
 *
 * In case there is not a full setup present and only the documentation should be build the following command sequence
 * can be used:
\verbatim  cmake -DWITH_GUDHI_THIRD_PARTY=OFF ..
make doxygen\endverbatim
 *
 * \subsection helloworld Hello world !
 * The <a target="_blank" href="https://github.com/GUDHI/hello-gudhi-world">Hello world for GUDHI</a>
 * project is an example to help developers to make their own C++ project on top of the GUDHI library.
 *
 * \section optionallibrary Optional third-party library
 * \subsection gmp GMP
 * The multi-field persistent homology algorithm requires GMP which is a free library for arbitrary-precision
 * arithmetic, operating on signed integers, rational numbers, and floating point numbers.
 * 
 * The following example requires the <a target="_blank" href="https://gmplib.org/">GNU Multiple Precision Arithmetic
 * Library</a> (GMP) and will not be built if GMP is not installed:
 * \li \gudhi_example_link{Persistent_cohomology,rips_multifield_persistence.cpp}
 *
 * Having GMP version 4.2 or higher installed is recommended.
 * 
 * \subsection cgal CGAL
 * Some GUDHI modules (cf. \ref main_page "modules list"), and few examples require CGAL, a C++ library that provides
 * easy access to efficient and reliable geometric algorithms.
 *
 * \note There is no need to install CGAL, you can just <CODE>cmake -DCMAKE_BUILD_TYPE=Release . && make</CODE> CGAL
 * (or even <CODE>cmake -DCMAKE_BUILD_TYPE=Release -DCGAL_HEADER_ONLY=ON .</CODE>), thereafter you will be able to
 * compile GUDHI by calling <CODE>cmake -DCMAKE_BUILD_TYPE=Release -DCGAL_DIR=/your/path/to/CGAL-X.Y .. && make</CODE>
 * 
 * The procedure to install this library according to
 * your operating system is detailed here http://doc.cgal.org/latest/Manual/installation.html
 * 
 * The following examples/utilities require the <a target="_blank" href="http://www.cgal.org/">Computational Geometry Algorithms
 * Library</a> (CGAL \cite cgal:eb-19b) and will not be built if CGAL version 4.11.0 or higher is not installed:
 * \li \gudhi_example_link{Simplex_tree,example_alpha_shapes_3_simplex_tree_from_off_file.cpp}
 * \li \gudhi_example_link{Witness_complex,strong_witness_persistence.cpp}
 * \li \gudhi_example_link{Witness_complex,weak_witness_persistence.cpp}
 * \li \gudhi_example_link{Witness_complex,example_strong_witness_complex_off.cpp}
 * \li \gudhi_example_link{Witness_complex,example_witness_complex_off.cpp}
 * \li \gudhi_example_link{Witness_complex,example_witness_complex_sphere.cpp}
 * \li \gudhi_example_link{Alpha_complex,Alpha_complex_from_off.cpp}
 * \li \gudhi_example_link{Alpha_complex,Alpha_complex_from_points.cpp}
 * \li \gudhi_example_link{Alpha_complex,alpha_complex_persistence.cpp}
 * \li \gudhi_example_link{Persistent_cohomology,custom_persistence_sort.cpp}
 * \li \gudhi_example_link{Bottleneck_distance,alpha_rips_persistence_bottleneck_distance.cpp}
 * \li \gudhi_example_link{Bottleneck_distance,bottleneck_basic_example.cpp}
 * \li \gudhi_example_link{Bottleneck_distance,bottleneck_distance.cpp}
 * \li \gudhi_example_link{Nerve_GIC,CoordGIC.cpp}
 * \li \gudhi_example_link{Nerve_GIC,FuncGIC.cpp}
 * \li \gudhi_example_link{Nerve_GIC,Nerve.cpp}
 * \li \gudhi_example_link{Nerve_GIC,VoronoiGIC.cpp}
 * \li \gudhi_example_link{Spatial_searching,example_spatial_searching.cpp}
 * \li \gudhi_example_link{Subsampling,example_choose_n_farthest_points.cpp}
 * \li \gudhi_example_link{Subsampling,example_pick_n_random_points.cpp}
 * \li \gudhi_example_link{Subsampling,example_sparsify_point_set.cpp}
 * \li \gudhi_example_link{Tangential_complex,example_basic.cpp}
 * \li \gudhi_example_link{Tangential_complex,example_with_perturb.cpp}
 * \li \gudhi_example_link{Alpha_complex,Weighted_alpha_complex_3d_from_points.cpp}
 * \li \gudhi_example_link{Alpha_complex,alpha_complex_3d_persistence.cpp}
 * \li \gudhi_example_link{Coxeter_triangulation,manifold_tracing_flat_torus_with_boundary.cpp}
 *
 * \subsection eigen Eigen
 * Some GUDHI modules (cf. \ref main_page "modules list"), and few examples require
 * <a target="_blank" href="https://eigen.tuxfamily.org">Eigen</a> is a C++ template library for linear algebra:
 * matrices, vectors, numerical solvers, and related algorithms.
 * 
 * The following examples/utilities require the <a target="_blank" href="https://eigen.tuxfamily.org">Eigen</a> and will not be
 * built if Eigen is not installed:
 * \li \gudhi_example_link{Alpha_complex,Alpha_complex_from_off.cpp}
 * \li \gudhi_example_link{Alpha_complex,Alpha_complex_from_points.cpp}
 * \li \gudhi_example_link{Alpha_complex,alpha_complex_persistence.cpp}
 * \li \gudhi_example_link{Alpha_complex,alpha_complex_3d_persistence.cpp}
 * \li \gudhi_example_link{Alpha_complex,Weighted_alpha_complex_3d_from_points.cpp}
 * \li \gudhi_example_link{Bottleneck_distance,alpha_rips_persistence_bottleneck_distance.cpp}
 * \li \gudhi_example_link{Persistent_cohomology,custom_persistence_sort.cpp}
 * \li \gudhi_example_link{Spatial_searching,example_spatial_searching.cpp}
 * \li \gudhi_example_link{Subsampling,example_choose_n_farthest_points.cpp}
 * \li \gudhi_example_link{Subsampling,example_pick_n_random_points.cpp}
 * \li \gudhi_example_link{Subsampling,example_sparsify_point_set.cpp}
 * \li \gudhi_example_link{Tangential_complex,example_basic.cpp}
 * \li \gudhi_example_link{Tangential_complex,example_with_perturb.cpp}
 * \li \gudhi_example_link{Witness_complex,strong_witness_persistence.cpp}
 * \li \gudhi_example_link{Witness_complex,weak_witness_persistence.cpp}
 * \li \gudhi_example_link{Witness_complex,example_strong_witness_complex_off.cpp}
 * \li \gudhi_example_link{Witness_complex,example_witness_complex_off.cpp}
 * \li \gudhi_example_link{Witness_complex,example_witness_complex_sphere.cpp}
 * \li \gudhi_example_link{Coxeter_triangulation,cell_complex_from_basic_circle_manifold.cpp}
 * \li \gudhi_example_link{Coxeter_triangulation,manifold_tracing_custom_function.cpp}
 * \li \gudhi_example_link{Coxeter_triangulation,manifold_tracing_flat_torus_with_boundary.cpp}
 *
 * \subsection tbb oneAPI Threading Building Blocks
 * <a target="_blank" href="https://github.com/oneapi-src/oneTBB">Intel&reg; oneAPI TBB</a> lets you easily write
 * parallel C++ programs that take full advantage of multicore performance, that are portable and composable, and that
 * have future-proof scalability.
 * 
 * Having Intel&reg; oneAPI TBB (version 20.19.7 or higher) installed is recommended to parallelize and accelerate some
 * GUDHI computations.
 * 
 * The following examples/utilities are using Intel&reg; oneAPI TBB if installed:
 * \li \gudhi_example_link{Alpha_complex,Alpha_complex_from_off.cpp}
 * \li \gudhi_example_link{Alpha_complex,Alpha_complex_from_points.cpp}
 * \li \gudhi_example_link{Alpha_complex,alpha_complex_3d_persistence.cpp}
 * \li \gudhi_example_link{Alpha_complex,alpha_complex_persistence.cpp}
 * \li \gudhi_example_link{Bitmap_cubical_complex,cubical_complex_persistence.cpp}
 * \li \gudhi_example_link{Bitmap_cubical_complex,periodic_cubical_complex_persistence.cpp}
 * \li \gudhi_example_link{Bitmap_cubical_complex,Random_bitmap_cubical_complex.cpp}
 * \li \gudhi_example_link{Nerve_GIC,CoordGIC.cpp}
 * \li \gudhi_example_link{Nerve_GIC,FuncGIC.cpp}
 * \li \gudhi_example_link{Nerve_GIC,Nerve.cpp}
 * \li \gudhi_example_link{Nerve_GIC,VoronoiGIC.cpp}
 * \li \gudhi_example_link{Simplex_tree,simple_simplex_tree.cpp}
 * \li \gudhi_example_link{Simplex_tree,example_alpha_shapes_3_simplex_tree_from_off_file.cpp}
 * \li \gudhi_example_link{Simplex_tree,simplex_tree_from_cliques_of_graph.cpp}
 * \li \gudhi_example_link{Simplex_tree,graph_expansion_with_blocker.cpp}
 * \li \gudhi_example_link{Persistent_cohomology,alpha_complex_3d_persistence.cpp}
 * \li \gudhi_example_link{Persistent_cohomology,alpha_complex_persistence.cpp}
 * \li \gudhi_example_link{Persistent_cohomology,rips_persistence_via_boundary_matrix.cpp}
 * \li \gudhi_example_link{Persistent_cohomology,persistence_from_file.cpp}
 * \li \gudhi_example_link{Persistent_cohomology,persistence_from_simple_simplex_tree.cpp}
 * \li \gudhi_example_link{Persistent_cohomology,plain_homology.cpp}
 * \li \gudhi_example_link{Persistent_cohomology,rips_multifield_persistence.cpp}
 * \li \gudhi_example_link{Persistent_cohomology,rips_persistence_step_by_step.cpp}
 * \li \gudhi_example_link{Persistent_cohomology,custom_persistence_sort.cpp}
 * \li \gudhi_example_link{Rips_complex,example_one_skeleton_rips_from_points.cpp}
 * \li \gudhi_example_link{Rips_complex,example_rips_complex_from_off_file.cpp}
 * \li \gudhi_example_link{Rips_complex,rips_distance_matrix_persistence.cpp}
 * \li \gudhi_example_link{Rips_complex,rips_persistence.cpp}
 * \li \gudhi_example_link{Witness_complex,strong_witness_persistence.cpp}
 * \li \gudhi_example_link{Witness_complex,weak_witness_persistence.cpp}
 * \li \gudhi_example_link{Witness_complex,example_nearest_landmark_table.cpp}
 *
 * \section Contributions Bug reports and contributions
 * Please help us improving the quality of the GUDHI library.
 * You may <a href="https://github.com/GUDHI/gudhi-devel/issues">report bugs</a> or
 * <a href="https://gudhi.inria.fr/contact/">contact us</a> for any suggestions.
 * 
 * GUDHI is open to external contributions. If you want to join our development team, please take some time to read our
 * <a href="https://github.com/GUDHI/gudhi-devel/blob/master/.github/CONTRIBUTING.md">contributing guide</a>.
 * 
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

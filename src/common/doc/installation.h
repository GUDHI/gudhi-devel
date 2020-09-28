/*! \page installation GUDHI installation
 *  \tableofcontents
 * As GUDHI is a header only library, there is no need to install the library.
 * 
 * Examples of GUDHI headers inclusion can be found in \ref utilities.
 * 
 * \section compiling Compiling
 * The library uses c++14 and requires <a target="_blank" href="http://www.boost.org/">Boost</a>  &ge; 1.56.0
 * and <a target="_blank" href="https://www.cmake.org/">CMake</a> &ge; 3.5.
 * It is a multi-platform library and compiles on Linux, Mac OSX and Visual Studio 2015.
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
 * \subsection testsuites Test suites
 * To test your build, run the following command in a terminal:
 * \verbatim make test \endverbatim
 * `make test` is using <a href="https://cmake.org/cmake/help/latest/manual/ctest.1.html">Ctest</a> (CMake test driver
 * program). If some of the tests are failing, please send us the result of the following command:
 * \verbatim ctest --output-on-failure \endverbatim
 * 
 * \subsection documentationgeneration Documentation
 * To generate the documentation, <a target="_blank" href="http://www.doxygen.org/">Doxygen</a> is required.
 * Run the following command in a terminal:
\verbatim
make doxygen
# Documentation will be generated in the folder YYYY-MM-DD-hh-mm-ss_GUDHI_X.Y.Z/doc/html/
# You can customize the directory name by calling `cmake -DUSER_VERSION_DIR=/my/custom/folder`
\endverbatim
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
 * The following example requires the <a target="_blank" href="http://gmplib.org/">GNU Multiple Precision Arithmetic
 * Library</a> (GMP) and will not be built if GMP is not installed:
 * \li <a href="_persistent_cohomology_2rips_multifield_persistence_8cpp-example.html">
 * Persistent_cohomology/rips_multifield_persistence.cpp</a>
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
 * \li <a href="_simplex_tree_2example_alpha_shapes_3_simplex_tree_from_off_file_8cpp-example.html">
 * Simplex_tree/example_alpha_shapes_3_simplex_tree_from_off_file.cpp</a>
 * \li <a href="_witness_complex_2strong_witness_persistence_8cpp-example.html">
 * Witness_complex/strong_witness_persistence.cpp</a>
 * \li <a href="_witness_complex_2weak_witness_persistence_8cpp-example.html">
 * Witness_complex/weak_witness_persistence.cpp</a>
 * \li <a href="_witness_complex_2example_strong_witness_complex_off_8cpp-example.html">
 * Witness_complex/example_strong_witness_complex_off.cpp</a>
 * \li <a href="_witness_complex_2example_witness_complex_off_8cpp-example.html">
 * Witness_complex/example_witness_complex_off.cpp</a>
 * \li <a href="_witness_complex_2example_witness_complex_sphere_8cpp-example.html">
 * Witness_complex/example_witness_complex_sphere.cpp</a>
 * \li <a href="_alpha_complex_2_alpha_complex_from_off_8cpp-example.html">
 * Alpha_complex/Alpha_complex_from_off.cpp</a>
 * \li <a href="_alpha_complex_2_alpha_complex_from_points_8cpp-example.html">
 * Alpha_complex/Alpha_complex_from_points.cpp</a>
 * \li <a href="_alpha_complex_2alpha_complex_persistence_8cpp-example.html">
 * Alpha_complex/alpha_complex_persistence.cpp</a>
 * \li <a href="_persistent_cohomology_2custom_persistence_sort_8cpp-example.html">
 * Persistent_cohomology/custom_persistence_sort.cpp</a>
 * \li <a href="_bottleneck_distance_2alpha_rips_persistence_bottleneck_distance_8cpp-example.html">
 * Bottleneck_distance/alpha_rips_persistence_bottleneck_distance.cpp.cpp</a>
 * \li <a href="_bottleneck_distance_2bottleneck_basic_example_8cpp-example.html">
 * Bottleneck_distance/bottleneck_basic_example.cpp</a>
 * \li <a href="_bottleneck_distance_2bottleneck_read_file_8cpp-example.html">
 * Bottleneck_distance/bottleneck_distance.cpp</a>
 * \li <a href="_nerve__g_i_c_2_coord_g_i_c_8cpp-example.html">
 * Nerve_GIC/CoordGIC.cpp</a>
 * \li <a href="_nerve__g_i_c_2_func_g_i_c_8cpp-example.html">
 * Nerve_GIC/FuncGIC.cpp</a>
 * \li <a href="_nerve__g_i_c_2_nerve_8cpp-example.html">
 * Nerve_GIC/Nerve.cpp</a>
 * \li <a href="_nerve__g_i_c_2_voronoi_g_i_c_8cpp-example.html">
 * Nerve_GIC/VoronoiGIC.cpp</a>
 * \li <a href="_spatial_searching_2example_spatial_searching_8cpp-example.html">
 * Spatial_searching/example_spatial_searching.cpp</a>
 * \li <a href="_subsampling_2example_choose_n_farthest_points_8cpp-example.html">
 * Subsampling/example_choose_n_farthest_points.cpp</a>
 * \li <a href="_subsampling_2example_custom_kernel_8cpp-example.html">
 * Subsampling/example_custom_kernel.cpp</a>
 * \li <a href="_subsampling_2example_pick_n_random_points_8cpp-example.html">
 * Subsampling/example_pick_n_random_points.cpp</a>
 * \li <a href="_subsampling_2example_sparsify_point_set_8cpp-example.html">
 * Subsampling/example_sparsify_point_set.cpp</a>
 * \li <a href="_tangential_complex_2example_basic_8cpp-example.html">
 * Tangential_complex/example_basic.cpp</a>
 * \li <a href="_tangential_complex_2example_with_perturb_8cpp-example.html">
 * Tangential_complex/example_with_perturb.cpp</a>
 * \li <a href="_alpha_complex_2_weighted_alpha_complex_3d_from_points_8cpp-example.html">
 * Alpha_complex/Weighted_alpha_complex_3d_from_points.cpp</a>
 * \li <a href="_alpha_complex_2alpha_complex_3d_persistence_8cpp-example.html">
 * Alpha_complex/alpha_complex_3d_persistence.cpp</a>
 *
 * \subsection eigen Eigen
 * Some GUDHI modules (cf. \ref main_page "modules list"), and few examples require
 * <a target="_blank" href="http://eigen.tuxfamily.org/">Eigen</a> is a C++ template library for linear algebra:
 * matrices, vectors, numerical solvers, and related algorithms.
 * 
 * The following examples/utilities require the <a target="_blank" href="http://eigen.tuxfamily.org/">Eigen</a> and will not be
 * built if Eigen is not installed:
 * \li <a href="_alpha_complex_2_alpha_complex_from_off_8cpp-example.html">
 * Alpha_complex/Alpha_complex_from_off.cpp</a>
 * \li <a href="_alpha_complex_2_alpha_complex_from_points_8cpp-example.html">
 * Alpha_complex/Alpha_complex_from_points.cpp</a>
 * \li <a href="_alpha_complex_2alpha_complex_persistence_8cpp-example.html">
 * Alpha_complex/alpha_complex_persistence.cpp</a>
 * \li <a href="_alpha_complex_2alpha_complex_3d_persistence_8cpp-example.html">
 * Alpha_complex/alpha_complex_3d_persistence.cpp</a>
 * \li <a href="_alpha_complex_2_weighted_alpha_complex_3d_from_points_8cpp-example.html">
 * Alpha_complex/Weighted_alpha_complex_3d_from_points.cpp</a>
 * \li <a href="_bottleneck_distance_2alpha_rips_persistence_bottleneck_distance_8cpp-example.html">
 * Bottleneck_distance/alpha_rips_persistence_bottleneck_distance.cpp.cpp</a>
 * \li <a href="_persistent_cohomology_2custom_persistence_sort_8cpp-example.html">
 * Persistent_cohomology/custom_persistence_sort.cpp</a>
 * \li <a href="_spatial_searching_2example_spatial_searching_8cpp-example.html">
 * Spatial_searching/example_spatial_searching.cpp</a>
 * \li <a href="_subsampling_2example_choose_n_farthest_points_8cpp-example.html">
 * Subsampling/example_choose_n_farthest_points.cpp</a>
 * \li <a href="_subsampling_2example_custom_kernel_8cpp-example.html">
 * Subsampling/example_custom_kernel.cpp</a>
 * \li <a href="_subsampling_2example_pick_n_random_points_8cpp-example.html">
 * Subsampling/example_pick_n_random_points.cpp</a>
 * \li <a href="_subsampling_2example_sparsify_point_set_8cpp-example.html">
 * Subsampling/example_sparsify_point_set.cpp</a>
 * \li <a href="_tangential_complex_2example_basic_8cpp-example.html">
 * Tangential_complex/example_basic.cpp</a>
 * \li <a href="_tangential_complex_2example_with_perturb_8cpp-example.html">
 * Tangential_complex/example_with_perturb.cpp</a>
 * \li <a href="_witness_complex_2strong_witness_persistence_8cpp-example.html">
 * Witness_complex/strong_witness_persistence.cpp</a>
 * \li <a href="_witness_complex_2weak_witness_persistence_8cpp-example.html">
 * Witness_complex/weak_witness_persistence.cpp</a>
 * \li <a href="_witness_complex_2example_strong_witness_complex_off_8cpp-example.html">
 * Witness_complex/example_strong_witness_complex_off.cpp</a>
 * \li <a href="_witness_complex_2example_witness_complex_off_8cpp-example.html">
 * Witness_complex/example_witness_complex_off.cpp</a>
 * \li <a href="_witness_complex_2example_witness_complex_sphere_8cpp-example.html">
 * Witness_complex/example_witness_complex_sphere.cpp</a>
 *
 * \subsection tbb Threading Building Blocks
 * <a target="_blank" href="https://www.threadingbuildingblocks.org/">Intel&reg; TBB</a> lets you easily write parallel
 * C++ programs that take full advantage of multicore performance, that are portable and composable, and that have
 * future-proof scalability.
 * 
 * Having Intel&reg; TBB installed is recommended to parallelize and accelerate some GUDHI computations.
 * 
 * The following examples/utilities are using Intel&reg; TBB if installed:
 * \li <a href="_alpha_complex_2_alpha_complex_from_off_8cpp-example.html">
 * Alpha_complex/Alpha_complex_from_off.cpp</a>
 * \li <a href="_alpha_complex_2_alpha_complex_from_points_8cpp-example.html">
 * Alpha_complex/Alpha_complex_from_points.cpp</a>
 * \li <a href="_alpha_complex_2alpha_complex_3d_persistence_8cpp-example.html">
 * Alpha_complex/alpha_complex_3d_persistence.cpp</a>
 * \li <a href="_alpha_complex_2alpha_complex_persistence_8cpp-example.html">
 * Alpha_complex/alpha_complex_persistence.cpp</a>
 * \li <a href="_bitmap_cubical_complex_2_bitmap_cubical_complex_8cpp-example.html">
 * Bitmap_cubical_complex/cubical_complex_persistence.cpp</a>
 * \li <a href="_bitmap_cubical_complex_2_bitmap_cubical_complex_periodic_boundary_conditions_8cpp-example.html">
 * Bitmap_cubical_complex/periodic_cubical_complex_persistence.cpp</a>
 * \li <a href="_bitmap_cubical_complex_2_random_bitmap_cubical_complex_8cpp-example.html">
 * Bitmap_cubical_complex/Random_bitmap_cubical_complex.cpp</a>
 * \li <a href="_nerve__g_i_c_2_coord_g_i_c_8cpp-example.html">
 * Nerve_GIC/CoordGIC.cpp</a>
 * \li <a href="_nerve__g_i_c_2_func_g_i_c_8cpp-example.html">
 * Nerve_GIC/FuncGIC.cpp</a>
 * \li <a href="_nerve__g_i_c_2_nerve_8cpp-example.html">
 * Nerve_GIC/Nerve.cpp</a>
 * \li <a href="_nerve__g_i_c_2_voronoi_g_i_c_8cpp-example.html">
 * Nerve_GIC/VoronoiGIC.cpp</a>
 * \li <a href="_simplex_tree_2simple_simplex_tree_8cpp-example.html">
 * Simplex_tree/simple_simplex_tree.cpp</a>
 * \li <a href="_simplex_tree_2example_alpha_shapes_3_simplex_tree_from_off_file_8cpp-example.html">
 * Simplex_tree/example_alpha_shapes_3_simplex_tree_from_off_file.cpp</a>
 * \li <a href="_simplex_tree_2simplex_tree_from_cliques_of_graph_8cpp-example.html">
 * Simplex_tree/simplex_tree_from_cliques_of_graph.cpp</a>
 * \li <a href="_simplex_tree_2graph_expansion_with_blocker_8cpp-example.html">
 * Simplex_tree/graph_expansion_with_blocker.cpp</a>
 * \li <a href="_persistent_cohomology_2alpha_complex_3d_persistence_8cpp-example.html">
 * Persistent_cohomology/alpha_complex_3d_persistence.cpp</a>
 * \li <a href="_persistent_cohomology_2alpha_complex_persistence_8cpp-example.html">
 * Persistent_cohomology/alpha_complex_persistence.cpp</a>
 * \li <a href="_persistent_cohomology_2rips_persistence_via_boundary_matrix_8cpp-example.html">
 * Persistent_cohomology/rips_persistence_via_boundary_matrix.cpp</a>
 * \li <a href="_persistent_cohomology_2persistence_from_file_8cpp-example.html">
 * Persistent_cohomology/persistence_from_file.cpp</a>
 * \li <a href="_persistent_cohomology_2persistence_from_simple_simplex_tree_8cpp-example.html">
 * Persistent_cohomology/persistence_from_simple_simplex_tree.cpp</a>
 * \li <a href="_persistent_cohomology_2plain_homology_8cpp-example.html">
 * Persistent_cohomology/plain_homology.cpp</a>
 * \li <a href="_persistent_cohomology_2rips_multifield_persistence_8cpp-example.html">
 * Persistent_cohomology/rips_multifield_persistence.cpp</a>
 * \li <a href="_persistent_cohomology_2rips_persistence_step_by_step_8cpp-example.html">
 * Persistent_cohomology/rips_persistence_step_by_step.cpp</a>
 * \li <a href="_persistent_cohomology_2custom_persistence_sort_8cpp-example.html">
 * Persistent_cohomology/custom_persistence_sort.cpp</a>
 * \li <a href="_rips_complex_2example_one_skeleton_rips_from_points_8cpp-example.html">
 * Rips_complex/example_one_skeleton_rips_from_points.cpp</a>
 * \li <a href="_rips_complex_2example_rips_complex_from_off_file_8cpp-example.html">
 * Rips_complex/example_rips_complex_from_off_file.cpp</a>
 * \li <a href="_rips_complex_2rips_distance_matrix_persistence_8cpp-example.html">
 * Rips_complex/rips_distance_matrix_persistence.cpp</a>
 * \li <a href="_rips_complex_2rips_persistence_8cpp-example.html">
 * Rips_complex/rips_persistence.cpp</a>
 * \li <a href="_witness_complex_2strong_witness_persistence_8cpp-example.html">
 * Witness_complex/strong_witness_persistence.cpp</a>
 * \li <a href="_witness_complex_2weak_witness_persistence_8cpp-example.html">
 * Witness_complex/weak_witness_persistence.cpp</a>
 * \li <a href="_witness_complex_2example_nearest_landmark_table_8cpp-example.html">
 * Witness_complex/example_nearest_landmark_table.cpp</a>
 *
 * \section Contributions Bug reports and contributions
 * Please help us improving the quality of the GUDHI library. You may report bugs or suggestions to:
 * \verbatim  Contact: gudhi-users@lists.gforge.inria.fr \endverbatim
 * 
 * GUDHI is open to external contributions. If you want to join our development team, please contact us.
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

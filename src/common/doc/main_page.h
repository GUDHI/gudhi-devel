/*! \mainpage
 *  \tableofcontents
 * \image html "Gudhi_banner.png" "" width=20cm
 * 
 * \section Introduction Introduction
 * The Gudhi library (Geometry Understanding in Higher Dimensions) is a generic open source C++ library for
 * Computational Topology and Topological Data Analysis
 * (<a class="el" target="_blank" href="https://en.wikipedia.org/wiki/Topological_data_analysis">TDA</a>).
 * The GUDHI library intends  to help the development of new algorithmic solutions in TDA and their transfer to
 * applications. It provides robust, efficient, flexible and easy to use implementations of state-of-the-art
 * algorithms and data structures.
 * 
 * The current release of the GUDHI library includes:
 * 
 * \li Data structures to represent, construct and manipulate simplicial complexes.
 * \li Simplification of simplicial complexes by edge contraction.
 * \li Algorithms to compute persistent homology persistent homology.
 * 
 * All data-structures are generic and several of their aspects can be parameterized via template classes.
 * We refer to \cite gudhilibrary_ICMS14 for a detailed description of the design of the library.
 *
 \section DataStructures Data structures
 \subsection AlphaComplexDataStructure Alpha complex
 \image html "alpha_complex_representation.png" "Alpha complex representation"
<table border="0">
  <tr>
    <td width="25%">
      <b>Author:</b> Vincent Rouvreau<br>
      <b>Introduced in:</b> GUDHI 1.3.0<br>
      <b>Copyright:</b> GPL v3<br>
    </td>
    <td width="75%">
    Alpha_complex is a simplicial complex constructed from the finite cells of a Delaunay Triangulation.<br>
    The filtration value of each simplex is computed as the square of the circumradius of the simplex if the
    circumsphere is empty (the simplex is then said to be Gabriel), and as the minimum of the filtration
    values of the codimension 1 cofaces that make it not Gabriel otherwise.
    All simplices that have a filtration value strictly greater than a given alpha squared value are not inserted into
    the complex.<br>
    <b>User manual:</b> \ref alpha_complex - <b>Reference manual:</b> Gudhi::alpha_complex::Alpha_complex
    </td>
 </tr>
</table>
 \subsection CubicalComplexDataStructure Cubical complex
 \image html "Cubical_complex_representation.png" "Cubical complex representation"
<table border="0">
  <tr>
    <td width="25%">
      <b>Author:</b> Pawel Dlotko<br>
      <b>Introduced in:</b> GUDHI 1.3.0<br>
      <b>Copyright:</b> GPL v3<br>
    </td>
    <td width="75%">
    The cubical complex is an example of a structured complex useful in computational mathematics (specially
    rigorous numerics) and image analysis.<br>
    <b>User manual:</b> \ref cubical_complex - <b>Reference manual:</b> Gudhi::cubical_complex::Bitmap_cubical_complex
    </td>
 </tr>
 \subsection RipsComplexDataStructure Rips complex
 \image html "rips_complex_representation.png" "Rips complex representation"
<table border="0">
  <tr>
    <td width="25%">
      <b>Author:</b> Cl&eacute;ment Maria, Pawel Dlotko, Vincent Rouvreau<br>
      <b>Introduced in:</b> GUDHI 1.4.0<br>
      <b>Copyright:</b> GPL v3<br>
    </td>
    <td width="75%">
    Rips_complex is a simplicial complex constructed from a one skeleton graph.<br>
    The filtration value of each edge is computed from a user-given distance function and is inserted until a
    user-given threshold value.<br>
    This complex can be built from a point cloud and a distance function, or from a distance matrix.<br>
    <b>User manual:</b> \ref rips_complex - <b>Reference manual:</b> Gudhi::rips_complex::Rips_complex
    </td>
 </tr>
</table>
</table>
 \subsection SimplexTreeDataStructure Simplex tree
 \image html "Simplex_tree_representation.png" "Simplex tree representation"
<table border="0">
  <tr>
    <td width="25%">
      <b>Author:</b> Cl&eacute;ment Maria<br>
      <b>Introduced in:</b> GUDHI 1.0.0<br>
      <b>Copyright:</b> GPL v3<br>
    </td>
    <td width="75%">
    The simplex tree is an efficient and flexible
 data structure for representing general (filtered) simplicial complexes. The data structure
 is described in \cite boissonnatmariasimplextreealgorithmica .<br>
    <b>User manual:</b> \ref simplex_tree - <b>Reference manual:</b> Gudhi::Simplex_tree
    </td>
 </tr>
</table>
 \subsection SkeletonBlockerDataStructure Skeleton blocker
 \image html "ds_representation.png" "Skeleton blocker representation"
<table border="0">
  <tr>
    <td width="25%">
      <b>Author:</b> David Salinas<br>
      <b>Introduced in:</b> GUDHI 1.1.0<br>
      <b>Copyright:</b> GPL v3<br>
    </td>
    <td width="75%">
    The Skeleton-Blocker data-structure proposes a light encoding for simplicial complexes by storing only an *implicit*
    representation of its simplices \cite socg_blockers_2011,\cite blockers2012. Intuitively, it just stores the
    1-skeleton of a simplicial complex with a graph and the set of its "missing faces" that is very small in practice.
    This data-structure handles all simplicial complexes operations such as simplex enumeration or simplex removal but
    operations that are particularly efficient are operations that do not require simplex enumeration such as edge
    iteration, link computation or simplex contraction.<br>
    <b>User manual:</b> \ref skbl - <b>Reference manual:</b> Gudhi::skeleton_blocker::Skeleton_blocker_complex
    </td>
 </tr>
</table>
 \subsection TangentialComplexDataStructure Tangential complex
 \image html "tc_examples.png" "Tangential complex representation"
<table border="0">
  <tr>
    <td width="25%">
      <b>Author:</b> Cl&eacute;ment Jamin<br>
      <b>Introduced in:</b> GUDHI 1.4.0<br>
      <b>Copyright:</b> GPL v3<br>
    </td>
    <td width="75%">
    A Tangential Delaunay complex is a <a target="_blank" href="https://en.wikipedia.org/wiki/Simplicial_complex">simplicial complex</a>
    designed to reconstruct a \f$ k \f$-dimensional manifold embedded in \f$ d \f$-dimensional Euclidean space. 
    The input is a point sample coming from an unknown manifold.
    The running time depends only linearly on the extrinsic dimension \f$ d \f$
    and exponentially on the intrinsic dimension \f$ k \f$.<br>
    <b>User manual:</b> \ref tangential_complex - <b>Reference manual:</b> Gudhi::tangential_complex::Tangential_complex
    </td>
 </tr>
</table>
 \subsection WitnessComplexDataStructure Witness complex
 \image html "Witness_complex_representation.png" "Witness complex representation"
<table border="0">
  <tr>
    <td width="25%">
      <b>Author:</b> Siargey Kachanovich<br>
      <b>Introduced in:</b> GUDHI 1.3.0<br>
      <b>Copyright:</b> GPL v3<br>
    </td>
    <td width="75%">
    Witness complex \f$ Wit(W,L) \f$  is a simplicial complex defined on two sets of points in \f$\mathbb{R}^D\f$.
    The data structure is described in \cite boissonnatmariasimplextreealgorithmica .<br>
    <b>User manual:</b> \ref witness_complex - <b>Reference manual:</b> Gudhi::witness_complex::SimplicialComplexForWitness
    </td>
 </tr>
</table>
 
 \section Toolbox Toolbox
 \subsection ContractionToolbox Contraction
 \image html "sphere_contraction_representation.png" "Sphere contraction example"
<table border="0">
  <tr>
    <td width="25%">
      <b>Author:</b> David Salinas<br>
      <b>Introduced in:</b> GUDHI 1.1.0<br>
      <b>Copyright:</b> GPL v3<br>
    </td>
    <td width="75%">
    The purpose of this package is to offer a user-friendly interface for edge contraction simplification of huge
    simplicial complexes. It uses the \ref skbl data-structure whose size remains small  during simplification of most
    used geometrical complexes of topological data analysis such as the Rips or the Delaunay complexes. In practice,
    the size of this data-structure is even much lower than the total number of simplices.<br>
    <b>User manual:</b> \ref contr
    </td>
 </tr>
</table>
 \subsection PersistentCohomologyToolbox Persistent Cohomology
 \image html "3DTorus_poch.png" "Rips Persistent Cohomology on a 3D Torus"
<table border="0">
  <tr>
    <td width="25%">
      <b>Author:</b> Cl&eacute;ment Maria<br>
      <b>Introduced in:</b> GUDHI 1.0.0<br>
      <b>Copyright:</b> GPL v3<br>
    </td>
    <td width="75%">
    The theory of homology consists in attaching to a topological space a sequence of (homology) groups, capturing
    global topological features like connected components, holes, cavities, etc. Persistent homology studies the
    evolution -- birth, life and death -- of these features when the topological space is changing. Consequently, the
    theory is essentially composed of three elements: topological spaces, their homology groups and an evolution
    scheme.
    Computation of persistent cohomology using the algorithm of \cite DBLP:journals/dcg/SilvaMV11 and
    \cite DBLP:journals/corr/abs-1208-5018 and the Compressed Annotation Matrix implementation of
    \cite DBLP:conf/esa/BoissonnatDM13 .<br>
    <b>User manual:</b> \ref persistent_cohomology - <b>Reference manual:</b> Gudhi::persistent_cohomology::Persistent_cohomology
    </td>
 </tr>
</table>
*/

/*! \page installation Gudhi installation
 *  \tableofcontents
 * As Gudhi is a header only library, there is no need to install the library.
 * 
 * Examples of Gudhi headers inclusion can be found in \ref demos.
 * 
 * \section compiling Compiling
 * The library uses c++11 and requires <a target="_blank" href="http://www.boost.org/">Boost</a> with version 1.48.0 or
 * more recent. It is a multi-platform library and compiles on Linux, Mac OSX and Visual Studio 2015.
 * 
 * \subsection demos Demos and examples
 * To build the demos and examples, run the following commands in a terminal:
\verbatim  cd /path-to-gudhi/
mkdir build
cd build/
cmake ..
make \endverbatim
 * A list of examples is available <a href="examples.html">here</a>.
 * 
 * \subsection testsuites Test suites
 * To test your build, run the following command in a terminal:
 * \verbatim  make test \endverbatim
 * 
 * \section optionallibrary Optional third-party library
 * \subsection gmp GMP:
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
 * \subsection cgal CGAL:
 * The \ref alpha_complex data structure and few examples requires CGAL, which is a C++ library which provides easy
 * access to efficient and reliable geometric algorithms.
 * 
 * Having CGAL version 4.4 or higher installed is recommended. The procedure to install this library according to
 * your operating system is detailed here http://doc.cgal.org/latest/Manual/installation.html
 * 
 * The following examples require the <a target="_blank" href="http://www.cgal.org/">Computational Geometry Algorithms
 * Library</a> (CGAL \cite cgal:eb-15b) and will not be built if CGAL is not installed:
 * \li <a href="_persistent_cohomology_2alpha_complex_3d_persistence_8cpp-example.html">
 * Persistent_cohomology/alpha_complex_3d_persistence.cpp</a>
 * \li <a href="_simplex_tree_2example_alpha_shapes_3_simplex_tree_from_off_file_8cpp-example.html">
 * Simplex_tree/example_alpha_shapes_3_simplex_tree_from_off_file.cpp</a>
 * 
 * The following example requires CGAL version &ge; 4.6:
 * \li <a href="_witness_complex_2witness_complex_sphere_8cpp-example.html">
 * Witness_complex/witness_complex_sphere.cpp</a>
 * 
 * The following example requires CGAL version &ge; 4.7:
 * \li <a href="_alpha_complex_2_alpha_complex_from_off_8cpp-example.html">
 * Alpha_complex/Alpha_complex_from_off.cpp</a>
 * \li <a href="_alpha_complex_2_alpha_complex_from_points_8cpp-example.html">
 * Alpha_complex/Alpha_complex_from_points.cpp</a>
 * \li <a href="_persistent_cohomology_2alpha_complex_persistence_8cpp-example.html">
 * Persistent_cohomology/alpha_complex_persistence.cpp</a>
 * \li <a href="_persistent_cohomology_2periodic_alpha_complex_3d_persistence_8cpp-example.html">
 * Persistent_cohomology/periodic_alpha_complex_3d_persistence.cpp</a>
 * \li <a href="_persistent_cohomology_2custom_persistence_sort_8cpp-example.html">
 * Persistent_cohomology/custom_persistence_sort.cpp</a>
 * 
 * \subsection eigen3 Eigen3:
 * The \ref alpha_complex data structure and few examples requires
 * <a target="_blank" href="http://eigen.tuxfamily.org/">Eigen3</a> is a C++ template library for linear algebra:
 * matrices, vectors, numerical solvers, and related algorithms.
 * 
 * The following example requires the <a target="_blank" href="http://eigen.tuxfamily.org/">Eigen3</a> and will not be
 * built if Eigen3 is not installed:
 * \li <a href="_alpha_complex_2_alpha_complex_from_off_8cpp-example.html">
 * Alpha_complex/Alpha_complex_from_off.cpp</a> (requires also Eigen3)
 * \li <a href="_alpha_complex_2_alpha_complex_from_points_8cpp-example.html">
 * Alpha_complex/Alpha_complex_from_points.cpp</a> (requires also Eigen3)
 * \li <a href="_persistent_cohomology_2alpha_complex_persistence_8cpp-example.html">
 * Persistent_cohomology/alpha_complex_persistence.cpp</a>
 * \li <a href="_persistent_cohomology_2periodic_alpha_complex_3d_persistence_8cpp-example.html">
 * Persistent_cohomology/periodic_alpha_complex_3d_persistence.cpp</a>
 * \li <a href="_persistent_cohomology_2custom_persistence_sort_8cpp-example.html">
 * Persistent_cohomology/custom_persistence_sort.cpp</a>
 * 
 * \subsection tbb Threading Building Blocks:
 * <a target="_blank" href="https://www.threadingbuildingblocks.org/">Intel&reg; TBB</a> lets you easily write parallel
 * C++ programs that take full advantage of multicore performance, that are portable and composable, and that have
 * future-proof scalability.
 * 
 * Having Intel&reg; TBB installed is recommended to parallelize and accelerate some GUDHI computations.
 * 
 * The following examples are using Intel&reg; TBB if installed:
 * \li <a href="_alpha_complex_2_alpha_complex_from_off_8cpp-example.html">
 * Alpha_complex/Alpha_complex_from_off.cpp</a>
 * \li <a href="_alpha_complex_2_alpha_complex_from_points_8cpp-example.html">
 * Alpha_complex/Alpha_complex_from_points.cpp</a>
 * \li <a href="_bitmap_cubical_complex_2_bitmap_cubical_complex_8cpp-example.html">
 * Bitmap_cubical_complex/Bitmap_cubical_complex.cpp</a>
 * \li <a href="_bitmap_cubical_complex_2_bitmap_cubical_complex_periodic_boundary_conditions_8cpp-example.html">
 * Bitmap_cubical_complex/Bitmap_cubical_complex_periodic_boundary_conditions.cpp</a>
 * \li <a href="_bitmap_cubical_complex_2_random_bitmap_cubical_complex_8cpp-example.html">
 * Bitmap_cubical_complex/Random_bitmap_cubical_complex.cpp</a>
 * \li <a href="_persistent_cohomology_2alpha_complex_3d_persistence_8cpp-example.html">
 * Persistent_cohomology/alpha_complex_3d_persistence.cpp</a>
 * \li <a href="_persistent_cohomology_2alpha_complex_persistence_8cpp-example.html">
 * Persistent_cohomology/alpha_complex_persistence.cpp</a>
 * \li <a href="_simplex_tree_2simple_simplex_tree_8cpp-example.html">
 * Simplex_tree/simple_simplex_tree.cpp</a>
 * \li <a href="_simplex_tree_2example_alpha_shapes_3_simplex_tree_from_off_file_8cpp-example.html">
 * Simplex_tree/example_alpha_shapes_3_simplex_tree_from_off_file.cpp</a>
 * \li <a href="_simplex_tree_2simplex_tree_from_cliques_of_graph_8cpp-example.html">
 * Simplex_tree/simplex_tree_from_cliques_of_graph.cpp</a>
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
 * \li <a href="_persistent_cohomology_2rips_persistence_8cpp-example.html">
 * Persistent_cohomology/rips_persistence.cpp</a>
 * \li <a href="_persistent_cohomology_2periodic_alpha_complex_3d_persistence_8cpp-example.html">
 * Persistent_cohomology/periodic_alpha_complex_3d_persistence.cpp</a>
 * \li <a href="_persistent_cohomology_2custom_persistence_sort_8cpp-example.html">
 * Persistent_cohomology/custom_persistence_sort.cpp</a>
 * \li <a href="_rips_complex_2example_one_skeleton_rips_from_points_8cpp-example.html">
 * Rips_complex/example_one_skeleton_rips_from_points.cpp</a>
 * \li <a href="_rips_complex_2example_rips_complex_from_off_file_8cpp-example.html">
 * Rips_complex/example_rips_complex_from_off_file.cpp</a>
 * 
 * \section Contributions Bug reports and contributions
 * Please help us improving the quality of the GUDHI library. You may report bugs or suggestions to:
 * \verbatim  Contact: gudhi-users@lists.gforge.inria.fr \endverbatim
 * 
 * Gudhi is open to external contributions. If you want to join our development team, please contact us.
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

// List of Gudhi examples - Doxygen needs at least a file tag to analyse comments
/*! @file Examples
 * @example Alpha_complex/Alpha_complex_from_off.cpp
 * @example Alpha_complex/Alpha_complex_from_points.cpp
 * @example Bitmap_cubical_complex/Bitmap_cubical_complex.cpp
 * @example Bitmap_cubical_complex/Bitmap_cubical_complex_periodic_boundary_conditions.cpp
 * @example Bitmap_cubical_complex/Random_bitmap_cubical_complex.cpp
 * @example common/example_CGAL_3D_points_off_reader.cpp
 * @example common/example_CGAL_points_off_reader.cpp
 * @example Contraction/Garland_heckbert.cpp
 * @example Contraction/Rips_contraction.cpp
 * @example Persistent_cohomology/alpha_complex_3d_persistence.cpp
 * @example Persistent_cohomology/alpha_complex_persistence.cpp
 * @example Persistent_cohomology/rips_persistence_via_boundary_matrix.cpp
 * @example Persistent_cohomology/periodic_alpha_complex_3d_persistence.cpp
 * @example Persistent_cohomology/persistence_from_file.cpp
 * @example Persistent_cohomology/persistence_from_simple_simplex_tree.cpp
 * @example Persistent_cohomology/plain_homology.cpp
 * @example Persistent_cohomology/rips_multifield_persistence.cpp
 * @example Persistent_cohomology/rips_persistence.cpp
 * @example Persistent_cohomology/custom_persistence_sort.cpp
 * @example Rips_complex/example_one_skeleton_rips_from_points.cpp
 * @example Rips_complex/example_rips_complex_from_off_file.cpp
 * @example Simplex_tree/mini_simplex_tree.cpp
 * @example Simplex_tree/simple_simplex_tree.cpp
 * @example Simplex_tree/example_alpha_shapes_3_simplex_tree_from_off_file.cpp
 * @example Simplex_tree/simplex_tree_from_cliques_of_graph.cpp
 * @example Skeleton_blocker/Skeleton_blocker_from_simplices.cpp
 * @example Skeleton_blocker/Skeleton_blocker_iteration.cpp
 * @example Skeleton_blocker/Skeleton_blocker_link.cpp
 * @example Tangential_complex/example_basic.cpp
 * @example Tangential_complex/example_with_perturb.cpp
 * @example Witness_complex/witness_complex_from_file.cpp
 * @example Witness_complex/witness_complex_sphere.cpp
 */


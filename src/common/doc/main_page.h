/*! \mainpage The C++ library
 *  \tableofcontents
 * \image html "Gudhi_banner.png" "" width=20cm
 * 
 * \section Introduction Introduction
 * The GUDHI library (Geometry Understanding in Higher Dimensions) is a generic open source
 * <a class="el" target="_blank" href="http://gudhi.gforge.inria.fr/doc/latest/">C++ library</a> for
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
 * \li Algorithms to compute persistent homology and bottleneck distance.
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
      <b>Requires:</b> \ref eigen3 and<br>
      \ref cgal &ge; 4.7.0 for Alpha_complex<br>
      \ref cgal &ge; 4.11.0 for Alpha_complex_3d
    </td>
    <td width="75%">
    Alpha_complex is a simplicial complex constructed from the finite cells of a Delaunay Triangulation.<br>
    The filtration value of each simplex is computed as the square of the circumradius of the simplex if the
    circumsphere is empty (the simplex is then said to be Gabriel), and as the minimum of the filtration
    values of the codimension 1 cofaces that make it not Gabriel otherwise.
    All simplices that have a filtration value strictly greater than a given alpha squared value are not inserted into
    the complex.<br>
    <b>User manual:</b> \ref alpha_complex - <b>Reference manual:</b> Gudhi::alpha_complex::Alpha_complex and
    Gudhi::alpha_complex::Alpha_complex_3d
    </td>
 </tr>
</table>
 \subsection CechComplexDataStructure Čech complex
 \image html "cech_complex_representation.png" "Čech complex representation"
<table border="0">
  <tr>
    <td width="25%">
      <b>Author:</b> Vincent Rouvreau<br>
      <b>Introduced in:</b> GUDHI 2.2.0<br>
      <b>Copyright:</b> GPL v3<br>
    </td>
    <td width="75%">
    The Čech complex is a simplicial complex constructed from a proximity graph.<br>
    The set of all simplices is filtered by the radius of their minimal enclosing ball.<br>
    <b>User manual:</b> \ref cech_complex - <b>Reference manual:</b> Gudhi::cech_complex::Cech_complex
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
</table>
 \subsection RipsComplexDataStructure Rips complex
 \image html "rips_complex_representation.png" "Rips complex representation"
<table border="0">
  <tr>
    <td width="25%">
      <b>Author:</b> Cl&eacute;ment Maria, Pawel Dlotko, Vincent Rouvreau, Marc Glisse<br>
      <b>Introduced in:</b> GUDHI 2.0.0<br>
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
 \subsection CoverComplexDataStructure Cover Complexes
 \image html "gicvisu.jpg" "Graph Induced Complex of a point cloud."
<table border="0">
  <tr>
    <td width="25%">
      <b>Author:</b> Mathieu Carri&egrave;re<br>
      <b>Introduced in:</b> GUDHI 2.1.0<br>
      <b>Copyright:</b> GPL v3<br>
      <b>Requires:</b> \ref cgal &ge; 4.8.1
    </td>
    <td width="75%">
    Nerves and Graph Induced Complexes are cover complexes, i.e. simplicial complexes that provably contain
    topological information about the input data. They can be computed with a cover of the
    data, that comes i.e. from the preimage of a family of intervals covering the image
    of a scalar-valued function defined on the data. <br>
    <b>User manual:</b> \ref cover_complex - <b>Reference manual:</b> Gudhi::cover_complex::Cover_complex
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
 \subsection StrongCollapseDataStructure Strong collapse
 \image html "strong_collapse_representation.png" "Strong collapse representation"
<table border="0">
  <tr>
    <td width="25%">
      <b>Author:</b> Siddharth Pritam<br>
      <b>Introduced in:</b> GUDHI 2.4.0<br>
      <b>Copyright:</b> GPL v3<br>
    </td>
    <td width="75%">
    An elementary strong collapse of a simplicial complex $K$ is removal of a vertex $v$ whose link $lk_K(v)$ is a simplicial cone $v^\prime*L$.
    A series of elementary strong collapses is called a strong collapse. Strong collapse is a special type of classical simple collapse. 
    To compute strong collapse of a general simplicial complex we only need a data-structure which consists of maximal simplices and vertices. 
    And for a flag complex, strong collapse can be computed using only the 1-skeleton of the flag complex. 
    Moreover unlike simple collapse, a strong collapse can be associated to a simplical map from the complex $K$ to an unique core $K^0$. 
	
	Computation of strong collapse and persistent homology of a filtered flag complex via strong collapse as described in \cite strongcollapsesocg2019.
     <br>
    <b>User manual:</b> \ref strong_collapse - <b>Reference manual:</b> Gudhi::strong_collapse::Flag_complex_strong_collapse
    </td>
 </tr>
</table>
 \subsection TangentialComplexDataStructure Tangential complex
 \image html "tc_examples.png" "Tangential complex representation"
<table border="0">
  <tr>
    <td width="25%">
      <b>Author:</b> Cl&eacute;ment Jamin<br>
      <b>Introduced in:</b> GUDHI 2.0.0<br>
      <b>Copyright:</b> GPL v3<br>
      <b>Requires:</b> \ref cgal &ge; 4.8.1 and \ref eigen3
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
      <b>Euclidean version requires:</b> \ref cgal &ge; 4.6.0 and \ref eigen3
    </td>
    <td width="75%">
    Witness complex \f$ Wit(W,L) \f$  is a simplicial complex defined on two sets of points in \f$\mathbb{R}^D\f$.
    The data structure is described in \cite boissonnatmariasimplextreealgorithmica .<br>
    <b>User manual:</b> \ref witness_complex - <b>Reference manual:</b> Gudhi::witness_complex::SimplicialComplexForWitness
    </td>
 </tr>
</table>
 
 \section Toolbox Toolbox
 
 \subsection BottleneckDistanceToolbox Bottleneck distance
 \image html "perturb_pd.png" "Bottleneck distance is the length of the longest edge"
<table border="0">
  <tr>
    <td width="25%">
      <b>Author:</b> Fran&ccedil;ois Godi<br>
      <b>Introduced in:</b> GUDHI 2.0.0<br>
      <b>Copyright:</b> GPL v3<br>
      <b>Requires:</b> \ref cgal &ge; 4.8.1
    </td>
    <td width="75%">
    Bottleneck distance measures the similarity between two persistence diagrams. 
    It's the shortest distance b for which there exists a perfect matching between 
    the points of the two diagrams (+ all the diagonal points) such that 
    any couple of matched points are at distance at most b.
    <br>
    <b>User manual:</b> \ref bottleneck_distance
    </td>
 </tr>
</table>
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
 \subsection PersistenceRepresentationsToolbox Persistence representations
 \image html "average_landscape.png" "Persistence representations"
<table border="0">
  <tr>
    <td width="25%">
      <b>Author:</b> Pawel Dlotko<br>
      <b>Introduced in:</b> GUDHI 2.1.0<br>
      <b>Copyright:</b> GPL v3<br>      
    </td>
    <td width="75%">
    It contains implementation of various representations of persistence diagrams; diagrams themselves, persistence
    landscapes (rigorous and grid version), persistence heath maps, vectors and others. It implements basic
    functionalities which are neccessary to use persistence in statistics and machine learning.<br>
    <b>User manual:</b> \ref Persistence_representations
    </td>
 </tr>
</table>
 
*/

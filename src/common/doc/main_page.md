[TOC]

# The C++ library {#main_page}
\image html "Gudhi_banner.png"
<br><br><br><br>

## Complexes {#Complexes}
### Cubical complex
 
<table>
  <tr>
    <td width="35%" rowspan=2>
      \image html "Cubical_complex_representation.png"
    </td>
    <td width="50%">
    The cubical complex is an example of a structured complex useful in computational mathematics (specially
    rigorous numerics) and image analysis.<br>
    </td>
    <td width="15%">
      <b>Author:</b> Pawel Dlotko<br>
      <b>Introduced in:</b> GUDHI 1.3.0<br>
      <b>Copyright:</b> MIT<br>
    </td>
 </tr>
 <tr>
    <td colspan=2 height="25">
    <b>User manual:</b> \ref cubical_complex
    </td>
 </tr>
</table>

### Simplicial complex

#### Alpha complex

<table>
  <tr>
    <td width="35%" rowspan=2>
      \image html "alpha_complex_representation.png"
    </td>
    <td width="50%">
    Alpha complex is a simplicial complex constructed from the finite cells of a Delaunay Triangulation.<br>
    The filtration value of each simplex is computed as the square of the circumradius of the simplex if the
    circumsphere is empty (the simplex is then said to be Gabriel), and as the minimum of the filtration
    values of the codimension 1 cofaces that make it not Gabriel otherwise.
    All simplices that have a filtration value strictly greater than a given alpha squared value are not inserted into
    the complex.<br>
    </td>
    <td width="15%">
      <b>Author:</b> Vincent Rouvreau<br>
      <b>Introduced in:</b> GUDHI 1.3.0<br>
      <b>Copyright:</b> MIT [(GPL v3)](../../licensing/)<br>
      <b>Requires:</b> \ref eigen &ge; 3.1.0 and \ref cgal &ge; 4.11.0
    </td>
 </tr>
 <tr>
    <td colspan=2 height="25">
    <b>User manual:</b> \ref alpha_complex
    </td>
 </tr>
</table>

#### Čech complex

<table>
 <tr>
    <td width="35%" rowspan=2>
      \image html "cech_complex_representation.png"
    </td>
    <td width="50%">
    The Čech complex is a simplicial complex constructed from a proximity graph.
    The set of all simplices is filtered by the radius of their minimal enclosing ball.
    </td>
    <td width="15%">
      <b>Author:</b> Vincent Rouvreau<br>
      <b>Introduced in:</b> GUDHI 2.2.0<br>
      <b>Copyright:</b> MIT [(GPL v3)](../../licensing/)<br>
      <b>Includes:</b> [Miniball](https://people.inf.ethz.ch/gaertner/subdir/software/miniball.html)<br>
    </td>
 </tr>
 <tr>
    <td colspan=2 height="25">
    <b>User manual:</b> \ref cech_complex
    </td>
 </tr>
</table>

#### Rips complex

<table>
  <tr>
    <td width="35%" rowspan=2>
      \image html "rips_complex_representation.png"
    </td>
    <td width="50%">
    Rips complex is a simplicial complex constructed from a one skeleton graph.<br>
    The filtration value of each edge is computed from a user-given distance function and is inserted until a
    user-given threshold value.<br>
    This complex can be built from a point cloud and a distance function, or from a distance matrix.
    </td>
    <td width="15%">
      <b>Author:</b> Cl&eacute;ment Maria, Pawel Dlotko, Vincent Rouvreau, Marc Glisse<br>
      <b>Introduced in:</b> GUDHI 2.0.0<br>
      <b>Copyright:</b> MIT<br>
    </td>
 </tr>
 <tr>
    <td colspan=2 height="25">
    <b>User manual:</b> \ref rips_complex
    </td>
 </tr>
</table>

#### Witness complex

<table>
  <tr>
    <td width="35%" rowspan=2>
      \image html "Witness_complex_representation.png"
    </td>
    <td width="50%">
    Witness complex \f$ Wit(W,L) \f$  is a simplicial complex defined on two sets of points in \f$\mathbb{R}^D\f$.
    The data structure is described in \cite boissonnatmariasimplextreealgorithmica .
    </td>
    <td width="15%">
      <b>Author:</b> Siargey Kachanovich<br>
      <b>Introduced in:</b> GUDHI 1.3.0<br>
      <b>Copyright:</b> MIT ([GPL v3](../../licensing/) for Euclidean version)<br>
      <b>Euclidean version requires:</b> \ref eigen &ge; 3.1.0 and \ref cgal &ge; 4.11.0
    </td>
 </tr>
 <tr>
    <td colspan=2 height="25">
    <b>User manual:</b> \ref witness_complex
    </td>
 </tr>
</table>

### Cover Complexes
<table>
  <tr>
    <td width="35%" rowspan=2>
      \image html "gicvisu.jpg"
    </td>
    <td width="50%">
    Nerves and Graph Induced Complexes are cover complexes, i.e. simplicial complexes that provably contain
    topological information about the input data. They can be computed with a cover of the
    data, that comes i.e. from the preimage of a family of intervals covering the image
    of a scalar-valued function defined on the data. <br>
    </td>
    <td width="15%">
      <b>Author:</b> Mathieu Carri&egrave;re<br>
      <b>Introduced in:</b> GUDHI 2.1.0<br>
      <b>Copyright:</b> MIT [(GPL v3)](../../licensing/)<br>
      <b>Requires:</b> \ref cgal &ge; 4.11.0
    </td>
 </tr>
 <tr>
    <td colspan=2 height="25">
    <b>User manual:</b> \ref cover_complex
    </td>
 </tr>
</table>

## Data structures and basic operations {#DataStructuresAndBasicOperations}

### Data structures

#### Simplex tree
<table>
  <tr>
    <td width="35%" rowspan=2>
      \image html "Simplex_tree_representation.png"
    </td>
    <td width="50%">
    The simplex tree is an efficient and flexible
 data structure for representing general (filtered) simplicial complexes. The data structure
 is described in \cite boissonnatmariasimplextreealgorithmica .
    </td>
    <td width="15%">
      <b>Author:</b> Cl&eacute;ment Maria<br>
      <b>Introduced in:</b> GUDHI 1.0.0<br>
      <b>Copyright:</b> MIT<br>
    </td>
 </tr>
 <tr>
    <td colspan=2 height="25">
    <b>User manual:</b> \ref simplex_tree
    </td>
 </tr>
</table>

#### Skeleton blocker

<table>
  <tr>
    <td width="35%" rowspan=2>
      \image html "ds_representation.png"
    </td>
    <td width="50%">
    The Skeleton-Blocker data-structure proposes a light encoding for simplicial complexes by storing only an *implicit*
    representation of its simplices \cite socg_blockers_2011,\cite blockers2012. Intuitively, it just stores the
    1-skeleton of a simplicial complex with a graph and the set of its "missing faces" that is very small in practice.
    This data-structure handles all simplicial complexes operations such as simplex enumeration or simplex removal but
    operations that are particularly efficient are operations that do not require simplex enumeration such as edge
    iteration, link computation or simplex contraction.
    </td>
    <td width="15%">
      <b>Author:</b> David Salinas<br>
      <b>Introduced in:</b> GUDHI 1.1.0<br>
      <b>Copyright:</b> MIT<br>
    </td>
 </tr>
 <tr>
    <td colspan=2 height="25">
    <b>User manual:</b> \ref skbl
    </td>
 </tr>
</table>

#### Toplex Map

<table>
  <tr>
    <td width="35%" rowspan=2>
      \image html "map.png"
    </td>
    <td width="50%">
    The Toplex map data structure is composed firstly of a raw storage of toplices (the maximal simplices)
    and secondly of a map which associate any vertex to a set of pointers toward all toplices
    containing this vertex.
    </td>
    <td width="15%">
      <b>Author:</b> Fran&ccedil;ois Godi<br>
      <b>Introduced in:</b> GUDHI 2.1.0<br>
      <b>Copyright:</b> MIT<br>
    </td>
 </tr>
 <tr>
    <td colspan=2 height="25">
    <b>User manual:</b> \ref toplex_map
    </td>
 </tr>
</table>

### Basic operations

#### Contraction

<table>
  <tr>
    <td width="35%" rowspan=2>
      \image html "sphere_contraction_representation.png"
    </td>
    <td width="50%">
    The purpose of this package is to offer a user-friendly interface for edge contraction simplification of huge
    simplicial complexes. It uses the \ref skbl data-structure whose size remains small  during simplification of most
    used geometrical complexes of topological data analysis such as the Rips or the Delaunay complexes. In practice,
    the size of this data-structure is even much lower than the total number of simplices.
    </td>
    <td width="15%">
      <b>Author:</b> David Salinas<br>
      <b>Introduced in:</b> GUDHI 1.1.0<br>
      <b>Copyright:</b> MIT [(LGPL v3)](../../licensing/)<br>
      <b>Requires:</b> \ref cgal &ge; 4.11.0
    </td>
 </tr>
 <tr>
    <td colspan=2 height="25">
    <b>User manual:</b> \ref contr
    </td>
 </tr>
</table>

## Topological descriptors computation {#TopologicalDescriptorsComputation}

### Persistent Cohomology

<table>
  <tr>
    <td width="35%" rowspan=2>
      \image html "3DTorus_poch.png"
    </td>
    <td width="50%">
    The theory of homology consists in attaching to a topological space a sequence of (homology) groups, capturing
    global topological features like connected components, holes, cavities, etc. Persistent homology studies the
    evolution -- birth, life and death -- of these features when the topological space is changing. Consequently, the
    theory is essentially composed of three elements: topological spaces, their homology groups and an evolution
    scheme.
    Computation of persistent cohomology using the algorithm of \cite DBLP:journals/dcg/SilvaMV11 and
    \cite DBLP:journals/corr/abs-1208-5018 and the Compressed Annotation Matrix implementation of
    \cite DBLP:conf/esa/BoissonnatDM13 .
    </td>
    <td width="15%">
      <b>Author:</b> Cl&eacute;ment Maria<br>
      <b>Introduced in:</b> GUDHI 1.0.0<br>
      <b>Copyright:</b> MIT<br>
    </td>
 </tr>
 <tr>
    <td colspan=2 height="25">
    <b>User manual:</b> \ref persistent_cohomology
    </td>
 </tr>
</table>

## Manifold reconstruction {#ManifoldReconstruction}

### Tangential complex

<table>
  <tr>
    <td width="35%" rowspan=2>
      \image html "tc_examples.png"
    </td>
    <td width="50%">
    A Tangential Delaunay complex is a <a target="_blank" href="https://en.wikipedia.org/wiki/Simplicial_complex">simplicial complex</a>
    designed to reconstruct a \f$ k \f$-dimensional manifold embedded in \f$ d \f$-dimensional Euclidean space. 
    The input is a point sample coming from an unknown manifold.
    The running time depends only linearly on the extrinsic dimension \f$ d \f$
    and exponentially on the intrinsic dimension \f$ k \f$.
    </td>
    <td width="15%">
      <b>Author:</b> Cl&eacute;ment Jamin<br>
      <b>Introduced in:</b> GUDHI 2.0.0<br>
      <b>Copyright:</b> MIT [(GPL v3)](../../licensing/)<br>
      <b>Requires:</b> \ref eigen &ge; 3.1.0 and \ref cgal &ge; 4.11.0
    </td>
 </tr>
 <tr>
    <td colspan=2 height="25">
    <b>User manual:</b> \ref tangential_complex
    </td>
 </tr>
</table>

## Topological descriptors tools {#TopologicalDescriptorsTools}

### Bottleneck distance

<table>
  <tr>
    <td width="35%" rowspan=2>
      \image html "perturb_pd.png"
    </td>
    <td width="50%">
    Bottleneck distance measures the similarity between two persistence diagrams. 
    It's the shortest distance b for which there exists a perfect matching between 
    the points of the two diagrams (+ all the diagonal points) such that 
    any couple of matched points are at distance at most b,
    where the distance between points is the sup norm in \f$\mathbb{R}^2\f$
    (not the Euclidean distance).
    </td>
    <td width="15%">
      <b>Author:</b> Fran&ccedil;ois Godi<br>
      <b>Introduced in:</b> GUDHI 2.0.0<br>
      <b>Copyright:</b> MIT [(GPL v3)](../../licensing/)<br>
      <b>Requires:</b> \ref cgal &ge; 4.11.0
    </td>
 </tr>
 <tr>
    <td colspan=2 height="25">
    <b>User manual:</b> \ref bottleneck_distance
    </td>
 </tr>
</table>

### Persistence representations

<table>
  <tr>
    <td width="35%" rowspan=2>
      \image html "average_landscape.png"
    </td>
    <td width="50%">
    It contains implementation of various representations of persistence diagrams; diagrams themselves, persistence
    landscapes (rigorous and grid version), persistence heat maps, vectors and others. It implements basic
    functionalities which are necessary to use persistence in statistics and machine learning.
    </td>
    <td width="15%">
      <b>Author:</b> Pawel Dlotko<br>
      <b>Introduced in:</b> GUDHI 2.1.0<br>
      <b>Copyright:</b> MIT<br>      
    </td>
 </tr>
 <tr>
    <td colspan=2 height="25">
    <b>User manual:</b> \ref Persistence_representations
    </td>
 </tr>
</table>

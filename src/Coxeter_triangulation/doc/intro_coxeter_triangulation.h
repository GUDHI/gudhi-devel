/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef DOC_COXETER_TRIANGULATION_INTRO_COXETER_TRIANGULATION_H_
#define DOC_COXETER_TRIANGULATION_INTRO_COXETER_TRIANGULATION_H_

// needs namespaces for Doxygen to link on classes
namespace Gudhi {
namespace coxeter_triangulation {

/**  \defgroup coxeter_triangulation Coxeter triangulation

\author    Siargey Kachanovich

@{

\section overview Module overview

Coxeter triangulation module is designed to provide tools for constructing a piecewise-linear approximation of an
\f$m\f$-dimensional smooth manifold embedded in \f$ \mathbb{R}^d \f$ using an ambient triangulation.
For a more detailed description of the module see \cite KachanovichThesis.

\section manifoldtracing Manifold tracing algorithm
The central piece of the module is the manifold tracing algorithm represented by the class
\ref Gudhi::coxeter_triangulation::Manifold_tracing "Manifold_tracing".
The manifold tracing algorithm takes as input a manifold of some dimension \f$m\f$ embedded in \f$\mathbb{R}^d\f$
represented by an intersection oracle (see Section \ref intersectionoracle "Intersection oracle"), a point on the
manifold and an ambient triangulation (see Section \ref ambienttriangulations "Ambient triangulations").
The output consists of one map (or two maps in the case of manifolds with boundary) from the \f$(d-m)\f$-dimensional
(and \f$(d-m+1)\f$-dimensional in the case of manifolds with boundary) simplices in the ambient triangulation that
intersect the manifold to their intersection points.
From this output, it is possible to construct the cell complex of the piecewise-linear approximation of the input
manifold.

There are two methods that execute the manifold tracing algorithm: the method
\ref Gudhi::coxeter_triangulation::Manifold_tracing::manifold_tracing_algorithm() "Manifold_tracing::manifold_tracing_algorithm(seed_points, triangulation, oracle, out_simplex_map)"
for manifolds without boundary and
\ref Gudhi::coxeter_triangulation::Manifold_tracing::manifold_tracing_algorithm() "Manifold_tracing::manifold_tracing_algorithm(seed_points, triangulation, oracle, interior_simplex_map,boundary_simplex_map)"
for manifolds with boundary. The algorithm functions as follows. It starts at the specified seed points and inserts a
\f$(d-m)\f$-dimensional simplices nearby each seed point that intersect the manifold into the output. Starting from
this simplex, the algorithm propagates the search for other \f$(d-m)\f$-dimensional simplices that intersect the
manifold by marching from a simplex to neighbouring simplices via their common cofaces.

This class \ref Gudhi::coxeter_triangulation::Manifold_tracing "Manifold_tracing" has one template parameter
`Triangulation_` which specifies the ambient triangulation which is used by the algorithm.
The template type `Triangulation_` has to be a model of the concept
\ref Gudhi::coxeter_triangulation::TriangulationForManifoldTracing "TriangulationForManifoldTracing".

The module also provides two static methods:
\ref Gudhi::coxeter_triangulation::manifold_tracing_algorithm() "manifold_tracing_algorithm(seed_points, triangulation, oracle, out_simplex_map)"
for manifolds without boundary and
\ref manifold_tracing_algorithm() "manifold_tracing_algorithm(seed_points, triangulation, oracle, interior_simplex_map, boundary_simplex_map)"
for manifolds with boundary. For these static methods it is not necessary to specify any template arguments.

\section ambienttriangulations Ambient triangulations

The ambient triangulations supported by the manifold tracing algorithm have to be models of the concept
\ref Gudhi::coxeter_triangulation::TriangulationForManifoldTracing "TriangulationForManifoldTracing".
This module offers two such models: the class
\ref Gudhi::coxeter_triangulation::Freudenthal_triangulation "Freudenthal_triangulation" and the derived class
\ref Gudhi::coxeter_triangulation::Coxeter_triangulation "Coxeter_triangulation".

Both these classes encode affine transformations of the so-called Freudenthal-Kuhn triangulation of \f$\mathbb{R}^d\f$.
The Freudenthal-Kuhn triangulation of \f$\mathbb{R}^d\f$ is defined as the simplicial subdivision of the unit cubic
partition of \f$\mathbb{R}^d\f$.
Each simplex is encoded using the permutahedral representation, which consists of an integer-valued vector \f$y\f$ that
positions the simplex in a specific cube in the cubical partition and an ordered partition \f$\omega\f$ of the set
\f$\{1,\ldots,d+1\}\f$, which positions the simplex in the simplicial subdivision of the cube.
The default constructor
\ref Gudhi::coxeter_triangulation::Freudenthal_triangulation::Freudenthal_triangulation(std::size_t)
"Freudenthal_triangulation(d)" the Freudenthal-Kuhn triangulation of \f$\mathbb{R}^d\f$. The class
\ref Gudhi::coxeter_triangulation::Freudenthal_triangulation "Freudenthal_triangulation" can also encode any affine
transformation of the Freudenthal-Kuhn triangulation of \f$\mathbb{R}^d\f$ using an invertible matrix \f$\Lambda\f$ and
an offset vector \f$b\f$ that can be specified in the constructor and which can be changed using the methods
change_matrix and change_offset. The class
\ref Gudhi::coxeter_triangulation::Coxeter_triangulation "Coxeter_triangulation" is derived from
\ref Gudhi::coxeter_triangulation::Freudenthal_triangulation "Freudenthal_triangulation" and its default constructor
\ref Gudhi::coxeter_triangulation::Coxeter_triangulation::Coxeter_triangulation(std::size_t) "Coxeter_triangulation(d)"
builds a Coxeter triangulation of type \f$\tilde{A}_d\f$, which has the best simplex quality of all linear
transformations of the Freudenthal-Kuhn triangulation of \f$\mathbb{R}^d\f$.

\image html two_triangulations.png "Coxeter (on the left) and Freudenthal-Kuhn triangulation (on the right)"


\section intersectionoracle Intersection oracle

The input \f$m\f$-dimensional manifold in \f$\mathbb{R}^d\f$ needs to be given via the intersection oracle that answers
the following query: given a \f$(d-m)\f$-dimensional simplex, does it intersect the manifold?
The concept \ref Gudhi::coxeter_triangulation::IntersectionOracle "IntersectionOracle" describes all requirements for
an intersection oracle class to be compatible with the class
\ref Gudhi::coxeter_triangulation::Manifold_tracing "Manifold_tracing".
This module offers one model of the concept
\ref Gudhi::coxeter_triangulation::IntersectionOracle "IntersectionOracle", which is the class
\ref Gudhi::coxeter_triangulation::Implicit_manifold_intersection_oracle "Implicit_manifold_intersection_oracle".
This class represents a manifold given as the zero-set of a specified function
\f$F: \mathbb{R}^d \rightarrow \mathbb{R}^{d-m}\f$.
The function \f$F\f$ is given by a class which is a model of the concept
\ref Gudhi::coxeter_triangulation::FunctionForImplicitManifold "FunctionForImplicitManifold".
There are multiple function classes that are already implemented in this module.

\li \ref Gudhi::coxeter_triangulation::Constant_function(std::size_t, std::size_t, Eigen::VectorXd)
"Constant_function(d,k,v)" defines a constant function \f$F\f$ such that for all \f$x \in \mathbb{R}^d\f$, we have
  \f$F(x) = v \in \mathbb{R}^k\f$.
  The class Constant_function does not define an implicit manifold, but is useful as the domain function when defining
  boundaryless implicit manifolds.
\li \ref Gudhi::coxeter_triangulation::Function_affine_plane_in_Rd(N,b) "Function_affine_plane_in_Rd(N,b)" defines an
  \f$m\f$-dimensional implicit affine plane in the \f$d\f$-dimensional Euclidean space given by a normal matrix \f$N\f$
  and an offset vector \f$b\f$.
\li \ref Gudhi::coxeter_triangulation::Function_Sm_in_Rd(r,m,d,center) "Function_Sm_in_Rd(r,m,d,center)" defines an
  \f$m\f$-dimensional implicit sphere embedded in the \f$d\f$-dimensional Euclidean space of radius \f$r\f$ centered at
  the point 'center'.
\li \ref Gudhi::coxeter_triangulation::Function_moment_curve_in_Rd(r,d) "Function_moment_curve(r,d)" defines the moment
  curve in the \f$d\f$-dimensional Euclidean space of radius \f$r\f$ given as the parameterized curve (but implemented
  as an implicit curve):
  \f[ (r, rt, \ldots, rt^{d-1}) \in \mathbb{R}^d,\text{ for $t \in \mathbb{R}$.} \f]
\li \ref Gudhi::coxeter_triangulation::Function_torus_in_R3(R, r) "Function_torus_in_R3(R, r)" defines a torus in
  \f$\mathbb{R}^3\f$ with the outer radius \f$R\f$ and the inner radius, given by the equation:
  \f[ z^2 + (\sqrt{x^2 + y^2} - r)^2 - R^2 = 0. \f]
\li \ref Gudhi::coxeter_triangulation::Function_chair_in_R3(a, b, k) "Function_chair_in_R3(a, b, k)" defines the
  \"Chair\" surface in \f$\mathbb{R}^3\f$ defined by the equation:
  \f[ (x^2 + y^2 + z^2 - ak^2)^2 - b((z-k)^2 - 2x^2)((z+k)^2 - 2y^2) = 0. \f]
\li \ref Gudhi::coxeter_triangulation::Function_iron_in_R3() "Function_iron_in_R3()" defines the \"Iron\" surface in
  \f$\mathbb{R}^3\f$ defined by the equation:
  \f[ \frac{-x^6-y^6-z^6}{300} + \frac{xy^2z}{2.1} + y^2 + (z-2)^2 = 1. \f]
\li \ref Gudhi::coxeter_triangulation::Function_lemniscate_revolution_in_R3(a) "Function_lemniscate_revolution_in_R3(a)"
  defines a revolution surface in \f$\mathbb{R}^3\f$ obtained from the lemniscate of Bernoulli defined by the equation:
  \f[ (x^2 + y^2 + z^2)^2 - 2a^2(x^2 - y^2 - z^2) = 0. \f]
\li \ref Gudhi::coxeter_triangulation::Function_whitney_umbrella_in_R3() "Function_whitney_umbrella_in_R3()" defines
  the Whitney umbrella surface in \f$\mathbb{R}^3\f$ defined by the equation:
  \f[ x^2 - y^2z = 0. \f]

The base function classes above can be composed or modified into new functions using the following classes and methods:

\li \ref Gudhi::coxeter_triangulation::Cartesian_product "Cartesian_product(functions...)" expresses the Cartesian
  product \f$F_1^{-1}(0) \times \ldots \times F_k^{-1}(0)\f$ of multiple implicit manifolds as an implicit manifold.
  For convenience, a static function
  \ref Gudhi::coxeter_triangulation::make_product_function() "make_product_function(functions...)" is provided that
  takes a pack of function-typed objects as the argument.
\li \ref Gudhi::coxeter_triangulation::Embed_in_Rd "Embed_in_Rd(F, d)" expresses an implicit manifold given as the
  zero-set of a function \f$F\f$ embedded in a higher-dimensional Euclidean space \f$\mathbb{R}^d\f$.
  For convenience, a static function \ref Gudhi::coxeter_triangulation::make_embedding() "make_embedding(F, d)" is
  provided.
\li \ref Gudhi::coxeter_triangulation::Linear_transformation "Linear_transformation(F, M)" applies a linear
  transformation given by a matrix \f$M\f$ on an implicit manifold given as the zero-set of the function \f$F\f$.
  For convenience, a static function
  \ref Gudhi::coxeter_triangulation::make_linear_transformation() "make_linear_transformation(F, M)" is provided.
\li \ref Gudhi::coxeter_triangulation::Translate "Translate(F, v)" translates an implicit manifold given as the
  zero-set of ththe function \f$F\f$ by a vector \f$v\f$.
  For convenience, a static function \ref Gudhi::coxeter_triangulation::translate() "translate(F, v)" is provided.
\li \ref Gudhi::coxeter_triangulation::Negation() "Negation(F)" defines the negative of the given function \f$F\f$.
  This class is useful to define the complementary of a given domain, when defining a manifold with boundary.
  For convenience, a static function \ref Gudhi::coxeter_triangulation::negation() "negation(F)" is provided.
\li \ref Gudhi::coxeter_triangulation::PL_approximation "PL_approximation(F, T)" defines a piecewise-linear
  approximation of a given function \f$F\f$ induced by an ambient triangulation \f$T\f$.
  The purpose of this class is to define a piecewise-linear function that is compatible with the requirements for the
  domain function \f$D\f$ when defining a manifold with boundary.
  For convenience, a static function
  \ref Gudhi::coxeter_triangulation::make_pl_approximation() "make_pl_approximation(F, T)" is provided.
  The type of \f$T\f$ is required to be a model of the concept
  \ref Gudhi::coxeter_triangulation::TriangulationForManifoldTracing "TriangulationForManifoldTracing".

\section cellcomplex Cell complex construction

The output of the manifold tracing algorithm can be transformed into the Hasse diagram of a cell complex that
approximates the input manifold using the class \ref Gudhi::coxeter_triangulation::Cell_complex "Cell_complex".
The type of the cells in the Hasse diagram is
\ref Gudhi::Hasse_diagram::Hasse_diagram_cell "Hasse_cell<int, double, bool>" provided by the module Hasse diagram.
The cells in the cell complex given by an object of the class
\ref Gudhi::coxeter_triangulation::Cell_complex "Cell_complex" are accessed through several maps that are accessed
through the following methods.

\li The method
\ref Gudhi::coxeter_triangulation::Cell_complex::interior_simplex_cell_maps() "interior_simplex_cell_maps()"
returns a vector of maps from the cells of various dimensions in the interior of the cell complex to the permutahedral
representations of the corresponding simplices in the ambient triangulation.
Each individual map for cells of a specific dimension \f$l\f$ can be accessed using the method
\ref Gudhi::coxeter_triangulation::Cell_complex::interior_simplex_cell_map() "interior_simplex_cell_map(l)".
\li The method
\ref Gudhi::coxeter_triangulation::Cell_complex::boundary_simplex_cell_maps() "boundary_simplex_cell_maps()"
returns a vector of maps from the cells of various dimensions on the boundary of the cell complex to the permutahedral
representations of the corresponding simplices in the ambient triangulation.
Each individual map for cells of a specific dimension \f$l\f$ can be accessed using the method
\ref Gudhi::coxeter_triangulation::Cell_complex::boundary_simplex_cell_map() "boundary_simplex_cell_map(l)".
\li The method \ref Gudhi::coxeter_triangulation::Cell_complex::cell_simplex_map() "cell_simplex_map()" returns a map
from the cells in the cell complex to the permutahedral representations of the corresponding simplices in the ambient
triangulation.
\li The method \ref Gudhi::coxeter_triangulation::Cell_complex::cell_point_map() "cell_point_map()" returns a map from
the vertex cells in the cell complex to their Cartesian coordinates.

\section example Examples

\subsection examplewithoutboundaries Basic example without boundaries
\include Coxeter_triangulation/cell_complex_from_basic_circle_manifold.cpp

The program output is:

\include Coxeter_triangulation/cell_complex_from_basic_circle_manifold_for_doc.txt

\subsection exampleswithboundaries Example with boundaries

Here is an example of constructing a piecewise-linear approximation of a flat torus embedded in \f$\mathbb{R}^4\f$,
rotated by a random rotation in \f$\mathbb{R}^4\f$ and cut by a hyperplane.

\include Coxeter_triangulation/manifold_tracing_flat_torus_with_boundary.cpp

The output in <a target="_blank" href="https://www.ljll.math.upmc.fr/frey/software.html">medit</a> is:

\image html "flat_torus_with_boundary.png" "Output from the example of a flat torus with boundary"

\subsection exampleswithcustomfunction Example with a custom function

In the following more complex example, we define a custom function for the implicit manifold.

\include Coxeter_triangulation/manifold_tracing_custom_function.cpp

The output in <a target="_blank" href="https://www.ljll.math.upmc.fr/frey/software.html">medit</a> looks as follows:

\image html "custom_function.png" "Output from the example with a custom function"


 */
/** @} */  // end defgroup coxeter_triangulation

}  // namespace coxeter_triangulation

}  // namespace Gudhi

#endif  // DOC_COXETER_TRIANGULATION_INTRO_COXETER_TRIANGULATION_H_

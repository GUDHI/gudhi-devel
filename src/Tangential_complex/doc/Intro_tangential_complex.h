/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clement Jamin
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef DOC_TANGENTIAL_COMPLEX_INTRO_TANGENTIAL_COMPLEX_H_
#define DOC_TANGENTIAL_COMPLEX_INTRO_TANGENTIAL_COMPLEX_H_

// needs namespaces for Doxygen to link on classes
namespace Gudhi {
namespace tangential_complex {

/**  \defgroup tangential_complex Tangential complex

\author    Cl&eacute;ment Jamin

@{

\section tangentialdefinition Definition

A Tangential Delaunay complex is a
<a target="_blank" href="https://en.wikipedia.org/wiki/Simplicial_complex">simplicial complex</a>
designed to reconstruct a \f$k\f$-dimensional smooth manifold embedded in \f$d\f$-dimensional Euclidean space. 
The input is a point sample coming from an unknown manifold, which means that the points lie close to a structure of
"small" intrinsic dimension.
The running time depends only linearly on the extrinsic dimension \f$ d \f$
and exponentially on the intrinsic dimension \f$ k \f$.

An extensive description of the Tangential complex can be found in \cite tangentialcomplex2014.

\subsection whatisthetc What is a Tangential Complex?

Let us start with the description of the Tangential complex of a simple example, with \f$ k=1 \f$ and \f$ d=2 \f$.
The point set \f$ \mathscr P \f$ is located on a closed curve embedded in 2D.
Only 4 points will be displayed (more are required for PCA) to simplify the figures.
\image html "tc_example_01.png" "The input"
For each point \f$ P \f$, estimate its tangent subspace \f$ T_P \f$ using PCA.
\image html "tc_example_02.png" "The estimated normals"
Let us add the Voronoi diagram of the points in orange. For each point \f$ P \f$, construct its star in the Delaunay
triangulation of \f$ \mathscr P \f$ restricted to \f$ T_P \f$.
\image html "tc_example_03.png" "The Voronoi diagram"
The Tangential Delaunay complex is the union of those stars.

In practice, neither the ambient Voronoi diagram nor the ambient Delaunay triangulation is computed.
Instead, local \f$ k \f$-dimensional regular triangulations are computed with a limited number of points as we only
need the star of each point. More details can be found in \cite tangentialcomplex2014.

\subsection inconsistencies Inconsistencies

Inconsistencies between the stars can occur.
An inconsistency occurs when a simplex is not in the star of all its vertices.

Let us take the same example.
\image html "tc_example_07_before.png" "Before"
Let us slightly move the tangent subspace \f$ T_Q \f$
\image html "tc_example_07_after.png" "After"
Now, the star of \f$ Q \f$ contains \f$ QP \f$, but the star of \f$ P \f$ does not contain \f$ QP \f$. We have an inconsistency.
\image html "tc_example_08.png" "After"

One way to solve inconsistencies is to randomly perturb the positions of the points involved in an inconsistency.
In the current implementation, this perturbation is done in the tangent subspace of each point. 
The maximum perturbation radius is given as a parameter to the constructor.

In most cases, we recommend to provide a point set where the minimum distance between any two points
is not too small. This can be achieved using the functions provided by the Subsampling module. Then, a good value to start with for
the maximum perturbation radius would be around half the minimum distance between any two points. 
The \ref example_with_perturb below shows an example of such a process.

In most cases, this process is able to dramatically reduce the number of inconsistencies, but is not guaranteed to succeed.

\subsection output Output

The result of the computation is exported as a `Simplex_tree`. It is the union of the stars of all the input points.
A vertex in the Simplex Tree is the index of the point in the range provided by the user. 
The point corresponding to a vertex can also be obtained through the `Tangential_complex::get_point` function.
Note that even if the positions of the points are perturbed, their original positions are kept (e.g. `Tangential_complex::get_point` returns the original position of the point).

The result can be obtained after the computation of the Tangential complex itself and/or after the perturbation process.

\section simple_example Simple example

This example builds the Tangential complex of point set.
Note that the dimension of the kernel here is dynamic, which is slower, but more flexible:
the intrinsic and ambient dimensions does not have to be known at compile-time.

\include Tangential_complex/example_basic.cpp

\section example_with_perturb Example with perturbation

This example builds the Tangential complex of a point set, then tries to solve inconsistencies
by perturbing the positions of points involved in inconsistent simplices.
Note that the dimension of the kernel here is static, which is the best choice when the
dimensions are known at compile-time.

\include Tangential_complex/example_with_perturb.cpp

 */
/** @} */  // end defgroup tangential_complex

}  // namespace tangential_complex

}  // namespace Gudhi

#endif  // DOC_TANGENTIAL_COMPLEX_INTRO_TANGENTIAL_COMPLEX_H_

/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef DOC_TANGENTIAL_COMPLEX_INTRO_COXETER_TRIANGULATION_H_
#define DOC_TANGENTIAL_COMPLEX_INTRO_COXETER_TRIANGULATION_H_

// needs namespaces for Doxygen to link on classes
namespace Gudhi {
namespace tangential_complex {

/**  \defgroup coxter_triangulation Coxeter triangulation

\author    Siargey Kachanovich

@{

\section Module overview

Coxeter triangulation module is designed to provide tools for constructing a piecewise-linear approximation of an \f$m\f$-dimensional smooth manifold embedded in \f$\mathbb{R}^d\f$ using an ambient triangulation.

\section manifoldtracing Manifold tracing algorithm
The central piece of the module is the manifold tracing algorithm represented by the class Gudhi::Manifold_tracing.

\section ambienttriangulations Ambient triangulations

The ambient triangulations available in the module are \f$\mathcal{T}\f$ which is a linear transformation of the Freudenthal-Kuhn triangulation of \f$\mathbb{R}^d\f$.

\section intersectionoracle Implicit manifold intersection oracle


\section cellcomplex Cell complex construction

 */
/** @} */  // end defgroup coxeter_triangulation

}  // namespace coxeter_triangulation

}  // namespace Gudhi

#endif  // DOC_TANGENTIAL_COMPLEX_INTRO_TANGENTIAL_COMPLEX_H_

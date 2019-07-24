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

#ifndef CONCEPT_COXETER_TRIANGULATION_FUNCTION_FOR_IMPLICIT_MANIFOLD_H_
#define CONCEPT_COXETER_TRIANGULATION_FUNCTION_FOR_IMPLICIT_MANIFOLD_H_

namespace Gudhi {

namespace coxeter_triangulation {

/** \brief The concept FunctionForImplicitManifold describes the requirements 
 * for a type to implement an implicit function class used for example in Manifold_tracing.
 */
struct FunctionForImplicitManifold {

  /** \brief Value of the function at a specified point 'p'. */
  void evaluate(const Eigen::VectorXd& p, Eigen::VectorXd& result) const;

  /** \brief Returns the domain (ambient) dimension. */
  std::size_t amb_d() const;
  
  /** \brief Returns the codomain dimension. */
  std::size_t cod_d() const;

  /** \brief Returns a point on the zero-set of the function. */
  void seed(Eigen::VectorXd& result) const;
};

}  // namespace coxeter_triangulation

}  // namespace Gudhi


#endif

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

#ifndef CONCEPT_COXETER_TRIANGULATION_DOMAIN_FOR_MANIFOLD_TRACING_H_
#define CONCEPT_COXETER_TRIANGULATION_DOMAIN_FOR_MANIFOLD_TRACING_H_

namespace Gudhi {

namespace coxeter_triangulation {

/** \brief The concept DomainForManifoldTracing describes the requirements 
 * for a type to implement a domain class used for example in Implicit_manifold_intersection_oracle.
 */
struct DomainForManifoldTracing {
  
  /** \brief Returns true if the specified point 'p' lies inside the domain. */
  bool operator()(const Eigen::VectorXd& p) const;

};

}  // namespace coxeter_triangulation

}  // namespace Gudhi


#endif

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

namespace Gudhi {

namespace coxeter_triangulation {

/** \class Query_result
 *  \brief The result of a query by an oracle such as Implicit_manifold_intersection_oracle. 
 *
 *  \tparam Simplex_handle The class of the query simplex.
 *
 *  \ingroup coxeter_triangulation
 */
template <class Simplex_handle>
struct Query_result {
  /** \brief The potentially lower-dimensional face of the query simplex
   *   that contains the intersection point. */
  Simplex_handle face;
  /** \brief The intersection point. */
  Eigen::VectorXd intersection;
  /** \brief True if the query simplex intersects the manifold. */
  bool success;
}

} // namespace coxeter_triangulation 

} // namespace Gudhi

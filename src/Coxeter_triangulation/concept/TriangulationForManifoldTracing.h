/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CONCEPT_COXETER_TRIANGULATION_TRIANGULATION_FOR_MANIFOLD_TRACING_H_
#define CONCEPT_COXETER_TRIANGULATION_TRIANGULATION_FOR_MANIFOLD_TRACING_H_

#include <Eigen/Dense>

namespace Gudhi {

namespace coxeter_triangulation {

/** \brief The concept TriangulationForManifoldTracing describes the requirements
 * for a type to implement a triangulation class used for example in Manifold_tracing.
 */
struct TriangulationForManifoldTracing {
  /** \brief Type of the simplices in the triangulation.
   *  Needs to be a model of the concept SimplexInCoxeterTriangulation. */
  typedef Simplex_handle;

  /** \brief Type of the vertices in the triangulation.
   *  Needs to be a random-access range of integer values. */
  typedef Vertex_handle;

  /** \brief Returns the permutahedral representation of the simplex in the
   *  triangulation that contains a given query point 'p'.
   * \tparam Point_d A class that represents a point in d-dimensional Euclidean space.
   * The coordinates should be random-accessible. Needs to provide the method size().
   * @param[in] point The query point.
   */
  template <class Point_d>
  Simplex_handle locate_point(const Point_d& point) const;

  /** \brief Returns the Cartesian coordinates of the given vertex 'v'.
   * @param[in] v The input vertex.
   */
  Eigen::VectorXd cartesian_coordinates(const Vertex_handle& v) const;

  /** \brief Returns the Cartesian coordinates of the barycenter of a given simplex 's'.
   * @param[in] s The input simplex given by permutahedral representation.
   */
  Eigen::VectorXd barycenter(const Simplex_handle& s) const;
};

}  // namespace coxeter_triangulation

}  // namespace Gudhi

#endif

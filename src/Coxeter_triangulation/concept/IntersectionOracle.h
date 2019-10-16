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

#ifndef CONCEPT_COXETER_TRIANGULATION_INTERSECTION_ORACLE_H_
#define CONCEPT_COXETER_TRIANGULATION_INTERSECTION_ORACLE_H_

namespace Gudhi {

namespace coxeter_triangulation {

/** \brief The concept IntersectionOracle describes the requirements 
 * for a type to implement an intersection oracle class used for example in Manifold_tracing.
 *
 */
struct IntersectionOracle {

  /** \brief Returns the domain (ambient) dimension of the underlying manifold. */
  std::size_t amb_d() const;
  
  /** \brief Returns the codomain dimension of the underlying manifold. */
  std::size_t cod_d() const;

  /** \brief Intersection query with the relative interior of the manifold.
   *
   *  \details The returned structure Query_result contains the boolean value
   *   that is true only if the intersection point of the query simplex and
   *   the relative interior of the manifold exists, the intersection point
   *   and the face of the query simplex that contains 
   *   the intersection point.
   *   
   *  \tparam Simplex_handle The class of the query simplex.
   *   Needs to be a model of the concept SimplexInCoxeterTriangulation.
   *  \tparam Triangulation The class of the triangulation.
   *   Needs to be a model of the concept TriangulationForManifoldTracing.
   *  
   *  @param[in] simplex The query simplex. The dimension of the simplex
   *   should be the same as the codimension of the manifold 
   *   (the codomain dimension of the function).
   *  @param[in] triangulation The ambient triangulation. The dimension of 
   *   the triangulation should be the same as the ambient dimension of the manifold 
   *   (the domain dimension of the function).
   */
  template <class Simplex_handle,
	    class Triangulation>
  Query_result<Simplex_handle> intersects(const Simplex_handle& simplex,
					  const Triangulation& triangulation) const;

  /** \brief Intersection query with the boundary of the manifold.
   *  
   *  \details The returned structure Query_result contains the boolean value
   *   that is true only if the intersection point of the query simplex and
   *   the boundary of the manifold exists, the intersection point
   *   and the face of the query simplex that contains 
   *   the intersection point.
   *   
   *  \tparam Simplex_handle The class of the query simplex.
   *   Needs to be a model of the concept SimplexInCoxeterTriangulation.
   *  \tparam Triangulation The class of the triangulation.
   *   Needs to be a model of the concept TriangulationForManifoldTracing.
   *
   *  @param[in] simplex The query simplex. The dimension of the simplex
   *   should be the same as the codimension of the boundary of the manifold 
   *   (the codomain dimension of the function + 1).
   *  @param[in] triangulation The ambient triangulation. The dimension of 
   *   the triangulation should be the same as the ambient dimension of the manifold 
   *   (the domain dimension of the function).
   */
  template <class Simplex_handle,
	    class Triangulation>
  Query_result<Simplex_handle> intersects_boundary(const Simplex_handle& simplex,
						   const Triangulation& triangulation) const;

  /** \brief Returns true if the input point lies inside the piecewise-linear
   *   domain induced by the given ambient triangulation that defines the relative
   *   interior of the piecewise-linear approximation of the manifold.
   *
   * @param p The input point. Needs to have the same dimension as the ambient
   *  dimension of the manifold (the domain dimension of the function).
   * @param triangulation The ambient triangulation. Needs to have the same
   *  dimension as the ambient dimension of the manifold 
   *  (the domain dimension of the function).
   */
  template <class Triangulation>
  bool lies_in_domain(const Eigen::VectorXd& p,
		      const Triangulation& triangulation) const {
    Eigen::VectorXd pl_p = make_pl_approximation(domain_fun_, triangulation)(p);
    return pl_p(0) < 0;
  }

  /** \brief Returns the function that defines the interior of the manifold */
  const Function_& function() const;

};

}  // namespace coxeter_triangulation

}  // namespace Gudhi


#endif

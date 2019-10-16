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

#ifndef IMPLICIT_MANIFOLD_INTERSECTION_ORACLE_H_
#define IMPLICIT_MANIFOLD_INTERSECTION_ORACLE_H_

#include <Eigen/Dense>

#include <gudhi/Permutahedral_representation/face_from_indices.h>
#include <gudhi/Functions/Constant_function.h>
#include <gudhi/Functions/PL_approximation.h>
#include <gudhi/Query_result.h>

#include <vector>

namespace Gudhi {

namespace coxeter_triangulation {



/** \class Implicit_manifold_intersection_oracle
 *  \brief An oracle that supports the intersection query on an implicit manifold.
 *
 *  \tparam Function_ The function template parameter. Should be a model of 
 *   the concept FunctionForImplicitManifold.
 *  \tparam Domain_function_ The domain function template parameter. Should be a model of
 *   the concept FunctionForImplicitManifold.
 *
 *  \ingroup coxeter_triangulation
 */
template<class Function_,
	 class Domain_function_ = Constant_function>
class Implicit_manifold_intersection_oracle {

  /* Computes the affine coordinates of the intersection point of the implicit manifold
   * and the affine hull of the simplex. */
  template <class Simplex_handle, 
	    class Triangulation>
  Eigen::VectorXd compute_lambda(const Simplex_handle& simplex,
				 const Triangulation& triangulation) const {
    std::size_t cod_d = this->cod_d();
    Eigen::MatrixXd matrix(cod_d + 1, cod_d + 1);
    for (std::size_t i = 0; i < cod_d + 1; ++i)
      matrix(0, i) = 1;
    std::size_t j = 0;
    for (auto v: simplex.vertex_range()) {
      Eigen::VectorXd v_coords = fun_(triangulation.cartesian_coordinates(v));
      for (std::size_t i = 1; i < cod_d + 1; ++i)
	matrix(i, j) = v_coords(i-1);
      j++;
    }
    Eigen::VectorXd z(cod_d + 1);
    z(0) = 1;
    for (std::size_t i = 1; i < cod_d + 1; ++i)
      z(i) = 0;
    Eigen::VectorXd lambda = matrix.colPivHouseholderQr().solve(z);
    return lambda;
  }

  /* Computes the affine coordinates of the intersection point of the boundary
   * of the implicit manifold and the affine hull of the simplex. */
  template <class Simplex_handle, 
	    class Triangulation>
  Eigen::VectorXd compute_boundary_lambda(const Simplex_handle& simplex,
					  const Triangulation& triangulation) const {
    std::size_t cod_d = this->cod_d();
    Eigen::MatrixXd matrix(cod_d + 2, cod_d + 2);
    for (std::size_t i = 0; i < cod_d + 2; ++i)
      matrix(0, i) = 1;
    std::size_t j = 0;
    for (auto v: simplex.vertex_range()) {
      Eigen::VectorXd v_coords = fun_(triangulation.cartesian_coordinates(v));
      for (std::size_t i = 1; i < cod_d + 1; ++i)
	matrix(i, j) = v_coords(i-1);
      Eigen::VectorXd bv_coords = domain_fun_(triangulation.cartesian_coordinates(v));
      matrix(cod_d + 1, j) = bv_coords(0);
      j++;
    }
    Eigen::VectorXd z(cod_d + 2);
    z(0) = 1;
    for (std::size_t i = 1; i < cod_d + 2; ++i)
      z(i) = 0;
    Eigen::VectorXd lambda = matrix.colPivHouseholderQr().solve(z);
    return lambda;
  }

  /* Computes the intersection result for a given simplex in a triangulation. */
  template <class Simplex_handle,
	    class Triangulation>
  Query_result<Simplex_handle> intersection_result(const Eigen::VectorXd& lambda,
						   const Simplex_handle& simplex,
						   const Triangulation& triangulation) const {
    using QR = Query_result<Simplex_handle>;
    std::size_t amb_d = triangulation.dimension();
    std::size_t cod_d = simplex.dimension();

    for (std::size_t i = 0; i < (std::size_t)lambda.size(); ++i)
      if (lambda(i) < 0 || lambda(i) > 1)
	return QR({Eigen::VectorXd(), false});

    Eigen::MatrixXd vertex_matrix(cod_d + 1, amb_d);
    auto v_range = simplex.vertex_range();
    auto v_it = v_range.begin();
    for (std::size_t i = 0; i < cod_d + 1 && v_it != v_range.end(); ++v_it, ++i) {
      Eigen::VectorXd v_coords = triangulation.cartesian_coordinates(*v_it);
      for (std::size_t j = 0; j < amb_d; ++j)
	vertex_matrix(i, j) = v_coords(j);
    }
    Eigen::VectorXd intersection = lambda.transpose()*vertex_matrix;
    return QR({intersection, true});
  }
  
public:

  /** \brief Ambient dimension of the implicit manifold. */
  std::size_t amb_d() const {
    return fun_.amb_d();
  }
  
  /** \brief Codimension of the implicit manifold. */
  std::size_t cod_d() const {
    return fun_.cod_d();
  }

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
					  const Triangulation& triangulation) const {
    Eigen::VectorXd lambda = compute_lambda(simplex, triangulation);
    return intersection_result(lambda, simplex, triangulation);
  }

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
						   const Triangulation& triangulation) const {
    Eigen::VectorXd lambda = compute_boundary_lambda(simplex, triangulation);
    return intersection_result(lambda, simplex, triangulation);
  }

  
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
  const Function_& function() const {
    return fun_;
  }
  
  /** \brief Constructs an intersection oracle for an implicit manifold potentially 
   *   with boundary from given function and domain.
   *
   *  @param function The input function that represents the implicit manifold
   *   before the restriction with the domain.
   *  @param domain_function The input domain function that can be used to define an implicit
   *   manifold with boundary.
   */
  Implicit_manifold_intersection_oracle(const Function_& function,
					const Domain_function_& domain_function)
    : fun_(function), domain_fun_(domain_function) {}

  /** \brief Constructs an intersection oracle for an implicit manifold 
   *   without boundary from a given function.
   *
   *   \details To use this constructor, the template Domain_function_ needs to be left 
   *   at its default value (Gudhi::coxeter_triangulation::Constant_function).
   *
   *  @param function The input function that represents the implicit manifold
   *   without boundary.
   */
  Implicit_manifold_intersection_oracle(const Function_& function)
    : fun_(function),
      domain_fun_(function.amb_d(), 1, Eigen::VectorXd::Constant(1,-1)) {}
  
private:
  Function_ fun_;
  Domain_function_ domain_fun_;
};

/** \brief Static constructor of an intersection oracle from a function with a domain.
 *
 *  @param function The input function that represents the implicit manifold
 *   before the restriction with the domain.
 *  @param domain_function The input domain function that can be used to define an implicit
 *   manifold with boundary.
 *
 *  \ingroup coxeter_triangulation
 */
template<class Function_,
	 class Domain_function_>
Implicit_manifold_intersection_oracle<Function_, Domain_function_>
make_oracle(const Function_& function,
	    const Domain_function_& domain_function){
  return Implicit_manifold_intersection_oracle<Function_, Domain_function_>(function,
									    domain_function);
}


/** \brief Static constructor of an intersection oracle from a function without a domain.
 *
 *  @param function The input function that represents the implicit manifold
 *   without boundary.
 *
 *  \ingroup coxeter_triangulation
 */
template<class Function_>
Implicit_manifold_intersection_oracle<Function_> make_oracle(const Function_& function){
  return Implicit_manifold_intersection_oracle<Function_>(function);
}

} // namespace coxeter_triangulation 

} // namespace Gudhi

#endif

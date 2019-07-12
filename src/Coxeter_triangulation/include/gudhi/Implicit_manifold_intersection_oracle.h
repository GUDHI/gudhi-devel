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
#include <gudhi/Functions/Domain_from_function.h>
#include <gudhi/Functions/Constant_function.h>
#include <gudhi/Query_result.h>

#include <vector>

namespace Gudhi {

namespace coxeter_triangulation {

using Infinite_domain = Domain_from_function<Constant_function>;

/** \class Implicit_manifold_intersection_oracle
 *  \brief An oracle that supports the intersection query on an implicit manifold.
 *
 *  \tparam Function_ The function template parameter. Should be a model of 
 *   the concept FunctionForImplicitManifold.
 *  \tparam Domain_ The domain template parameter. Should be a model of
 *   the concept DomainForManifoldTracing.
 *
 *  \ingroup coxeter_triangulation
 */
template<class Function_,
	 class Domain_ = Infinite_domain>
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
    return matrix.colPivHouseholderQr().solve(z);
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
      matrix(cod_d + 1, j) = domain_.function_(triangulation.cartesian_coordinates(v))(0);
      j++;
    }
    Eigen::VectorXd z(cod_d + 2);
    z(0) = 1;
    for (std::size_t i = 1; i < cod_d + 2; ++i)
      z(i) = 0;
    return matrix.colPivHouseholderQr().solve(z);
  }

  /* Computes the intersection result for a given simplex in a triangulation. */
  template <class Simplex_handle,
	    class Triangulation>
  Query_result<Simplex_handle> intersection_result(const Eigen::VectorXd& lambda,
						   const Simplex_handle& simplex,
						   const Triangulation& triangulation) const {
    using QR = Query_result<Simplex_handle>;
    std::size_t amb_d = triangulation.dimension();
    std::size_t cod_d = fun_->cod_d();

    std::vector<std::size_t> snapping_indices;
    for (std::size_t i = 0; i < lambda.size(); ++i) {
      if (lambda(i) < -threshold_ || lambda(i) > 1 + threshold_)
	return QR({Simplex_handle(), Eigen::VectorXd(), false});
      if (lambda(i) >= threshold_ && lambda(i) <= 1 - threshold_)
	snapping_indices.push_back(i);
    }

    std::size_t snap_d = snapping_indices.size();
    std::size_t i = 0;
    std::size_t num_line = 0;
    Eigen::MatrixXd vertex_matrix(snap_d, amb_d);
    Eigen::VectorXd reduced_lambda(snap_d);
    auto v_range = simplex.vertex_range();
    auto v_it = v_range.begin();
    for (; num_line < snap_d && v_it != v_range.end(); ++v_it, ++i) {
      if (i == snapping_indices[num_line]) {
	Eigen::VectorXd v_coords = triangulation.cartesian_coordinates(*v_it);
	for (std::size_t j = 0; j < amb_d; ++j)
	  vertex_matrix(num_line, j) = v_coords(j);
	reduced_lambda(num_line) = lambda(i);
	num_line++;
      }
    }
    reduced_lambda /= reduced_lambda.sum();
    Eigen::VectorXd intersection = reduced_lambda.transpose()*vertex_matrix;
    return QR({face_from_indices(simplex, snapping_indices), intersection, true});
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
					  const Triangulation& triangulation) {
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
						   const Triangulation& triangulation) {
    Eigen::VectorXd lambda = compute_boundary_lambda(simplex, triangulation);
    return intersection_result(lambda, simplex, triangulation);
  }

  bool lies_in_domain(const Eigen::VectorXd& p) {
    return true;
  }
  
  /** \brief Constructs an intersection oracle for an implicit manifold potentially 
   *   with boundary from given function and domain.
   *
   *  @param function The input function that represents the implicit manifold
   *   before the restriction with the domain.
   *  @param domain The input domain that can be used to define an implicit
   *   manifold with boundary.
   *  @param threshold The input parameter that defines the distance in terms
   *   of affine coordinates to a lower-dimensional face for the intersection 
   *   point to be considered lying at the lower-dimensional face.
   */
  Implicit_manifold_intersection_oracle(const Function_& function,
					const Domain_& domain,
					double threshold = 0)
    : fun_(function), domain_(domain), threshold_(threshold) {}

  /** \brief Constructs an intersection oracle for an implicit manifold 
   *   without boundary from a given function.
   *
   *  @param function The input function that represents the implicit manifold
   *   without boundary.
   *  @param threshold The input parameter that defines the distance in terms
   *   of affine coordinates to a lower-dimensional face for the intersection 
   *   point to be considered lying at the lower-dimensional face.
   */
  Implicit_manifold_intersection_oracle(const Function_& function, double threshold = 0)
    : fun_(function),
      domain_(Constant_function(function.amb_d(), 1, Eigen::VectorXd::Constant(1,-1))),
      threshold_(threshold) {}
  
private:
  Function_ fun_;
  Domain_ domain_;
  double threshold_ = 0;
};

template<class Function_,
	 class Domain_>
Implicit_manifold_intersection_oracle<Function_, Domain_> make_oracle(const Function_& f,
								     const Domain_& dom,
								     double threshold = 0){
  return Implicit_manifold_intersection_oracle<Function_, Domain_>(f, dom, threshold);
}


template<class Function_>
Implicit_manifold_intersection_oracle<Function_> make_oracle(const Function_& f){
  return Implicit_manifold_intersection_oracle<Function_>(f);
}

} // namespace coxeter_triangulation 

} // namespace Gudhi

#endif

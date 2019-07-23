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

#ifndef FUNCTIONS_PL_APPROXIMATION_H_
#define FUNCTIONS_PL_APPROXIMATION_H_

#include <cstdlib>

#include <gudhi/Functions/Function.h>
#include <Eigen/Dense>

namespace Gudhi {

namespace coxeter_triangulation {

/* \class PL_approximation
 * \brief Constructs a piecewise-linear approximation of a function induced by
 * an ambient triangulation.
 *
 * \tparam Function The function template parameter. Should be a model of 
 * the concept FunctionForImplicitManifold.
 * \tparam Triangulation The triangulation template parameter. Should be a model of
 * the concept TriangulationForManifoldTracing.
 *
 * \ingroup coxeter_triangulation
 */
template <class Function,
	  class Triangulation>
struct PL_approximation : public Function {
  
  /** 
   * \brief Value of the function at a specified point.
   * @param[in] p The input point. The dimension needs to coincide with the ambient dimension.
   */
  void evaluate(const Eigen::VectorXd& p, Eigen::VectorXd& result) const {
    std::size_t cod_d = this->cod_d();
    std::size_t amb_d = this->amb_d();
    auto s = tr_.locate_point(p);
    Eigen::MatrixXd matrix(cod_d, s.dimension() + 1);
    Eigen::MatrixXd vertex_matrix(amb_d + 1, s.dimension() + 1);
    for (std::size_t i = 0; i < s.dimension() + 1; ++i)
      vertex_matrix(0, i) = 1;
    std::size_t j = 0;
    for (auto v: s.vertex_range()) {
      Eigen::VectorXd pt_v = tr_.cartesian_coordinates(v);
      Eigen::VectorXd fun_v;
      fun_.evaluate(pt_v, fun_v);
      for (std::size_t i = 1; i < amb_d + 1; ++i)
	vertex_matrix(i, j) = pt_v(i-1); 
      for (std::size_t i = 0; i < cod_d; ++i)
	matrix(i, j) = fun_v(i);
      j++;
    }
    assert(j == s.dimension()+1);
    Eigen::VectorXd z(amb_d + 1);
    z(0) = 1;
    for (std::size_t i = 1; i < amb_d + 1; ++i)
      z(i) = p(i-1);
    Eigen::VectorXd lambda = vertex_matrix.colPivHouseholderQr().solve(z);    
    result = matrix * lambda;
  }

  /** \brief Returns the domain (ambient) dimension. */
  std::size_t amb_d() const {return fun_.amb_d();}

  /** \brief Returns the codomain dimension. */
  std::size_t cod_d() const {return fun_.cod_d();}

  /** \brief Returns a point on the zero-set. */
  void seed(Eigen::VectorXd& result) const {
    // Eigen::VectorXd seed_k(fun_.amb_d());
    // auto c = tr_.locate_point(fun_.seed());
    // for (auto f: c.face_range(fun_.cod_d())) {      
    // }
    // TODO: not finished. Should use an oracle.
    // Also need to make sure, s is of the highest dimension.
    // return seed_k;
  }

  /** 
   * \brief Constructor of the piecewise-linear approximation of a function 
   * induced by an ambient triangulation.
   *
   * @param[in] function The function.
   * @param[in] triangulation The ambient triangulation.
   */
  PL_approximation(const Function& function, const Triangulation& triangulation) :
    fun_(function), tr_(triangulation) {}

private:  
  Function fun_;
  Triangulation tr_;
};


/** 
 * \brief Static constructor of the piecewise-linear approximation of a function 
 * induced by an ambient triangulation.
 *
 * @param[in] function The function.
 * @param[in] triangulation The ambient triangulation.
 *
 * \tparam Function The function template parameter. Should be a model of 
 * the concept FunctionForImplicitManifold.
 */
template <class Function,
	  class Triangulation>
PL_approximation<Function, Triangulation> make_pl_approximation(const Function& function,
								const Triangulation& triangulation) {
  return PL_approximation<Function, Triangulation>(function, triangulation);
}

} // namespace coxeter_triangulation

} // namespace Gudhi

#endif

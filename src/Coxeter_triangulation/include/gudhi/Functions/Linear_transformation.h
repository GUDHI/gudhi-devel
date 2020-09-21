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

#ifndef FUNCTIONS_LINEAR_TRANSFORMATION_H_
#define FUNCTIONS_LINEAR_TRANSFORMATION_H_

#include <cstdlib>

#include <gudhi/Functions/Function.h>
#include <Eigen/Dense>

namespace Gudhi {

namespace coxeter_triangulation {

/** \class Linear_transformation
 * \brief Transforms the zero-set of the function by a given linear transformation.
 * The underlying function corresponds to f(M*x), where M is the transformation matrix.
 *
 * \tparam Function_ The function template parameter. Should be a model of 
 * the concept FunctionForImplicitManifold.
 *
 * \ingroup coxeter_triangulation
 */
template <class Function_>
struct Linear_transformation : public Function {
  
  /** 
   * \brief Value of the function at a specified point.
   * @param[in] p The input point. The dimension needs to coincide with the ambient dimension.
   */
  Eigen::VectorXd operator()(const Eigen::VectorXd& p) const {
    Eigen::VectorXd result = fun_(matrix_.householderQr().solve(p));
    return result;
  }

  /** \brief Returns the domain (ambient) dimension. */
  std::size_t amb_d() const {return fun_.amb_d();}

  /** \brief Returns the codomain dimension. */
  std::size_t cod_d() const {return fun_.cod_d();}

  /** \brief Returns a point on the zero-set. */
  Eigen::VectorXd seed() const {
    Eigen::VectorXd result = fun_.seed();
    result = matrix_ * result;
    return result;
  }

  /** 
   * \brief Constructor of a linearly transformed function.
   *
   * @param[in] function The function to be linearly transformed.
   * @param[in] matrix The transformation matrix. Its dimension should be d*d,
   * where d is the domain (ambient) dimension of 'function'.
   */
  Linear_transformation(const Function_& function, const Eigen::MatrixXd& matrix) :
    fun_(function), matrix_(matrix) {
  }

private:
  Function_ fun_;
  Eigen::MatrixXd matrix_;
};


/** 
 * \brief Static constructor of a linearly transformed function.
 *
 * @param[in] function The function to be linearly transformed.
 * @param[in] matrix The transformation matrix. Its dimension should be d*d,
 * where d is the domain (ambient) dimension of 'function'.
 *
 * \tparam Function_ The function template parameter. Should be a model of 
 * the concept FunctionForImplicitManifold.
 */
template <class Function_>
Linear_transformation<Function_> make_linear_transformation(const Function_& function,
						       const Eigen::MatrixXd& matrix) {
  return Linear_transformation<Function_>(function, matrix); 
}

} // namespace coxeter_triangulation

} // namespace Gudhi

#endif

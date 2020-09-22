/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef FUNCTIONS_FUNCTION_WHITNEY_UMBRELLA_IN_R3_H_
#define FUNCTIONS_FUNCTION_WHITNEY_UMBRELLA_IN_R3_H_

#include <cstdlib>  // for std::size_t

#include <gudhi/Functions/Function.h>

#include <Eigen/Dense>

namespace Gudhi {

namespace coxeter_triangulation {

/**
 * \class Function_whitney_umbrella_in_R3
 * \brief A class that encodes the function, the zero-set of which is the Whitney umbrella
 * surface embedded in R^3.
 *
 * \ingroup coxeter_triangulation
 */
struct Function_whitney_umbrella_in_R3 : public Function {
  /**
   * \brief Value of the function at a specified point.
   * @param[in] p The input point. The dimension needs to coincide with the ambient dimension.
   */
  virtual Eigen::VectorXd operator()(const Eigen::VectorXd& p) const override {
    double x = p(0) - off_[0], y = p(1) - off_[1], z = p(2) - off_[2];
    Eigen::VectorXd result(cod_d());
    result(0) = x * x - y * y * z;
    return result;
  }

  /** \brief Returns the (ambient) domain dimension.*/
  virtual std::size_t amb_d() const override { return 3; };

  /** \brief Returns the codomain dimension. */
  virtual std::size_t cod_d() const override { return 1; };

  /** \brief Returns a point on the surface. This seed point is only one of
   * two necessary seed points for the manifold tracing algorithm.
   * See the method seed2() for the other point.
   */
  virtual Eigen::VectorXd seed() const override {
    Eigen::Vector3d result(1 + off_[0], 1 + off_[1], 1 + off_[2]);
    return result;
  }

  /** \brief Returns a point on the surface. This seed point is only one of
   * two necessary seed points for the manifold tracing algorithm.
   * See the method seed() for the other point.
   */
  Eigen::VectorXd seed2() const {
    Eigen::Vector3d result(-1 + off_[0], -1 + off_[1], 1 + off_[2]);
    return result;
  }

  /**
   * \brief Constructor of the function that defines the Whitney umbrella in R^3.
   *
   * @param[in] off Offset vector.
   */
  Function_whitney_umbrella_in_R3(Eigen::Vector3d off = Eigen::Vector3d::Zero()) : off_(off) {}

 private:
  Eigen::Vector3d off_;
};

}  // namespace coxeter_triangulation

}  // namespace Gudhi

#endif

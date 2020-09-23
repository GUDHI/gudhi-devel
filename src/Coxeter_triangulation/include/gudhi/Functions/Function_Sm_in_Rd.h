/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef FUNCTIONS_FUNCTION_SM_IN_RD_H_
#define FUNCTIONS_FUNCTION_SM_IN_RD_H_

#include <cstdlib>  // for std::size_t

#include <gudhi/Functions/Function.h>

#include <Eigen/Dense>

namespace Gudhi {

namespace coxeter_triangulation {

/**
 * \class Function_Sm_in_Rd
 * \brief A class for the function that defines an m-dimensional implicit sphere embedded
 * in the d-dimensional Euclidean space.
 */
struct Function_Sm_in_Rd : public Function {
  /** \brief Value of the function at a specified point.
   * @param[in] p The input point. The dimension needs to coincide with the ambient dimension.
   */
  virtual Eigen::VectorXd operator()(const Eigen::VectorXd& p) const override {
    Eigen::VectorXd x = p;
    for (std::size_t i = 0; i < d_; ++i) x(i) -= center_[i];
    Eigen::VectorXd result = Eigen::VectorXd::Zero(k_);
    for (std::size_t i = 0; i < m_ + 1; ++i) result(0) += x(i) * x(i);
    result(0) -= r_ * r_;
    for (std::size_t j = 1; j < k_; ++j) result(j) = x(m_ + j);
    return result;
  }

  /** \brief Returns the domain dimension. Same as the ambient dimension of the sphere. */
  virtual std::size_t amb_d() const override { return d_; };

  /** \brief Returns the codomain dimension. Same as the codimension of the sphere. */
  virtual std::size_t cod_d() const override { return k_; };

  /** \brief Returns a point on the sphere. */
  virtual Eigen::VectorXd seed() const override {
    Eigen::VectorXd result = Eigen::VectorXd::Zero(d_);
    result(0) += r_;
    for (std::size_t i = 0; i < d_; ++i) result(i) += center_[i];
    return result;
  }

  /**
   * \brief Constructor of the function that defines an m-dimensional implicit sphere embedded
   * in the d-dimensional Euclidean space.
   *
   * @param[in] r The radius of the sphere.
   * @param[in] m The dimension of the sphere.
   * @param[in] d The ambient dimension of the sphere.
   * @param[in] center The center of the sphere.
   */
  Function_Sm_in_Rd(double r, std::size_t m, std::size_t d, Eigen::VectorXd center)
      : m_(m), k_(d - m), d_(d), r_(r), center_(center) {}

  /**
   * \brief Constructor of the function that defines an m-dimensional implicit sphere embedded
   * in the d-dimensional Euclidean space centered at the origin.
   *
   * @param[in] r The radius of the sphere.
   * @param[in] m The dimension of the sphere.
   * @param[in] d The ambient dimension of the sphere.
   */
  Function_Sm_in_Rd(double r, std::size_t m, std::size_t d)
      : m_(m), k_(d - m), d_(d), r_(r), center_(Eigen::VectorXd::Zero(d_)) {}

  /**
   * \brief Constructor of the function that defines an m-dimensional implicit sphere embedded
   * in the (m+1)-dimensional Euclidean space.
   *
   * @param[in] r The radius of the sphere.
   * @param[in] m The dimension of the sphere.
   * @param[in] center The center of the sphere.
   */
  Function_Sm_in_Rd(double r, std::size_t m, Eigen::VectorXd center)
      : m_(m), k_(1), d_(m_ + 1), r_(r), center_(center) {}

  /**
   * \brief Constructor of the function that defines an m-dimensional implicit sphere embedded
   * in the (m+1)-dimensional Euclidean space centered at the origin.
   *
   * @param[in] r The radius of the sphere.
   * @param[in] m The dimension of the sphere.
   */
  Function_Sm_in_Rd(double r, std::size_t m) : m_(m), k_(1), d_(m_ + 1), r_(r), center_(Eigen::VectorXd::Zero(d_)) {}

  Function_Sm_in_Rd(const Function_Sm_in_Rd& rhs) : Function_Sm_in_Rd(rhs.r_, rhs.m_, rhs.d_, rhs.center_) {}

 private:
  std::size_t m_, k_, d_;
  double r_;
  Eigen::VectorXd center_;
};

}  // namespace coxeter_triangulation

}  // namespace Gudhi

#endif

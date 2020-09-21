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

#ifndef FUNCTIONS_FUNCTION_MOMENT_CURVE_IN_RD_H_
#define FUNCTIONS_FUNCTION_MOMENT_CURVE_IN_RD_H_

#include <gudhi/Functions/Function.h>
#include <Eigen/Dense>

namespace Gudhi {

namespace coxeter_triangulation {

/** 
 * \class Function_moment_curve_in_Rd 
 * \brief A class for the function that defines an implicit moment curve
 * in the d-dimensional Euclidean space.
 *
 * \ingroup coxeter_triangulation
 */
struct Function_moment_curve_in_Rd : public Function {

/** \brief Value of the function at a specified point.
 * @param[in] p The input point. The dimension needs to coincide with the ambient dimension.
 */
  virtual Eigen::VectorXd operator()(const Eigen::VectorXd& p) const override {
    Eigen::VectorXd result(k_);
    for (std::size_t i = 1; i < d_; ++i)
      result(i-1) = p(i) - p(0) * p(i-1);
    return result;
  }
  
  /** \brief Returns the domain (ambient) dimension.. */
  virtual std::size_t amb_d() const override {return d_;};

  /** \brief Returns the codomain dimension. */
  virtual std::size_t cod_d() const override {return k_;};

  /** \brief Returns a point on the moment curve. */
  virtual Eigen::VectorXd seed() const override {
    Eigen::VectorXd result = Eigen::VectorXd::Zero(d_);
    return result;
  }
  
  /** 
   * \brief Constructor of the function that defines an implicit moment curve
   * in the d-dimensional Euclidean space.
   *
   * @param[in] r Numerical parameter.
   * @param[in] d The ambient dimension.
   */
  Function_moment_curve_in_Rd(double r,
			      std::size_t d)
    : m_(1), k_(d-1), d_(d), r_(r) {}

  /** 
   * \brief Constructor of the function that defines an implicit moment curve
   * in the d-dimensional Euclidean space.
   *
   * @param[in] r Numerical parameter.
   * @param[in] d The ambient dimension.
   * @param[in] offset The offset of the moment curve.
   */
  Function_moment_curve_in_Rd(double r,
			      std::size_t d,
			      Eigen::VectorXd& offset)
    : m_(1), k_(d-1), d_(d), r_(r), off_(offset) {}
  
protected:
  std::size_t m_, k_, d_;
  double r_;
  Eigen::VectorXd off_;
};

} // namespace coxeter_triangulation

} // namespace Gudhi


#endif

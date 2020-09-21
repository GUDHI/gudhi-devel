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

#ifndef FUNCTIONS_FUNCTION_CHAIR_IN_R3_H_
#define FUNCTIONS_FUNCTION_CHAIR_IN_R3_H_

#include <cstdlib>

#include <gudhi/Functions/Function.h>
#include <Eigen/Dense>

namespace Gudhi {

namespace coxeter_triangulation {

/** 
 * \class Function_chair_in_R3 
 * \brief A class that encodes the function, the zero-set of which is a so-called
 * "chair" surface embedded in R^3.
 *
 * \ingroup coxeter_triangulation
 */
struct Function_chair_in_R3 : public Function {

  /** 
   * \brief Value of the function at a specified point.
   * @param[in] p The input point. The dimension needs to coincide with the ambient dimension.
   */
  virtual Eigen::VectorXd operator()(const Eigen::VectorXd& p) const override {
    double x = p(0)-off_[0], y = p(1)-off_[1], z = p(2)-off_[2];
    Eigen::VectorXd result(cod_d());
    result(0) = std::pow(x*x + y*y + z*z - a_*k_*k_, 2) - b_*((z-k_)*(z-k_) - 2*x*x)*((z+k_)*(z+k_) - 2*y*y);
    return result;
  }

  /** \brief Returns the domain (ambient) dimension. */
  virtual std::size_t amb_d() const override {return 3;}

  /** \brief Returns the codomain dimension. */
  virtual std::size_t cod_d() const override {return 1;}

  /** \brief Returns a point on the surface. */
  virtual Eigen::VectorXd seed() const override {
    double t1 = a_-b_;
    double discr = t1*t1 - (1.0 - b_)*(a_*a_ - b_);
    double z0 = k_*std::sqrt((t1+std::sqrt(discr))/(1-b_));
    Eigen::Vector3d result(off_[0], off_[1], z0+off_[2]);
    return result;
  }

  /** 
   * \brief Constructor of the function that defines the 'chair' surface
   * embedded in R^3.
   *
   * @param[in] a A numerical parameter.
   * @param[in] b A numerical parameter.
   * @param[in] k A numerical parameter.
   * @param[in] off Offset vector.
   */
  Function_chair_in_R3(double a = 0.8,
		       double b = 0.4,
		       double k = 1.0,
		       Eigen::Vector3d off = Eigen::Vector3d::Zero()) :
    a_(a), b_(b), k_(k), off_(off) {}
  
protected:
  double a_, b_, k_;
  Eigen::Vector3d off_;
};

} // namespace coxeter_triangulation

} // namespace Gudhi

#endif

// (x^2 + y^2 + z^2 - a*k^2)^2 - b*((z-k)^2 - 2*x^2)*((z+k)^2 - 2*y^2)
// sqrt(k/(1-b))*sqrt(a-b + sqrt((a-b)^2 - (1-b)*(a^2 - b)*k^2))

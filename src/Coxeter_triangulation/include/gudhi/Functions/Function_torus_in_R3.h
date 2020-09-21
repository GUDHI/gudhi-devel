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

#ifndef FUNCTIONS_FUNCTION_TORUS_IN_R3_H_
#define FUNCTIONS_FUNCTION_TORUS_IN_R3_H_

#include <gudhi/Functions/Function.h>
#include <Eigen/Dense>

namespace Gudhi {

namespace coxeter_triangulation {

/** 
 * \class Function_torus_in_R3 
 * \brief A class that encodes the function, the zero-set of which is a torus
 * surface embedded in R^3.
 *
 * \ingroup coxeter_triangulation
 */
struct Function_torus_in_R3 : public Function {

  /** 
   * \brief Value of the function at a specified point.
   * @param[in] p The input point. The dimension needs to coincide with the ambient dimension.
   */
  virtual Eigen::VectorXd operator()(const Eigen::VectorXd& p) const override {
    double x = p(0)-off_[0], y = p(1)-off_[1], z = p(2)-off_[2];
    Eigen::VectorXd result(cod_d());
    result(0) = (z*z + (std::sqrt(x*x + y*y) - r_)*(std::sqrt(x*x + y*y) - r_) - R_*R_);
    return result;
  }

  /** \brief Returns the domain (ambient) dimension. */
  virtual std::size_t amb_d() const override {return 3;};

  /** \brief Returns the codomain dimension. */
  virtual std::size_t cod_d() const override {return 1;};

  /** \brief Returns a point on the surface. */
  virtual Eigen::VectorXd seed() const override {
    Eigen::Vector3d result(R_ + r_ +off_[0], off_[1], off_[2]);
    return result;
  }
  
  /** 
   * \brief Constructor of the function that defines a torus embedded in R^3.
   *
   * @param[in] R The outer radius of the torus.
   * @param[in] r The inner radius of the torus.
   * @param[in] off Offset vector.
   */
  Function_torus_in_R3(double R = 1,
		       double r = 0.5,
		       Eigen::Vector3d off = Eigen::Vector3d::Zero()) :
    R_(R), r_(r), off_(off) {}
  
protected:
  double R_, r_;
  Eigen::Vector3d off_;
};

} // namespace coxeter_triangulation

} // namespace Gudhi


#endif

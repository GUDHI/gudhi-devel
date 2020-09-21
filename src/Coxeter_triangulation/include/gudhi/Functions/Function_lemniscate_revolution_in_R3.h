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

#ifndef FUNCTIONS_FUNCTION_LEMNISCATE_REVOLUTION_IN_R3_H_
#define FUNCTIONS_FUNCTION_LEMNISCATE_REVOLUTION_IN_R3_H_

#include <gudhi/Functions/Function.h>
#include <Eigen/Dense>

namespace Gudhi {

namespace coxeter_triangulation {

/** 
 * \class Function_lemniscate_revolution_in_R3 
 * \brief A class that encodes the function, the zero-set of which is a surface of revolution
 * around the x axis based on the lemniscate of Bernoulli embedded in R^3.
 *
 * \ingroup coxeter_triangulation
 */
struct Function_lemniscate_revolution_in_R3 : public Function {

  /** 
   * \brief Value of the function at a specified point.
   * @param[in] p The input point. The dimension needs to coincide with the ambient dimension.
   */
  Eigen::VectorXd operator()(const Eigen::VectorXd& p) const {
    double x = p(0)-off_[0], y = p(1)-off_[1], z = p(2)-off_[2];
    Eigen::VectorXd result(cod_d());
    double x2 = x*x, y2 = y*y, z2 = z*z, a2 = a_*a_;
    double t1 = x2 + y2 + z2;
    result(0) = t1*t1 - 2*a2*(x2 - y2 - z2);
    return result;
  }

  /** \brief Returns the (ambient) domain dimension.*/
  std::size_t amb_d() const {return 3;};

  /** \brief Returns the codomain dimension. */
  std::size_t cod_d() const {return 1;};

  /** \brief Returns a point on the surface. This seed point is only one of 
   * two necessary seed points for the manifold tracing algorithm.
   * See the method seed2() for the other point. 
   */
  Eigen::VectorXd seed() const {
    Eigen::Vector3d result(sqrt(2*a_)+off_[0], off_[1], off_[2]);
    return result;
  }

  /** \brief Returns a point on the surface. This seed point is only one of 
   * two necessary seed points for the manifold tracing algorithm. 
   * See the method seed() for the other point. 
   */
  Eigen::VectorXd seed2() const {
    Eigen::Vector3d result(-sqrt(2*a_)+off_[0], off_[1], off_[2]);
    return result;
  }
  
  /** 
   * \brief Constructor of the function that defines a surface of revolution
   * around the x axis based on the lemniscate of Bernoulli embedded in R^3.
   *
   * @param[in] a A numerical parameter.
   * @param[in] off Offset vector.
   */
  Function_lemniscate_revolution_in_R3(double a = 1,
				       Eigen::Vector3d off = Eigen::Vector3d::Zero())
    : a_(a), off_(off) {}

protected:
  double a_;
  Eigen::Vector3d off_;
};

} // namespace coxeter_triangulation

} // namespace Gudhi


#endif

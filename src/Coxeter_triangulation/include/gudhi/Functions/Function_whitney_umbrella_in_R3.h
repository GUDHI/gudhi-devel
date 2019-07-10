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

#ifndef FUNCTIONS_FUNCTION_WHITNEY_UMBRELLA_IN_R3_H_
#define FUNCTIONS_FUNCTION_WHITNEY_UMBRELLA_IN_R3_H_

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
struct Function_whitney_umbrella_in_R3 {

  /** 
   * \brief Value of the function at a specified point.
   * @param[in] p The input point. The dimension needs to coincide with the ambient dimension.
   */
  Eigen::VectorXd operator()(const Eigen::VectorXd& p) const {
    double x = p(0)-off_[0], y = p(1)-off_[1], z = p(2)-off_[2];
    Eigen::VectorXd coords(cod_d());
    coords(0) = x*x - y*y*z;
    return coords;
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
    return Eigen::Vector3d(1+off_[0], 1+off_[1], 1+off_[2]);
  }

  /** \brief Returns a point on the surface. This seed point is only one of 
   * two necessary seed points for the manifold tracing algorithm. 
   * See the method seed() for the other point. 
   */
  Eigen::VectorXd seed2() const {
    return Eigen::Vector3d(-1+off_[0], -1+off_[1], 1+off_[2]);
  }
  
  /** 
   * \brief Constructor of the function that defines the Whitney umbrella in R^3.
   *
   * @param[in] off Offset vector.
   */
  Function_whitney_umbrella_in_R3(Eigen::Vector3d off = Eigen::Vector3d::Zero())
    : off_(off) {}
  Eigen::Vector3d off_;
};

} // namespace coxeter_triangulation

} // namespace Gudhi


#endif

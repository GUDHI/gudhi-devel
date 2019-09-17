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

#ifndef FUNCTIONS_FUNCTION_RP2_IN_R4_H_
#define FUNCTIONS_FUNCTION_RP2_IN_R4_H_

#include <gudhi/Functions/Function.h>
#include <Eigen/Dense>

namespace Gudhi {

namespace coxeter_triangulation {

/** 
 * \class Function_RP2_in_R4 
 * \brief A class that encodes the function, the zero-set of which is a real projective plane
 * surface embedded in R^4.
 *
 * \ingroup coxeter_triangulation
 */
struct Function_RP2_in_R4 : public Function {

  /** 
   * \brief Value of the function at a specified point.
   * @param[in] p The input point. The dimension needs to coincide with the ambient dimension.
   */
  void evaluate(const Eigen::VectorXd& p, Eigen::VectorXd& result) const {
    double x = p(0)-off_[0], y = p(1)-off_[1], z = p(2)-off_[2], t = p(3)-off_[3];
    result.resize(cod_d());
    double x2 = x*x, y2 = y*y, z2 = z*z;
    result(0) = x2*y2 + y2*z2 + x2*z2 - x*y*z;
    result(1) = t - R_*R_;
  }

  /** \brief Returns the domain (ambient) dimension. */
  std::size_t amb_d() const {return 4;};

  /** \brief Returns the codomain dimension. */
  std::size_t cod_d() const {return 2;};

  /** \brief Returns a point on the surface. */
  void seed(Eigen::VectorXd& result) const {
    result = Eigen::VectorXd(4);
    for (std::size_t i = 0; i < 4; ++i)
      result(i) = off_[i];
    result(3) += R_*R_;
  }

  /** 
   * \brief Constructor of the function that defines a real projective plane embedded in R^4.
   *
   * @param[in] R The radius of the reference two-dimensional sphere.
   */
  Function_RP2_in_R4(double R = 1) :
    R_(R), off_(Eigen::VectorXd::Zero(4)) {}

  
  /** 
   * \brief Constructor of the function that defines a real projective plane embedded in R^4.
   *
   * @param[in] R The radius of the reference two-dimensional sphere.
   * @param[in] off Offset vector.
   */
  Function_RP2_in_R4(double R,
		     Eigen::VectorXd off) :
    R_(R), off_(off) {}
  
protected:
  double R_;
  Eigen::VectorXd off_;
};

} // namespace coxeter_triangulation

} // namespace Gudhi


#endif

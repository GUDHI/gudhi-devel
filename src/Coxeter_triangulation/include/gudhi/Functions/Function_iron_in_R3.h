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

#ifndef FUNCTIONS_FUNCTION_IRON_IN_R3_H_
#define FUNCTIONS_FUNCTION_IRON_IN_R3_H_

#include <gudhi/Functions/Function.h>
#include <Eigen/Dense>

namespace Gudhi {

namespace coxeter_triangulation {

/** 
 * \class Function_iron_in_R3 
 * \brief A class that encodes the function, the zero-set of which is a surface
 * embedded in R^3 that ressembles an iron.
 *
 * \ingroup coxeter_triangulation
 */
struct Function_iron_in_R3 : public Function {

  /** 
   * \brief Value of the function at a specified point.
   * @param[in] p The input point. The dimension needs to coincide with the ambient dimension.
   */
  void evaluate(const Eigen::VectorXd& p, Eigen::VectorXd& result) const {
    double x = p(0), y = p(1), z = p(2);
    result.resize(cod_d());
    result(0) = - std::pow(x,6)/300. - std::pow(y,6)/300. - std::pow(z,6)/300. + x*y*y*z/2.1 + y*y + std::pow(z-2, 4) - 1;
  }

  /** \brief Returns the domain (ambient) dimension. */
  std::size_t amb_d() const {return 3;};

  /** \brief Returns the codomain dimension. */
  std::size_t cod_d() const {return 1;};

  /** \brief Returns a point on the surface. */
  void seed(Eigen::VectorXd& result) const {
    result = Eigen::Vector3d(std::pow(4500, 1./6), 0, 0);
  }
  
  /** 
   * \brief Constructor of the function that defines a surface embedded in R^3
   * that ressembles an iron.
   *
   * @param[in] off Offset vector.
   */
  Function_iron_in_R3(Eigen::Vector3d off = Eigen::Vector3d::Zero()) :
    off_(off) {}
  
protected:
  Eigen::Vector3d off_;
};

} // namespace coxeter_triangulation

} // namespace Gudhi


#endif

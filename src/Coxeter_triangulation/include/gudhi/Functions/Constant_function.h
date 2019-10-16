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

#ifndef FUNCTIONS_CONSTANT_FUNCTION_H_
#define FUNCTIONS_CONSTANT_FUNCTION_H_

#include <gudhi/Functions/Function.h>
#include <Eigen/Dense>

namespace Gudhi {

namespace coxeter_triangulation {

/** 
 * \class Constant_function 
 * \brief A class that encodes a constant function from R^d to R^k.
 * This class does not have any implicit manifold in correspondence.
 *
 * \ingroup coxeter_triangulation
 */
struct Constant_function : public Function {

  /** \brief Value of the function at a specified point. The value is constant.
   * @param[in] p The input point. The dimension needs to coincide with the ambient dimension.
   */
  Eigen::VectorXd operator()(const Eigen::VectorXd& p) const {
    Eigen::VectorXd result = value_;
    return result;
  }
  
  /** \brief Returns the domain dimension. Same as the ambient dimension of the sphere. */
  std::size_t amb_d() const {return d_;};

  /** \brief Returns the codomain dimension. Same as the codimension of the sphere. */
  std::size_t cod_d() const {return k_;};

  /** \brief No seed point is available. Throws an exception on evocation. */
  Eigen::VectorXd seed() const {
    throw "Seed invoked on a constant function.\n";
  }

  Constant_function() {}
  
  /** 
   * \brief Constructor of a constant function from R^d to R^m.
   *
   * @param[in] d The domain dimension.
   * @param[in] k The codomain dimension.
   * @param[in] value The constant value of the function.
   */
  Constant_function(std::size_t d, std::size_t k, const Eigen::VectorXd& value)
    : d_(d), k_(k), value_(value) {}
  
  std::size_t d_, k_;
  Eigen::VectorXd value_;
};

} // namespace coxeter_triangulation

} // namespace Gudhi

#endif

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

#ifndef FUNCTIONS_FUNCTION_H_
#define FUNCTIONS_FUNCTION_H_

#include <Eigen/Dense>

namespace Gudhi {

namespace coxeter_triangulation {

/** 
 * \class Function
 * \brief The parent class for all functions implemented in the module.
 *  Contains virtual methods needed to be a model of the concept FunctionForImplicitManifold.
 *
 * \ingroup coxeter_triangulation
 */
struct Function {
  
  /** \brief Virtual method for the value of the function at a specified point.
   * @param[in] p The input point.
   */
  virtual Eigen::VectorXd operator()(const Eigen::VectorXd& p) const {
    return Eigen::VectorXd();
  }
  
  /** \brief Virtual method for the domain dimension. */
  virtual std::size_t amb_d() const {return 0;};

  /** \brief Virtual method for the codomain dimension. */
  virtual std::size_t cod_d() const {return 0;};

  /** \brief Virtual method for the seed point. */
  virtual Eigen::VectorXd seed() const {
    return Eigen::VectorXd();
  }

  virtual Function* clone() const = 0;

  /** \brief Virtual destructor. */
  virtual ~Function() {}
};

} // namespace coxeter_triangulation

} // namespace Gudhi


#endif

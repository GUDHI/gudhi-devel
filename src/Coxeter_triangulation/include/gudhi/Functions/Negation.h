/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef FUNCTIONS_NEGATION_H_
#define FUNCTIONS_NEGATION_H_

#include <cstdlib>  // for std::size_t

#include <gudhi/Functions/Function.h>

#include <Eigen/Dense>

namespace Gudhi {

namespace coxeter_triangulation {

/** 
 *\class Negation
 * \brief Constructs the "minus" function. The zero-set is the same, but
 * the values at other points are the negative of their original value.
 *
 * \tparam Function_ The function template parameter. Should be a model of 
 * the concept FunctionForImplicitManifold.
 *
 * \ingroup coxeter_triangulation
 */
template <class Function_>
struct Negation : public Function {
  
  /** 
   * \brief Value of the function at a specified point.
   * @param[in] p The input point. The dimension needs to coincide with the ambient dimension.
   */
  virtual Eigen::VectorXd operator()(const Eigen::VectorXd& p) const override {
    Eigen::VectorXd result = -fun_(p);
    return result;
  }

  /** \brief Returns the domain (ambient) dimension. */
  virtual std::size_t amb_d() const override {return fun_.amb_d();}

  /** \brief Returns the codomain dimension. */
  virtual std::size_t cod_d() const override {return fun_.cod_d();}

  /** \brief Returns a point on the zero-set. */
  virtual Eigen::VectorXd seed() const override {
    Eigen::VectorXd result = fun_.seed();
    return result;
  }

  /** 
   * \brief Constructor of the negative function.
   *
   * @param[in] function The function to be negated.
   */
  Negation(const Function_& function) :
    fun_(function) {
  }

 private:
  Function_ fun_;
};


/** 
 * \brief Static constructor of the negative function.
 *
 * @param[in] function The function to be translated.
 * domain (ambient) dimension of 'function'.
 *
 * \tparam Function_ The function template parameter. Should be a model of 
 * the concept FunctionForImplicitManifold.
 *
 * \ingroup coxeter_triangulation
 */
template <class Function_>
Negation<Function_> negation(const Function_& function) {
  return Negation<Function_>(function);
}

} // namespace coxeter_triangulation

} // namespace Gudhi

#endif

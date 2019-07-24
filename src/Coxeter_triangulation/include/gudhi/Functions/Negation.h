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

#ifndef FUNCTIONS_NEGATION_H_
#define FUNCTIONS_NEGATION_H_

#include <cstdlib>

#include <gudhi/Functions/Function.h>
#include <Eigen/Dense>

namespace Gudhi {

namespace coxeter_triangulation {

/* \class Negation
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
  void operator()(const Eigen::VectorXd& p, Eigen::VectorXd& result) const {
    fun_.evaluate(p);
    result = -result;
  }

  /** \brief Returns the domain (ambient) dimension. */
  std::size_t amb_d() const {return fun_.amb_d();}

  /** \brief Returns the codomain dimension. */
  std::size_t cod_d() const {return fun_.cod_d();}

  /** \brief Returns a point on the zero-set. */
  void seed(Eigen::VectorXd& result) const {
    fun_.seed(result);
  }

  /** 
   * \brief Constructor of the negative function.
   *
   * @param[in] function The function to be negated.
   */
  Negation(const Function_& function) :
    fun_(function) {
  }
  Function_ fun_;
};


/** 
 * \brief Static constructor of the negative function.
 *
 * @param[in] function The function to be translated.
 * @param[in] off The offset vector. The dimension should correspond to the 
 * domain (ambient) dimension of 'function'.
 *
 * \tparam Function_ The function template parameter. Should be a model of 
 * the concept FunctionForImplicitManifold.
 */
template <class Function_>
Negation<Function_>  negation(const Function_& function) {
  return Negation<Function_>(function);
}

} // namespace coxeter_triangulation

} // namespace Gudhi

#endif

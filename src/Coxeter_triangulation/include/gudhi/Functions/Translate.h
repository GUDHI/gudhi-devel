/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef FUNCTIONS_TRANSLATE_H_
#define FUNCTIONS_TRANSLATE_H_

#include <cstdlib>  // for std::size_t

#include <Eigen/Dense>

namespace Gudhi {

namespace coxeter_triangulation {

/**
 * \class Translate
 * \brief Translates the zero-set of the function by a vector.
 * The underlying function corresponds to f(x-off), where off is the offset vector.
 *
 * \tparam Function_ The function template parameter. Should be a model of
 * the concept FunctionForImplicitManifold.
 */
template <class Function_>
struct Translate {
  /**
   * \brief Value of the function at a specified point.
   * @param[in] p The input point. The dimension needs to coincide with the ambient dimension.
   */
  Eigen::VectorXd operator()(const Eigen::VectorXd& p) const {
    Eigen::VectorXd result = fun_(p - off_);
    return result;
  }

  /** \brief Returns the domain (ambient) dimension. */
  std::size_t amb_d() const { return fun_.amb_d(); }

  /** \brief Returns the codomain dimension. */
  std::size_t cod_d() const { return fun_.cod_d(); }

  /** \brief Returns a point on the zero-set. */
  Eigen::VectorXd seed() const {
    Eigen::VectorXd result = fun_.seed();
    result += off_;
    return result;
  }

  /**
   * \brief Constructor of the translated function.
   *
   * @param[in] function The function to be translated.
   * @param[in] off The offset vector. The dimension should correspond to the
   * domain (ambient) dimension of 'function'.
   */
  Translate(const Function_& function, const Eigen::VectorXd& off) : fun_(function), off_(off) {}

 private:
  Function_ fun_;
  Eigen::VectorXd off_;
};

/**
 * \brief Static constructor of a translated function.
 *
 * @param[in] function The function to be translated.
 * @param[in] off The offset vector. The dimension should correspond to the
 * domain (ambient) dimension of 'function'.
 *
 * \tparam Function_ The function template parameter. Should be a model of
 * the concept FunctionForImplicitManifold.
 *
 * \ingroup coxeter_triangulation
 */
template <class Function_>
Translate<Function_> translate(const Function_& function, Eigen::VectorXd off) {
  return Translate<Function_>(function, off);
}

}  // namespace coxeter_triangulation

}  // namespace Gudhi

#endif

/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef FUNCTIONS_EMBED_IN_RD_H_
#define FUNCTIONS_EMBED_IN_RD_H_

#include <cstdlib>  // for std::size_t

#include <gudhi/Functions/Function.h>

#include <Eigen/Dense>


namespace Gudhi {

namespace coxeter_triangulation {

/**
 * \class Embed_in_Rd
 * \brief Embedding of an implicit manifold in a higher dimension.
 *
 * \tparam Function_ The function template parameter. Should be a model of 
 * the concept FunctionForImplicitManifold.
 *
 * \ingroup coxeter_triangulation
 */
template <class Function_>
struct Embed_in_Rd : public Function {
  
  /** 
   * \brief Value of the function at a specified point.
   * @param[in] p The input point. The dimension needs to coincide with the ambient dimension.
   */
  virtual Eigen::VectorXd operator()(const Eigen::VectorXd& p) const override {
    Eigen::VectorXd x = p;
    Eigen::VectorXd x_k(fun_.amb_d()), x_rest(d_ - fun_.amb_d());
    for (std::size_t i = 0; i < fun_.amb_d(); ++i)
      x_k(i) = x(i);
    for (std::size_t i = fun_.amb_d(); i < d_; ++i)
      x_rest(i - fun_.amb_d()) = x(i);
    Eigen::VectorXd result = fun_(x_k);
    result.conservativeResize(this->cod_d());
    for (std::size_t i = fun_.cod_d(); i < this->cod_d(); ++i)
      result(i) = x_rest(i - fun_.cod_d());
    return result;
  }

  /** \brief Returns the domain (ambient) dimension. */
  virtual std::size_t amb_d() const override {return d_;}

  /** \brief Returns the codomain dimension. */
  virtual std::size_t cod_d() const override {return d_-(fun_.amb_d() - fun_.cod_d());}

  /** \brief Returns a point on the zero-set of the embedded function. */
  virtual Eigen::VectorXd seed() const override {
    Eigen::VectorXd result = fun_.seed();
    result.conservativeResize(d_);
    for (std::size_t l = fun_.amb_d(); l < d_; ++l)
      result(l) = 0;
    return result;
  }

  /** 
   * \brief Constructor of the embedding function.
   *
   * @param[in] function The function to be embedded in higher dimension.
   * @param[in] d Embedding dimension.
   */
  Embed_in_Rd(const Function_& function, std::size_t d) :
    fun_(function), d_(d) {
  }

 private:
  Function_ fun_;
  std::size_t d_;
};


/** 
 * \brief Static constructor of an embedding function.
 *
 * @param[in] function The function to be embedded in higher dimension.
 * @param[in] d Embedding dimension.
 *
 * \tparam Function_ The function template parameter. Should be a model of 
 * the concept FunctionForImplicitManifold.
 *
 * \ingroup coxeter_triangulation
 */
template <class Function_>
Embed_in_Rd<Function_> make_embedding(const Function_& function, std::size_t d) {
  return Embed_in_Rd<Function_>(function, d);
}

} // namespace coxeter_triangulation

} // namespace Gudhi

#endif

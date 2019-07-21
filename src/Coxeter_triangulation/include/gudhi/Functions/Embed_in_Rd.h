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

#ifndef FUNCTIONS_EMBED_IN_RD_H_
#define FUNCTIONS_EMBED_IN_RD_H_

#include <cstdlib>
#include <random>

#include <gudhi/Functions/Function.h>
#include <Eigen/Dense>


namespace Gudhi {

namespace coxeter_triangulation {

/* \class Embed_in_Rd
 * \brief Embedding of an implicit manifold in a higher dimension.
 *
 * \tparam Function The function template parameter. Should be a model of 
 * the concept FunctionForImplicitManifold.
 *
 * \ingroup coxeter_triangulation
 */
template <class Function>
struct Embed_in_Rd : public Function {
  
  /** 
   * \brief Value of the function at a specified point.
   * @param[in] p The input point. The dimension needs to coincide with the ambient dimension.
   */
  Eigen::VectorXd operator()(const Eigen::VectorXd& p) const {
    Eigen::VectorXd x = p;
    Eigen::VectorXd x_k(fun_.amb_d()), x_rest(d_ - fun_.amb_d());
    for (std::size_t i = 0; i < fun_.amb_d(); ++i)
      x_k(i) = x(i);
    for (std::size_t i = fun_.amb_d(); i < d_; ++i)
      x_rest(i - fun_.amb_d()) = x(i);
    Eigen::VectorXd res = fun_(x_k);
    res.conservativeResize(this->cod_d());
    for (std::size_t i = fun_.cod_d(); i < this->cod_d(); ++i)
      res(i) = x_rest(i - fun_.cod_d());
    return res;
  }

  /** \brief Returns the domain (ambient) dimension. */
  std::size_t amb_d() const {return d_;}

  /** \brief Returns the codomain dimension. */
  std::size_t cod_d() const {return d_-(fun_.amb_d() - fun_.cod_d());}

  /** \brief Returns a point on the zero-set of the embedded function. */
  Eigen::VectorXd seed() const {
    Eigen::VectorXd seed_k = fun_.seed();
    seed_k.conservativeResize(d_);
    for (std::size_t l = fun_.seed().size(); l < d_; ++l)
      seed_k(l) = 0;
    return seed_k;
  }

  /** 
   * \brief Constructor of the embedding function.
   *
   * @param[in] function The function to be embedded in higher dimension.
   * @param[in] d Embedding dimension.
   */
  Embed_in_Rd(const Function& function, std::size_t d) :
    fun_(function), d_(d) {
  }
  Function fun_;
  std::size_t d_;
};


/** 
 * \brief Static constructor of an embedding function.
 *
 * @param[in] function The function to be embedded in higher dimension.
 * @param[in] d Embedding dimension.
 *
 * \tparam Function The function template parameter. Should be a model of 
 * the concept FunctionForImplicitManifold.
 */
template <class Function>
Embed_in_Rd<Function> make_embedding(const Function& function, std::size_t d) {
  return Embed_in_Rd<Function>(function, d);
}

} // namespace coxeter_triangulation

} // namespace Gudhi

#endif

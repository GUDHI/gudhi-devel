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

#ifndef FUNCTIONS_CARTESIAN_PRODUCT_H_
#define FUNCTIONS_CARTESIAN_PRODUCT_H_

#include <cstdlib>
#include <tuple>
#include <utility>

#include <gudhi/Functions/Function.h>
#include <Eigen/Dense>

namespace Gudhi {

namespace coxeter_triangulation {

/* Get the domain dimension of the tuple of functions.
 */
template <std::size_t I = 0, typename... T>
inline typename std::enable_if<I == sizeof... (T), std::size_t>::type
get_amb_d (const std::tuple<T...>& tuple) {
  return 0;
}

template <std::size_t I = 0, typename... T>
inline typename std::enable_if<I != sizeof... (T), std::size_t>::type
get_amb_d (const std::tuple<T...>& tuple) {
  return std::get<I>(tuple).amb_d() + get_amb_d<I+1, T...>(tuple);
}

/* Get the codomain dimension of the tuple of functions.
 */
template <std::size_t I = 0, typename... T>
inline typename std::enable_if<I == sizeof... (T), std::size_t>::type
get_cod_d (const std::tuple<T...>& tuple) {
  return 0;
}

template <std::size_t I = 0, typename... T>
inline typename std::enable_if<I != sizeof... (T), std::size_t>::type
get_cod_d (const std::tuple<T...>& tuple) {
  return std::get<I>(tuple).cod_d() + get_cod_d<I+1, T...>(tuple);
}

/* Get the seed of the tuple of functions.
 */
template <std::size_t I = 0, typename... T>
inline typename std::enable_if<I == sizeof... (T), void>::type
get_seed (const std::tuple<T...>& tuple, Eigen::VectorXd& point, std::size_t i = 0) {
}

template <std::size_t I = 0, typename... T>
inline typename std::enable_if<I != sizeof... (T), void>::type
get_seed (const std::tuple<T...>& tuple, Eigen::VectorXd& point, std::size_t i = 0) {
  const auto& f = std::get<I>(tuple);
  std::size_t n = f.amb_d();
  Eigen::VectorXd seed = f.seed();
  for (std::size_t j = 0; j < n; ++j)
    point(i+j) = seed(j);
  get_seed<I+1, T...>(tuple, point, i+n);  
}

/* Get the seed of the tuple of functions.
 */
template <std::size_t I = 0, typename... T>
inline typename std::enable_if<I == sizeof... (T), void>::type
get_value (const std::tuple<T...>& tuple, const Eigen::VectorXd& x, Eigen::VectorXd& point, std::size_t i = 0, std::size_t j = 0) {
}

template <std::size_t I = 0, typename... T>
inline typename std::enable_if<I != sizeof... (T), void>::type
get_value (const std::tuple<T...>& tuple, const Eigen::VectorXd& x, Eigen::VectorXd& point, std::size_t i = 0, std::size_t j = 0) {
  const auto& f = std::get<I>(tuple);
  std::size_t n = f.amb_d();
  std::size_t k = f.cod_d();
  Eigen::VectorXd x_i(n);
  for (std::size_t l = 0; l < n; ++l)
    x_i(l) = x(i+l);
  Eigen::VectorXd res = f(x_i);
  for (std::size_t l = 0; l < k; ++l)
    point(j+l) = res(l);
  get_value<I+1, T...>(tuple, x, point, i+n, j+k);  
}


/* \class Cartesian_product
 * \brief Constructs the function the zero-set of which is the Cartesian product
 * of the zero-sets of some given functions.
 *
 * \tparam Functions A pack template parameter for functions. All functions should be models of 
 * the concept FunctionForImplicitManifold.
 *
 * \ingroup coxeter_triangulation
 */
template <class... Functions>
struct Cartesian_product : public Function {
  
  /** 
   * \brief Value of the function at a specified point.
   * @param[in] p The input point. The dimension needs to coincide with the ambient dimension.
   */
  Eigen::VectorXd operator()(const Eigen::VectorXd& p) const {
    Eigen::VectorXd result(cod_d_);
    get_value(function_tuple_, p, result, 0, 0);
    return result;
  }

  /** \brief Returns the domain (ambient) dimension. */
  std::size_t amb_d() const {return amb_d_;}

  /** \brief Returns the codomain dimension. */
  std::size_t cod_d() const {return cod_d_;}

  /** \brief Returns a point on the zero-set. */
  Eigen::VectorXd seed() const {
    Eigen::VectorXd result(amb_d_);
    get_seed(function_tuple_, result, 0);
    return result;
  }
  
  /** 
   * \brief Constructor of the Cartesian product function.
   *
   * @param[in] functions The functions the zero-sets of which are factors in the
   * Cartesian product of the resulting function.
   */
  Cartesian_product(const Functions&... functions)
    : function_tuple_(std::make_tuple(functions...)) {
    amb_d_ = get_amb_d(function_tuple_);
    cod_d_ = get_cod_d(function_tuple_);
  }

private:
  std::tuple<Functions...> function_tuple_;
  std::size_t amb_d_, cod_d_;
};


/** 
 * \brief Static constructor of a Cartesian product function.
 *
 * @param[in] functions The functions the zero-sets of which are factors in the
 * Cartesian product of the resulting function.
 *
 * \tparam Functions A pack template parameter for functions. All functions should be models of 
 * the concept FunctionForImplicitManifold.
 */
template <typename... Functions>
Cartesian_product<Functions...> make_product_function(const Functions&... functions) {
  return Cartesian_product<Functions...>(functions...);
}

template <typename... Functions>
Cartesian_product<Functions...>* new_product_function(const Functions&... functions) {
  return new Cartesian_product<Functions...>(functions...);
}


} // namespace coxeter_triangulation

} // namespace Gudhi

#endif

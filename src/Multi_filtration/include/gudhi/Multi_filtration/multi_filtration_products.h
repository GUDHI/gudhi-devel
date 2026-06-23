/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber, David Loiseaux
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file multi_filtration_products.h
 * @author Hannah Schreiber, David Loiseaux
 */

#ifndef MF_PRODUCTS_H_
#define MF_PRODUCTS_H_

#include <algorithm>
#include <cstddef>
#include <cmath>
#include <initializer_list>
#include <vector>

#ifdef GUDHI_USE_TBB
#include <oneapi/tbb/parallel_for.h>
#endif

#include <gudhi/Multi_filtration/multi_filtration_utils.h>

namespace Gudhi {

namespace multi_filtration {

/**
 * @ingroup multi_filtration
 *
 * @private
 */
template <typename U>
inline bool _is_less(U a, U b)
{
  // workaround -Ofast optimization which is default on Windows, where x < NaN / NaN < x is not well defined
  if (_is_nan(a)) return !_is_nan(b);
  if (_is_nan(b)) return false;
  return a < b;
};

/**
 * @ingroup multi_filtration
 *
 * @brief For the first argument, computes the smallest (resp. the greatest if `MultiFiltrationValue::Co` is true)
 * scalar product of all its generators with the given range.
 *
 * @tparam U Arithmetic type of the result.
 * @tparam MultiFiltrationValue Multi filtration value type with methods num_parameters(), num_generators() and
 * operator(g, p), as well as has_negative_cones() and type definition `size_type`. For example:
 * @ref Gudhi::multi_filtration::Multi_parameter_filtration,
 * @ref Gudhi::multi_filtration::Dynamic_multi_parameter_filtration or
 * @ref Gudhi::multi_filtration::Degree_rips_bifiltration.
 * @tparam CoefficientRange Range of values convertible into `U` with a begin() and size() method.
 * @param f Multi filtration value.
 * @param x Vector of coefficients.
 * @return Scalar product of @p f with @p x of type `U`.
 */
template <typename U, class MultiFiltrationValue, class CoefficientRange = std::initializer_list<U> >
inline U compute_linear_projection(const MultiFiltrationValue &f, const CoefficientRange &x)
{
  using size_type = typename MultiFiltrationValue::size_type;

  auto project_generator = [&](size_type g) -> U {
    U projection = 0;
    std::size_t size = std::min(x.size(), f.num_parameters());
    auto it = x.begin();
    for (std::size_t i = 0; i < size; i++) {
      // workaround -Ofast optimization which is default on Windows
      // otherwise 0 += NaN can be equal to 0 with -Ofast which we don't want
      if (_is_nan(f(g, i))) return f(g, i);
      projection += static_cast<U>(*it) * static_cast<U>(f(g, i));
      ++it;
    }
    return projection;
  };

  if (f.num_generators() == 1) return project_generator(0);

  auto is_less = [](U a, U b) { return _is_less(a, b); };

#ifdef GUDHI_USE_TBB
  std::vector<U> projections(f.num_generators());
  tbb::parallel_for(size_type{0}, f.num_generators(), [&](size_type g) { projections[g] = project_generator(g); });
  if constexpr (MultiFiltrationValue::has_negative_cones()) {
    return *std::max_element(projections.begin(), projections.end(), is_less);
  } else {
    return *std::min_element(projections.begin(), projections.end(), is_less);
  }
#else
  if constexpr (MultiFiltrationValue::has_negative_cones()) {
    U projection = std::numeric_limits<U>::lowest();
    for (size_type g = 0; g < f.num_generators(); ++g) {
      // Order in the max important to spread possible NaNs
      projection = std::max(project_generator(g), projection, is_less);
    }
    return projection;
  } else {
    U projection = std::numeric_limits<U>::max();
    for (size_type g = 0; g < f.num_generators(); ++g) {
      // Order in the min important to spread possible NaNs
      projection = std::min(project_generator(g), projection, is_less);
    }
    return projection;
  }
#endif
}

// enables compute_linear_projection<U>(...) as well as compute_linear_projection(...) with default value
// MultiFiltrationValue::value_type and not a fixed type like double
/**
 * @ingroup multi_filtration
 *
 * @brief For the first argument, computes the smallest (resp. the greatest if `MultiFiltrationValue::Co` is true)
 * scalar product of all its generators with the given range.
 *
 * @tparam MultiFiltrationValue Multi filtration value type with methods num_parameters(), num_generators() and
 * operator(g, p), as well as has_negative_cones() and type definition `size_type`. For example:
 * @ref Gudhi::multi_filtration::Multi_parameter_filtration,
 * @ref Gudhi::multi_filtration::Dynamic_multi_parameter_filtration or
 * @ref Gudhi::multi_filtration::Degree_rips_bifiltration.
 * @tparam CoefficientRange Range of values convertible into `U` with a begin() and size() method.
 * @param f Multi filtration value.
 * @param x Vector of coefficients.
 * @return Scalar product of @p f with @p x of type `MultiFiltrationValue::value_type`.
 */
template <class MultiFiltrationValue,
          class CoefficientRange = std::initializer_list<typename MultiFiltrationValue::value_type> >
inline auto compute_linear_projection(const MultiFiltrationValue &f, const CoefficientRange &x)
{
  return compute_linear_projection<typename MultiFiltrationValue::value_type, MultiFiltrationValue, CoefficientRange>(
      f, x);
}

/**
 * @ingroup multi_filtration
 *
 * @private
 */
template <typename T>
inline T _sqrt(T v)
{
  if constexpr (std::is_integral_v<T>) {
    // to avoid Windows issue that don't know how to cast integers for cmath methods
    return std::sqrt(static_cast<double>(v));
  } else {
    return std::sqrt(v);
  }
}

/**
 * @ingroup multi_filtration
 *
 * @private
 */
template <typename T, class F, typename size_type>
inline T _compute_frobenius_norm(size_type number_of_elements, F &&norm)
{
  if (number_of_elements == 1) return std::forward<F>(norm)(0);

  T out = 0;
  for (size_type p = 0; p < number_of_elements; ++p) {
    T v = std::forward<F>(norm)(p);
    // workaround -Ofast optimization which is default on Windows
    // otherwise 0 += NaN can be equal to 0 with -Ofast which we don't want
    if (_is_nan(v)) return v;
    out += v * v;
  }
  return out;
}

/**
 * @ingroup multi_filtration
 *
 * @brief Computes the euclidean distance from the first parameter to the second parameter as the minimum of
 * all Euclidean distances between a generator of @p f1 and a generator of @p f2.
 *
 * @tparam U Arithmetic type of the result.
 * @tparam MultiFiltrationValue Multi filtration value type with methods num_parameters(), num_generators() and
 * operator(g, p), as well as ensures_1_criticality() and type definitions `size_type` and `value_type`. For example:
 * @ref Gudhi::multi_filtration::Multi_parameter_filtration,
 * @ref Gudhi::multi_filtration::Dynamic_multi_parameter_filtration or
 * @ref Gudhi::multi_filtration::Degree_rips_bifiltration.
 * @param f1 Source filtration value.
 * @param f2 Target filtration value.
 * @return Euclidean distance between @p f1 and @p f2 with type `U`.
 */
template <typename U, class MultiFiltrationValue>
inline U compute_euclidean_distance_to(const MultiFiltrationValue &f1, const MultiFiltrationValue &f2)
{
  using size_type = typename MultiFiltrationValue::size_type;
  using T = typename MultiFiltrationValue::value_type;

  GUDHI_CHECK(f1.num_parameters() == f2.num_parameters(),
              std::invalid_argument("We cannot compute the distance between two points of different dimensions."));

  // TODO: verify if this really makes a differences in the 1-critical case, otherwise just keep the general case
  if constexpr (MultiFiltrationValue::ensures_1_criticality()) {
    return _sqrt<U>(
        _compute_frobenius_norm<T>(f1.num_parameters(), [&](size_type p) -> T { return f1(0, p) - f2(0, p); }));
  } else {
    auto is_less = [](U a, U b) { return _is_less(a, b); };
    U res = std::numeric_limits<U>::max();
    for (size_type g1 = 0; g1 < f1.num_generators(); ++g1) {
      for (size_type g2 = 0; g2 < f2.num_generators(); ++g2) {
        // Euclidean distance as a Frobenius norm with matrix 1 x p and values 'f(g1, p) - other(g2, p)'
        // Order in the min important to spread possible NaNs
        res = std::min(static_cast<U>(_compute_frobenius_norm<T>(
                           f1.num_parameters(), [&](size_type p) -> T { return f1(g1, p) - f2(g2, p); })),
                       res,
                       is_less);
      }
    }
    return _sqrt(res);
  }
}

// enables compute_euclidean_distance_to<U>(...) as well as compute_euclidean_distance_to(...) with default value
// MultiFiltrationValue::value_type and not a fixed type like double
/**
 * @ingroup multi_filtration
 *
 * @brief Computes the euclidean distance from the first parameter to the second parameter as the minimum of
 * all Euclidean distances between a generator of @p f1 and a generator of @p f2.
 *
 * @tparam MultiFiltrationValue Multi filtration value type with methods num_parameters(), num_generators() and
 * operator(g, p), as well as ensures_1_criticality() and type definitions `size_type` and `value_type`. For example:
 * @ref Gudhi::multi_filtration::Multi_parameter_filtration,
 * @ref Gudhi::multi_filtration::Dynamic_multi_parameter_filtration or
 * @ref Gudhi::multi_filtration::Degree_rips_bifiltration.
 * @param f1 Source filtration value.
 * @param f2 Target filtration value.
 * @return Euclidean distance between @p f1 and @p f2 with type `MultiFiltrationValue::value_type`.
 */
template <class MultiFiltrationValue>
inline auto compute_euclidean_distance_to(const MultiFiltrationValue &f1, const MultiFiltrationValue &f2)
{
  return compute_euclidean_distance_to<typename MultiFiltrationValue::value_type, MultiFiltrationValue>(f1, f2);
}

/**
 * @ingroup multi_filtration
 *
 * @brief Computes the norm of the given filtration value.
 *
 * The filtration value is seen as a \f$ num_generators x num_parameters \f$ matrix and a standard Frobenius norm
 * is computed from it: the square root of the sum of the squares of all elements in the matrix.
 *
 * @tparam U Arithmetic type of the result.
 * @tparam MultiFiltrationValue Multi filtration value type with methods num_parameters(), num_generators() and
 * operator(g, p), as well as type definitions `size_type` and `value_type`. For example:
 * @ref Gudhi::multi_filtration::Multi_parameter_filtration,
 * @ref Gudhi::multi_filtration::Dynamic_multi_parameter_filtration or
 * @ref Gudhi::multi_filtration::Degree_rips_bifiltration.
 * @param f Multi filtration value.
 * @return The norm of @p f of type `U`.
 */
template <typename U, class MultiFiltrationValue>
inline U compute_norm(const MultiFiltrationValue &f)
{
  using size_type = typename MultiFiltrationValue::size_type;
  using T = typename MultiFiltrationValue::value_type;

  // Frobenius norm with matrix g x p based on Euclidean norm
  U out = 0;
  for (size_type g = 0; g < f.num_generators(); ++g) {
    auto v = _compute_frobenius_norm<T>(f.num_parameters(), [&](size_type p) -> T { return f(g, p); });
    // workaround -Ofast optimization which is default on Windows
    // otherwise 0 += NaN can be equal to 0 with -Ofast which we don't want
    if (_is_nan(v)) return v;
    out += v;
  }
  return _sqrt(out);
}

// enables compute_norm<U>(...) as well as compute_norm(...) with default value
// MultiFiltrationValue::value_type and not a fixed type like double
/**
 * @ingroup multi_filtration
 *
 * @brief Computes the norm of the given filtration value.
 *
 * The filtration value is seen as a \f$ num_generators x num_parameters \f$ matrix and a standard Frobenius norm
 * is computed from it: the square root of the sum of the squares of all elements in the matrix.
 *
 * @tparam MultiFiltrationValue Multi filtration value type with methods num_parameters(), num_generators() and
 * operator(g, p), as well as type definitions `size_type` and `value_type`. For example:
 * @ref Gudhi::multi_filtration::Multi_parameter_filtration,
 * @ref Gudhi::multi_filtration::Dynamic_multi_parameter_filtration or
 * @ref Gudhi::multi_filtration::Degree_rips_bifiltration.
 * @param f Multi filtration value.
 * @return The norm of @p f of type `MultiFiltrationValue::value_type`.
 */
template <class MultiFiltrationValue>
inline auto compute_norm(const MultiFiltrationValue &f)
{
  return compute_norm<typename MultiFiltrationValue::value_type, MultiFiltrationValue>(f);
}

}  // namespace multi_filtration

}  // namespace Gudhi

#endif  // MF_PRODUCTS_H_

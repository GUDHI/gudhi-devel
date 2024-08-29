/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Loiseaux
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - 2024/08 Hannah Schreiber: Generalization to all signed arithmetic types for T
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef ONE_CRITICAL_FILTRATIONS_H_
#define ONE_CRITICAL_FILTRATIONS_H_

#include <algorithm>  //std::lower_bound
#include <cmath>      //std::isnan
#include <cstddef>    //std::size_t
#include <cstdint>    //std::int32_t
#include <ostream>    //std::ostream
#include <limits>     //std::numerical_limits
#include <stdexcept>  //std::logic_error
#include <vector>

#include <gudhi/Debug_utils.h>

namespace Gudhi::multi_filtration {

/**
 * \brief Vector-like filtration value, for multiparameter persistence, with
 * numpy-like methods.
 *
 * \ingroup multi_filtration
 *
 * \details Child of `std::vector<T>` that has numpy-like pointwise operators.
 * One critical simplicial filtrations are filtrations such that the lifetime
 * each simplex is a positive cone, e.g.
 *  - \f$ \{ x \in  \mathbb R^2 : x>=(1,2)\} \f$ is valid, while
 *  - \f$ \{ x \in  \mathbb R^2 : x>=(1,2)\} \cap \{x \in  \mathbb R^2 :
 * x>=(2,1)\} \f$ is not- Finitely critical filtrations are filtrations such
 * that the lifetime of each simplex is a union of such cones, e.g.,
 *  - \f$ \{ x \in  \mathbb R^2 : x>=(1,2)\} \cap \{ x \in  \mathbb R^2 :
 * x>=(2,1)\} \f$ is finitely critical, and more particularly 2-critical, while
 *  - \f$ \{ x \in  \mathbb R^2 : x>= \mathrm{epigraph}(y\mapsto e^{-y})\} \f$
 * is not.
 * Use \ref Multi_critical_filtration for multi-critical filtrations.
 * \tparam T value type of the vector-like. Has to implement std::isnan(), std::numeric_limits<T>::has_quiet_NaN,
 * std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::has_infinity, std::numeric_limits<T>::infinity(),
 * std::numeric_limits<T>::max().
 */
template <typename T>
class One_critical_filtration : public std::vector<T>
{
 private:
  using Base = std::vector<T>;

 public:
  template <typename U = T>
  using Single_point = One_critical_filtration<U>;

  // CONSTRUCTORS

  One_critical_filtration() : Base() {};
  // warning: can be problematic if the user never updates the values and let it like that, {-inf, -inf, ...} is not
  // considered as -inf.
  One_critical_filtration(int n) : Base(n, -T_inf) {};  // minus infinity by default
  One_critical_filtration(int n, T value) : Base(n, value) {};
  One_critical_filtration(std::initializer_list<T> init) : Base(init) {};
  One_critical_filtration(const std::vector<T> &v) : Base(v) {};
  One_critical_filtration(std::vector<T> &&v) : Base(std::move(v)) {};
  One_critical_filtration(typename std::vector<T>::iterator it_begin, typename std::vector<T>::iterator it_end)
      : Base(it_begin, it_end) {};
  One_critical_filtration(typename std::vector<T>::const_iterator it_begin,
                          typename std::vector<T>::const_iterator it_end)
      : Base(it_begin, it_end) {};

  One_critical_filtration &operator=(const One_critical_filtration &a)
  {
    Base::operator=(a);
    return *this;
  }

  // HERITAGE

  using std::vector<T>::operator[];
  using value_type = T;

  // CONVERTERS

  operator std::vector<T> &() const { return *this; }

  operator std::vector<T>() const { return static_cast<std::vector<T> >(*this); }

  // like numpy
  template <typename U>
  One_critical_filtration<U> as_type() const
  {
    One_critical_filtration<U> out;
    out.reserve(Base::size());
    for (std::size_t i = 0u; i < Base::size(); i++) out.push_back(static_cast<U>(Base::operator[](i)));
    return out;
  }

  // ACCESS

  std::size_t num_parameters() const { return Base::size(); }

  constexpr static One_critical_filtration inf() { return {T_inf}; }

  constexpr static One_critical_filtration minus_inf() { return {-T_inf}; }

  constexpr static One_critical_filtration nan()
  {
    if constexpr (std::numeric_limits<T>::has_quiet_NaN) {
      return {std::numeric_limits<T>::quiet_NaN()};
    } else {
      return {};  // to differentiate it from 0, an empty filtration value can't do much anyway.
    }
  }

  // DESCRIPTORS

  bool is_inf() const
  {
    if (Base::size() != 1) return false;
    return (Base::operator[](0) == T_inf);
  }

  bool is_minus_inf() const
  {
    if constexpr (std::is_same<T, bool>::value) {
      return false;  // suppresses a warning
    } else {
      if (Base::size() != 1) return false;
      return (Base::operator[](0) == -T_inf);
    }
  }

  bool is_nan() const
  {
    if constexpr (std::numeric_limits<T>::has_quiet_NaN) {
      if (Base::size() != 1) return false;
      return std::isnan(Base::operator[](0));
    } else {
      return Base::empty();
    }
  }

  bool is_finite() const
  {
    if (Base::size() > 1) return true;
    if (Base::size() == 0) return false;
    auto first_value = Base::operator[](0);  // TODO : Maybe check all entries ?
    if (std::isnan(first_value) || first_value == -T_inf || first_value == T_inf) return false;
    return true;
  }

  // COMPARAISON OPERATORS

  friend bool operator<(const One_critical_filtration &a, const One_critical_filtration &b)
  {
    if (a.is_inf() || a.is_nan() || b.is_nan() || b.is_minus_inf()) return false;
    if (b.is_inf() || a.is_minus_inf()) return true;
    bool isSame = true;
    auto n = a.size();
    GUDHI_CHECK(a.size() == b.size(), "Two filtration points with different number of parameters are not comparable.");
    for (auto i = 0u; i < n; ++i) {
      if (a[i] > b[i]) return false;
      if (isSame && a[i] != b[i]) isSame = false;
    }
    return !isSame;
  }

  friend bool operator<=(const One_critical_filtration &a, const One_critical_filtration &b)
  {
    if (a.is_nan() || b.is_nan()) return false;
    if (b.is_inf() || a.is_minus_inf()) return true;
    if (a.is_inf() || b.is_minus_inf()) return false;
    auto n = a.size();
    GUDHI_CHECK(a.size() == b.size(), "Two filtration points with different number of parameters are not comparable.");
    for (std::size_t i = 0u; i < n; ++i) {
      if (a[i] > b[i]) return false;
    }
    return true;
  }

  friend bool operator>(const One_critical_filtration &a, const One_critical_filtration &b) { return b < a; }

  friend bool operator>=(const One_critical_filtration &a, const One_critical_filtration &b) { return b <= a; }

  friend bool operator==(const One_critical_filtration &a, const One_critical_filtration &b)
  {
    if (a.num_parameters() != b.num_parameters()) return false;
    for (auto i = 0u; i < a.num_parameters(); i++) {
      if (a[i] != b[i]) return false;
    }
    return true;
  }

  friend bool operator!=(const One_critical_filtration &a, const One_critical_filtration &b) { return !(a == b); }

  // ARITHMETIC OPERATORS

  // opposite
  friend One_critical_filtration operator-(const One_critical_filtration &f)
  {
    One_critical_filtration result;
    result.reserve(f.size());
    for (auto val : f) {
      result.push_back(-val);
    }
    return result;
  }

  // One_critical_filtration &operator-()
  // {
  //   for (auto &val : *this) {
  //     val = -val;
  //   }
  // }

  // subtraction
  friend One_critical_filtration operator-(One_critical_filtration result, const One_critical_filtration &to_subtract)
  {
    result -= to_subtract;
    return result;
  }

  friend One_critical_filtration operator-(One_critical_filtration result, const T &to_subtract)
  {
    result -= to_subtract;
    return result;
  }

  friend One_critical_filtration operator-(const T &value, One_critical_filtration result)
  {
    // TODO: in one go
    result = -result;
    result += value;
    return result;
  }

  friend One_critical_filtration &operator-=(One_critical_filtration &result,
                                             const One_critical_filtration &to_subtract)
  {
    if (result.empty()) return result;

    if (result.is_nan() || to_subtract.is_nan() || (result.is_inf() && to_subtract.is_inf()) ||
        (result.is_minus_inf() && to_subtract.is_minus_inf())) {
      result = nan();
      return result;
    }
    if (result.is_inf() || to_subtract.is_minus_inf()) {
      result = inf();
      return result;
    }
    if (result.is_minus_inf() || to_subtract.is_inf()) {
      result = minus_inf();
      return result;
    }

    GUDHI_CHECK(result.size() == to_subtract.size(),
                "Two filtration points with different number of parameters cannot be subtracted.");

    return apply_operation_with_finite_values_(result, to_subtract, subtract_);
  }

  friend One_critical_filtration &operator-=(One_critical_filtration &result, const T &to_subtract)
  {
    if (result.empty()) return result;

    if (result.is_nan() || std::isnan(to_subtract) || (result.is_inf() && to_subtract == T_inf) ||
        (result.is_minus_inf() && to_subtract == -T_inf)) {
      result = nan();
      return result;
    }
    if (result.is_inf() || to_subtract == -T_inf) {
      result = inf();
      return result;
    }
    if (result.is_minus_inf() || to_subtract == T_inf) {
      result = minus_inf();
      return result;
    }

    return apply_scalar_operation_on_finite_value_(result, to_subtract, subtract_);
  }

  // addition
  friend One_critical_filtration operator+(One_critical_filtration result, const One_critical_filtration &to_add)
  {
    result += to_add;
    return result;
  }

  friend One_critical_filtration operator+(One_critical_filtration result, const T &to_add)
  {
    result += to_add;
    return result;
  }

  friend One_critical_filtration operator+(const T &to_add, One_critical_filtration result)
  {
    result += to_add;
    return result;
  }

  friend One_critical_filtration &operator+=(One_critical_filtration &result, const One_critical_filtration &to_add)
  {
    if (result.empty()) return result;

    if (result.is_nan() || to_add.is_nan() || (result.is_inf() && to_add.is_minus_inf()) ||
        (result.is_minus_inf() && to_add.is_inf())) {
      result = nan();
      return result;
    }
    if (result.is_inf() || to_add.is_inf()) {
      result = inf();
      return result;
    }
    if (result.is_minus_inf() || to_add.is_minus_inf()) {
      result = minus_inf();
      return result;
    }

    GUDHI_CHECK(result.size() == to_add.size(),
                "Two filtration points with different number of parameters cannot be added.");

    return apply_operation_with_finite_values_(result, to_add, add_);
  }

  friend One_critical_filtration &operator+=(One_critical_filtration &result, const T &to_add)
  {
    if (result.empty()) return result;

    if (result.is_nan() || std::isnan(to_add) || (result.is_inf() && to_add == -T_inf) ||
        (result.is_minus_inf() && to_add == T_inf)) {
      result = nan();
      return result;
    }
    if (result.is_inf() || to_add == T_inf) {
      result = inf();
      return result;
    }
    if (result.is_minus_inf() || to_add == -T_inf) {
      result = minus_inf();
      return result;
    }

    return apply_scalar_operation_on_finite_value_(result, to_add, add_);
  }

  // multiplication
  friend One_critical_filtration operator*(One_critical_filtration result, const One_critical_filtration &to_mul)
  {
    result *= to_mul;
    return result;
  }

  friend One_critical_filtration operator*(One_critical_filtration result, const T &to_mul)
  {
    result *= to_mul;
    return result;
  }

  friend One_critical_filtration operator*(const T &to_mul, One_critical_filtration result)
  {
    result *= to_mul;
    return result;
  }

  friend One_critical_filtration &operator*=(One_critical_filtration &result, const One_critical_filtration &to_mul)
  {
    if (result.empty()) return result;

    if (result.is_nan() || to_mul.is_nan()) {
      result = nan();
      return result;
    }

    bool res_is_infinite = result.is_inf() || result.is_minus_inf();
    bool to_mul_is_infinite = to_mul.is_inf() || to_mul.is_minus_inf();

    if (res_is_infinite && to_mul_is_infinite) {
      if (to_mul.is_minus_inf()) {
        result[0] = -result[0];
      }
      return result;
    }

    if (res_is_infinite || to_mul_is_infinite) {
      const One_critical_filtration &finite = res_is_infinite ? to_mul : result;
      const T infinite = res_is_infinite ? result[0] : to_mul[0];
      result = finite;
      return apply_scalar_operation_on_finite_value_(result, infinite, multiply_);
    }

    GUDHI_CHECK(result.size() == to_mul.size(),
                "Two filtration points with different number of parameters cannot be multiplied.");

    return apply_operation_with_finite_values_(result, to_mul, multiply_);
  }

  friend One_critical_filtration &operator*=(One_critical_filtration &result, const T &to_mul)
  {
    if (result.empty()) return result;

    if (result.is_nan() || std::isnan(to_mul)) {
      result = nan();
      return result;
    }

    if (result.is_inf() || result.is_minus_inf()) {
      if (to_mul == 0)
        result = nan();
      else if (to_mul < 0)
        result[0] = -result[0];
      return result;
    }

    return apply_scalar_operation_on_finite_value_(result, to_mul, multiply_);
  }

  // division
  friend One_critical_filtration operator/(One_critical_filtration result, const One_critical_filtration &to_div)
  {
    result /= to_div;
    return result;
  }

  friend One_critical_filtration operator/(One_critical_filtration result, const T &to_div)
  {
    result /= to_div;
    return result;
  }

  friend One_critical_filtration operator/(const T &value, const One_critical_filtration &f)
  {
    if (f.empty()) return f;
    if (std::isnan(value) || f.is_nan()) return nan();

    One_critical_filtration result(f.size(), value);
    result /= f;
    return result;
  }

  friend One_critical_filtration &operator/=(One_critical_filtration &result, const One_critical_filtration &to_div)
  {
    if (result.empty()) return result;

    bool res_is_infinite = result.is_inf() || result.is_minus_inf();
    bool to_div_is_infinite = to_div.is_inf() || to_div.is_minus_inf();

    if (result.is_nan() || to_div.is_nan() || (res_is_infinite && to_div_is_infinite)) {
      result = nan();
      return result;
    }

    if (to_div_is_infinite) {
      return apply_scalar_operation_on_finite_value_with_all_nan_possible_(result, to_div[0], divide_);
    }

    GUDHI_CHECK(res_is_infinite || result.size() == to_div.size(),
                "Two filtration points with different number of parameters cannot be divided.");

    if (res_is_infinite) {
      result.resize(to_div.size(), result[0]);
    }

    return apply_operation_with_finite_values_(result, to_div, divide_);
  }

  friend One_critical_filtration &operator/=(One_critical_filtration &result, const T &to_div)
  {
    if (result.empty()) return result;

    bool res_is_infinite = result.is_inf() || result.is_minus_inf();
    bool to_div_is_infinite = to_div == T_inf || to_div == -T_inf;

    if (to_div == 0 || std::isnan(to_div) || result.is_nan() || (res_is_infinite && to_div_is_infinite)) {
      result = nan();
      return result;
    }

    if (res_is_infinite) {
      if (to_div < 0) result[0] = -result[0];
      return result;
    }

    return apply_scalar_operation_on_finite_value_with_all_nan_possible_(result, to_div, divide_);
  }

  // MODIFIERS

  /** \brief This functions take the filtration value `this` and pushes it to
   * the cone \f$ \{ y\in \mathbb R^n : y>=x \} \f$. After calling this method,
   * the value of this is updated to \f$ \mathrm{this} = \min \{ y\in \mathbb
   * R^n : y>=this \}\cap \{ y\in \mathbb R^n : y>=x \}
   * @param[in] x The target filtration value on which to push `this`.
   */
  void push_to(const One_critical_filtration &x)
  {
    if (this->is_inf() || this->is_nan() || x.is_nan() || x.is_minus_inf()) return;
    if (x.is_inf() || this->is_minus_inf()) {
      *this = x;
      return;
    }

    GUDHI_CHECK(this->num_parameters() == x.num_parameters(),
                "A filtration value cannot be pushed to another one with different numbers of parameters.");

    for (std::size_t i = 0; i < x.num_parameters(); i++)
      Base::operator[](i) = Base::operator[](i) > x[i] ? Base::operator[](i) : x[i];
  }

  /** \brief This functions take the filtration value `this` and pulls it to the
   * cone \f$ \{ y\in \mathbb R^n : y<=x \} \f$. After calling this method, the
   * value of this is updated to \f$ \mathrm{this} = \max \{ y\in \mathbb R^n :
   * y<=this \}\cap \{ y\in \mathbb R^n : y<=x \}
   * @param[in] x The target filtration value on which to push `this`.
   */
  void pull_to(const One_critical_filtration &x)
  {
    if (x.is_inf() || this->is_nan() || x.is_nan() || this->is_minus_inf()) return;
    if (this->is_inf() || x.is_minus_inf()) {
      *this = x;
      return;
    }

    GUDHI_CHECK(this->num_parameters() == x.num_parameters(),
                "A filtration value cannot be pulled to another one with different numbers of parameters.");

    for (std::size_t i = 0u; i < x.num_parameters(); i++)
      Base::operator[](i) = Base::operator[](i) > x[i] ? x[i] : Base::operator[](i);
  }

  /*
   * Same as `compute_coordinates_in_grid` but does the operation in-place
   */
  template <typename oned_array>
  void project_onto_grid(const std::vector<oned_array> &grid, bool coordinate = true)
  {
    GUDHI_CHECK(grid.size() >= Base::size(),
                "The grid should not be smaller than the number of parameters in the filtration value.");
    for (std::size_t parameter = 0u; parameter < Base::size(); ++parameter) {
      const auto &filtration = grid[parameter];
      auto d =
          std::distance(filtration.begin(),
                        std::lower_bound(filtration.begin(),
                                         filtration.end(),
                                         static_cast<typename oned_array::value_type>(Base::operator[](parameter))));
      Base::operator[](parameter) = coordinate ? static_cast<T>(d) : static_cast<T>(filtration[d]);
    }
  }

  // FONCTIONNALITIES

  // scalar product of a filtration value with x.
  template <typename U = T>
  friend U compute_linear_projection(const One_critical_filtration &f, const std::vector<U> &x)
  {
    U projection = 0;
    std::size_t size = std::min(x.size(), f.size());
    for (std::size_t i = 0u; i < size; i++) projection += x[i] * static_cast<U>(f[i]);
    return projection;
  }

  friend T compute_norm(const One_critical_filtration &f)
  {
    T out = 0;
    for (auto &val : f) out += (val * val);
    return std::sqrt(out);
  }

  friend T compute_euclidean_distance_to(const One_critical_filtration &f, const One_critical_filtration &other)
  {
    T out = 0;
    for (std::size_t i = 0u; i < other.size(); i++) {
      out += (f[i] - other[i]) * (f[i] - other[i]);
    }
    return std::sqrt(out);
  }

  /**
   * Given a grid in an array of shape (num_parameters, filtration_values of this parameter),
   * projects itself into this grid, and returns the coordinates of this projected points
   * in the given grid
   */
  template <typename out_type = std::int32_t, typename U = T>
  friend One_critical_filtration<out_type> compute_coordinates_in_grid(const One_critical_filtration &f,
                                                                       const std::vector<std::vector<U> > &grid)
  {
    One_critical_filtration<out_type> coords = f.as_type<out_type>();
    coords.project_onto_grid(grid);
    return coords;
  }

  /**
   * Given a grid in an array of shape (num_parameters, filtration_values of this parameter),
   * and assuming that `this` correspond to the coordinates in this grid,
   * returns the points evaluated in this grid
   */
  template <typename U>
  friend One_critical_filtration<U> evaluate_coordinates_in_grid(const One_critical_filtration &f,
                                                                 const std::vector<std::vector<U> > &grid)
  {
    One_critical_filtration<U> pushed_value(f.size());

    GUDHI_CHECK(grid.size() == f.size(),
                "The size of the grid should correspond to the number of parameters in the filtration value.");

    U grid_inf = One_critical_filtration<U>::T_inf;

    for (std::size_t parameter = 0u; parameter < grid.size(); ++parameter) {
      const auto &filtration = grid[parameter];
      const auto &c = f[parameter];
      pushed_value[parameter] = c == f.T_inf ? grid_inf : filtration[c];
    }
    return pushed_value;
  }

  // UTILITIES

  // easy debug
  friend std::ostream &operator<<(std::ostream &stream, const One_critical_filtration &f)
  {
    if (f.is_inf()) {
      stream << "[inf, ..., inf]";
      return stream;
    }
    if (f.is_minus_inf()) {
      stream << "[-inf, ..., -inf]";
      return stream;
    }
    if (f.is_nan()) {
      stream << "[NaN]";
      return stream;
    }
    if (f.empty()) {
      stream << "[]";
      return stream;
    }
    stream << "[";
    for (std::size_t i = 0; i < f.size() - 1; i++) {
      stream << f[i] << ", ";
    }
    if (!f.empty()) stream << f.back();
    stream << "]";
    return stream;
  }

 public:
  // TODO : maybe add the {inf}, minus inf, nan  there as static members? this
  // would make comparisons faster (just compare the ptr)
  // TODO : I'm not sure why constexpr doesn't work anymore
  constexpr static const T T_inf =
      std::numeric_limits<T>::has_infinity ? std::numeric_limits<T>::infinity() : std::numeric_limits<T>::max();

  // for compiler
  constexpr static bool is_multi_critical = false;

 private:
  constexpr static bool subtract_(T &v1, T v2) { return add_(v1, -v2); }

  constexpr static bool add_(T &v1, T v2)
  {
    if (std::isnan(v1) || std::isnan(v2) || (v1 == T_inf && v2 == -T_inf) || (v1 == -T_inf && v2 == T_inf)) {
      v1 = std::numeric_limits<T>::quiet_NaN();
      return false;
    }
    if (v1 == T_inf || v1 == -T_inf) {
      return true;
    }
    if (v2 == T_inf || v2 == -T_inf) {
      v1 = v2;
      return true;
    }

    v1 += v2;
    return true;
  }

  constexpr static bool multiply_(T &v1, T v2)
  {
    bool v1_is_infinite = v1 == T_inf || v1 == -T_inf;
    bool v2_is_infinite = v2 == T_inf || v2 == -T_inf;

    if (std::isnan(v1) || std::isnan(v2) || (v1_is_infinite && v2 == 0) || (v1 == 0 && v2_is_infinite)) {
      v1 = std::numeric_limits<T>::quiet_NaN();
      return false;
    }

    if ((v1 == T_inf && v2 > 0) || (v1 == -T_inf && v2 < 0) || (v1 < 0 && v2 == -T_inf) || (v1 > 0 && v2 == T_inf)) {
      v1 = T_inf;
      return true;
    }

    if ((v1 == T_inf && v2 < 0) || (v1 == -T_inf && v2 > 0) || (v1 > 0 && v2 == -T_inf) || (v1 < 0 && v2 == T_inf)) {
      v1 = -T_inf;
      return true;
    }

    v1 *= v2;
    return true;
  }

  constexpr static bool divide_(T &v1, T v2)
  {
    bool v1_is_infinite = v1 == T_inf || v1 == -T_inf;
    bool v2_is_infinite = v2 == T_inf || v2 == -T_inf;

    if (std::isnan(v1) || std::isnan(v2) || v2 == 0 || (v1_is_infinite && v2_is_infinite)) {
      v1 = std::numeric_limits<T>::quiet_NaN();
      return false;
    }

    if (v1 == 0 || (v1_is_infinite && v2 > 0)) return true;

    if (v1_is_infinite && v2 < 0) {
      v1 = -v1;
      return true;
    }

    if (v2_is_infinite) {
      v1 = 0;
      return true;
    }

    v1 /= v2;
    return true;
  }

  constexpr static bool update_sign_(T toComp, int &sign)
  {
    if (toComp == T_inf) {
      if (sign == 0)
        sign = 1;
      else if (sign == -1)
        return false;
    } else if (toComp == -T_inf) {
      if (sign == 0)
        sign = -1;
      else if (sign == 1)
        return false;
    } else {
      return false;
    }

    return true;
  }

  template <typename F>
  static One_critical_filtration &apply_operation_with_finite_values_(One_critical_filtration &result,
                                                                      const One_critical_filtration &to_operate,
                                                                      F &&operate)
  {
    bool allSameInf = true;
    bool allNaN = true;
    int sign = 0;
    for (auto i = 0u; i < result.size(); ++i) {
      if (operate(result[i], to_operate[i])) {
        allNaN = false;
      } else {
        if constexpr (!std::numeric_limits<T>::has_quiet_NaN) {
          result = nan();
          return result;
        }
      }
      if (allSameInf) allSameInf = update_sign_(result[i], sign);
    }

    if (allSameInf) result = (sign == 1 ? inf() : minus_inf());
    if (allNaN) result = nan();

    return result;
  }

  template <typename F>
  static One_critical_filtration &apply_scalar_operation_on_finite_value_(One_critical_filtration &result,
                                                                          const T &to_operate,
                                                                          F &&operate)
  {
    for (auto &val : result) {
      if constexpr (std::numeric_limits<T>::has_quiet_NaN) {
        operate(val, to_operate);
      } else {
        if (!operate(val, to_operate)) {
          result = nan();
          return result;
        }
      }
    }

    return result;
  }

  template <typename F>
  static One_critical_filtration &apply_scalar_operation_on_finite_value_with_all_nan_possible_(
      One_critical_filtration &result,
      const T &to_operate,
      F &&operate)
  {
    bool allNaN = true;

    for (auto &val : result) {
      if (operate(val, to_operate)) {
        allNaN = false;
      } else {
        if constexpr (!std::numeric_limits<T>::has_quiet_NaN) {
          result = nan();
          return result;
        }
      }
    }
    if (allNaN) result = nan();

    return result;
  }
};

}  // namespace Gudhi::multi_filtration

namespace std {

template <typename T>
class numeric_limits<Gudhi::multi_filtration::One_critical_filtration<T> >
{
 public:
  static constexpr bool has_infinity = true;

  static constexpr Gudhi::multi_filtration::One_critical_filtration<T> infinity() noexcept
  {
    return Gudhi::multi_filtration::One_critical_filtration<T>::inf();
  };

  // non-standard
  static constexpr Gudhi::multi_filtration::One_critical_filtration<T> minus_infinity() noexcept
  {
    return Gudhi::multi_filtration::One_critical_filtration<T>::minus_inf();
  };

  static constexpr Gudhi::multi_filtration::One_critical_filtration<T> max() noexcept(false)
  {
    throw std::logic_error(
        "The maximal value cannot be represented with no finite numbers of parameters."
        "Use `max(number_of_parameters)` instead");
  };

  // non-standard, so I don't want to define a default value.
  static constexpr Gudhi::multi_filtration::One_critical_filtration<T> max(unsigned int n) noexcept
  {
    return Gudhi::multi_filtration::One_critical_filtration<T>(n, std::numeric_limits<T>::max());
  };

  static constexpr Gudhi::multi_filtration::One_critical_filtration<T> quiet_NaN() noexcept
  {
    return Gudhi::multi_filtration::One_critical_filtration<T>::nan();
  };
};

}  // namespace std

#endif  // ONE_CRITICAL_FILTRATIONS_H_

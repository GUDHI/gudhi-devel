/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Loiseaux
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef ONE_CRITICAL_FILTRATIONS_H_
#define ONE_CRITICAL_FILTRATIONS_H_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <vector>

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
 * Use \ref Multi_critical_filtration for multicritical filtrations.
 * \tparam T value type of the vector-like.
 */
template <typename T>
class One_critical_filtration : public std::vector<T> {
public:
  One_critical_filtration() : std::vector<T>(){};
  One_critical_filtration(int n)
      : std::vector<T>(n, -T_inf){}; // minus infinity by default
  One_critical_filtration(int n, T value)
      : std::vector<T>(n, value){};
  One_critical_filtration(std::initializer_list<T> init)
      : std::vector<T>(init){};
  One_critical_filtration(const std::vector<T> &v)
      : std::vector<T>(v){};
  One_critical_filtration(std::vector<T> &&v)
      : std::vector<T>(std::move(v)){};
  One_critical_filtration(typename std::vector<T>::iterator it_begin,
                                     typename std::vector<T>::iterator it_end)
      : std::vector<T>(it_begin, it_end){};
  One_critical_filtration(
      typename std::vector<T>::const_iterator it_begin,
      typename std::vector<T>::const_iterator it_end)
      : std::vector<T>(it_begin, it_end){};

  using std::vector<T>::operator[];
  using value_type = T;
  using OneCritical = One_critical_filtration<T>;
  template <typename value_type>
  using Base = One_critical_filtration<value_type>;
  operator std::vector<T> &() const { return *this; }
  std::vector<T> get_vector() const {
    return static_cast<std::vector<T>>(*this);
  }
  std::size_t num_parameters() const { return this->size(); }
  // TODO : REMOVE SIZE ( confusing when kcritical) 
  // static std::size_t num_generators() {return 1;}

  operator std::vector<T>() const { return static_cast<std::vector<T>>(*this); }

  inline const std::vector<T> & as_vector() const { return *this; }

  inline bool is_inf() const {
    if (this->size() != 1)
      return false;
    return (this->operator[](0) == T_inf);
  }
  inline bool is_minus_inf() const {
    if constexpr (std::is_same<T, bool>::value) {
      return false; // suppresses a warning
    } else {
      if (this->size() != 1)
        return false;
      return (this->operator[](0) == -T_inf);
    }
  }
  inline bool is_nan() const {
    if (this->size() != 1)
      return false;
    return std::isnan(*(this->begin()));
  }
  inline bool is_finite() const {
    if (this->size() > 1)
      return true;
    if (this->size() == 0)
      return false;
    auto first_value = *(this->begin()); // TODO : Maybe check all entries ?
    if (std::isnan(first_value) || first_value == -T_inf ||
        first_value == T_inf)
      return false;
    return true;
  }

  inline friend bool operator<(const One_critical_filtration<T> &a,
                               const One_critical_filtration<T> &b) {
    if (a.is_inf() || a.is_nan() || b.is_nan() || b.is_minus_inf())
      return false;
    if (b.is_inf() || a.is_minus_inf())
      return true;
    bool isSame = true;
    auto n = a.size();
    assert(a.size() == b.size());
    for (auto i = 0u; i < n; ++i) {
      if (a[i] > b[i])
        return false;
      if (isSame && a[i] != b[i])
        isSame = false;
    }
    if (isSame)
      return false;
    return true;
  }
  inline friend bool operator<=(const One_critical_filtration<T> &a,
                                const One_critical_filtration<T> &b) {
    if (a.is_nan() || b.is_nan())
      return false;
    if (b.is_inf() || a.is_minus_inf())
      return true;
    if (a.is_inf() || b.is_minus_inf())
      return false;
    auto n = a.size();
    assert(a.size() == b.size());
    for (std::size_t i = 0u; i < n; ++i) {
      if (a[i] > b[i])
        return false;
    }
    return true;
  }

  // GREATER THAN OPERATORS
  inline friend bool operator>(const One_critical_filtration<T> &a,
                               const One_critical_filtration<T> &b) {
    return b < a;
  }
  inline friend bool operator>=(const One_critical_filtration<T> &a,
                                const One_critical_filtration<T> &b) {
    return b <= a;
  }

  inline One_critical_filtration<T> &
  operator=(const One_critical_filtration<T> &a) {
    std::vector<T>::operator=(a);
    return *this;
  }
  inline One_critical_filtration<T>& operator-() {
    for (auto& truc : *this) {
      truc = -truc;
    }
  }
  inline One_critical_filtration<T> copy() const { return *this; }

  std::vector<T> &_convert_back() { return *this; }
  constexpr static One_critical_filtration<T> inf() {
    return {T_inf};
  }
  constexpr static One_critical_filtration<T> minus_inf() {
    return {-T_inf};
  }
  constexpr static One_critical_filtration<T> nan() {
    return {std::numeric_limits<T>::quiet_NaN()};
  }

  // operators *=
  inline friend One_critical_filtration<T> &
  operator-=(One_critical_filtration<T> &result,
             const One_critical_filtration<T> &to_substract) {
    if (result.is_nan() || to_substract.is_nan() ||
        (result.is_inf() && to_substract.is_inf()) ||
        (result.is_minus_inf() && to_substract.is_minus_inf())) [[unlikely]] {
      result = std::numeric_limits<T>::quiet_NaN();
      return result;
    }
    if (result.is_inf() || to_substract.is_minus_inf()) [[unlikely]] {
      result = inf();
      return result;
    }
    if (result.is_minus_inf() || to_substract.is_inf()) [[unlikely]] {
      result = minus_inf();
      return result;
    }
    std::transform(result.begin(), result.end(), to_substract.begin(),
                   result.begin(), std::minus<T>());
    return result;
  }
  inline friend One_critical_filtration<T> &
  operator+=(One_critical_filtration<T> &result,
             const One_critical_filtration<T> &to_add) {

    if (result.is_nan() || to_add.is_nan() ||
        (result.is_inf() && to_add.is_minus_inf()) ||
        (result.is_minus_inf() && to_add.is_inf())) [[unlikely]] {
      result = nan();
      return result;
    }
    if (result.is_inf() || to_add.is_inf()) [[unlikely]] {
      result = inf();
      return result;
    }
    if (result.is_minus_inf() || to_add.is_minus_inf()) [[unlikely]] {
      result = minus_inf();
      return result;
    }

    std::transform(result.begin(), result.end(), to_add.begin(), result.begin(),
                   std::plus<T>());
    return result;
  }
  inline friend One_critical_filtration<T> &
  operator*=(One_critical_filtration<T> &result,
             const One_critical_filtration<T> &to_add) {

    assert(result.is_finite() || to_add.is_finite());
    std::transform(result.begin(), result.end(), to_add.begin(), result.begin(),
                   std::multiplies<T>());
    return result;
  }

  inline friend One_critical_filtration<T> &
  operator-=(One_critical_filtration<T> &result,
             const T &to_substract) {
    for (auto &truc : result) {
      truc -= to_substract;
    }
    return result;
  }
  inline friend One_critical_filtration<T> &
  operator+=(One_critical_filtration<T> &result, const T &to_add) {
    for (auto &truc : result) {
      truc += to_add;
    }
    return result;
  }
  inline friend One_critical_filtration<T> &
  operator*=(One_critical_filtration<T> &result, const T &to_add) {
    if (to_add == T_inf || to_add == -T_inf) [[unlikely]] {
      for (auto &truc : result) {
        if (truc > 0)
          truc = to_add;
        else if (truc < 0)
          truc = -to_add;
        // 0*inf = 0
      }
      return result;
    }
    for (auto &truc : result) {
      truc *= to_add;
    }
    return result;
  }
  inline friend One_critical_filtration<T> &
  operator/=(One_critical_filtration<T> &result, const T &to_add) {
    for (auto &truc : result) {
      truc /= to_add;
    }
    return result;
  }

  /// OPERATORS *
  inline friend One_critical_filtration<T>
  operator+(One_critical_filtration<T> result, const T &to_add) {
    result += to_add;
    return result;
  }

  inline friend One_critical_filtration<T>
  operator+(One_critical_filtration<T> result,
            const One_critical_filtration<T> &to_add) {
    result += to_add;
    return result;
  }
  inline friend One_critical_filtration<T>
  operator-(One_critical_filtration<T> result,
            const One_critical_filtration<T> &to_add) {
    result -= to_add;
    return result;
  }
  inline friend One_critical_filtration<T>
  operator-(One_critical_filtration<T> result, const T &to_add) {
    result -= to_add;
    return result;
  }
  inline friend One_critical_filtration<T>
  operator*(One_critical_filtration<T> result, const T &to_add) {
    result *= to_add;
    return result;
  }

  inline friend One_critical_filtration<T>
  operator+(const T &to_add,
            const One_critical_filtration<T> &result) {
    return result + to_add;
  }
  inline friend One_critical_filtration<T>
  operator-(const T &to_add,
            const One_critical_filtration<T> &result) {
    return result - to_add;
  }
  inline friend One_critical_filtration<T>
  operator*(const T &to_add,
            const One_critical_filtration<T> &result) {
    return result * to_add;
  }

  // template<class array_like>
  inline friend bool
  operator==(const One_critical_filtration<T> &self,
             const One_critical_filtration<T> &to_compare) {
    if (self.num_parameters() != to_compare.num_parameters())
      return false;
    auto it = to_compare.begin();
    for (auto i = 0u; i < self.num_parameters(); i++) {
      if (self.at(i) != *(it++))
        return false;
    }
    return true;
  }

  inline static std::vector<std::vector<T>> to_python(
      const std::vector<One_critical_filtration<T>> &to_convert) {
    return std::vector<std::vector<T>>(to_convert.begin(), to_convert.end());
  }

  inline static std::vector<One_critical_filtration<T>>
  from_python(const std::vector<std::vector<T>> &to_convert) {
    return std::vector<One_critical_filtration<T>>(
        to_convert.begin(), to_convert.end());
  }

  /** \brief This functions take the filtration value `this` and pushes it to
   * the cone \f$ \{ y\in \mathbb R^n : y>=x \} \f$. After calling this method,
   * the value of this is updated to \f$ \mathrm{this} = \min \{ y\in \mathbb
   * R^n : y>=this \}\cap \{ y\in \mathbb R^n : y>=x \}
   * @param[in] x The target filtration value on which to push `this`.
   */
  inline void push_to(const One_critical_filtration<T> &x) {
    if (this->is_inf() || this->is_nan() || x.is_nan() || x.is_minus_inf())
      return;
    if (x.is_inf() || this->is_minus_inf()) {
      *this = x;
      return;
    }
    if (this->num_parameters() != x.num_parameters()) {
      std::cerr << "Does only work with 1-critical filtrations ! Sizes "
                << this->num_parameters() << " and " << x.num_parameters() << "are different !"
                << std::endl;
      std::cerr << "This : " << *this << std::endl;
      std::cerr << "arg : " << x << std::endl;
      throw std::logic_error("Bad sizes");
    }
    for (std::size_t i = 0; i < x.num_parameters(); i++)
      this->operator[](i) =
          this->operator[](i) > x[i] ? this->operator[](i) : x[i];
  }

  /** \brief This functions take the filtration value `this` and pulls it to the
   * cone \f$ \{ y\in \mathbb R^n : y<=x \} \f$. After calling this method, the
   * value of this is updated to \f$ \mathrm{this} = \max \{ y\in \mathbb R^n :
   * y<=this \}\cap \{ y\in \mathbb R^n : y<=x \}
   * @param[in] x The target filtration value on which to push `this`.
   */
  inline void pull_to(const One_critical_filtration<T> &x) {
    if (x.is_inf() || this->is_nan() || x.is_nan() || this->is_minus_inf())
      return;
    if (this->is_inf() || x.is_minus_inf()) {
      *this = x;
      return;
    }
    if (this->num_parameters() != x.num_parameters()) {
      std::cerr << "Does only work with 1-critical filtrations ! Sizes "
                << this->num_parameters() << " and " << x.num_parameters() << "are different !"
                << std::endl;
      std::cerr << "This : " << *this << std::endl;
      std::cerr << "arg : " << x << std::endl;
      throw std::logic_error("Bad sizes");
    }
    for (std::size_t i = 0u; i < x.num_parameters(); i++)
      this->operator[](i) =
          this->operator[](i) > x[i] ? x[i] : this->operator[](i);
  }
  // Warning, this function  assumes that the comparisons checks have already
  // been made !
  inline void insert_new(One_critical_filtration<T> to_concatenate) {
    this->insert(this->end(), std::move_iterator(to_concatenate.begin()),
                 std::move_iterator(to_concatenate.end()));
  }

  // scalar product of a filtration value with x.
  template <typename U=T>
  inline U linear_projection(const std::vector<U> &x) {
    U projection = 0;
    std::size_t size = std::min(x.size(), this->size());
    for (std::size_t i = 0u; i < size; i++)
      projection += x[i] * static_cast<U>(this->operator[](i));
    return projection;
  }

  // easy debug
  inline friend std::ostream &
  operator<<(std::ostream &stream,
             const One_critical_filtration<T> &truc) {
    if (truc.is_inf()) {
      stream << "[inf, ..., inf]";
      return stream;
    }
    if (truc.is_minus_inf()) {
      stream << "[-inf, ..., -inf]";
      return stream;
    }
    if (truc.is_nan()) {
      stream << "[NaN]";
      return stream;
    }
    if (truc.empty()) {
      stream << "[]";
      return stream;
    }
    stream << "[";
    for (std::size_t i = 0; i < truc.size() - 1; i++) {
      stream << truc[i] << ", ";
    }
    if (!truc.empty())
      stream << truc.back();
    stream << "]";
    return stream;
  }
  inline T norm() const {
    T out = 0;
    for (auto &stuff : *this)
      out += std::pow(stuff, 2);
    return sqrt(out);
  }
  inline T distance(const One_critical_filtration<T> &other) const {
    T out = 0;
    for (std::size_t i = 0u; i < other.size(); i++) {
      out += std::pow((*this)[i] - other[i], 2);
    }
    return sqrt(out);
  }
  /**
   * Given a grid in an array of shape (num_parameters, filtration_values of this parameter),
   * projects itself into this grid, and returns the coordinates of this projected points
   * in the given grid
   */
  template <typename out_type=std::int32_t,typename U>
  inline One_critical_filtration<out_type>
  coordinates_in_grid(const std::vector<std::vector<U>>& grid) const {
    One_critical_filtration<out_type> coords(this->size());
    assert(grid.size() >= this->size());
    for (std::size_t parameter = 0u; parameter < grid.size(); ++parameter) {
      const auto &filtration = grid[parameter];
      auto C = std::distance(
          filtration.begin(),
          std::lower_bound(filtration.begin(), filtration.end(),
                           static_cast<U>(this->operator[](parameter))));
      coords[parameter] = static_cast<out_type>(C);
    }
    return coords;
  }
  /**
   * Given a grid in an array of shape (num_parameters, filtration_values of this parameter),
   * and assuming that `this` correspond to the coordinates in this grid, 
   * returns the points evaluated in this grid
   */
  template <typename U>
  inline One_critical_filtration<U>
  evaluate_in_grid(const std::vector<std::vector<U>>& grid) const {
    One_critical_filtration<U> pushed_value(this->size());
    assert(grid.size() == this->size());
    U grid_inf = std::numeric_limits<U>::has_infinity
                     ? std::numeric_limits<U>::infinity()
                     : std::numeric_limits<U>::max();
    for (std::size_t parameter = 0u; parameter < grid.size(); ++parameter) {
      const auto &filtration = grid[parameter];
      const auto& c = this->operator[](parameter);
      pushed_value[parameter] = c == T_inf ? grid_inf : filtration[c];
    }
    return pushed_value;
  }

  /*
   * Same as `evaluate_in_grid` but does the operation in-place
   */
  template <typename oned_array>
  inline void coordinates_in_grid_inplace(const std::vector<oned_array>& grid,
                                          bool coordinate = true) {
    One_critical_filtration<std::size_t> coords(this->size());
    assert(grid.size() >= this->size());
    for (std::size_t parameter = 0u; parameter < grid.size(); ++parameter) {
      const auto &filtration = grid[parameter];
      auto d =
          std::distance(filtration.begin(),
                        std::lower_bound(filtration.begin(), filtration.end(),
                                         this->operator[](parameter)));
      this->operator[](parameter) =
          coordinate ? static_cast<T>(d) : static_cast<T>(grid[parameter][d]);
      ;
    }
  }
  // like numpy
  template <typename U> inline One_critical_filtration<U> astype() const {
    One_critical_filtration<U> out(this->size());
    for (std::size_t i = 0u; i < this->size(); i++)
      out[i] = static_cast<U>(this->operator[](i));
    return out;
  }

public:
  // TODO : maybe add the {inf}, minus inf, nan  there as static members? this
  // would make comparisons faster (just compare the ptr)
  // TODO : I'm not sure why constexpr doens't work anymore
  constexpr static const T T_inf = std::numeric_limits<T>::has_infinity
                                       ? std::numeric_limits<T>::infinity()
                                       : std::numeric_limits<T>::max();

  // for compiler
  constexpr static bool is_multi_critical = false;
};

} // namespace Gudhi::multi_filtration

namespace std {

template <typename T>
class numeric_limits<Gudhi::multi_filtration::One_critical_filtration<T>> {
public:
  static constexpr bool has_infinity = std::numeric_limits<T>::has_infinity;

  static Gudhi::multi_filtration::One_critical_filtration<T> infinity() throw() {
    return Gudhi::multi_filtration::One_critical_filtration<T>(1, std::numeric_limits<T>::infinity());
  };
  static Gudhi::multi_filtration::One_critical_filtration<T> minus_infinity() throw() {
    return Gudhi::multi_filtration::One_critical_filtration<T>(1, -std::numeric_limits<T>::infinity());
  };
  static Gudhi::multi_filtration::One_critical_filtration<T> max() throw() {
    return Gudhi::multi_filtration::One_critical_filtration<T>(1, std::numeric_limits<T>::max());
  };
  static Gudhi::multi_filtration::One_critical_filtration<T>quiet_NaN() throw() {
    return Gudhi::multi_filtration::One_critical_filtration<T>(1, numeric_limits<T>::quiet_NaN());
  };
};

} // namespace std

#endif // ONE_CRITICAL_FILTRATIONS_H_

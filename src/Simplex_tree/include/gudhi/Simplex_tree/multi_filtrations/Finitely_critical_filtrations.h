/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which
 * is released under MIT. See file LICENSE or go to
 * https://gudhi.inria.fr/licensing/ for full license details. Author(s): David

 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef FINITELY_CRITICAL_FILTRATIONS_H_
#define FINITELY_CRITICAL_FILTRATIONS_H_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <vector>

namespace Gudhi::multiparameter::multi_filtrations {

/** 
 * \brief Vector-like filtration value, for multiparameter persistence, with
 * numpy-like methods.
 *
 * \ingroup multiparameter
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

/*
 * Multi-critical filtration extension to \ref One_critical_filtration . 
 * If the `co` parameter is set to true, it reverses the poset order,
 * i.e., the order \f$\le$\f  in \f$\mathbb R^n\f$ becomes \f$\ge$\f. 
 *
 * if `multi_filtration_` contains the points \f$(a_1, a_2, \ldots, a_k) \in (\mathbb R^n)^k\f$,
 * then a new point \f$x$\f will be in this filtration if there exists an \f$a_i\f$ such that 
 * \f$a_i \le x\f$.
 *
 */
template <typename T, bool co=false> class Multi_critical_filtration {
public:
  using value_type = T;
  using OneCritical = One_critical_filtration<T>; // Type of the One critical subtype
  template <typename value_type>
  using Base = Multi_critical_filtration<value_type>;
  std::size_t num_parameters() const { return multi_filtration_.size() == 0 ? 0 : multi_filtration_[0].size(); }
  Multi_critical_filtration() : multi_filtration_(1)
  {
  multi_filtration_[0] = co ? One_critical_filtration<T>::inf() : One_critical_filtration<T>::minus_inf();
  }; // initialize at 1, with 0 parameter (consistent with below)
  Multi_critical_filtration(int n)
      : multi_filtration_({One_critical_filtration<T>(
            n)}){}; // minus infinity by default
  Multi_critical_filtration(int n, T value)
      : multi_filtration_({One_critical_filtration<T>(
            n, value)}){}; // std::vector<T>(n, value){};
  Multi_critical_filtration(std::initializer_list<T> init)
      : multi_filtration_({One_critical_filtration<T>(
            init)}){}; // : std::vector<T>({std::vector<T>(init)};
  Multi_critical_filtration(const std::vector<T> &v) : multi_filtration_({v}){};
  Multi_critical_filtration(std::vector<T> &&v) : multi_filtration_({std::move(v)}){};
  Multi_critical_filtration(const std::vector<One_critical_filtration<T>> &v) : multi_filtration_(v){};
  Multi_critical_filtration(std::vector<One_critical_filtration<T>> &&v) : multi_filtration_(std::move(v)){};
  Multi_critical_filtration(typename std::vector<T>::iterator it_begin,
                      typename std::vector<T>::iterator it_end)
      : multi_filtration_(
            {One_critical_filtration<T>(it_begin, it_end)}){};
  Multi_critical_filtration(typename std::vector<T>::const_iterator it_begin,
                      typename std::vector<T>::const_iterator it_end)
      : multi_filtration_({One_critical_filtration<T>(
            it_begin, it_end)}){}; // : std::vector<T>(it_begin, it_end){};
  Multi_critical_filtration(const Multi_critical_filtration<T> &other)
      : multi_filtration_(other.multi_filtration_){};

  Multi_critical_filtration &operator=(const Multi_critical_filtration<T> &other) {
    this->multi_filtration_ = other.multi_filtration_;
    return *this;
  }
  inline friend bool
  operator==(const Multi_critical_filtration<T,co> &self,
             const Multi_critical_filtration<T,co> &to_compare) {
    if (self.num_generators() != to_compare.num_generators())
      return false;
    auto it = to_compare.begin();
    for (auto i = 0u; i < self.num_generators(); i++) {
      if (self[i] != *(it++))
        return false;
    }
    return true;
  }
  void reserve(std::size_t n) { this->multi_filtration_.reserve(n); }
  void set_num_generators(std::size_t n) {this->multi_filtration_.resize(n);}

  inline bool is_inf() const {
    return this->multi_filtration_.size() == 1 &&
           this->multi_filtration_[0].is_inf();
  }
  inline bool is_minus_inf() const {
    return this->multi_filtration_.size() == 1 &&
           this->multi_filtration_[0].is_minus_inf();
  }
  inline bool is_nan() const {
    return this->multi_filtration_.size() == 1 &&
           this->multi_filtration_[0].is_nan();
  }
  inline bool is_finite() const {
    if (this->empty())
      return false;
    for (auto &stuff : *this) {
      if (!stuff.is_finite())
        return false;
    }
    return true;
  }

  operator OneCritical() const {
    assert(this->num_generators() == 1);
    return this->multi_filtration_[0];
  }

  OneCritical &operator[](std::size_t i) { return this->multi_filtration_[i]; }
  inline OneCritical factorize_below() const {
    if (this->num_generators() == 0) [[unlikely]]
      return OneCritical();
    OneCritical result(multi_filtration_[0].num_parameters(), OneCritical::T_inf);
    for (const auto &stuff : this->multi_filtration_) {
      if (stuff.is_nan() || stuff.is_minus_inf())
        return stuff;
      if (stuff.is_inf())
        continue;
      for (std::size_t i = 0; i < stuff.num_parameters(); ++i) {
        result[i] = std::min(result[i], stuff[i]);
      }
    }
    return result;
  }
  /*
   * returns the smallest value for the poset order that is bigger than all of the values
   * in this multi filtration
   */
  inline OneCritical factorize_above() const {
    if (this->num_generators() == 0) [[unlikely]]
      return OneCritical();
    OneCritical result(multi_filtration_[0].num_parameters(), -OneCritical::T_inf);
    for (auto &stuff : this->multi_filtration_) {
      if (stuff.is_nan() || stuff.is_inf())
        return stuff;
      if (stuff.is_minus_inf())
        continue;
      for (std::size_t i = 0; i < stuff.num_parameters(); ++i) {
        result[i] = std::max(result[i], stuff[i]);
      }
    }
    return result;
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
    for (auto &stuff : *this) {
      stuff.push_to(x);
    }
  }
  // TODO : this is not well defined for kcritical <-> kcritical

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
    for (auto &stuff : *this) {
      stuff.pull_to(x);
    }
  } 
  // cannot be const, as gudhi iterators are not const
  template <typename U> inline U linear_projection(const std::vector<U> &x) {
    if constexpr (co) {
      U projection = std::numeric_limits<U>::lowest();
      for (const auto &y : *this) {
        projection = std::max(projection, y.linear_projection(x));
      }
      return projection;
    } else {
      U projection = std::numeric_limits<U>::max();
      for (auto &y : *this) { // cannot be const (Gudhi)
        projection = std::min(projection, y.linear_projection(x));
      }
      return projection;
    }
  }

  /*
  * Checks if b is cleanable with respect to a
  */
  static inline bool dominates(const OneCritical&a, const OneCritical&b, value_type max_error) {
    if constexpr (co)
      return a - max_error <= b;
    else {
      return a + max_error >= b;
    }
  }

  static inline bool dominates(const OneCritical&a, const OneCritical&b) {
    if constexpr (co)
      return a <= b;
    else {
      return a >= b;
    }
  }
  static inline bool strictly_dominates(const OneCritical&a, const OneCritical&b) {
    if constexpr (co)
      return a < b;
    else {
      return a > b;
    }
  }

  /*
  * Same method as the one in OneCriticalFiltration. 
  * Given a grid, and assuming that `this` are the coordianates
  * in this grid, evaluate `this` into this grid
  */
  template <typename U>
  inline Multi_critical_filtration<U>
  evaluate_in_grid(const std::vector<std::vector<U>>& grid) const {
    Multi_critical_filtration<U> out(this->num_generators());
    for (std::size_t i = 0; i < this->num_generators(); ++i){
      out[i] = this->operator[](i).evaluate_in_grid(grid);
    }
    return out;
  }
  
  /*
  * Remove redundant values 
  */
  inline void clean(value_type max_error = 0) {
    // A bit ugly, todo : erase+removeif ?
    for (std::size_t i = 0; i < multi_filtration_.size(); i++) {
      for (std::size_t j = 0; j < multi_filtration_.size(); j++) {
        if (i == j)
          continue;
        if (dominates(multi_filtration_[j], multi_filtration_[i], max_error)) {
          multi_filtration_[i].clear();
        }
        while (j < multi_filtration_.size() && dominates(multi_filtration_[i], multi_filtration_[j], max_error)) {
          multi_filtration_[j].clear();
          i--;
        }
      }
    }
    multi_filtration_.erase(
        std::remove_if(multi_filtration_.begin(), multi_filtration_.end(),
                       [](const One_critical_filtration<T> &a) {
                         return a.empty();
                       }),
        multi_filtration_.end());
  }


  /*
  * Adds a birth point to the list of births, 
  * if it is useful (according to the method `dominate`)
  */
  inline void add_point(const One_critical_filtration<T> &x) {
    const bool verbose = false;
    if (multi_filtration_.empty()) {
      if constexpr (verbose)
        std::cout << "Adding x=" << x << " (currently empty)" << std::endl;
      multi_filtration_.push_back(x);
      return;
    }
    for (const auto &y : multi_filtration_){
      if (dominates(x,y)){
        return;
      }
    }
    if constexpr (verbose)
      std::cout << "x: " << x << " is useful, removing unnecessary entries"
                << std::endl;
    multi_filtration_.erase(
        std::remove_if(multi_filtration_.begin(), multi_filtration_.end(),
                       [&x](const One_critical_filtration<T> &y) {
                         if constexpr (verbose) {
                           if (dominates(y, x)) {
                             std::cout << "Removing y=" << y << std::endl;
                           } else {
                             std::cout << "Keeping y=" << y << std::endl;
                           }
                         }
                         return dominates(y, x);
                       }),
        multi_filtration_.end());

    multi_filtration_.push_back(x);
  }
  inline void re_clean(){
    // Ensures all points are useful again. Can be useful if points are added manually.
    // TODO : maybe optimize
    Multi_critical_filtration<value_type> out; // should be inf
    out.multi_filtration_.reserve(this->multi_filtration_.size());
    for (const auto& x : multi_filtration_){
      out.add_point(x);
    }
    std::swap(multi_filtration_, out);
  }

  // easy debug
  inline friend std::ostream &operator<<(std::ostream &stream,
                                         const Multi_critical_filtration<T,co> &truc) {
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
    stream << "(k=" << truc.multi_filtration_.size() << ")[";
    for (const auto &machin : truc) {
      stream << machin << "; ";
    }
    if (truc.multi_filtration_.size() > 0) {
      stream << "\b"
             << "\b";
    }
    stream << "]";
    return stream;
  }
  inline void clear() { multi_filtration_.clear(); }
  inline bool empty() const { return multi_filtration_.empty(); }
  inline const OneCritical& operator[] (std::size_t i) const { return multi_filtration_[i]; }

  inline typename std::vector<One_critical_filtration<T>>::iterator
  begin() {
    return multi_filtration_.begin();
  }
  inline typename std::vector<One_critical_filtration<T>>::iterator
  end() {
    return multi_filtration_.end();
  }
  inline typename std::vector<
      One_critical_filtration<T>>::const_iterator
  begin() const {
    return multi_filtration_.begin();
  }
  inline typename std::vector<
      One_critical_filtration<T>>::const_iterator
  end() const {
    return multi_filtration_.end();
  }


  /*
  * Same as its one critical counterpart.
  * If a grid is given, projects multifiltration_ on this grid and returns
  * multi critical filtration composed of the coordinates in the given grid
  *
  */
  template <typename oned_array>
  inline Multi_critical_filtration<std::int32_t>
  coordinates_in_grid(const std::vector<oned_array>& grid) const {
    assert(grid.size() >= this->num_generators());
    Multi_critical_filtration<std::int32_t> out(this->num_generators());
    for (std::size_t i = 0u; i < this->num_generators(); ++i) {
      out[i] = this->multi_filtration_[i].coordinates_in_grid(grid);
    }
    return out;
  }
  /*
  * Same as `coordinates_in_grid`, but does the operation in-place.
  *
  */
  template <typename oned_array>
  inline void coordinates_in_grid_inplace(const std::vector<oned_array>& grid,
                                          bool coordinate = true) {
    assert(grid.size() >= this->num_generators());
    for (auto &x : this->multi_filtration_) {
      x.coordinates_in_grid_inplace(grid, coordinate);
    }
  }
  template <typename U> inline Multi_critical_filtration<U> astype() const {
    std::vector<One_critical_filtration<U>> out(this->num_generators());
    for (std::size_t i = 0u; i < this->num_generators(); ++i) {
      out[i] = this->multi_filtration_[i].template astype<U>();
    }
    return Multi_critical_filtration<U>(std::move(out));
  }
  inline void push_back(const OneCritical& x){
    multi_filtration_.push_back(x);
  }
  inline const std::vector<One_critical_filtration<T>>& as_vector() const {
    return multi_filtration_;
  }

  inline std::vector<std::vector<T>> as_VECTOR() const {
    std::vector<std::vector<T>> out(this->num_generators(),
                                    std::vector<T>(this->num_parameters()));
    for (std::size_t i = 0; i < this->num_generators(); ++i) {
      for (std::size_t j = 0; j < this->num_parameters(); ++j) {
        out[i][j] = multi_filtration_[i][j];
      }
    }
    return out;
  }

  inline void _clean(bool keep_inf=true) {
    multi_filtration_.erase(std::remove_if(multi_filtration_.begin(), multi_filtration_.end(),
                              [keep_inf](const OneCritical &a) {
                                return a.empty() ||
                                       ((!keep_inf) &&
                                        (a.is_inf() || a.is_minus_inf()));
                              }),
                            multi_filtration_.end());
  }
  inline std::size_t num_generators() const {
    return multi_filtration_.size();
  }



  // TODO : this costs a lot... optimize / cheat in some way for python ?
  /*
  * Checks if `this`, seen as a birth curve is under the `other` birth curve,
  *
  */
  inline bool operator<(const Multi_critical_filtration<T,co> &other) const {
    //check if this curves is below other's curve 
    // ie for each guy in this, check if there is a guy in other that dominates him
    for (std::size_t i = 0u; i < multi_filtration_.size(); ++i) {
      for (std::size_t j = 0u; j < other.multi_filtration_.size(); ++j) {
        // i<j
        if (strictly_dominates(other.multi_filtration_[j],multi_filtration_[i]))
          continue;
      }
      return false;
    }
    return true;
  }
  /*
  * Checks if `this`, seen as a birth curve is over the `other` birth curve,
  */
  inline bool operator>(const Multi_critical_filtration<T,co> &other) const {
    return other < *this;
  }

  /*
  * Checks if `this`, seen as a birth curve is under the `other` birth curve,
  */
  inline bool operator<=(const Multi_critical_filtration<T,co> &other) const {
    //check if this curves is below other's curve 
    // ie for each guy in this, check if there is a guy in other that dominates him
    for (std::size_t i = 0u; i < multi_filtration_.size(); ++i) {
      for (std::size_t j = 0u; j < other.multi_filtration_.size(); ++j) {
        // i <= j 
        if (dominates(other.multi_filtration_[j], multi_filtration_[i]))
          continue;
      }
      return false;
    }
    return true;
  }
  /*
  * Checks if `this`, seen as a birth curve is over the `other` birth curve,
  */
  inline bool operator>=(const Multi_critical_filtration<T,co> &other) const {
    return other <= *this;
  }

public:
  // for compiler
  constexpr static const T T_inf = One_critical_filtration<T>::T_inf;
  constexpr static const bool is_multi_critical = true;

private:
  std::vector<One_critical_filtration<T>> multi_filtration_;
};

} // namespace Gudhi::multiparameter::multi_filtrations

namespace std {

template <typename T>
class numeric_limits<Gudhi::multiparameter::multi_filtrations::
                         One_critical_filtration<T>> {
public:
  static constexpr bool has_infinity = std::numeric_limits<T>::has_infinity;

  static Gudhi::multiparameter::multi_filtrations::
      One_critical_filtration<T>
      infinity() throw() {
    return Gudhi::multiparameter::multi_filtrations::
        One_critical_filtration<T>(
            1, std::numeric_limits<T>::infinity());
  };
  static Gudhi::multiparameter::multi_filtrations::
      One_critical_filtration<T>
      minus_infinity() throw() {
    return Gudhi::multiparameter::multi_filtrations::
        One_critical_filtration<T>(
            1, -std::numeric_limits<T>::infinity());
  };
  static Gudhi::multiparameter::multi_filtrations::
      One_critical_filtration<T>
      max() throw() {
    return Gudhi::multiparameter::multi_filtrations::
        One_critical_filtration<T>(1, std::numeric_limits<T>::max());
  };
  static Gudhi::multiparameter::multi_filtrations::
      One_critical_filtration<T>
      quiet_NaN() throw() {
    return Gudhi::multiparameter::multi_filtrations::
        One_critical_filtration<T>(1,
                                              numeric_limits<T>::quiet_NaN());
  };
};

} // namespace std

#endif // FINITELY_CRITICAL_FILTRATIONS_H_

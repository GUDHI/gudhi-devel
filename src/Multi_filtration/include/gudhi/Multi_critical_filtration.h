/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Loiseaux
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef MULTI_CRITICAL_FILTRATIONS_H_
#define MULTI_CRITICAL_FILTRATIONS_H_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <vector>

#include <gudhi/One_critical_filtration.h>

namespace Gudhi::multi_filtration {

/**
 * Multi-critical filtration extension to \ref One_critical_filtration .
 * If the `co` parameter is set to true, it reverses the poset order,
 * i.e., the order \f$\le\f$  in \f$\mathbb R^n\f$ becomes \f$\ge\f$.
 *
 * if `multi_filtration_` contains the points \f$(a_1, a_2, \ldots, a_k) \in (\mathbb R^n)^k\f$,
 * then a new point \f$x\f$ will be in this filtration if there exists an \f$a_i\f$ such that
 * \f$a_i \le x\f$.
 *
 * \ingroup multi_filtration

 */
template <typename T, bool co = false>
class Multi_critical_filtration {
 public:
  using value_type = T;
  using OneCritical = One_critical_filtration<T>;  // Type of the One critical subtype
  template <typename value_type>
  using Base = Multi_critical_filtration<value_type>;
  std::size_t num_parameters() const { return multi_filtration_.size() == 0 ? 0 : multi_filtration_[0].size(); }
  Multi_critical_filtration() : multi_filtration_(1) {
    multi_filtration_[0] = co ? One_critical_filtration<T>::inf() : One_critical_filtration<T>::minus_inf();
  };  // initialize at 1, with 0 parameter (consistent with below)
  Multi_critical_filtration(int n)
      : multi_filtration_({One_critical_filtration<T>(n)}) {};  // minus infinity by default
  Multi_critical_filtration(int n, T value)
      : multi_filtration_({One_critical_filtration<T>(n, value)}) {};  // std::vector<T>(n, value){};
  Multi_critical_filtration(std::initializer_list<T> init)
      : multi_filtration_({One_critical_filtration<T>(init)}) {};  // : std::vector<T>({std::vector<T>(init)};
  Multi_critical_filtration(const std::vector<T> &v) : multi_filtration_({v}) {};
  Multi_critical_filtration(std::vector<T> &&v) : multi_filtration_({std::move(v)}) {};
  Multi_critical_filtration(const std::vector<One_critical_filtration<T>> &v) : multi_filtration_(v) {};
  Multi_critical_filtration(std::vector<One_critical_filtration<T>> &&v) : multi_filtration_(std::move(v)) {};
  Multi_critical_filtration(typename std::vector<T>::iterator it_begin, typename std::vector<T>::iterator it_end)
      : multi_filtration_({One_critical_filtration<T>(it_begin, it_end)}) {};
  Multi_critical_filtration(typename std::vector<T>::const_iterator it_begin,
                            typename std::vector<T>::const_iterator it_end)
      : multi_filtration_({One_critical_filtration<T>(it_begin, it_end)}) {};  // : std::vector<T>(it_begin, it_end){};
  Multi_critical_filtration(const Multi_critical_filtration<T> &other) : multi_filtration_(other.multi_filtration_) {};

  Multi_critical_filtration &operator=(const Multi_critical_filtration<T> &other) {
    this->multi_filtration_ = other.multi_filtration_;
    return *this;
  }
  inline friend bool operator==(const Multi_critical_filtration<T, co> &self,
                                const Multi_critical_filtration<T, co> &to_compare) {
    if (self.num_generators() != to_compare.num_generators()) return false;
    auto it = to_compare.begin();
    for (auto i = 0u; i < self.num_generators(); i++) {
      if (self[i] != *(it++)) return false;
    }
    return true;
  }
  void reserve(std::size_t n) { this->multi_filtration_.reserve(n); }
  void set_num_generators(std::size_t n) { this->multi_filtration_.resize(n); }

  inline bool is_inf() const { return this->multi_filtration_.size() == 1 && this->multi_filtration_[0].is_inf(); }
  inline bool is_minus_inf() const {
    return this->multi_filtration_.size() == 1 && this->multi_filtration_[0].is_minus_inf();
  }
  inline bool is_nan() const { return this->multi_filtration_.size() == 1 && this->multi_filtration_[0].is_nan(); }
  inline bool is_finite() const {
    if (this->empty()) return false;
    for (auto &stuff : *this) {
      if (!stuff.is_finite()) return false;
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
      if (stuff.is_nan() || stuff.is_minus_inf()) return stuff;
      if (stuff.is_inf()) continue;
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
      if (stuff.is_nan() || stuff.is_inf()) return stuff;
      if (stuff.is_minus_inf()) continue;
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
    if (this->is_inf() || this->is_nan() || x.is_nan() || x.is_minus_inf()) return;
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
    if (x.is_inf() || this->is_nan() || x.is_nan() || this->is_minus_inf()) return;
    if (this->is_inf() || x.is_minus_inf()) {
      *this = x;
      return;
    }
    for (auto &stuff : *this) {
      stuff.pull_to(x);
    }
  }
  // cannot be const, as gudhi iterators are not const
  template <typename U>
  inline U linear_projection(const std::vector<U> &x) {
    if constexpr (co) {
      U projection = std::numeric_limits<U>::lowest();
      for (const auto &y : *this) {
        projection = std::max(projection, y.linear_projection(x));
      }
      return projection;
    } else {
      U projection = std::numeric_limits<U>::max();
      for (auto &y : *this) {  // cannot be const (Gudhi)
        projection = std::min(projection, y.linear_projection(x));
      }
      return projection;
    }
  }

  /*
   * Checks if b is cleanable with respect to a
   */
  static inline bool dominates(const OneCritical &a, const OneCritical &b, value_type max_error) {
    if constexpr (co)
      return a - max_error <= b;
    else {
      return a + max_error >= b;
    }
  }

  static inline bool dominates(const OneCritical &a, const OneCritical &b) {
    if constexpr (co)
      return a <= b;
    else {
      return a >= b;
    }
  }
  static inline bool strictly_dominates(const OneCritical &a, const OneCritical &b) {
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
  inline Multi_critical_filtration<U> evaluate_in_grid(const std::vector<std::vector<U>> &grid) const {
    Multi_critical_filtration<U> out(this->num_generators());
    for (std::size_t i = 0; i < this->num_generators(); ++i) {
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
        if (i == j) continue;
        if (dominates(multi_filtration_[j], multi_filtration_[i], max_error)) {
          multi_filtration_[i].clear();
        }
        while (j < multi_filtration_.size() && dominates(multi_filtration_[i], multi_filtration_[j], max_error)) {
          multi_filtration_[j].clear();
          i--;
        }
      }
    }
    multi_filtration_.erase(std::remove_if(multi_filtration_.begin(), multi_filtration_.end(),
                                           [](const One_critical_filtration<T> &a) { return a.empty(); }),
                            multi_filtration_.end());
  }

  /*
   * Adds a birth point to the list of births,
   * if it is useful (according to the method `dominate`)
   */
  inline void add_point(const One_critical_filtration<T> &x) {
    const bool verbose = false;
    if (multi_filtration_.empty()) {
      if constexpr (verbose) std::cout << "Adding x=" << x << " (currently empty)" << std::endl;
      multi_filtration_.push_back(x);
      return;
    }
    for (const auto &y : multi_filtration_) {
      if (dominates(x, y)) {
        return;
      }
    }
    if constexpr (verbose) std::cout << "x: " << x << " is useful, removing unnecessary entries" << std::endl;
    multi_filtration_.erase(std::remove_if(multi_filtration_.begin(), multi_filtration_.end(),
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
  inline void re_clean() {
    // Ensures all points are useful again. Can be useful if points are added manually.
    // TODO : maybe optimize
    Multi_critical_filtration<value_type> out;  // should be inf
    out.multi_filtration_.reserve(this->multi_filtration_.size());
    for (const auto &x : multi_filtration_) {
      out.add_point(x);
    }
    std::swap(multi_filtration_, out);
  }

  // easy debug
  inline friend std::ostream &operator<<(std::ostream &stream, const Multi_critical_filtration<T, co> &truc) {
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
  inline const OneCritical &operator[](std::size_t i) const { return multi_filtration_[i]; }

  inline typename std::vector<One_critical_filtration<T>>::iterator begin() { return multi_filtration_.begin(); }
  inline typename std::vector<One_critical_filtration<T>>::iterator end() { return multi_filtration_.end(); }
  inline typename std::vector<One_critical_filtration<T>>::const_iterator begin() const {
    return multi_filtration_.begin();
  }
  inline typename std::vector<One_critical_filtration<T>>::const_iterator end() const {
    return multi_filtration_.end();
  }

  /*
   * Same as its one critical counterpart.
   * If a grid is given, projects multifiltration_ on this grid and returns
   * multi critical filtration composed of the coordinates in the given grid
   *
   */
  template <typename oned_array>
  inline Multi_critical_filtration<std::int32_t> coordinates_in_grid(const std::vector<oned_array> &grid) const {
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
  inline void coordinates_in_grid_inplace(const std::vector<oned_array> &grid, bool coordinate = true) {
    assert(grid.size() >= this->num_generators());
    for (auto &x : this->multi_filtration_) {
      x.coordinates_in_grid_inplace(grid, coordinate);
    }
  }
  template <typename U>
  inline Multi_critical_filtration<U> astype() const {
    std::vector<One_critical_filtration<U>> out(this->num_generators());
    for (std::size_t i = 0u; i < this->num_generators(); ++i) {
      out[i] = this->multi_filtration_[i].template astype<U>();
    }
    return Multi_critical_filtration<U>(std::move(out));
  }
  inline void push_back(const OneCritical &x) { multi_filtration_.push_back(x); }
  inline const std::vector<One_critical_filtration<T>> &as_vector() const { return multi_filtration_; }

  inline std::vector<std::vector<T>> as_VECTOR() const {
    std::vector<std::vector<T>> out(this->num_generators(), std::vector<T>(this->num_parameters()));
    for (std::size_t i = 0; i < this->num_generators(); ++i) {
      for (std::size_t j = 0; j < this->num_parameters(); ++j) {
        out[i][j] = multi_filtration_[i][j];
      }
    }
    return out;
  }

  inline void _clean(bool keep_inf = true) {
    multi_filtration_.erase(std::remove_if(multi_filtration_.begin(), multi_filtration_.end(),
                                           [keep_inf](const OneCritical &a) {
                                             return a.empty() || ((!keep_inf) && (a.is_inf() || a.is_minus_inf()));
                                           }),
                            multi_filtration_.end());
  }
  inline std::size_t num_generators() const { return multi_filtration_.size(); }

  // TODO : this costs a lot... optimize / cheat in some way for python ?
  /*
   * Checks if `this`, seen as a birth curve is under the `other` birth curve,
   *
   */
  inline bool operator<(const Multi_critical_filtration<T, co> &other) const {
    // check if this curves is below other's curve
    //  ie for each guy in this, check if there is a guy in other that dominates him
    for (std::size_t i = 0u; i < multi_filtration_.size(); ++i) {
      for (std::size_t j = 0u; j < other.multi_filtration_.size(); ++j) {
        // i<j
        if (strictly_dominates(other.multi_filtration_[j], multi_filtration_[i])) continue;
      }
      return false;
    }
    return true;
  }
  /*
   * Checks if `this`, seen as a birth curve is over the `other` birth curve,
   */
  inline bool operator>(const Multi_critical_filtration<T, co> &other) const { return other < *this; }

  /*
   * Checks if `this`, seen as a birth curve is under the `other` birth curve,
   */
  inline bool operator<=(const Multi_critical_filtration<T, co> &other) const {
    // check if this curves is below other's curve
    //  ie for each guy in this, check if there is a guy in other that dominates him
    for (std::size_t i = 0u; i < multi_filtration_.size(); ++i) {
      for (std::size_t j = 0u; j < other.multi_filtration_.size(); ++j) {
        // i <= j
        if (dominates(other.multi_filtration_[j], multi_filtration_[i])) continue;
      }
      return false;
    }
    return true;
  }
  /*
   * Checks if `this`, seen as a birth curve is over the `other` birth curve,
   */
  inline bool operator>=(const Multi_critical_filtration<T, co> &other) const { return other <= *this; }

 public:
  // for compiler
  constexpr static const T T_inf = One_critical_filtration<T>::T_inf;
  constexpr static const bool is_multi_critical = true;

 private:
  std::vector<One_critical_filtration<T>> multi_filtration_;
};

}  // namespace Gudhi::multi_filtration

#endif  // MULTI_CRITICAL_FILTRATIONS_H_

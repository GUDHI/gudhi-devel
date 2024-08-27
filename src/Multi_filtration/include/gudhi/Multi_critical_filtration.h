/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Loiseaux
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - 2024/08 Hannah Schreiber: Optimization and correction + numeric_limits
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef MULTI_CRITICAL_FILTRATIONS_H_
#define MULTI_CRITICAL_FILTRATIONS_H_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>
#include "gudhi/Debug_utils.h"

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
class Multi_critical_filtration
{
 public:
  using Single_point = One_critical_filtration<T>;  // Type of the One critical subtype
  // TODO: not the best name...
  using Points = std::vector<Single_point>;
  using iterator = typename Points::iterator;
  using const_iterator = typename Points::const_iterator;

  // CONSTRUCTORS

  // initialize at 1, with 1 parameter (consistent with below)
  Multi_critical_filtration() : multi_filtration_(1, co ? Single_point::inf() : Single_point::minus_inf()) {};
  // warning: can be problematic if the user never updates the values and let it like that, {-inf, -inf, ...} is not
  // considered as -inf.
  Multi_critical_filtration(int n) : multi_filtration_({Single_point(n)}) {};
  // std::vector<T>(n, value){};
  Multi_critical_filtration(int n, T value) : multi_filtration_({Single_point(n, value)}) {};
  // : std::vector<T>({std::vector<T>(init)};
  Multi_critical_filtration(std::initializer_list<T> init) : multi_filtration_({Single_point(init)}) {};
  Multi_critical_filtration(const std::vector<T> &v) : multi_filtration_({v}) {};
  Multi_critical_filtration(std::vector<T> &&v) : multi_filtration_({std::move(v)}) {};
  // assumes all of the same size and that v is minimal (add simplify?)
  Multi_critical_filtration(const std::vector<Single_point> &v) : multi_filtration_(v) {};
  // assumes all of the same size and that v is minimal (add simplify?)
  Multi_critical_filtration(std::vector<Single_point> &&v) : multi_filtration_(std::move(v)) {};
  Multi_critical_filtration(typename std::vector<T>::iterator it_begin, typename std::vector<T>::iterator it_end)
      : multi_filtration_({Single_point(it_begin, it_end)}) {};
  // : std::vector<T>(it_begin, it_end){};
  Multi_critical_filtration(typename std::vector<T>::const_iterator it_begin,
                            typename std::vector<T>::const_iterator it_end)
      : multi_filtration_({Single_point(it_begin, it_end)}) {};
  Multi_critical_filtration(const Multi_critical_filtration<T> &other) : multi_filtration_(other.multi_filtration_) {};

  Multi_critical_filtration &operator=(const Multi_critical_filtration<T> &other)
  {
    multi_filtration_ = other.multi_filtration_;
    return *this;
  }

  // VECTOR-LIKE

  using value_type = T;

  Single_point &operator[](std::size_t i) { return multi_filtration_[i]; }

  const Single_point &operator[](std::size_t i) const { return multi_filtration_[i]; }

  iterator begin() { return multi_filtration_.begin(); }

  iterator end() { return multi_filtration_.end(); }

  const_iterator begin() const { return multi_filtration_.begin(); }

  const_iterator end() const { return multi_filtration_.end(); }

  bool empty() const { return multi_filtration_.empty(); }

  void reserve(std::size_t n) { multi_filtration_.reserve(n); }

  void clear() { multi_filtration_.clear(); }

  // CONVERTERS

  operator Single_point() const
  {
    GUDHI_CHECK(num_generators() == 1,
                "Casting a " + std::to_string(num_generators()) +
                    "-critical filtration value into an 1-critical filtration value.");
    return multi_filtration_[0];
  }

  // like numpy
  template <typename U>
  Multi_critical_filtration<U> as_type() const
  {
    std::vector<One_critical_filtration<U>> out(num_generators());
    for (std::size_t i = 0u; i < num_generators(); ++i) {
      out[i] = multi_filtration_[i].template as_type<U>();
    }
    return Multi_critical_filtration<U>(std::move(out));
  }

  // only needed once in multipers and does not need to be a friend, so could just be defined
  // at the place it is needed and removed from here.
  friend std::vector<std::vector<T>> get_content(const Multi_critical_filtration &f)
  {
    std::vector<std::vector<T>> out(f.num_generators(), std::vector<T>(f.num_parameters()));
    const auto &cont = f.get_underlying_container();
    for (std::size_t i = 0; i < f.num_generators(); ++i) {
      for (std::size_t j = 0; j < f.num_parameters(); ++j) {
        out[i][j] = cont[i][j];
      }
    }
    return out;
  }

  // ACCESS

  const Points &get_underlying_container() const { return multi_filtration_; }

  std::size_t num_parameters() const { return multi_filtration_.empty() ? 0 : multi_filtration_[0].num_parameters(); }

  std::size_t num_generators() const { return multi_filtration_.size(); }

  constexpr static Multi_critical_filtration inf() { return Multi_critical_filtration(Single_point::inf()); }

  constexpr static Multi_critical_filtration minus_inf()
  {
    return Multi_critical_filtration(Single_point::minus_inf());
  }

  constexpr static Multi_critical_filtration nan() { return Multi_critical_filtration(Single_point::nan()); }

  // DESCRIPTORS

  bool is_inf() const { return multi_filtration_.size() == 1 && multi_filtration_[0].is_inf(); }

  bool is_minus_inf() const { return multi_filtration_.size() == 1 && multi_filtration_[0].is_minus_inf(); }

  bool is_nan() const { return multi_filtration_.size() == 1 && multi_filtration_[0].is_nan(); }

  bool is_finite() const
  {
    if (empty()) return false;
    if (multi_filtration_.size() > 1) return true;
    return multi_filtration_[0].is_finite();
  }

  // COMPARAISON OPERATORS

  // TODO : this costs a lot... optimize / cheat in some way for python ?
  /*
   * Checks if `this`, seen as a birth curve is under the `other` birth curve,
   *
   */
  bool operator<(const Multi_critical_filtration &other) const
  {
    // check if this curves is below other's curve
    //  ie for each guy in this, check if there is a guy in other that dominates him
    for (std::size_t i = 0u; i < other.multi_filtration_.size(); ++i) {
      //for each point in other, verify if it is strictly in the cone of at least one point of this
      bool isContained = false;
      for (std::size_t j = 0u; j < multi_filtration_.size() && !isContained; ++j) {
        // i<j
        isContained = _strictly_contains(multi_filtration_[j], other.multi_filtration_[i]);
      }
      if (!isContained) return false;
    }
    return true;
  }

  /*
   * Checks if `this`, seen as a birth curve is over the `other` birth curve,
   */
  bool operator>(const Multi_critical_filtration &other) const { return other < *this; }

  /*
   * Checks if `this`, seen as a birth curve is under the `other` birth curve,
   */
  bool operator<=(const Multi_critical_filtration &other) const
  {
    // check if this curves is below other's curve
    //  ie for each guy in this, check if there is a guy in other that dominates him
    for (std::size_t i = 0u; i < other.multi_filtration_.size(); ++i) {
      //for each point in other, verify if it is in the cone of at least one point of this
      bool isContained = false;
      for (std::size_t j = 0u; j < multi_filtration_.size() && !isContained; ++j) {
        // i<j
        isContained = _contains(multi_filtration_[j], other.multi_filtration_[i]);
      }
      if (!isContained) return false;
    }
    return true;
  }

  /*
   * Checks if `this`, seen as a birth curve is over the `other` birth curve,
   */
  bool operator>=(const Multi_critical_filtration &other) const { return other <= *this; }

  friend bool operator==(const Multi_critical_filtration &self, const Multi_critical_filtration &to_compare)
  {
    if (self.num_generators() != to_compare.num_generators()) return false;
    for (auto i = 0u; i < self.num_generators(); i++) {
      if (self[i] != to_compare[i]) return false;
    }
    return true;
  }

  friend bool operator!=(const Multi_critical_filtration &self, const Multi_critical_filtration &to_compare)
  {
    return !(self == to_compare);
  }

  // MODIFIERS

  void set_num_generators(std::size_t n) { multi_filtration_.resize(n); }

  /** \brief This functions take the filtration value `this` and pushes it to
   * the cone \f$ \{ y\in \mathbb R^n : y>=x \} \f$. After calling this method,
   * the value of this is updated to \f$ \mathrm{this} = \min \{ y\in \mathbb
   * R^n : y>=this \}\cap \{ y\in \mathbb R^n : y>=x \}
   * @param[in] x The target filtration value on which to push `this`.
   */
  void push_to(const Single_point &x)
  {
    if (this->is_inf() || this->is_nan() || x.is_nan() || x.is_minus_inf()) return;

    GUDHI_CHECK(x.is_inf() || x.num_parameters() == multi_filtration_[0].num_parameters() || !is_finite(),
                "Pushing to a filtration value with different number of parameters.");

    if (x.is_inf() || this->is_minus_inf()) {
      multi_filtration_ = {x};
      return;
    }
    for (auto &fil : *this) {
      fil.push_to(x);
    }

    simplify();
  }

  // TODO: this is not well defined for k-critical <-> k-critical

  /** \brief This functions take the filtration value `this` and pulls it to the
   * cone \f$ \{ y\in \mathbb R^n : y<=x \} \f$. After calling this method, the
   * value of this is updated to \f$ \mathrm{this} = \max \{ y\in \mathbb R^n :
   * y<=this \}\cap \{ y\in \mathbb R^n : y<=x \}
   * @param[in] x The target filtration value on which to push `this`.
   */
  void pull_to(const Single_point &x)
  {
    if (x.is_inf() || this->is_nan() || x.is_nan() || this->is_minus_inf()) return;

    GUDHI_CHECK(x.is_minus_inf() || x.num_parameters() == multi_filtration_[0].num_parameters() || !is_finite(),
                "Pulling to a filtration value with different number of parameters.");

    if (this->is_inf() || x.is_minus_inf()) {
      multi_filtration_ = {x};
      return;
    }
    for (auto &fil : *this) {
      fil.pull_to(x);
    }

    simplify();
  }

  /*
   * Adds a birth point to the list of births,
   * if it is useful (according to the method `dominate`)
   */
  bool add_point(const Single_point &x)
  {
    if (multi_filtration_.empty()) {
      multi_filtration_.push_back(x);
      return true;
    }

    GUDHI_CHECK(x.num_parameters() == multi_filtration_[0].num_parameters() || !is_finite() || !x.is_finite(),
                "Cannot add a point with different number of parameters.");

    std::size_t end = multi_filtration_.size();

    if (_point_can_be_added(x, 0, end)) {
      multi_filtration_.resize(end);
      multi_filtration_.push_back(x);
      return true;
    }

    return false;
  }

  // No security
  void add_guaranteed_point(const Single_point &x) { multi_filtration_.push_back(x); }

  /*
   * Same as `compute_coordinates_in_grid`, but does the operation in-place.
   */
  template <typename oned_array>
  void project_onto_grid(const std::vector<oned_array> &grid, bool coordinate = true)
  {
    // not sure why this is needed? Are the tests with the numbers of parameters not enough?
    GUDHI_CHECK(grid.size() >= num_generators(),
                "The grid should not be smaller than the number of generators in the filtration value.");

    for (auto &x : multi_filtration_) {
      x.project_onto_grid(grid, coordinate);
    }

    if (!coordinate) simplify();
  }

  // no NaN?
  void remove_empty_points(bool include_infinities = false)
  {
    multi_filtration_.erase(std::remove_if(multi_filtration_.begin(),
                                           multi_filtration_.end(),
                                           [include_infinities](const Single_point &a) {
                                             return a.empty() ||
                                                    ((include_infinities) && (a.is_inf() || a.is_minus_inf()));
                                           }),
                            multi_filtration_.end());
  }

  void simplify()
  {
    std::size_t end = 0;

    for (std::size_t curr = 0; curr < multi_filtration_.size(); ++curr) {
      if (_point_can_be_added(multi_filtration_[curr], 0, end)) {
        std::swap(multi_filtration_[end], multi_filtration_[curr]);
        ++end;
      }
    }

    multi_filtration_.resize(end);
  }

  // FONCTIONNALITIES

  friend Single_point factorize_below(const Multi_critical_filtration &f)
  {
    if (f.num_generators() == 0) return Single_point();
    Single_point result(f.num_parameters(), Single_point::T_inf);
    for (const auto &fil : f) {
      if (fil.is_nan() || fil.is_minus_inf()) return fil;
      if (fil.is_inf()) continue;
      for (std::size_t i = 0; i < f.num_parameters(); ++i) {
        result[i] = std::min(result[i], fil[i]);
      }
    }
    return result;
  }

  /*
   * returns the smallest value for the poset order that is bigger than all of the values
   * in this multi filtration
   */
  friend Single_point factorize_above(const Multi_critical_filtration &f)
  {
    if (f.num_generators() == 0) return Single_point();
    Single_point result(f.num_parameters(), -Single_point::T_inf);
    for (auto &fil : f) {
      if (fil.is_nan() || fil.is_inf()) return fil;
      if (fil.is_minus_inf()) continue;
      for (std::size_t i = 0; i < fil.num_parameters(); ++i) {
        result[i] = std::max(result[i], fil[i]);
      }
    }
    return result;
  }

  template <typename U = T>
  friend U compute_linear_projection(const Multi_critical_filtration &f, const std::vector<U> &x)
  {
    if constexpr (co) {
      U projection = std::numeric_limits<U>::lowest();
      for (const auto &y : f) {
        projection = std::max(projection, compute_linear_projection(y, x));
      }
      return projection;
    } else {
      U projection = std::numeric_limits<U>::max();
      for (const auto &y : f) {
        projection = std::min(projection, compute_linear_projection(y, x));
      }
      return projection;
    }
  }

  /*
   * Same as its one critical counterpart.
   * If a grid is given, projects multi_filtration_ on this grid and returns
   * multi critical filtration composed of the coordinates in the given grid
   *
   */
  template <typename out_type = std::int32_t, typename U = T>
  friend Multi_critical_filtration<out_type> compute_coordinates_in_grid(const Multi_critical_filtration &f,
                                                                         const std::vector<std::vector<U>> &grid)
  {
    Multi_critical_filtration<out_type> coords = f.as_type<out_type>();
    coords.project_onto_grid(grid);
    return coords;
  }

  /*
   * Same method as the one in OneCriticalFiltration.
   * Given a grid, and assuming that `this` are the coordinates
   * in this grid, evaluate `this` into this grid
   */
  template <typename U>
  friend Multi_critical_filtration<U> evaluate_coordinates_in_grid(const Multi_critical_filtration &f,
                                                                   const std::vector<std::vector<U>> &grid)
  {
    Multi_critical_filtration<U> out;
    out.set_num_generators(f.num_generators());
    for (std::size_t i = 0; i < f.num_generators(); ++i) {
      out[i] = evaluate_coordinates_in_grid(f[i], grid);
    }
    out.simplify();
    return out;
  }

  // UTILITIES

  // easy debug
  friend std::ostream &operator<<(std::ostream &stream, const Multi_critical_filtration &f)
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
    stream << "(k=" << f.multi_filtration_.size() << ")[";
    for (const auto &val : f) {
      stream << val << "; ";
    }
    if (f.multi_filtration_.size() > 0) {
      stream << "\b"
             << "\b";
    }
    stream << "]";
    return stream;
  }

 public:
  // for compiler
  constexpr static const bool is_multi_critical = true;

 private:
  Points multi_filtration_;

  /*
   * Checks if b is cleanable with respect to a
   */
  static bool _strictly_contains(const Single_point &a, const Single_point &b)
  {
    if constexpr (co)
      return a > b;
    else {
      return a < b;
    }
  }
  static bool _contains(const Single_point &a, const Single_point &b)
  {
    if constexpr (co)
      return a >= b;
    else {
      return a <= b;
    }
  }

  // 0 == equal
  // 1 == a dom b
  // 2 == b dom a
  // 3 == none
  static int _get_domination_relation(const Single_point &a, const Single_point &b)
  {
    if (a.is_nan() || b.is_nan()) return 3;

    GUDHI_CHECK(a.size() == b.size(),
                "Two points in the same k-critical value have to have the same numbers of parameters.");

    bool equal = true;
    bool allGreater = true;
    bool allSmaller = true;
    for (unsigned int i = 0; i < a.size(); ++i) {
      if (a[i] < b[i]) {
        if (!allSmaller) return 3;
        equal = false;
        allGreater = false;
      } else if (a[i] > b[i]) {
        if (!allGreater) return 3;
        equal = false;
        allSmaller = false;
      }
    }
    if (equal) return 0;

    if constexpr (co) {
      if (allSmaller) return 1;
      return 2;
    } else {
      if (allGreater) return 1;
      return 2;
    }
  }

  // assumes between 'curr' and 'end' everything is simplified:
  // no nan values and if there is an inf/-inf, then 'end - curr == 1'
  bool _point_can_be_added(const Single_point &x, std::size_t curr, std::size_t &end)
  {
    if (x.empty() || x.is_nan() || (x.is_inf() && end - curr != 0)) return false;

    if (x.is_minus_inf()) {
      if (end - curr == 1 && multi_filtration_[curr].is_minus_inf()) return false;
      // assumes that everything between curr and end is already simplified
      // so, if end - curr != 1, there can be no minus_inf anymore.
      end = curr;
      return true;
    }

    while (curr != end) {
      int res = _get_domination_relation(multi_filtration_[curr], x);
      if (res == 2 || res == 0) return false;  // x dominates or is equal
      if (res == 1) {                          // x is dominated
        --end;
        std::swap(multi_filtration_[curr], multi_filtration_[end]);
      } else {
        ++curr;
      }
    }
    return true;
  }
};

}  // namespace Gudhi::multi_filtration

namespace std {

template <typename T>
class numeric_limits<Gudhi::multi_filtration::Multi_critical_filtration<T>>
{
 public:
  static constexpr bool has_infinity = true;

  static constexpr Gudhi::multi_filtration::Multi_critical_filtration<T> infinity() noexcept
  {
    return Gudhi::multi_filtration::Multi_critical_filtration<T>::inf();
  };

  // non-standard
  static constexpr Gudhi::multi_filtration::Multi_critical_filtration<T> minus_infinity() noexcept
  {
    return Gudhi::multi_filtration::Multi_critical_filtration<T>::minus_inf();
  };

  static constexpr Gudhi::multi_filtration::Multi_critical_filtration<T> max() noexcept(false)
  {
    throw std::logic_error(
        "The maximal value cannot be represented with no finite numbers of generators."
        "Use `max(number_of_generators, number_of_parameters)` instead");
  };

  // non-standard, so I don't want to define default values.
  static constexpr Gudhi::multi_filtration::Multi_critical_filtration<T> max(unsigned int g, unsigned int n) noexcept
  {
    std::vector<typename Gudhi::multi_filtration::Multi_critical_filtration<T>::Single_point> v(
        g, std::vector<T>(n, std::numeric_limits<T>::max()));
    return Gudhi::multi_filtration::Multi_critical_filtration<T>(std::move(v));
  };

  static constexpr Gudhi::multi_filtration::Multi_critical_filtration<T> quiet_NaN() noexcept
  {
    return Gudhi::multi_filtration::Multi_critical_filtration<T>::nan();
  };
};

}  // namespace std

#endif  // MULTI_CRITICAL_FILTRATIONS_H_

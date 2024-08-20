/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Loiseaux
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef LINE_FILTRATION_TRANSLATION_H_INCLUDED
#define LINE_FILTRATION_TRANSLATION_H_INCLUDED

#include <cstddef>
#include <stdexcept>

#include <gudhi/One_critical_filtration.h>
#include <gudhi/Multi_critical_filtration.h>
#include <gudhi/Multi_persistence/Box.h>

namespace Gudhi::multi_persistence {

/* A line in \f$\mathbb R^n\f$, with some helpers to project points on it.
 * When the direction is not given, it is assumed to be diagonal.
 * As the line has a builtin parametrization, points in \f$\mathbb R^n\f$
 * that are on a line are given a time parameter in \f$\mathbb R\f$.
 * The method that end with a 2 returns the time t, while the other
 * ones return the full coordinates
 *
 * @ingroup multi_persistence
 */
template <typename T>
class Line {
 public:
  using point_type = Gudhi::multi_filtration::One_critical_filtration<T>;
  using kcritical_point_type = Gudhi::multi_filtration::Multi_critical_filtration<T>;
  /*
   * Checks that the argument define a correct, positively slopped line.
   */
  bool check_direction() const;
  Line();
  Line(const point_type &x);
  Line(point_type &&x);
  Line(const point_type &x, const point_type &v);
  /*
   * Returns the point whose intersection is \f$ \min\{ y\ge x \} \cap \mathrm{this}\f$
   */
  inline point_type push_forward(point_type x) const;
  /*
   * Retuns the time parameter of the coordinate given by push_forward.
   */
  template <typename U = T>
  inline U push_forward2(const point_type &x) const;
  /*
   * Retuns the time parameter of the coordinate given by push_forward.
   */
  template <typename U = T>
  inline U push_forward2(const kcritical_point_type &x) const;
  /*
   * Returns the point whose intersection is \f$ \max\{ y\le x \} \cap \mathrm{this}\f$
   */
  inline point_type push_back(point_type x) const;
  /*
   * Retuns the time parameter of the coordinate given by push_back.
   */
  template <typename U = T>
  inline U push_back2(const point_type &x) const;
  /*
   * Retuns the time parameter of the coordinate given by push_back.
   */
  template <typename U = T>
  inline U push_back2(const kcritical_point_type &x) const;
  inline int get_dim() const;
  /*
   * Given a box, returns the coordinates of the intersection of this box and `this` as a pair of points (low, high)
   * in this line, representing this interval.
   */
  std::pair<point_type, point_type> get_bounds(const Box<T> &box) const;
  /*
   * Retuns the times parameter of the coordinates in the pair given by get_bounds.
   */
  std::pair<T, T> get_bounds2(const Box<T> &box) const;

  // translation
  inline friend Line &operator+=(Line &to_translate, const point_type &x) {
    to_translate.basepoint_ -= x;
    return to_translate;
  }

  inline point_type &basepoint() { return basepoint_; }
  inline point_type &direction() { return direction_; }
  inline const point_type &basepoint() const { return basepoint_; }
  inline const point_type &direction() const { return direction_; }

 private:
  point_type basepoint_;  // any point on the line
  point_type direction_;  // direction of the line
};
template <typename T>
inline bool Line<T>::check_direction() const {
  bool is_trivial = true;
  for (const auto &stuff : basepoint_) {
    if (!stuff) {
      is_trivial = false;
    }
    if (stuff < 0) {
      throw std::invalid_argument("Direction should have positive entries.");
    }
  }
  if (is_trivial) {
    throw std::invalid_argument("Direction should have at least one non-trivial entry.");
  }
  if (direction_.size() && direction_.size() != basepoint_.size())
    throw std::invalid_argument("The dimensions of basepoint and direction are not equal.");
}
template <typename T>
Line<T>::Line() {}

template <typename T> Line<T>::Line(const point_type &x) : basepoint_(x) { check_direction();}

template <typename T>
Line<T>::Line(const point_type &x) : basepoint_(x) {
  check_direction();
}
template <typename T>
Line<T>::Line(point_type &&x) : basepoint_(std::move(x)) {
  check_direction();
}
template <typename T>
Line<T>::Line(const point_type &x, const point_type &v) : basepoint_(x), direction_(v) {
  check_direction();
}

template <typename T>
inline typename Line<T>::point_type Line<T>::push_forward(point_type x) const {  // TODO remove copy
  if (x.is_inf() || x.is_nan() || x.is_minus_inf()) return x;
  T t = this->push_forward2<T>(x);
  if (direction_.size() > 0) {
    for (std::size_t i = 0; i < x.size(); i++) x[i] = basepoint_[i] + t * direction_[i];
  } else {
    for (std::size_t i = 0; i < x.size(); i++) x[i] = basepoint_[i] + t;
  }
  return x;
}
template <typename T>
template <typename U>
inline U Line<T>::push_forward2(const point_type &x) const {
  constexpr const U inf =
      std::numeric_limits<U>::has_infinity ? std::numeric_limits<U>::infinity() : std::numeric_limits<U>::max();
  if (x.is_inf() || x.is_nan()) return inf;
  if (x.is_minus_inf()) return -inf;
  U t = -inf;
  if (direction_.size()) {
    for (std::size_t i = 0; i < x.size(); i++) {
      if (direction_[i] == 0) [[unlikely]] {
        if (x[i] < basepoint_[i])
          continue;
        else {
          return inf;
        }
      } else [[likely]] {
        t = std::max(t, (static_cast<U>(x[i]) - static_cast<U>(basepoint_[i])) / static_cast<U>((direction_[i])));
      }
    }
  } else {
    for (std::size_t i = 0; i < x.size(); i++) t = std::max(t, static_cast<U>(x[i]) - static_cast<U>(basepoint_[i]));
  }

  return t;
}
template <typename T>
template <typename U>
inline U Line<T>::push_forward2(const kcritical_point_type &x) const {
  constexpr const U inf =
      std::numeric_limits<U>::has_infinity ? std::numeric_limits<U>::infinity() : std::numeric_limits<U>::max();
  if (x.is_inf() || x.is_nan()) return inf;
  if (x.is_minus_inf()) return -inf;
  U t = inf;
  for (const auto &y : x) {
    t = std::min(t, this->push_forward2<U>(y));
  }
  return t;
}

template <typename T>
inline typename Line<T>::point_type Line<T>::push_back(point_type x) const {
  if (x.is_inf() || x.is_nan() || x.is_minus_inf()) return x;

  T t = this->push_back2(x);
  if (direction_.size() > 0) {
    for (std::size_t i = 0; i < x.size(); i++) x[i] = basepoint_[i] + t * direction_[i];
  } else
    for (std::size_t i = 0; i < x.size(); i++) x[i] = basepoint_[i] + t;

  return x;
}

template <typename T>
template <typename U>
inline U Line<T>::push_back2(const point_type &x) const {
  constexpr const U inf =
      std::numeric_limits<U>::has_infinity ? std::numeric_limits<U>::infinity() : std::numeric_limits<U>::max();
  if (x.is_inf()) return inf;
  if (x.is_minus_inf() || x.is_nan()) return -inf;
  U t = inf;

  if (direction_.size()) {
    for (std::size_t i = 0; i < x.size(); i++) {
      if (direction_[i] == 0) [[unlikely]] {
        if (x[i] > basepoint_[i])
          continue;
        else {
          return -inf;
        }
      } else [[likely]] {
        t = std::min(t, (static_cast<U>(x[i]) - static_cast<U>(basepoint_[i])) / static_cast<U>(direction_[i]));
      }
    }
  } else {
    for (std::size_t i = 0; i < x.size(); i++) t = std::min(t, static_cast<U>(x[i] - basepoint_[i]));
  }
  return t;
}

template <typename T>

template <typename U>
inline U Line<T>::push_back2(const kcritical_point_type &x) const {
  constexpr const U inf =
      std::numeric_limits<U>::has_infinity ? std::numeric_limits<U>::infinity() : std::numeric_limits<U>::max();
  if (x.is_inf()) return inf;
  if (x.is_minus_inf() || x.is_nan()) return -inf;
  U t = -inf;
  for (const auto &y : x) {
    t = std::max(t, this->push_back2<U>(y));
  }
  return t;
}

template <typename T>
inline int Line<T>::get_dim() const {
  return basepoint_.size();
}

template <typename T>
inline std::pair<T, T> Line<T>::get_bounds2(const Box<T> &box) const {
  return {this->push_forward2(box.get_bottom_corner()), this->push_back2(box.get_upper_corner())};
}

template <typename T>
inline std::pair<typename Line<T>::point_type, typename Line<T>::point_type> Line<T>::get_bounds(
    const Box<T> &box) const {
  return {this->push_forward(box.get_bottom_corner()), this->push_back(box.get_upper_corner())};
}
}  // namespace Gudhi::multi_persistence

#endif  // LINE_FILTRATION_TRANSLATION_H_INCLUDED

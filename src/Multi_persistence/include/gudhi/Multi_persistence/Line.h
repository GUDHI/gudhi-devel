/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Loiseaux
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - 2024/08 Hannah Schreiber: doc
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file Line.h
 * @author David Loiseaux
 * @brief Contains the @ref Gudhi::multi_persistence::Line class.
 */

#ifndef LINE_FILTRATION_TRANSLATION_H_INCLUDED
#define LINE_FILTRATION_TRANSLATION_H_INCLUDED

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <limits>

#include <gudhi/Debug_utils.h>
#include <gudhi/One_critical_filtration.h>
#include <gudhi/Multi_critical_filtration.h>
#include <gudhi/Multi_persistence/Box.h>

namespace Gudhi::multi_persistence {

/**
 * @class Line Line.h gudhi/Multi_persistence/Line.h
 * @ingroup multi_persistence
 *
 * @brief A line in \f$\mathbb R^n\f$, with some helpers to project points on it.
 *
 * @tparam T Type of the coordinate values. Has to follow the conditions of the template parameter of
 * @ref One_critical_filtration "".
 */
template <typename T>
class Line {
 public:
  /**
   * @brief Coordinates in \f$\mathbb R^n\f$.
   */
  using Point = Gudhi::multi_filtration::One_critical_filtration<T>;
  /**
   * @brief Set of coordinates in \f$\mathbb R^n\f$.
   */
  using K_critical_point = Gudhi::multi_filtration::Multi_critical_filtration<T>;

  /**
   * @brief Default constructor. Sets the number of coordinates to 0.
   */
  Line() : basePoint_(0), direction_(0) {}  // has to be explicitly set to 0, otherwise becomes -inf
  /**
   * @brief Constructs a line going through the given point with slope 1.
   *
   * @param x A point of the line.
   */
  Line(const Point &x) : basePoint_(x), direction_(0) {}  // default direction
  /**
   * @brief Constructs a line going through the given point with slope 1.
   *
   * @param x A point of the line. Will be moved.
   */
  Line(Point &&x) : basePoint_(std::move(x)), direction_(0) {}  // default direction
  /**
   * @brief Constructs a line going through the given point in the direction of the given vector.
   * If the vector has no coordinates, the slope is assumed to be 1.
   * Otherwise, the vector has to be non trivial and all its coordinates have to be positive.
   *
   * @param x A point of the line.
   * @param vector Direction of the line. Positive and non trivial.
   */
  Line(const Point &x, const Point &vector) : basePoint_(x), direction_(vector) { check_direction_(); }

  /**
   * @brief Returns the coordinates of the point on the line with "time" parameter `t`. That is, the point \f$ x \f$
   * such that \f$ x[i] = base\_point[i] + t \times direction[i] \f$ for all \f$ i \in [0, n - 1] \f$ with \f$ n \f$
   * the number of coordinates.
   */
  Point operator[](T t) const {
    GUDHI_CHECK(direction_.empty() || direction_.size() == basePoint_.size(),
                "Direction and base point do not have the same dimension.");

    if constexpr (std::numeric_limits<T>::has_quiet_NaN){   //to avoid windows error
      if (std::isnan(t)) return Point::nan();
    }
    if (t == Point::T_inf) return Point::inf();
    if (t == -Point::T_inf) return Point::minus_inf();

    Point x(basePoint_.size());

    if (direction_.size() > 0) {
      for (std::size_t i = 0; i < x.size(); i++) x[i] = basePoint_[i] + t * direction_[i];
    } else
      for (std::size_t i = 0; i < x.size(); i++) x[i] = basePoint_[i] + t;

    return x;
  }

  /**
   * @brief Translates the given line in the given direction.
   */
  friend Line &operator+=(Line &to_translate, const Point &v) {
    to_translate.basePoint_ += v;
    return to_translate;
  }

  /**
   * @brief Returns a reference to the current base point of the line.
   */
  Point &base_point() { return basePoint_; }
  /**
   * @brief Returns a const reference to the current base point of the line.
   */
  const Point &base_point() const { return basePoint_; }

  /**
   * @brief Returns a reference to the direction vector of the line.
   */
  Point &direction() { return direction_; }
  /**
   * @brief Returns a const reference to the direction vector of the line.
   */
  const Point &direction() const { return direction_; }

  // TODO: factorize forward and backward version by adding a `co` to One_critical_filtration?
  // Could make problems with One_critical_filtration being the type of basePoint_ and direction_

  /**
   * @brief Computes the "time" parameter \f$ t \f$ of the starting point \f$ p = base\_point + t \times direction \f$
   * of the intersection between the line and the closed positive cone originating at `x`.
   *
   * @tparam U Type of the time parameter.
   * @param x Origin of the closed positive cone.
   */
  template <typename U = T>
  U compute_forward_intersection(const Point &x) const {
    GUDHI_CHECK(direction_.empty() || direction_.size() == x.size(), "x has not as many parameters as the line.");

    constexpr const U inf =
        std::numeric_limits<U>::has_infinity ? std::numeric_limits<U>::infinity() : std::numeric_limits<U>::max();
    if (x.is_plus_inf() || x.is_nan()) return inf;
    if (x.is_minus_inf()) return -inf;
    U t = -inf;
    if (direction_.size()) {
      for (std::size_t i = 0; i < x.size(); i++) {
        if (direction_[i] == 0) {
          if (x[i] > basePoint_[i]) return inf;
        } else {
          t = std::max(t, (static_cast<U>(x[i]) - static_cast<U>(basePoint_[i])) / static_cast<U>((direction_[i])));
        }
      }
    } else {
      for (std::size_t i = 0; i < x.size(); i++) t = std::max(t, static_cast<U>(x[i]) - static_cast<U>(basePoint_[i]));
    }

    return t;
  }

  /**
   * @brief Computes the "time" parameter \f$ t \f$ of the starting point \f$ p = base\_point + t \times direction \f$
   * of the intersection between the line and the union of closed positive cones originating at the points in `x`.
   *
   * @tparam U Type of the time parameter.
   * @param x Set of origins for the closed positive cones.
   */
  template <typename U = T>
  U compute_forward_intersection(const K_critical_point &x) const {
    constexpr const U inf =
        std::numeric_limits<U>::has_infinity ? std::numeric_limits<U>::infinity() : std::numeric_limits<U>::max();
    if (x.is_plus_inf() || x.is_nan()) return inf;
    if (x.is_minus_inf()) return -inf;
    U t = inf;
    for (const auto &y : x) {
      t = std::min(t, compute_forward_intersection<U>(y));
    }
    return t;
  }

  /**
   * @brief Computes the "time" parameter \f$ t \f$ of the starting point \f$ p = base\_point + t \times direction \f$
   * of the intersection between the line and the open negative cone originating at `x`.
   *
   * @tparam U Type of the time parameter.
   * @param x Origin of the open negative cone.
   */
  template <typename U = T>
  U compute_backward_intersection(const Point &x) const {
    constexpr const U inf =
        std::numeric_limits<U>::has_infinity ? std::numeric_limits<U>::infinity() : std::numeric_limits<U>::max();
    if (x.is_plus_inf()) return inf;
    if (x.is_minus_inf() || x.is_nan()) return -inf;
    U t = inf;

    if (direction_.size()) {
      for (std::size_t i = 0; i < x.size(); i++) {
        if (direction_[i] == 0) {
          if (x[i] <= basePoint_[i]) return -inf;
        } else {
          t = std::min(t, (static_cast<U>(x[i]) - static_cast<U>(basePoint_[i])) / static_cast<U>(direction_[i]));
        }
      }
    } else {
      for (std::size_t i = 0; i < x.size(); i++) t = std::min(t, static_cast<U>(x[i] - basePoint_[i]));
    }
    return t;
  }

  /**
   * @brief Computes the "time" parameter \f$ t \f$ of the starting point \f$ p = base\_point + t \times direction \f$
   * of the intersection between the line and the union of open negative cones originating at the points in `x`.
   *
   * @tparam U Type of the time parameter.
   * @param x Set of origins for the open negative cones.
   */
  template <typename U = T>
  U compute_backward_intersection(const K_critical_point &x) const {
    constexpr const U inf =
        std::numeric_limits<U>::has_infinity ? std::numeric_limits<U>::infinity() : std::numeric_limits<U>::max();
    if (x.is_plus_inf()) return inf;
    if (x.is_minus_inf() || x.is_nan()) return -inf;
    U t = -inf;
    for (const auto &y : x) {
      t = std::max(t, compute_backward_intersection<U>(y));
    }
    return t;
  }

  /**
   * @brief Given a box, returns "time" parameter of the intersection of this box and the line.
   *
   * @param box Box to intersect.
   * @return A pair representing the two bounding points of the intersection, such that the first element is the
   * smallest of the two. If the box and the line do not intersect or the box is trivial, returns the pair {inf, -inf}.
   */
  std::pair<T, T> get_bounds(const Box<T> &box) const {
    if (box.is_trivial()) return {Point::T_inf, -Point::T_inf};

    T bottom = compute_forward_intersection(box.get_lower_corner());
    T top = compute_backward_intersection(box.get_upper_corner());

    if (bottom > top) return {Point::T_inf, -Point::T_inf};  // no intersection

    return {bottom, top};
  }

 private:
  Point basePoint_; /**< Any point on the line. */
  Point direction_; /**< Direction of the line. */

  /**
   * @brief Checks that the arguments define a correct and positively slopped line.
   */
  void check_direction_() const {
    if (direction_.size() == 0) return;  // default slope

    bool is_trivial = true;
    for (T v : direction_) {
      if (v) {
        is_trivial = false;
      }
      if (v < 0) {
        throw std::invalid_argument("Direction should have positive entries.");
      }
    }
    if (is_trivial) {
      throw std::invalid_argument("Direction should have at least one non-trivial entry.");
    }
    if (direction_.size() != basePoint_.size())
      throw std::invalid_argument("The dimensions of base point and direction are not equal.");
  }
};

}  // namespace Gudhi::multi_persistence

#endif  // LINE_FILTRATION_TRANSLATION_H_INCLUDED

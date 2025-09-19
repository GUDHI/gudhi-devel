/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Loiseaux
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - 2024/08 Hannah Schreiber: documentation
 *      - 2025/03 Hannah Schreiber: Change of point types.
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file Line.h
 * @author David Loiseaux
 * @brief Contains the @ref Gudhi::multi_persistence::Line class.
 */

#ifndef MP_LINE_FILTRATION_H_INCLUDED
#define MP_LINE_FILTRATION_H_INCLUDED

#include <cmath>
#include <cstddef>
#include <stdexcept>

#include <gudhi/Debug_utils.h>
#include <gudhi/Multi_persistence/Box.h>
#include <gudhi/Multi_persistence/Point.h>
#include <gudhi/Multi_filtration/multi_filtration_utils.h>

namespace Gudhi {
namespace multi_persistence {

/**
 * @class Line Line.h gudhi/Multi_persistence/Line.h
 * @ingroup multi_persistence
 *
 * @brief A line in \f$\mathbb R^n\f$, with some helpers to project points on it.
 *
 * @tparam T Type of the coordinate values.
 */
template <typename T>
class Line
{
 public:
  /**
   * @brief Coordinates in \f$\mathbb R^n\f$.
   */
  using Point_t = Point<T>;

  /**
   * @brief Default constructor. Sets the number of coordinates to 0.
   */
  Line() = default;

  /**
   * @brief Constructs a line going through the given point with slope 1.
   *
   * @param x A point of the line.
   */
  Line(const Point_t &x) : basePoint_(x), direction_() {}  // default direction

  /**
   * @brief Constructs a line going through the given point with slope 1.
   *
   * @param x A point of the line. Will be moved.
   */
  Line(Point_t &&x) : basePoint_(std::move(x)), direction_() {}  // default direction

  /**
   * @brief Constructs a line going through the given point in the direction of the given vector.
   * If the vector has no coordinates, the slope is assumed to be 1.
   * Otherwise, the vector has to be non trivial and all its coordinates have to be positive.
   *
   * @param x A point of the line.
   * @param vector Direction of the line. Positive and non trivial.
   */
  Line(const Point_t &x, const Point_t &vector) : basePoint_(x), direction_(vector) { _check_direction(); }

  /**
   * @brief Returns the coordinates of the point on the line with "time" parameter `t`. That is, the point \f$ x \f$
   * such that \f$ x[i] = base\_point[i] + t \times direction[i] \f$ for all \f$ i \in [0, n - 1] \f$ with \f$ n \f$
   * the number of coordinates.
   */
  Point_t operator[](T t) const
  {
    GUDHI_CHECK(direction_.size() == 0 || direction_.size() == basePoint_.size(),
                "Direction and base point do not have the same dimension.");

    if (Gudhi::multi_filtration::_is_nan(t) || t == Point_t::T_inf || t == Point_t::T_m_inf)
      return Point_t(basePoint_.size(), t);

    Point_t x(basePoint_.size());

    if (direction_.size() > 0) {
      for (std::size_t i = 0; i < x.size(); i++) x[i] = basePoint_[i] + t * direction_[i];
    } else
      for (std::size_t i = 0; i < x.size(); i++) x[i] = basePoint_[i] + t;

    return x;
  }

  /**
   * @brief Translates the given line in the given direction.
   */
  friend Line &operator+=(Line &to_translate, const Point_t &v)
  {
    to_translate.basePoint_ += v;
    return to_translate;
  }

  /**
   * @brief Returns a reference to the current base point of the line.
   */
  Point_t &base_point() { return basePoint_; }

  /**
   * @brief Returns a const reference to the current base point of the line.
   */
  const Point_t &base_point() const { return basePoint_; }

  /**
   * @brief Returns a reference to the direction vector of the line.
   */
  Point_t &direction() { return direction_; }

  /**
   * @brief Returns a const reference to the direction vector of the line.
   */
  const Point_t &direction() const { return direction_; }

  /**
   * @brief Computes the "time" parameter \f$ t \f$ of the starting point \f$ p = base\_point + t \times direction \f$
   * of the intersection between the line and the closed positive cone originating at point `x`.
   *
   * @tparam U Type of the time parameter.
   * @param x Origin of the closed positive cone.
   */
  template <typename U = T>
  U compute_forward_intersection(const Point_t &x) const
  {
    GUDHI_CHECK(basePoint_.size() == x.size(), "x has not as many parameters as the line.");

    constexpr const U inf = Point<U>::T_inf;

    U t = Point<U>::T_m_inf;
    for (unsigned int p = 0; p < x.size(); ++p) {
      if (Gudhi::multi_filtration::_is_nan(x[p])) return inf;
      auto div = direction_.size() == 0 ? 1 : direction_[p];
      if (div == 0) {
        if (x[p] > basePoint_[p]) return inf;
      } else {
        t = std::max(t, (static_cast<U>(x[p]) - static_cast<U>(basePoint_[p])) / static_cast<U>(div));
      }
    }

    return t;
  }

  /**
   * @brief Computes the "time" parameter \f$ t \f$ of the starting point \f$ p = base\_point + t \times direction \f$
   * of the intersection between the line and the union of closed positive cones originating at
   * multi-parameter filtration value `x`. If `x` contains a NaN value, returns +infinity.
   *
   * @tparam U Type of the time parameter.
   * @tparam FiltrationValue Type of a multi-parameter filtration value. Has to implement the following methods:
   * `num_parameters`, `num_generators`, `operator()(generator_index, parameter_index)`.
   * See @ref Gudhi::multi_filtration::Multi_parameter_filtration for an example.
   * @param x Origin of the closed positive cones.
   */
  template <typename U = T, class FiltrationValue>
  U compute_forward_intersection(const FiltrationValue &x) const
  {
    GUDHI_CHECK(basePoint_.size() == x.num_parameters(), "x has not as many parameters as the line.");

    constexpr const U inf = Point<U>::T_inf;

    U t = inf;
    for (unsigned int g = 0; g < x.num_generators(); ++g) {
      U tmp = Point<U>::T_m_inf;
      for (unsigned int p = 0; p < x.num_parameters(); ++p) {
        if (Gudhi::multi_filtration::_is_nan(x(g, p))) return inf;
        auto div = direction_.size() == 0 ? 1 : direction_[p];
        if (div == 0) {
          if (x(g, p) > basePoint_[p]) tmp = inf;
        } else {
          tmp = std::max(tmp, (static_cast<U>(x(g, p)) - static_cast<U>(basePoint_[p])) / static_cast<U>(div));
        }
      }
      t = std::min(t, tmp);
    }

    return t;
  }

  /**
   * @brief Computes the "time" parameter \f$ t \f$ of the starting point \f$ p = base\_point + t \times direction \f$
   * of the intersection between the line and the open negative cone originating at point `x`.
   *
   * @tparam U Type of the time parameter.
   * @param x Origin of the open negative cone.
   */
  template <typename U = T>
  U compute_backward_intersection(const Point_t &x) const
  {
    GUDHI_CHECK(basePoint_.size() == x.size(), "x has not as many parameters as the line.");

    constexpr const U m_inf = Point<U>::T_m_inf;

    U t = Point<U>::T_inf;
    for (unsigned int p = 0; p < x.size(); ++p) {
      if (Gudhi::multi_filtration::_is_nan(x[p])) return m_inf;
      auto div = direction_.size() == 0 ? 1 : direction_[p];
      if (div == 0) {
        if (x[p] <= basePoint_[p]) return m_inf;
      } else {
        t = std::min(t, (static_cast<U>(x[p]) - static_cast<U>(basePoint_[p])) / static_cast<U>(div));
      }
    }

    return t;
  }

  /**
   * @brief Computes the "time" parameter \f$ t \f$ of the starting point \f$ p = base\_point + t \times direction \f$
   * of the intersection between the line and the union of open negative cones originating at
   * multi-parameter filtration value `x`. If `x` contains a NaN value, returns -infinity.
   *
   * @tparam U Type of the time parameter.
   * @tparam FiltrationValue Type of a multi-parameter filtration value. Has to implement the following methods:
   * `num_parameters`, `num_generators`, `operator()(generator_index, parameter_index)`.
   * @param x Origin of the open negative cones.
   */
  template <typename U = T, class FiltrationValue>
  U compute_backward_intersection(const FiltrationValue &x) const
  {
    GUDHI_CHECK(basePoint_.size() == x.num_parameters(), "x has not as many parameters as the line.");

    constexpr const U m_inf = Point<U>::T_m_inf;

    U t = m_inf;
    for (unsigned int g = 0; g < x.num_generators(); ++g) {
      U tmp = Point<U>::T_inf;
      for (unsigned int p = 0; p < x.num_parameters(); ++p) {
        if (Gudhi::multi_filtration::_is_nan(x(g, p))) return m_inf;
        auto div = direction_.size() == 0 ? 1 : direction_[p];
        if (div == 0) {
          if (x(g, p) <= basePoint_[p]) tmp = m_inf;
        } else {
          tmp = std::min(tmp, (static_cast<U>(x(g, p)) - static_cast<U>(basePoint_[p])) / static_cast<U>(div));
        }
      }
      t = std::max(t, tmp);
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
  std::pair<T, T> get_bounds(const Box<T> &box) const
  {
    if (box.is_trivial()) return {Point_t::T_inf, Point_t::T_m_inf};

    T bottom = compute_forward_intersection(box.get_lower_corner());
    T top = compute_backward_intersection(box.get_upper_corner());

    if (bottom > top) return {Point_t::T_inf, Point_t::T_m_inf};  // no intersection

    return {bottom, top};
  }

 private:
  Point_t basePoint_; /**< Any point on the line. */
  Point_t direction_; /**< Direction of the line. */

  /**
   * @brief Checks that the arguments define a correct and positively slopped line.
   */
  void _check_direction() const
  {
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

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // MP_LINE_FILTRATION_H_INCLUDED

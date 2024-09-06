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
 * @file Box.h
 * @author David Loiseaux
 * @brief Contains the @ref Gudhi::multi_persistence::Box class.
 */

#ifndef BOX_H_INCLUDED
#define BOX_H_INCLUDED

#include <ostream>  //std::ostream

#include <gudhi/Debug_utils.h>
#include <gudhi/One_critical_filtration.h>

namespace Gudhi::multi_persistence {

/**
 * @class Box Box.h gudhi/Multi_persistence/Box.h
 * @ingroup multi_persistence
 *
 * @brief Simple box in \f$\mathbb R^n\f$ defined by two diametrically opposite corners.
 *
 * @tparam T Type of the coordinates of the Box. Has to follow the conditions of the template parameter of
 * @ref One_critical_filtration "".
 */
template <typename T>
class Box {
 public:
  using Point = Gudhi::multi_filtration::One_critical_filtration<T>; /**< Type of a point in \f$\mathbb R^n\f$. */

  /**
   * @brief Default constructor. Constructs a trivial box with corners at minus infinity.
   */
  Box() {}

  /**
   * @brief Constructs a box from the two given corners. Assumes that \f$ lowerCorner \le @p upperCorner \f$ and
   * if both are finite values, they have the same dimension.
   *
   * @param lowerCorner First corner of the box. Has to be smaller than `upperCorner`.
   * @param upperCorner Second corner of the box. Has to be greater than `lowerCorner`.
   */
  Box(const Point &lowerCorner, const Point &upperCorner) : lowerCorner_(lowerCorner), upperCorner_(upperCorner) {
    GUDHI_CHECK(lowerCorner.size() == upperCorner.size() && lowerCorner <= upperCorner, "This box is trivial !");
  }

  /**
   * @brief Constructs a box from the two given corners. Assumes that \f$ box.first \le @p box.second \f$ and
   * if both are finite values, they have the same dimension.
   *
   * @param box Pair of corners defining the wished box.
   */
  Box(const std::pair<Point, Point> &box) : Box(box.first, box.second) {}

  /**
   * @brief Returns the lowest of both defining corners.
   */
  const Point &get_lower_corner() const { return lowerCorner_; }

  /**
   * @brief Returns the lowest of both defining corners.
   */
  Point &get_lower_corner() { return lowerCorner_; }

  /**
   * @brief Returns the greatest of both defining corners.
   */
  Point &get_upper_corner() { return upperCorner_; }

  /**
   * @brief Returns the greatest of both defining corners.
   */
  const Point &get_upper_corner() const { return upperCorner_; }

  /**
   * @brief Returns a pair of const references to both defining corners.
   */
  std::pair<const Point &, const Point &> get_bounding_corners() const { return {lowerCorner_, upperCorner_}; }

  /**
   * @brief Returns a pair of references to both defining corners.
   */
  std::pair<Point &, Point &> get_bounding_corners() { return {lowerCorner_, upperCorner_}; }

  /**
   * @brief Returns true if and only if one of the following is true:
   * - one of the corners is empty
   * - one of the corners has value NaN
   * - both corners have value infinity
   * - both corners have value minus infinity
   * - both corners are finite but don't have the same dimension.
   */
  bool is_trivial() const {
    return lowerCorner_.empty() || upperCorner_.empty() || lowerCorner_.is_nan() || upperCorner_.is_nan() ||
           (lowerCorner_.is_plus_inf() && upperCorner_.is_plus_inf()) ||
           (lowerCorner_.is_minus_inf() && upperCorner_.is_minus_inf()) ||
           (lowerCorner_.is_finite() && upperCorner_.is_finite() &&
            lowerCorner_.num_parameters() != upperCorner_.num_parameters());
  }

  /**
   * @brief Returns true if and only if the given point is inside the box.
   * If the box is not {-infinity, infinity} and the given point is finite, but has not the same dimension
   * than the box, the point is considered outside.
   */
  bool contains(const Point &point) const {
    if (point.is_nan() || is_trivial()) return false;
    if (point.is_plus_inf()) return upperCorner_.is_plus_inf();
    if (point.is_minus_inf()) return lowerCorner_.is_minus_inf();

    if ((lowerCorner_.is_finite() && point.size() != lowerCorner_.size()) ||
        (upperCorner_.is_finite() && point.size() != upperCorner_.size())) {
      // TODO: make it a warning, with future GUDHI_CHECK version?
      // std::cerr << "Box and point are not of the same dimension." << std::endl;
      return false;
    }

    return lowerCorner_ <= point && point <= upperCorner_;
  }

  /**
   * @brief Returns the dimension of the box. If the box is trivial or both corners are infinite, the dimension is 0.
   */
  std::size_t dimension() const {
    if (is_trivial()) return 0;
    if (lowerCorner_.is_minus_inf() && upperCorner_.is_plus_inf()) return 0;  // not so sure what we want to do here
    return lowerCorner_.is_finite() ? lowerCorner_.size() : upperCorner_.size();
  }

  /**
   * @brief Inflates the box by delta.
   *
   * @param delta Inflation coefficient.
   */
  void inflate(T delta) {
    lowerCorner_ -= delta;
    upperCorner_ += delta;
  }

  /**
   * @brief Outstream operator.
   */
  friend std::ostream &operator<<(std::ostream &os, const Box<T> &box) {
    os << "Box -- Bottom corner : ";
    os << box.get_lower_corner();
    os << ", Top corner : ";
    os << box.get_upper_corner();
    return os;
  }

 private:
  Point lowerCorner_; /**< Lowest of defining corners. */
  Point upperCorner_; /**< Greatest of defining corners. */
};

}  // namespace Gudhi::multi_persistence

#endif  // BOX_H_INCLUDED

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
 * @file Box.h
 * @author David Loiseaux
 * @brief Contains the @ref Gudhi::multi_persistence::Box class.
 */

#ifndef MP_BOX_H_INCLUDED
#define MP_BOX_H_INCLUDED

#include <ostream>  //std::ostream

#include <gudhi/Debug_utils.h>
#include <gudhi/Multi_persistence/Point.h>
#include <gudhi/Multi_filtration/multi_filtration_utils.h>

namespace Gudhi {
namespace multi_persistence {

/**
 * @class Box Box.h gudhi/Multi_persistence/Box.h
 * @ingroup multi_persistence
 *
 * @brief Simple box in \f$\mathbb R^n\f$ defined by two diametrically opposite corners.
 *
 * @tparam T Type of the coordinates of the Box.
 */
template <typename T>
class Box
{
 public:
  using Point_t = Point<T>; /**< Type of a point in \f$\mathbb R^n\f$. */

  /**
   * @brief Default constructor. Constructs a trivial box with corners at minus infinity.
   */
  Box() = default;

  /**
   * @brief Constructs a box from the two given corners. Assumes that \f$ lowerCorner \le @p upperCorner \f$ and
   * if both are finite values, they have the same dimension.
   *
   * @param lowerCorner First corner of the box. Has to be smaller than `upperCorner`.
   * @param upperCorner Second corner of the box. Has to be greater than `lowerCorner`.
   */
  Box(const Point_t &lowerCorner, const Point_t &upperCorner) : lowerCorner_(lowerCorner), upperCorner_(upperCorner)
  {
    GUDHI_CHECK(lowerCorner.size() == upperCorner.size(), "The two corners of the box don't have the same dimension.");
    GUDHI_CHECK(lowerCorner <= upperCorner, "The first corner is not smaller than the second.");
  }

  /**
   * @brief Constructs a box from the two given corners. Assumes that \f$ box.first \le @p box.second \f$ and
   * if both are finite values, they have the same dimension.
   *
   * @param box Pair of corners defining the wished box.
   */
  Box(const std::pair<Point_t, Point_t> &box) : Box(box.first, box.second) {}

  /**
   * @brief Returns the lowest of both defining corners.
   */
  const Point_t &get_lower_corner() const { return lowerCorner_; }

  /**
   * @brief Returns the lowest of both defining corners.
   */
  Point_t &get_lower_corner() { return lowerCorner_; }

  /**
   * @brief Returns the greatest of both defining corners.
   */
  Point_t &get_upper_corner() { return upperCorner_; }

  /**
   * @brief Returns the greatest of both defining corners.
   */
  const Point_t &get_upper_corner() const { return upperCorner_; }

  /**
   * @brief Returns a pair of const references to both defining corners.
   */
  std::pair<const Point_t &, const Point_t &> get_bounding_corners() const { return {lowerCorner_, upperCorner_}; }

  /**
   * @brief Returns a pair of references to both defining corners.
   */
  std::pair<Point_t &, Point_t &> get_bounding_corners() { return {lowerCorner_, upperCorner_}; }

  /**
   * @brief Returns true if and only if one of the following is true:
   * - one of the corners is empty
   * - one of the corners contains the value NaN
   * - both corners have value infinity
   * - both corners have value minus infinity
   * - both corners don't have the same dimension.
   */
  [[nodiscard]] bool is_trivial() const
  {
    if (lowerCorner_.size() == 0 || upperCorner_.size() == 0) return true;
    if (lowerCorner_.size() != upperCorner_.size()) return false;  // should not happen?

    T inf = Point_t::T_inf;

    bool lowerIsInf = true, lowerIsMinusInf = true;
    bool upperIsInf = true, upperIsMinusInf = true;
    for (unsigned int i = 0; i < lowerCorner_.size(); ++i) {
      T lc = lowerCorner_[i];
      T uc = upperCorner_[i];
      if (Gudhi::multi_filtration::_is_nan(lc) || Gudhi::multi_filtration::_is_nan(uc)) return true;
      if (lc != inf) lowerIsInf = false;
      if (lc != -inf) lowerIsMinusInf = false;
      if (uc != inf) upperIsInf = false;
      if (uc != -inf) upperIsMinusInf = false;
      if ((!lowerIsInf && !lowerIsMinusInf) || (!upperIsInf && !upperIsMinusInf)) return false;
    }
    return (lowerIsInf && upperIsInf) || (lowerIsMinusInf && upperIsMinusInf);
  }

  /**
   * @brief Returns true if and only if the given point is inside the box.
   */
  bool contains(const Point_t &point) const
  {
    GUDHI_CHECK(point.size() == lowerCorner_.size(), "Point should not have a different dimension than the box.");

    for (unsigned int i = 0; i < point.size(); ++i) {
      T lc = lowerCorner_[i];
      T uc = upperCorner_[i];
      T p = point[i];
      if (Gudhi::multi_filtration::_is_nan(p) || Gudhi::multi_filtration::_is_nan(lc) ||
          Gudhi::multi_filtration::_is_nan(uc))
        return false;
      if (lc > p || uc < p) return false;
    }
    return true;
  }

  /**
   * @brief Returns the dimension of the box.
   */
  [[nodiscard]] std::size_t dimension() const { return lowerCorner_.size(); }

  /**
   * @brief Inflates the box by delta.
   *
   * @param delta Inflation coefficient.
   */
  void inflate(T delta)
  {
    lowerCorner_ -= delta;
    upperCorner_ += delta;
  }

  /**
   * @brief Outstream operator.
   */
  friend std::ostream &operator<<(std::ostream &os, const Box<T> &box)
  {
    os << "Box -- Bottom corner : ";
    os << box.get_lower_corner();
    os << ", Top corner : ";
    os << box.get_upper_corner();
    return os;
  }

 private:
  Point_t lowerCorner_; /**< Lowest of defining corners. */
  Point_t upperCorner_; /**< Greatest of defining corners. */
};

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // MP_BOX_H_INCLUDED

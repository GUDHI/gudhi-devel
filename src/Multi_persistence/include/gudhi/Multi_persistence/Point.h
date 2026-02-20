/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file Point.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::multi_persistence::Point class.
 */

#ifndef MP_POINT_H_INCLUDED
#define MP_POINT_H_INCLUDED

#include <algorithm>  //std::for_each
#include <ostream>    //std::ostream
#include <vector>

#include <gudhi/Debug_utils.h>
#include <gudhi/Multi_filtration/multi_filtration_utils.h>

namespace Gudhi {
namespace multi_persistence {

/**
 * @class Point Point.h gudhi/Multi_persistence/Point.h
 * @ingroup multi_persistence
 *
 * @brief Simple point class with possibility to translate and multiply coordinates.
 *
 * @tparam T Type of the coordinates of the point.
 */
template <typename T>
class Point
{
 public:
  using Container = std::vector<T>; /**< Type of coordinate container. */

  using value_type = typename Container::value_type;                         /**< Type of coordinates. */
  using allocator_type = typename Container::allocator_type;                 /**< Allocator type. */
  using size_type = typename Container::size_type;                           /**< Size type. */
  using difference_type = typename Container::difference_type;               /**< Difference type. */
  using reference = typename Container::reference;                           /**< Coordinate reference type. */
  using const_reference = typename Container::const_reference;               /**< Coordinate const reference type. */
  using pointer = typename Container::pointer;                               /**< Coordinate pointer type. */
  using iterator = typename Container::iterator;                             /**< Coordinate iterator type. */
  using const_iterator = typename Container::const_iterator;                 /**< Coordinate const iterator type. */
  using reverse_iterator = typename Container::reverse_iterator;             /**< Coordinate reverse iterator type. */
  using const_reverse_iterator = typename Container::const_reverse_iterator; /**< Coordinate reverse iterator type. */

  /**
   * @brief Default constructor.
   */
  Point() = default;

  /**
   * @brief Constructs a new point with given number of coordinates. All values are default initialized.
   */
  explicit Point(size_type count) : coordinates_(count) {}

  /**
   * @brief Constructs a new point with given number of coordinates. All values are initialized with given value.
   */
  Point(size_type count, const T &value) : coordinates_(count, value) {}

  /**
   * @brief Constructs a new point from the given range.
   *
   * @tparam InputIt Iterator type that must follow the condition of the corresponding vector constructor.
   */
  template <class InputIt>
  Point(InputIt first, InputIt last) : coordinates_(first, last)
  {}

  /**
   * @brief Constructs a new point from the given range.
   */
  Point(std::initializer_list<T> init) : coordinates_(init.begin(), init.end()) {}

  /**
   * @brief Constructs a new point by copying the given container.
   */
  Point(const Container &init) : coordinates_(init) {}

  /**
   * @brief Constructs a new point by moving the given container.
   */
  Point(Container &&init) : coordinates_(std::move(init)) {}

  /**
   * @brief Assign copy operator.
   */
  Point &operator=(std::initializer_list<value_type> ilist)
  {
    coordinates_ = Container(ilist.begin(), ilist.end());
    return *this;
  }

  operator std::vector<T>() const { return coordinates_; }

  /**
   * @brief At operator.
   */
  reference at(size_type pos) { return coordinates_.at(pos); }

  /**
   * @brief At operator.
   */
  const_reference at(size_type pos) const { return coordinates_.at(pos); }

  /**
   * @brief operator[].
   */
  reference operator[](size_type pos) { return coordinates_[pos]; }

  /**
   * @brief operator[].
   */
  const_reference operator[](size_type pos) const { return coordinates_[pos]; }

  /**
   * @brief Front operator.
   */
  reference front() { return coordinates_.front(); }

  /**
   * @brief Front operator.
   */
  const_reference front() const { return coordinates_.front(); }

  /**
   * @brief Back operator.
   */
  reference back() { return coordinates_.back(); }

  /**
   * @brief Back operator.
   */
  const_reference back() const { return coordinates_.back(); }

  /**
   * @brief Data operator.
   */
  T *data() noexcept { return coordinates_.data(); }

  /**
   * @brief Data operator.
   */
  const T *data() const noexcept { return coordinates_.data(); }

  /**
   * @brief begin.
   */
  iterator begin() noexcept { return coordinates_.begin(); }

  /**
   * @brief begin.
   */
  const_iterator begin() const noexcept { return coordinates_.begin(); }

  /**
   * @brief cbegin.
   */
  const_iterator cbegin() const noexcept { return coordinates_.cbegin(); }

  /**
   * @brief end.
   */
  iterator end() noexcept { return coordinates_.end(); }

  /**
   * @brief end.
   */
  const_iterator end() const noexcept { return coordinates_.end(); }

  /**
   * @brief cend.
   */
  const_iterator cend() const noexcept { return coordinates_.cend(); }

  /**
   * @brief rbegin.
   */
  reverse_iterator rbegin() noexcept { return coordinates_.rbegin(); }

  /**
   * @brief rbegin.
   */
  const_reverse_iterator rbegin() const noexcept { return coordinates_.rbegin(); }

  /**
   * @brief crbegin.
   */
  const_reverse_iterator crbegin() const noexcept { return coordinates_.crbegin(); }

  /**
   * @brief rend.
   */
  reverse_iterator rend() noexcept { return coordinates_.rend(); }

  /**
   * @brief rend.
   */
  const_reverse_iterator rend() const noexcept { return coordinates_.rend(); }

  /**
   * @brief crend.
   */
  const_reverse_iterator crend() const noexcept { return coordinates_.crend(); }

  /**
   * @brief Number of coordinates.
   */
  size_type size() const noexcept { return coordinates_.size(); }

  /**
   * @brief Swap operator.
   */
  void swap(Point &other) noexcept(std::allocator_traits<allocator_type>::propagate_on_container_swap::value ||
                                   std::allocator_traits<allocator_type>::is_always_equal::value)
  {
    coordinates_.swap(other.coordinates_);
  }

  /**
   * @brief Swap operator.
   */
  friend void swap(Point &p1, Point &p2) noexcept { p1.coordinates_.swap(p2.coordinates_); }

  /**
   * @brief Outstream operator.
   */
  friend std::ostream &operator<<(std::ostream &os, const Point<T> &point)
  {
    os << "[ ";
    for (const T &p : point) os << p << " ";
    os << " ]";
    return os;
  }

  /**
   * @brief Returns true if and only if all coordinates of the first argument are strictly smaller than the ones of
   * the second argument.
   *
   * Note that this order is not total.
   */
  friend bool operator<(const Point &a, const Point &b)
  {
    if (&a == &b) return false;
    GUDHI_CHECK(a.size() == b.size(), "Cannot compare two points with different number of coordinates.");
    bool isSame = true;
    for (size_type i = 0U; i < a.size(); ++i) {
      if (a[i] > b[i] || Gudhi::multi_filtration::_is_nan(a[i]) || Gudhi::multi_filtration::_is_nan(b[i])) return false;
      if (isSame && a[i] != b[i]) isSame = false;
    }
    return !isSame;
  }

  /**
   * @brief Returns true if and only if all coordinates of the first argument are smaller or equal to the ones of
   * the second argument.
   *
   * Note that this order is not total.
   */
  friend bool operator<=(const Point &a, const Point &b)
  {
    if (&a == &b) return true;
    GUDHI_CHECK(a.size() == b.size(), "Cannot compare two points with different number of coordinates.");
    for (size_type i = 0U; i < a.size(); ++i) {
      if (a[i] > b[i] || Gudhi::multi_filtration::_is_nan(a[i]) || Gudhi::multi_filtration::_is_nan(b[i])) return false;
    }
    return true;
  }

  /**
   * @brief Returns true if and only if all coordinates of the first argument are strictly greater than the ones of
   * the second argument.
   *
   * Note that this order is not total.
   */
  friend bool operator>(const Point &a, const Point &b) { return b < a; }

  /**
   * @brief Returns true if and only if all coordinates of the first argument are greater or equal to the ones of
   * the second argument.
   *
   * Note that this order is not total.
   */
  friend bool operator>=(const Point &a, const Point &b) { return b <= a; }

  /**
   * @brief Returns `true` if and only if all coordinates of both arguments are equal.
   */
  friend bool operator==(const Point &a, const Point &b)
  {
    if (&a == &b) return true;
    if (a.size() != b.size()) return false;
    return a.coordinates_ == b.coordinates_;
  }

  /**
   * @brief Returns `true` if and only if \f$ a == b \f$ returns `false`.
   */
  friend bool operator!=(const Point &a, const Point &b) { return !(a == b); }

  // ARITHMETIC OPERATORS

  // opposite
  /**
   * @brief Returns a filtration value such that an entry at index \f$ p \f$ is equal to \f$ -f(p) \f$.
   *
   * Used conventions:
   * - \f$ -NaN = NaN \f$.
   *
   * @param f Value to opposite.
   * @return The opposite of @p f.
   */
  friend Point operator-(const Point &f)
  {
    Point result(f.coordinates_);
    std::for_each(result.begin(), result.end(), [](T &v) { v = -v; });
    return result;
  }

  // subtraction
  /**
   * @brief Returns a filtration value such that an entry at index \f$ (p) \f$ is equal to \f$ f(p) - val(p) \f$.
   * Both points have to have the same dimension.
   *
   * Used conventions:
   * - \f$ inf - inf = NaN \f$,
   * - \f$ -inf - (-inf) = NaN \f$,
   * - \f$ NaN - b = NaN \f$,
   * - \f$ a - NaN = NaN \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @param f First element of the subtraction.
   * @param val Second element of the subtraction.
   */
  template <typename U = T>
  friend Point operator-(Point f, const Point<U> &val)
  {
    GUDHI_CHECK(f.size() == val.size(), "Cannot translate point by a direction with different dimension.");
    f -= val;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (p) \f$ is equal to \f$ f(p) - val \f$.
   *
   * Used conventions:
   * - \f$ inf - inf = NaN \f$,
   * - \f$ -inf - (-inf) = NaN \f$,
   * - \f$ NaN - b = NaN \f$,
   * - \f$ a - NaN = NaN \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @param f First element of the subtraction.
   * @param val Second element of the subtraction.
   */
  friend Point operator-(Point f, const T &val)
  {
    f -= val;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (p) \f$ is equal to \f$ val - f(p) \f$.
   *
   * Used conventions:
   * - \f$ inf - inf = NaN \f$,
   * - \f$ -inf - (-inf) = NaN \f$,
   * - \f$ NaN - b = NaN \f$,
   * - \f$ a - NaN = NaN \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @param val First element of the subtraction.
   * @param f Second element of the subtraction.
   */
  friend Point operator-(const T &val, Point f)
  {
    f._apply_operation(val, [](T &valF, const T &valR) {
      valF = -valF;
      Gudhi::multi_filtration::_add(valF, valR);
    });
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ (p) \f$ is equal to \f$ f(p) - val(p) \f$.
   * Both points have to have the same dimension.
   *
   * Used conventions:
   * - \f$ inf - inf = NaN \f$,
   * - \f$ -inf - (-inf) = NaN \f$,
   * - \f$ NaN - b = NaN \f$,
   * - \f$ a - NaN = NaN \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @param f First element of the subtraction.
   * @param val Second element of the subtraction.
   */
  template <typename U = T>
  friend Point &operator-=(Point &f, const Point<U> &val)
  {
    GUDHI_CHECK(f.size() == val.size(), "Cannot translate point by a direction with different dimension.");
    f._apply_operation(val, [](T &valF, const T &valR) { Gudhi::multi_filtration::_subtract(valF, valR); });
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ (p) \f$ is equal to \f$ f(p) - val \f$.
   *
   * Used conventions:
   * - \f$ inf - inf = NaN \f$,
   * - \f$ -inf - (-inf) = NaN \f$,
   * - \f$ NaN - b = NaN \f$,
   * - \f$ a - NaN = NaN \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @param f First element of the subtraction.
   * @param val Second element of the subtraction.
   */
  friend Point &operator-=(Point &f, const T &val)
  {
    f._apply_operation(val, [](T &valF, const T &valR) { Gudhi::multi_filtration::_subtract(valF, valR); });
    return f;
  }

  // addition
  /**
   * @brief Returns a filtration value such that an entry at index \f$ (p) \f$ is equal to \f$ f(p) + val(p) \f$.
   * Both points have to have the same dimension.
   *
   * Used conventions:
   * - \f$ inf + (-inf) = NaN \f$,
   * - \f$ -inf + inf = NaN \f$,
   * - \f$ NaN + b = NaN \f$,
   * - \f$ a + NaN = NaN \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @param f First element of the addition.
   * @param val Second element of the addition.
   */
  template <typename U = T>
  friend Point operator+(Point f, const Point<U> &val)
  {
    GUDHI_CHECK(f.size() == val.size(), "Cannot translate point by a direction with different dimension.");
    f += val;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (p) \f$ is equal to \f$ f(p) + val \f$.
   *
   * Used conventions:
   * - \f$ inf + (-inf) = NaN \f$,
   * - \f$ -inf + inf = NaN \f$,
   * - \f$ NaN + b = NaN \f$,
   * - \f$ a + NaN = NaN \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @param f First element of the addition.
   * @param val Second element of the addition.
   */
  friend Point operator+(Point f, const T &val)
  {
    f += val;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (p) \f$ is equal to \f$ val + f(p) \f$.
   *
   * Used conventions:
   * - \f$ inf + (-inf) = NaN \f$,
   * - \f$ -inf + inf = NaN \f$,
   * - \f$ NaN + b = NaN \f$,
   * - \f$ a + NaN = NaN \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @param val First element of the addition.
   * @param f Second element of the addition.
   */
  friend Point operator+(const T &val, Point f)
  {
    f += val;
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ (p) \f$ is equal to \f$ f(p) + val(p) \f$.
   * Both points have to have the same dimension.
   *
   * Used conventions:
   * - \f$ inf + (-inf) = NaN \f$,
   * - \f$ -inf + inf = NaN \f$,
   * - \f$ NaN + b = NaN \f$,
   * - \f$ a + NaN = NaN \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @param f First element of the addition.
   * @param val Second element of the addition.
   */
  template <typename U = T>
  friend Point &operator+=(Point &f, const Point<U> &val)
  {
    GUDHI_CHECK(f.size() == val.size(), "Cannot translate point by a direction with different dimension.");
    f._apply_operation(val, [](T &valF, const T &valR) { Gudhi::multi_filtration::_add(valF, valR); });
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ (p) \f$ is equal to \f$ f(p) + val \f$.
   *
   * Used conventions:
   * - \f$ inf + (-inf) = NaN \f$,
   * - \f$ -inf + inf = NaN \f$,
   * - \f$ NaN + b = NaN \f$,
   * - \f$ a + NaN = NaN \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @param f First element of the addition.
   * @param val Second element of the addition.
   */
  friend Point &operator+=(Point &f, const T &val)
  {
    f._apply_operation(val, [](T &valF, const T &valR) { Gudhi::multi_filtration::_add(valF, valR); });
    return f;
  }

  // multiplication
  /**
   * @brief Returns a filtration value such that an entry at index \f$ (p) \f$ is equal to \f$ f(p) * val(p) \f$.
   * Both points have to have the same dimension.
   *
   * Used conventions:
   * - \f$ inf * 0 = NaN \f$,
   * - \f$ 0 * inf = NaN \f$,
   * - \f$ -inf * 0 = NaN \f$,
   * - \f$ 0 * (-inf) = NaN \f$,
   * - \f$ NaN * b = NaN \f$,
   * - \f$ a * NaN = NaN \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @param f First element of the multiplication.
   * @param val Second element of the multiplication.
   */
  template <typename U = T>
  friend Point operator*(Point f, const Point<U> &val)
  {
    GUDHI_CHECK(f.size() == val.size(), "Cannot translate point by a direction with different dimension.");
    f *= val;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (p) \f$ is equal to \f$ f(p) * val \f$.
   *
   * Used conventions:
   * - \f$ inf * 0 = NaN \f$,
   * - \f$ 0 * inf = NaN \f$,
   * - \f$ -inf * 0 = NaN \f$,
   * - \f$ 0 * (-inf) = NaN \f$,
   * - \f$ NaN * b = NaN \f$,
   * - \f$ a * NaN = NaN \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @param f First element of the multiplication.
   * @param val Second element of the multiplication.
   */
  friend Point operator*(Point f, const T &val)
  {
    f *= val;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (p) \f$ is equal to \f$ val * f(p) \f$.
   *
   * Used conventions:
   * - \f$ inf * 0 = NaN \f$,
   * - \f$ 0 * inf = NaN \f$,
   * - \f$ -inf * 0 = NaN \f$,
   * - \f$ 0 * (-inf) = NaN \f$,
   * - \f$ NaN * b = NaN \f$,
   * - \f$ a * NaN = NaN \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @param val First element of the multiplication.
   * @param f Second element of the multiplication.
   */
  friend Point operator*(const T &val, Point f)
  {
    f *= val;
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ (p) \f$ is equal to \f$ f(p) * val(p) \f$.
   * Both points have to have the same dimension.
   *
   * Used conventions:
   * - \f$ inf * 0 = NaN \f$,
   * - \f$ 0 * inf = NaN \f$,
   * - \f$ -inf * 0 = NaN \f$,
   * - \f$ 0 * (-inf) = NaN \f$,
   * - \f$ NaN * b = NaN \f$,
   * - \f$ a * NaN = NaN \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @param f First element of the multiplication.
   * @param val Second element of the multiplication.
   */
  template <typename U = T>
  friend Point &operator*=(Point &f, const Point<U> &val)
  {
    GUDHI_CHECK(f.size() == val.size(), "Cannot translate point by a direction with different dimension.");
    f._apply_operation(val, [](T &valF, const T &valR) { Gudhi::multi_filtration::_multiply(valF, valR); });
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ (p) \f$ is equal to \f$ f(p) * val \f$.
   *
   * Used conventions:
   * - \f$ inf * 0 = NaN \f$,
   * - \f$ 0 * inf = NaN \f$,
   * - \f$ -inf * 0 = NaN \f$,
   * - \f$ 0 * (-inf) = NaN \f$,
   * - \f$ NaN * b = NaN \f$,
   * - \f$ a * NaN = NaN \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @param f First element of the multiplication.
   * @param val Second element of the multiplication.
   */
  friend Point &operator*=(Point &f, const T &val)
  {
    f._apply_operation(val, [](T &valF, const T &valR) { Gudhi::multi_filtration::_multiply(valF, valR); });
    return f;
  }

  // division
  /**
   * @brief Returns a filtration value such that an entry at index \f$ (p) \f$ is equal to \f$ f(p) / val(p) \f$.
   * Both points have to have the same dimension.
   *
   * Used conventions:
   * - \f$ a / 0 = NaN \f$,
   * - \f$ inf / inf = NaN \f$,
   * - \f$ -inf / inf = NaN \f$,
   * - \f$ inf / -inf = NaN \f$,
   * - \f$ -inf / -inf = NaN \f$,
   * - \f$ NaN / b = NaN \f$,
   * - \f$ a / NaN = NaN \f$,
   * - \f$ a / inf = 0 \f$,
   * - \f$ a / -inf = 0 \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @param f First element of the division.
   * @param val Second element of the division.
   */
  template <typename U = T>
  friend Point operator/(Point f, const Point<U> &val)
  {
    GUDHI_CHECK(f.size() == val.size(), "Cannot translate point by a direction with different dimension.");
    f /= val;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (p) \f$ is equal to \f$ f(p) / val \f$.
   *
   * Used conventions:
   * - \f$ a / 0 = NaN \f$,
   * - \f$ inf / inf = NaN \f$,
   * - \f$ -inf / inf = NaN \f$,
   * - \f$ inf / -inf = NaN \f$,
   * - \f$ -inf / -inf = NaN \f$,
   * - \f$ NaN / b = NaN \f$,
   * - \f$ a / NaN = NaN \f$,
   * - \f$ a / inf = 0 \f$,
   * - \f$ a / -inf = 0 \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @param f First element of the division.
   * @param val Second element of the division.
   */
  friend Point operator/(Point f, const T &val)
  {
    f /= val;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (p) \f$ is equal to \f$ val / f(p) \f$.
   *
   * Used conventions:
   * - \f$ a / 0 = NaN \f$,
   * - \f$ inf / inf = NaN \f$,
   * - \f$ -inf / inf = NaN \f$,
   * - \f$ inf / -inf = NaN \f$,
   * - \f$ -inf / -inf = NaN \f$,
   * - \f$ NaN / b = NaN \f$,
   * - \f$ a / NaN = NaN \f$,
   * - \f$ a / inf = 0 \f$,
   * - \f$ a / -inf = 0 \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @param val First element of the division.
   * @param f Second element of the division.
   */
  friend Point operator/(const T &val, Point f)
  {
    f._apply_operation(val, [](T &valF, const T &valR) {
      T tmp = valF;
      valF = valR;
      Gudhi::multi_filtration::_divide(valF, tmp);
    });
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ (p) \f$ is equal to \f$ f(p) / val(p) \f$.
   * Both points have to have the same dimension.
   *
   * Used conventions:
   * - \f$ a / 0 = NaN \f$,
   * - \f$ inf / inf = NaN \f$,
   * - \f$ -inf / inf = NaN \f$,
   * - \f$ inf / -inf = NaN \f$,
   * - \f$ -inf / -inf = NaN \f$,
   * - \f$ NaN / b = NaN \f$,
   * - \f$ a / NaN = NaN \f$,
   * - \f$ a / inf = 0 \f$,
   * - \f$ a / -inf = 0 \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @param f First element of the division.
   * @param val Second element of the division.
   */
  template <typename U = T>
  friend Point &operator/=(Point &f, const Point<U> &val)
  {
    GUDHI_CHECK(f.size() == val.size(), "Cannot translate point by a direction with different dimension.");
    f._apply_operation(val, [](T &valF, const T &valR) { Gudhi::multi_filtration::_divide(valF, valR); });
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ (p) \f$ is equal to \f$ f(p) / val \f$.
   *
   * Used conventions:
   * - \f$ a / 0 = NaN \f$,
   * - \f$ inf / inf = NaN \f$,
   * - \f$ -inf / inf = NaN \f$,
   * - \f$ inf / -inf = NaN \f$,
   * - \f$ -inf / -inf = NaN \f$,
   * - \f$ NaN / b = NaN \f$,
   * - \f$ a / NaN = NaN \f$,
   * - \f$ a / inf = 0 \f$,
   * - \f$ a / -inf = 0 \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @param f First element of the division.
   * @param val Second element of the division.
   */
  friend Point &operator/=(Point &f, const T &val)
  {
    f._apply_operation(val, [](T &valF, const T &valR) { Gudhi::multi_filtration::_divide(valF, valR); });
    return f;
  }

  /**
   * @brief Plus infinity value of an entry of the filtration value.
   */
  constexpr static const T T_inf = Gudhi::multi_filtration::MF_T_inf<T>;

  /**
   * @brief Minus infinity value of an entry of the filtration value.
   */
  constexpr static const T T_m_inf = Gudhi::multi_filtration::MF_T_m_inf<T>;

 private:
  Container coordinates_; /**< Coordinates of the point. */

  template <typename U = T, class F>
  void _apply_operation(const Point<U> &range, F &&operate)
  {
    for (unsigned int p = 0; p < coordinates_.size(); ++p) {
      std::forward<F>(operate)(coordinates_[p], range[p]);
    }
  }

  template <class F>
  void _apply_operation(const T &val, F &&operate)
  {
    for (unsigned int i = 0; i < coordinates_.size(); ++i) {
      std::forward<F>(operate)(coordinates_[i], val);
    }
  }
};

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // MP_POINT_H_INCLUDED

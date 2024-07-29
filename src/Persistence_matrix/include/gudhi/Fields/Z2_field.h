/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022-24 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file Z2_field.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_fields::Z2_field_element class.
 */

#ifndef MATRIX_FIELD_Z2_H_
#define MATRIX_FIELD_Z2_H_

#include <utility>

namespace Gudhi {
namespace persistence_fields {

/**
 * @class Z2_field_element Z2_field.h gudhi/Fields/Z2_field.h
 * @ingroup persistence_fields
 *
 * @brief Class representing an element of the \f$ \mathbb{F}_2 \f$ field.
 */
class Z2_field_element
{
 public:
  using Element = bool;  /**< Type for the elements in the field. */
  template <class T>
  using isInteger = std::enable_if_t<std::is_integral_v<T> >;

  /**
   * @brief Default constructor.
   */
  Z2_field_element();
  /**
   * @brief Constructor setting the element to the given value.
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic.
   * @param element Value of the element.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  Z2_field_element(Integer_type element);
  /**
   * @brief Copy constructor.
   * 
   * @param toCopy Element to copy.
   */
  Z2_field_element(const Z2_field_element& toCopy);
  /**
   * @brief Move constructor.
   * 
   * @param toMove Element to move.
   */
  Z2_field_element(Z2_field_element&& toMove) noexcept;

  /**
   * @brief operator+=
   */
  friend void operator+=(Z2_field_element& f1, Z2_field_element const& f2) {
    f1.element_ = (f1.element_ != f2.element_);
  }
  /**
   * @brief operator+
   */
  friend Z2_field_element operator+(Z2_field_element f1, Z2_field_element const& f2) {
    f1 += f2;
    return f1;
  }
  /**
   * @brief operator+=
   *
   * @tparam Integer_type A native integer type: int, unsigned int, long int, bool, etc.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend void operator+=(Z2_field_element& f, const Integer_type& v) { f.element_ = (f.element_ != _get_value(v)); }
  /**
   * @brief operator+
   *
   * @tparam Integer_type A native integer type: int, unsigned int, long int, bool, etc.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Z2_field_element operator+(Z2_field_element f, const Integer_type& v) {
    f += v;
    return f;
  }
  /**
   * @brief operator+
   *
   * @tparam Integer_type A native integer type: int, unsigned int, long int, bool, etc.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Integer_type operator+(const Integer_type& v, const Z2_field_element& f) {
    return f.element_ != _get_value(v);
  }

  /**
   * @brief operator-=
   */
  friend void operator-=(Z2_field_element& f1, Z2_field_element const& f2) {
    f1.element_ = (f1.element_ != f2.element_);
  }
  /**
   * @brief operator-
   */
  friend Z2_field_element operator-(Z2_field_element f1, Z2_field_element const& f2) {
    f1 -= f2;
    return f1;
  }
  /**
   * @brief operator-=
   *
   * @tparam Integer_type A native integer type: int, unsigned int, long int, bool, etc.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend void operator-=(Z2_field_element& f, const Integer_type& v) { f.element_ = (f.element_ != _get_value(v)); }
  /**
   * @brief operator-
   *
   * @tparam Integer_type A native integer type: int, unsigned int, long int, bool, etc.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Z2_field_element operator-(Z2_field_element f, const Integer_type& v) {
    f -= v;
    return f;
  }
  /**
   * @brief operator-
   *
   * @tparam Integer_type A native integer type: int, unsigned int, long int, bool, etc.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Integer_type operator-(const Integer_type v, Z2_field_element const& f) {
    return f.element_ != _get_value(v);
  }

  /**
   * @brief operator*=
   */
  friend void operator*=(Z2_field_element& f1, Z2_field_element const& f2) {
    f1.element_ = (f1.element_ && f2.element_);
  }
  /**
   * @brief operator*
   */
  friend Z2_field_element operator*(Z2_field_element f1, Z2_field_element const& f2) {
    f1 *= f2;
    return f1;
  }
  /**
   * @brief operator*=
   *
   * @tparam Integer_type A native integer type: int, unsigned int, long int, bool, etc.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend void operator*=(Z2_field_element& f, const Integer_type& v) { f.element_ = (f.element_ && _get_value(v)); }
  /**
   * @brief operator*
   *
   * @tparam Integer_type A native integer type: int, unsigned int, long int, bool, etc.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Z2_field_element operator*(Z2_field_element f, const Integer_type& v) {
    f *= v;
    return f;
  }
  /**
   * @brief operator*
   *
   * @tparam Integer_type A native integer type: int, unsigned int, long int, bool, etc.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Integer_type operator*(const Integer_type& v, Z2_field_element const& f) {
    return f.element_ && _get_value(v);
  }

  /**
   * @brief operator==
   */
  friend bool operator==(const Z2_field_element& f1, const Z2_field_element& f2) { return f1.element_ == f2.element_; }
  /**
   * @brief operator==
   * 
   * @tparam Integer_type A native integer type: int, unsigned int, long int, bool, etc.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator==(const Integer_type& v, const Z2_field_element& f) {
    return _get_value(v) == f.element_;
  }
  /**
   * @brief operator==
   * 
   * @tparam Integer_type A native integer type: int, unsigned int, long int, bool, etc.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator==(const Z2_field_element& f, const Integer_type& v) {
    return _get_value(v) == f.element_;
  }
  /**
   * @brief operator!=
   */
  friend bool operator!=(const Z2_field_element& f1, const Z2_field_element& f2) { return !(f1 == f2); }
  /**
   * @brief operator!=
   * 
   * @tparam Integer_type A native integer type: int, unsigned int, long int, bool, etc.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator!=(const Integer_type v, const Z2_field_element& f) {
    return !(v == f);
  }
  /**
   * @brief operator!=
   * 
   * @tparam Integer_type A native integer type: int, unsigned int, long int, bool, etc.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator!=(const Z2_field_element& f, const Integer_type v) {
    return !(v == f);
  }

  /**
   * @brief Assign operator.
   */
  Z2_field_element& operator=(Z2_field_element other);
  /**
   * @brief Assign operator.
   */
  Z2_field_element& operator=(const unsigned int value);
  /**
   * @brief Swap operator.
   */
  friend void swap(Z2_field_element& f1, Z2_field_element& f2) { std::swap(f1.element_, f2.element_); }

  /**
   * @brief Casts the element into an unsigned int.
   */
  operator unsigned int() const;

  /**
   * @brief Returns the inverse of the element.
   * 
   * @return The inverse.
   */
  Z2_field_element get_inverse() const;
  /**
   * @brief For interface purposes with multi-fields. Returns the inverse together with the second argument.
   * 
   * @param productOfCharacteristics Some value.
   * @return Pair whose first element is the inverse of @p e and the second element is @p productOfCharacteristics.
   */
  std::pair<Z2_field_element, unsigned int> get_partial_inverse(unsigned int productOfCharacteristics) const;

  /**
   * @brief Returns the additive identity of the field.
   * 
   * @return false.
   */
  static Z2_field_element get_additive_identity();
  /**
   * @brief Returns the multiplicative identity of the field.
   * 
   * @return true.
   */
  static Z2_field_element get_multiplicative_identity();
  /**
   * @brief For interface purposes with multi-fields. Returns the multiplicative identity of the field.
   * 
   * @param productOfCharacteristics Some value.
   * @return true.
   */
  static Z2_field_element get_partial_multiplicative_identity([[maybe_unused]] unsigned int productOfCharacteristics);
  /**
   * @brief Returns the characteristic of the field, that is `2`.
   * 
   * @return 2.
   */
  static constexpr unsigned int get_characteristic();

  /**
   * @brief Returns the value of the element.
   * 
   * @return Value of the element.
   */
  Element get_value() const;

  // static constexpr bool handles_only_z2() { return true; }

 private:
  Element element_;

  template <typename Integer_type, class = isInteger<Integer_type> >
  static constexpr Element _get_value(Integer_type e) {
    if constexpr (std::is_same_v<Integer_type, bool>) {
      return e;
    } else {
      return e < 2 && e >= 0 ? e : e % 2;  // returns bool, so %2 won't be negative and is optimized to &
    }
  }
};

inline Z2_field_element::Z2_field_element() : element_(false) {}

template <typename Integer_type, class>
inline Z2_field_element::Z2_field_element(Integer_type element) : element_(_get_value(element)) {}

inline Z2_field_element::Z2_field_element(const Z2_field_element& toCopy) : element_(toCopy.element_) {}

inline Z2_field_element::Z2_field_element(Z2_field_element&& toMove) noexcept
    : element_(std::exchange(toMove.element_, 0)) {}

inline Z2_field_element& Z2_field_element::operator=(Z2_field_element other) {
  std::swap(element_, other.element_);
  return *this;
}

inline Z2_field_element& Z2_field_element::operator=(unsigned int const value) {
  element_ = _get_value(value);
  return *this;
}

inline Z2_field_element::operator unsigned int() const { return element_; }

inline Z2_field_element Z2_field_element::get_inverse() const {
  return element_ ? Z2_field_element(1) : Z2_field_element();
}

inline std::pair<Z2_field_element, unsigned int> Z2_field_element::get_partial_inverse(
    unsigned int productOfCharacteristics) const {
  return {get_inverse(), productOfCharacteristics};
}

inline Z2_field_element Z2_field_element::get_additive_identity() { return Z2_field_element(); }

inline Z2_field_element Z2_field_element::get_multiplicative_identity() { return Z2_field_element(1); }

inline Z2_field_element Z2_field_element::get_partial_multiplicative_identity(
    [[maybe_unused]] unsigned int productOfCharacteristics) {
  return Z2_field_element(1);
}

inline constexpr unsigned int Z2_field_element::get_characteristic() { return 2; }

inline Z2_field_element::Element Z2_field_element::get_value() const { return element_; }

}  // namespace persistence_fields
}  // namespace Gudhi

#endif  // MATRIX_FIELD_Z2_H_

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
 * @file Zp_field.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_fields::Zp_field_element class.
 */

#ifndef MATRIX_FIELD_ZP_H_
#define MATRIX_FIELD_ZP_H_

#include <utility>
#include <array>
#include <limits.h>

namespace Gudhi {
namespace persistence_fields {

/**
 * @class Zp_field_element Zp_field.h gudhi/Fields/Zp_field.h
 * @ingroup persistence_fields
 *
 * @brief Class representing an element of the \f$ \mathbb{F}_p \f$ field for any prime number \f$ p \f$.
 *
 * @tparam characteristic Value of the characteristic of the field. Has to be a positive prime number.
 * @tparam Unsigned_integer_type A native unsigned integer type: unsigned int, long unsigned int, etc.
 * Will be used as the field element type.
 */
template <unsigned int characteristic, typename Unsigned_integer_type = unsigned int,
          class = std::enable_if_t<std::is_unsigned_v<Unsigned_integer_type> > >
class Zp_field_element {
 public:
  using Element = Unsigned_integer_type; /**< Type for the elements in the field. */
  using Characteristic = Element;   /**< Type for the field characteristic. */
  template <class T>
  using isInteger = std::enable_if_t<std::is_integral_v<T> >;

  /**
   * @brief Default constructor. Sets the element to 0.
   */
  Zp_field_element() : element_(0) { static_assert(_is_prime(), "Characteristic has to be a prime number."); }
  /**
   * @brief Constructor setting the element to the given value.
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   * @param element Value of the element.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  Zp_field_element(Integer_type element) : element_(_get_value(element)) {
    static_assert(_is_prime(), "Characteristic has to be a prime number.");
  }
  /**
   * @brief Copy constructor.
   * 
   * @param toCopy Element to copy.
   */
  Zp_field_element(const Zp_field_element& toCopy) : element_(toCopy.element_) {}
  /**
   * @brief Move constructor.
   * 
   * @param toMove Element to move.
   */
  Zp_field_element(Zp_field_element&& toMove) noexcept : element_(std::exchange(toMove.element_, 0)) {}

  /**
   * @brief operator+=
   */
  friend void operator+=(Zp_field_element& f1, const Zp_field_element& f2) {
    f1.element_ = Zp_field_element::_add(f1.element_, f2.element_);
  }
  /**
   * @brief operator+
   */
  friend Zp_field_element operator+(Zp_field_element f1, const Zp_field_element& f2) {
    f1 += f2;
    return f1;
  }
  /**
   * @brief operator+=
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend void operator+=(Zp_field_element& f, const Integer_type& v) {
    f.element_ = Zp_field_element::_add(f.element_, _get_value(v));
  }
  /**
   * @brief operator+
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Zp_field_element operator+(Zp_field_element f, const Integer_type& v) {
    f += v;
    return f;
  }
  /**
   * @brief operator+
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Integer_type operator+(const Integer_type& v, Zp_field_element f) {
    f += v;
    return f.element_;
  }

  /**
   * @brief operator-=
   */
  friend void operator-=(Zp_field_element& f1, const Zp_field_element& f2) {
    f1.element_ = Zp_field_element::_subtract(f1.element_, f2.element_);
  }
  /**
   * @brief operator-
   */
  friend Zp_field_element operator-(Zp_field_element f1, const Zp_field_element& f2) {
    f1 -= f2;
    return f1;
  }
  /**
   * @brief operator-=
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend void operator-=(Zp_field_element& f, const Integer_type& v) {
    f.element_ = Zp_field_element::_subtract(f.element_, _get_value(v));
  }
  /**
   * @brief operator-
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Zp_field_element operator-(Zp_field_element f, const Integer_type& v) {
    f -= v;
    return f;
  }
  /**
   * @brief operator-
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Integer_type operator-(const Integer_type& v, const Zp_field_element& f) {
    return Zp_field_element::_subtract(_get_value(v), f.element_);
  }

  /**
   * @brief operator*=
   */
  friend void operator*=(Zp_field_element& f1, const Zp_field_element& f2) {
    f1.element_ = Zp_field_element::_multiply(f1.element_, f2.element_);
  }
  /**
   * @brief operator*
   */
  friend Zp_field_element operator*(Zp_field_element f1, const Zp_field_element& f2) {
    f1 *= f2;
    return f1;
  }
  /**
   * @brief operator*=
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend void operator*=(Zp_field_element& f, const Integer_type& v) {
    f.element_ = Zp_field_element::_multiply(f.element_, _get_value(v));
  }
  /**
   * @brief operator*
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Zp_field_element operator*(Zp_field_element f, const Integer_type& v) {
    f *= v;
    return f;
  }
  /**
   * @brief operator*
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Integer_type operator*(const Integer_type& v, Zp_field_element f) {
    f *= v;
    return f.element_;
  }

  /**
   * @brief operator==
   */
  friend bool operator==(const Zp_field_element& f1, const Zp_field_element& f2) {
    return f1.element_ == f2.element_;
  }
  /**
   * @brief operator==
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator==(const Integer_type& v, const Zp_field_element& f) {
    return Zp_field_element::_get_value(v) == f.element_;
  }
  /**
   * @brief operator==
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator==(const Zp_field_element& f, const Integer_type& v) {
    return Zp_field_element::_get_value(v) == f.element_;
  }
  /**
   * @brief operator!=
   */
  friend bool operator!=(const Zp_field_element& f1, const Zp_field_element& f2) { return !(f1 == f2); }
  /**
   * @brief operator!=
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator!=(const Integer_type& v, const Zp_field_element& f) {
    return !(v == f);
  }
  /**
   * @brief operator!=
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator!=(const Zp_field_element& f, const Integer_type& v) {
    return !(v == f);
  }

  /**
   * @brief Assign operator.
   */
  Zp_field_element& operator=(Zp_field_element other) {
    std::swap(element_, other.element_);
    return *this;
  }
  /**
   * @brief Assign operator.
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  Zp_field_element& operator=(const Integer_type& value) {
    element_ = _get_value(value);
    return *this;
  }
  /**
   * @brief Swap operator.
   */
  friend void swap(Zp_field_element& f1, Zp_field_element& f2) { std::swap(f1.element_, f2.element_); }

  /**
   * @brief Casts the element into an unsigned int.
   */
  operator unsigned int() const { return element_; }

  /**
   * @brief Returns the inverse of the element in the field.
   * 
   * @return The inverse.
   */
  Zp_field_element get_inverse() const {
    if (element_ != 0 && inverse_[element_] == 0) {  // initialize everything at instantiation instead?
      inverse_[element_] = _get_inverse(element_);
    }

    return Zp_field_element<characteristic>(inverse_[element_]);
  }
  /**
   * @brief For interface purposes with multi-fields. Returns the inverse together with the argument.
   * 
   * @param productOfCharacteristics Some value.
   * @return Pair whose first element is the inverse and the second element is @p productOfCharacteristics.
   */
  std::pair<Zp_field_element, unsigned int> get_partial_inverse(unsigned int productOfCharacteristics) const {
    return {get_inverse(), productOfCharacteristics};
  }

  /**
   * @brief Returns the additive identity of the field.
   * 
   * @return 0.
   */
  static Zp_field_element get_additive_identity() { return Zp_field_element<characteristic>(); }
  /**
   * @brief Returns the multiplicative identity of the field.
   * 
   * @return 1.
   */
  static Zp_field_element get_multiplicative_identity() { return Zp_field_element<characteristic>(1); }
  /**
   * @brief For interface purposes with multi-fields. Returns the multiplicative identity of the field.
   * 
   * @param productOfCharacteristics Some value.
   * @return 1.
   */
  static Zp_field_element get_partial_multiplicative_identity([[maybe_unused]] unsigned int productOfCharacteristics) {
    return Zp_field_element<characteristic>(1);
  }
  /**
   * @brief Returns the current characteristic.
   * 
   * @return The value of the current characteristic.
   */
  static constexpr unsigned int get_characteristic() { return characteristic; }

  /**
   * @brief Returns the value of the element.
   * 
   * @return Value of the element.
   */
  Element get_value() const { return element_; }

  // static constexpr bool handles_only_z2() { return false; }

 private:
  Element element_;                                            /**< Field element. */
  static inline std::array<unsigned int, characteristic> inverse_;  /**< All inverse elements. */

  static Element _add(Element element, Element v) {
    if (UINT_MAX - element < v) {
      // automatic unsigned integer overflow behaviour will make it work
      element += v;
      element -= characteristic;
      return element;
    }

    element += v;
    if (element >= characteristic) element -= characteristic;

    return element;
  }
  static Element _subtract(Element element, Element v) {
    if (element < v) {
      element += characteristic;
    }
    element -= v;

    return element;
  }
  static Element _multiply(Element element, Element v) {
    Element a = element;
    element = 0;
    Element temp_b;

    while (a != 0) {
      if (a & 1) {
        if (v >= characteristic - element) element -= characteristic;
        element += v;
      }
      a >>= 1;

      temp_b = v;
      if (v >= characteristic - v) temp_b -= characteristic;
      v += temp_b;
    }

    return element;
  }
  static int _get_inverse(Element element) {
    // to solve: Ax + My = 1
    int M = characteristic;
    int A = element;
    int y = 0, x = 1;
    // extended euclidean division
    while (A > 1) {
      int quotient = A / M;
      int temp = M;

      M = A % M, A = temp;
      temp = y;

      y = x - quotient * y;
      x = temp;
    }

    if (x < 0) x += characteristic;

    return x;
  }

  template <typename Integer_type, class = isInteger<Integer_type> >
  static constexpr Element _get_value(Integer_type e) {
    if constexpr (std::is_signed_v<Integer_type>) {
      if (e < -static_cast<Integer_type>(characteristic)) e = e % characteristic;
      if (e < 0) return e += characteristic;
      return e < static_cast<Integer_type>(characteristic) ? e : e % characteristic;
    } else {
      return e < characteristic ? e : e % characteristic;
    }
  }

  static constexpr bool _is_prime() {
    if (characteristic <= 1) return false;
    if (characteristic <= 3) return true;
    if (characteristic % 2 == 0 || characteristic % 3 == 0) return false;

    for (long i = 5; i * i <= characteristic; i = i + 6)
      if (characteristic % i == 0 || characteristic % (i + 2) == 0) return false;

    return true;
  }
};

}  // namespace persistence_fields
}  // namespace Gudhi

#endif  // MATRIX_FIELD_ZP_H_

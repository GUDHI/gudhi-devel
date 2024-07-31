/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file Zp_field_shared.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_fields::Shared_Zp_field_element class.
 */

#ifndef MATRIX_FIELD_ZP_VAR_H_
#define MATRIX_FIELD_ZP_VAR_H_

#include <utility>
#include <vector>
#include <limits.h>
#include <stdexcept>

namespace Gudhi {
namespace persistence_fields {

/**
 * @class Shared_Zp_field_element Zp_field_shared.h gudhi/Fields/Zp_field_shared.h
 * @ingroup persistence_fields
 *
 * @brief Class representing an element of the \f$ \mathbb{F}_p \f$ field for any prime number \f$ p \f$.
 * If each instantiation of the class can represent another element, they all share the same characteristics.
 * That is if the characteristics are set for one, they will be set for all the others. The characteristics can
 * be set before instantiating the elements with the static @ref Shared_Zp_field_element::initialize method.
 *
 * @tparam Unsigned_integer_type A native unsigned integer type: unsigned int, long unsigned int, etc.
 * Will be used as the field element type.
 */
template <typename Unsigned_integer_type = unsigned int,
          class = std::enable_if_t<std::is_unsigned_v<Unsigned_integer_type> > >
class Shared_Zp_field_element {
 public:
  using Element = Unsigned_integer_type; /**< Type for the elements in the field. */
  using Characteristic = Element;   /**< Type for the field characteristic. */
  template <class T>
  using isInteger = std::enable_if_t<std::is_integral_v<T> >;

  /**
   * @brief Default constructor. Sets the element to 0.
   */
  Shared_Zp_field_element() : element_(0) {}
  /**
   * @brief Constructor setting the element to the given value.
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   * @param element Value of the element.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  Shared_Zp_field_element(Integer_type element) : element_(_get_value(element)) {}
  /**
   * @brief Copy constructor.
   * 
   * @param toCopy Element to copy.
   */
  Shared_Zp_field_element(const Shared_Zp_field_element& toCopy) : element_(toCopy.element_) {}
  /**
   * @brief Move constructor.
   * 
   * @param toMove Element to move.
   */
  Shared_Zp_field_element(Shared_Zp_field_element&& toMove) noexcept : element_(std::exchange(toMove.element_, 0)) {}

  /**
   * @brief Initialize the field to the given characteristic.
   * Should be called first before constructing the field elements.
   * 
   * @param characteristic Characteristic of the field. A positive prime number.
   */
  static void initialize(Characteristic characteristic) {
    if (characteristic <= 1)
      throw std::invalid_argument("Characteristic must be strictly positive and a prime number.");

    inverse_.resize(characteristic);
    inverse_[0] = 0;
    for (Element i = 1; i < characteristic; ++i) {
      Element inv = 1;
      Element mult = inv * i;
      while ((mult % characteristic) != 1) {
        ++inv;
        if (mult == characteristic) throw std::invalid_argument("Characteristic must be a prime number.");
        mult = inv * i;
      }
      inverse_[i] = inv;
    }

    characteristic_ = characteristic;
  }

  /**
   * @brief Returns the value of the element.
   * 
   * @return Value of the element.
   */
  Element get_value() const { return element_; }

  /**
   * @brief operator+=
   */
  friend void operator+=(Shared_Zp_field_element& f1, const Shared_Zp_field_element& f2) {
    f1.element_ = Shared_Zp_field_element::_add(f1.element_, f2.element_);
  }
  /**
   * @brief operator+
   */
  friend Shared_Zp_field_element operator+(Shared_Zp_field_element f1, const Shared_Zp_field_element& f2) {
    f1 += f2;
    return f1;
  }
  /**
   * @brief operator+=
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend void operator+=(Shared_Zp_field_element& f, const Integer_type& v) {
    f.element_ = Shared_Zp_field_element::_add(f.element_, _get_value(v));
  }
  /**
   * @brief operator+
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Shared_Zp_field_element operator+(Shared_Zp_field_element f, const Integer_type& v) {
    f += v;
    return f;
  }
  /**
   * @brief operator+
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Integer_type operator+(const Integer_type& v, Shared_Zp_field_element f) {
    f += v;
    return f.element_;
  }

  /**
   * @brief operator-=
   */
  friend void operator-=(Shared_Zp_field_element& f1, const Shared_Zp_field_element& f2) {
    f1.element_ = Shared_Zp_field_element::_subtract(f1.element_, f2.element_);
  }
  /**
   * @brief operator-
   */
  friend Shared_Zp_field_element operator-(Shared_Zp_field_element f1, const Shared_Zp_field_element& f2) {
    f1 -= f2;
    return f1;
  }
  /**
   * @brief operator-=
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend void operator-=(Shared_Zp_field_element& f, const Integer_type& v) {
    f.element_ = Shared_Zp_field_element::_subtract(f.element_, _get_value(v));
  }
  /**
   * @brief operator-
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Shared_Zp_field_element operator-(Shared_Zp_field_element f, const Integer_type& v) {
    f -= v;
    return f;
  }
  /**
   * @brief operator-
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Integer_type operator-(const Integer_type& v, const Shared_Zp_field_element& f) {
    return Shared_Zp_field_element::_subtract(_get_value(v), f.element_);
  }

  /**
   * @brief operator*=
   */
  friend void operator*=(Shared_Zp_field_element& f1, const Shared_Zp_field_element& f2) {
    f1.element_ = Shared_Zp_field_element::_multiply(f1.element_, f2.element_);
  }
  /**
   * @brief operator*
   */
  friend Shared_Zp_field_element operator*(Shared_Zp_field_element f1, const Shared_Zp_field_element& f2) {
    f1 *= f2;
    return f1;
  }
  /**
   * @brief operator*=
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend void operator*=(Shared_Zp_field_element& f, const Integer_type& v) {
    f.element_ = Shared_Zp_field_element::_multiply(f.element_, _get_value(v));
  }
  /**
   * @brief operator*
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Shared_Zp_field_element operator*(Shared_Zp_field_element f, const Integer_type& v) {
    f *= v;
    return f;
  }
  /**
   * @brief operator*
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic. if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Integer_type operator*(const Integer_type& v, Shared_Zp_field_element f) {
    f *= v;
    return f.element_;
  }

  /**
   * @brief operator==
   */
  friend bool operator==(const Shared_Zp_field_element& f1, const Shared_Zp_field_element& f2) {
    return f1.element_ == f2.element_;
  }
  /**
   * @brief operator==
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator==(const Integer_type& v, const Shared_Zp_field_element& f) {
    return Shared_Zp_field_element::_get_value(v) == f.element_;
  }
  /**
   * @brief operator==
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator==(const Shared_Zp_field_element& f, const Integer_type& v) {
    return Shared_Zp_field_element::_get_value(v) == f.element_;
  }
  /**
   * @brief operator!=
   */
  friend bool operator!=(const Shared_Zp_field_element& f1, const Shared_Zp_field_element& f2) { return !(f1 == f2); }
  /**
   * @brief operator!=
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator!=(const Integer_type& v, const Shared_Zp_field_element& f) {
    return !(v == f);
  }
  /**
   * @brief operator!=
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator!=(const Shared_Zp_field_element& f, const Integer_type& v) {
    return !(v == f);
  }

  /**
   * @brief Assign operator.
   */
  Shared_Zp_field_element& operator=(Shared_Zp_field_element other) {
    std::swap(element_, other.element_);
    return *this;
  }
  /**
   * @brief Assign operator.
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  Shared_Zp_field_element& operator=(const Integer_type& value) {
    element_ = Shared_Zp_field_element::_get_value(value);
    return *this;
  }
  /**
   * @brief Swap operator.
   */
  friend void swap(Shared_Zp_field_element& f1, Shared_Zp_field_element& f2) { std::swap(f1.element_, f2.element_); }

  /**
   * @brief Casts the element into an unsigned int.
   */
  operator unsigned int() const { return element_; }

  /**
   * @brief Returns the inverse of the element in the field.
   * 
   * @return The inverse.
   */
  Shared_Zp_field_element get_inverse() const { return Shared_Zp_field_element(inverse_[element_]); }
  /**
   * @brief For interface purposes with multi-fields. Returns the inverse together with the argument.
   * 
   * @param productOfCharacteristics Some value.
   * @return Pair whose first element is the inverse and the second element is @p productOfCharacteristics.
   */
  std::pair<Shared_Zp_field_element, Characteristic> get_partial_inverse(
      Characteristic productOfCharacteristics) const {
    return {get_inverse(), productOfCharacteristics};
  }

  /**
   * @brief Returns the additive identity of the field.
   * 
   * @return 0.
   */
  static Shared_Zp_field_element get_additive_identity() { return Shared_Zp_field_element(); }
  /**
   * @brief Returns the multiplicative identity of the field.
   * 
   * @return 1.
   */
  static Shared_Zp_field_element get_multiplicative_identity() { return Shared_Zp_field_element(1); }
  /**
   * @brief For interface purposes with multi-fields. Returns the multiplicative identity of the field.
   * 
   * @param productOfCharacteristics Some value.
   * @return 1.
   */
  static Shared_Zp_field_element get_partial_multiplicative_identity(
      [[maybe_unused]] Characteristic productOfCharacteristics) {
    return Shared_Zp_field_element(1);
  }
  /**
   * @brief Returns the current characteristic.
   * 
   * @return The value of the current characteristic.
   */
  static Characteristic get_characteristic() { return characteristic_; }

  // static constexpr bool handles_only_z2() { return false; }

 private:
  Element element_;                              /**< Field element. */
  static inline Characteristic characteristic_;  /**< Current characteristic of the field. */
  static inline std::vector<Element> inverse_;   /**< All inverse elements. */

  static Element _add(Element element, Element v) {
    if (UINT_MAX - element < v) {
      // automatic unsigned integer overflow behaviour will make it work
      element += v;
      element -= characteristic_;
      return element;
    }

    element += v;
    if (element >= characteristic_) element -= characteristic_;

    return element;
  }
  static Element _subtract(Element element, Element v) {
    if (element < v) {
      element += characteristic_;
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
        if (v >= characteristic_ - element) element -= characteristic_;
        element += v;
      }
      a >>= 1;

      temp_b = v;
      if (v >= characteristic_ - v) temp_b -= characteristic_;
      v += temp_b;
    }

    return element;
  }

  template <typename Integer_type, class = isInteger<Integer_type> >
  static constexpr Element _get_value(Integer_type e) {
    if constexpr (std::is_signed_v<Integer_type>){
      if (e < -static_cast<Integer_type>(characteristic_)) e = e % characteristic_;
      if (e < 0) return e += characteristic_;
      return e < static_cast<Integer_type>(characteristic_) ? e : e % characteristic_;
    } else {
      return e < characteristic_ ? e : e % characteristic_;
    }
  }
};

}  // namespace persistence_fields
}  // namespace Gudhi

#endif  // MATRIX_FIELD_ZP_VAR_H_

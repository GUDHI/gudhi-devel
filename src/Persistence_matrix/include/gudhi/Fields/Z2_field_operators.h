/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file Z2_field_operators.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_fields::Z2_field_operators class.
 */

#ifndef MATRIX_FIELD_Z2_OPERATORS_H_
#define MATRIX_FIELD_Z2_OPERATORS_H_

#include <utility>

namespace Gudhi {
namespace persistence_fields {

/**
 * @class Z2_field_operators Z2_field_operators.h gudhi/Fields/Z2_field_operators.h
 * @ingroup persistence_fields
 *
 * @brief Class defining operators for the \f$ \mathbb{F}_2 \f$ field.
 */
class Z2_field_operators
{
 public:
  using Element = bool;                /**< Type for the elements in the field. */
  using Characteristic = unsigned int; /**< Type for the field characteristic. */
  template <class T>
  using isUnsignedInteger = std::enable_if_t<std::is_unsigned_v<T> >;
  template <class T>
  using isInteger = std::enable_if_t<std::is_integral_v<T> >;

  /**
   * @brief Default constructor.
   */
  Z2_field_operators(){};

  /**
   * @brief Returns the characteristic of the field, that is `2`.
   * 
   * @return 2.
   */
  static constexpr Characteristic get_characteristic() { return 2; }

  /**
   * @brief Returns the value of an integer in the field.
   * That is the positive value of the integer modulo the current characteristic.
   * 
   * @tparam Integer_type A native integer type: int, unsigned int, long int, bool, etc.
   * @param e Integer to return the value from.
   * @return A boolean representing `e % 2`.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  static Element get_value(Integer_type e) {
    if constexpr (std::is_same_v<Integer_type, bool>) {
      return e;
    } else {
      return e < 2 && e >= 0 ? e : e % 2;  // returns bool, so %2 won't be negative and is optimized to &
    }
  }

  /**
   * @brief Returns the sum of two elements in the field.
   * 
   * @tparam Unsigned_integer_type A native unsigned integer type: unsigned int, bool, etc.
   * @param e1 First element.
   * @param e2 Second element.
   * @return `(e1 + e2) % 2` as a boolean.
   */
  template <typename Unsigned_integer_type, class = isUnsignedInteger<Unsigned_integer_type> >
  static Element add(Unsigned_integer_type e1, Unsigned_integer_type e2) {
    if constexpr (std::is_same_v<Unsigned_integer_type, bool>) {
      return e1 != e2;
    } else {
      return get_value(e1) != get_value(e2);
    }
  }

  /**
   * @brief Stores in the first element the sum of two given elements in the field, that is
   * `(e1 + e2) % 2`, such that the result is positive.
   * 
   * @tparam Unsigned_integer_type A native unsigned integer type: unsigned int, bool, etc.
   * @param e1 First element.
   * @param e2 Second element.
   */
  template <typename Unsigned_integer_type, class = isUnsignedInteger<Unsigned_integer_type> >
  static void add_inplace(Unsigned_integer_type& e1, Unsigned_integer_type e2) {
    if constexpr (std::is_same_v<Unsigned_integer_type, bool>) {
      e1 = e1 != e2;
    } else {
      e1 = get_value(e1) != get_value(e2);
    }
  }

  /**
   * @brief Returns the subtraction in the field of the first element by the second element.
   * 
   * @tparam Unsigned_integer_type A native unsigned integer type: unsigned int, bool, etc.
   * @param e1 First element.
   * @param e2 Second element.
   * @return `(e1 - e2) % 2` as a boolean.
   */
  template <typename Unsigned_integer_type, class = isUnsignedInteger<Unsigned_integer_type> >
  static Element subtract(Unsigned_integer_type e1, Unsigned_integer_type e2) {
    if constexpr (std::is_same_v<Unsigned_integer_type, bool>) {
      return e1 != e2;
    } else {
      return get_value(e1) != get_value(e2);
    }
  }

  /**
   * @brief Stores in the first element the subtraction in the field of the first element by the second element,
   * that is `(e1 - e2) % 2`, such that the result is positive.
   * 
   * @tparam Unsigned_integer_type A native unsigned integer type: unsigned int, bool, etc.
   * @param e1 First element.
   * @param e2 Second element.
   */
  template <typename Unsigned_integer_type, class = isUnsignedInteger<Unsigned_integer_type> >
  static void subtract_inplace_front(Unsigned_integer_type& e1, Unsigned_integer_type e2) {
    if constexpr (std::is_same_v<Unsigned_integer_type, bool>) {
      e1 = e1 != e2;
    } else {
      e1 = get_value(e1) != get_value(e2);
    }
  }
  /**
   * @brief Stores in the second element the subtraction in the field of the first element by the second element,
   * that is `(e1 - e2) % 2`, such that the result is positive.
   * 
   * @tparam Unsigned_integer_type A native unsigned integer type: unsigned int, bool, etc.
   * @param e1 First element.
   * @param e2 Second element.
   */
  template <typename Unsigned_integer_type, class = isUnsignedInteger<Unsigned_integer_type> >
  static void subtract_inplace_back(Unsigned_integer_type e1, Unsigned_integer_type& e2) {
    if constexpr (std::is_same_v<Unsigned_integer_type, bool>) {
      e2 = e1 != e2;
    } else {
      e2 = get_value(e1) != get_value(e2);
    }
  }

  /**
   * @brief Returns the multiplication of two elements in the field.
   * 
   * @tparam Unsigned_integer_type A native unsigned integer type: unsigned int, bool, etc.
   * @param e1 First element.
   * @param e2 Second element.
   * @return `(e1 * e2) % 2` as a boolean.
   */
  template <typename Unsigned_integer_type, class = isUnsignedInteger<Unsigned_integer_type> >
  static Element multiply(Unsigned_integer_type e1, Unsigned_integer_type e2) {
    if constexpr (std::is_same_v<Unsigned_integer_type, bool>) {
      return e1 && e2;
    } else {
      return get_value(e1) ? get_value(e2) : false;
    }
  }

  /**
   * @brief Stores in the first element the multiplication of two given elements in the field,
   * that is `(e1 * e2) % 2`, such that the result is positive.
   * 
   * @tparam Unsigned_integer_type A native unsigned integer type: unsigned int, bool, etc.
   * @param e1 First element.
   * @param e2 Second element.
   */
  template <typename Unsigned_integer_type, class = isUnsignedInteger<Unsigned_integer_type> >
  static void multiply_inplace(Unsigned_integer_type& e1, Unsigned_integer_type e2) {
    if constexpr (std::is_same_v<Unsigned_integer_type, bool>) {
      e1 = e1 && e2;
    } else {
      e1 = get_value(e1) ? get_value(e2) : false;
    }
  }

  /**
   * @brief Multiplies the first element with the second one and adds the third one. Returns the result in the field.
   * 
   * @tparam Unsigned_integer_type A native unsigned integer type: unsigned int, bool, etc.
   * @param e First element.
   * @param m Second element.
   * @param a Third element.
   * @return `(e * m + a) % 2` as a boolean.
   */
  template <typename Unsigned_integer_type, class = isUnsignedInteger<Unsigned_integer_type> >
  static Element multiply_and_add(Unsigned_integer_type e, Unsigned_integer_type m, Unsigned_integer_type a) {
    if constexpr (std::is_same_v<Unsigned_integer_type, bool>) {
      return (e && m) != a;
    } else {
      return multiply(e, m) != get_value(a);
    }
  }

  /**
   * @brief Multiplies the first element with the second one and adds the third one, that is
   * `(e * m + a) % 2`, such that the result is positive. Stores the result in the first element.
   * 
   * @tparam Unsigned_integer_type A native unsigned integer type: unsigned int, bool, etc.
   * @param e First element.
   * @param m Second element.
   * @param a Third element.
   */
  template <typename Unsigned_integer_type, class = isUnsignedInteger<Unsigned_integer_type> >
  static void multiply_and_add_inplace_front(Unsigned_integer_type& e, 
                                             Unsigned_integer_type m,
                                             Unsigned_integer_type a) {
    if constexpr (std::is_same_v<Unsigned_integer_type, bool>) {
      e = (e && m) != a;
    } else {
      e = multiply(e, m) != get_value(a);
    }
  }

  /**
   * @brief Multiplies the first element with the second one and adds the third one, that is
   * `(e * m + a) % 2`, such that the result is positive. Stores the result in the third element.
   * 
   * @tparam Unsigned_integer_type A native unsigned integer type: unsigned int, bool, etc.
   * @param e First element.
   * @param m Second element.
   * @param a Third element.
   */
  template <typename Unsigned_integer_type, class = isUnsignedInteger<Unsigned_integer_type> >
  static void multiply_and_add_inplace_back(Unsigned_integer_type e, 
                                            Unsigned_integer_type m,
                                            Unsigned_integer_type& a) {
    if constexpr (std::is_same_v<Unsigned_integer_type, bool>) {
      a = (e && m) != a;
    } else {
      a = multiply(e, m) != get_value(a);
    }
  }

  /**
   * @brief Adds the first element to the second one and multiplies the third one with it.
   * Returns the result in the field.
   * 
   * @tparam Unsigned_integer_type A native unsigned integer type: unsigned int, bool, etc.
   * @param e First element.
   * @param a Second element.
   * @param m Third element.
   * @return `((e + a) * m) % 2` as a boolean.
   */
  template <typename Unsigned_integer_type, class = isUnsignedInteger<Unsigned_integer_type> >
  static Element add_and_multiply(Unsigned_integer_type e, Unsigned_integer_type a, Unsigned_integer_type m) {
    if constexpr (std::is_same_v<Unsigned_integer_type, bool>) {
      return (e != a) && m;
    } else {
      return add(e, a) ? get_value(m) : false;
    }
  }

  /**
   * @brief Adds the first element to the second one and multiplies the third one with it, that is
   * `((e + a) * m) % 2`, such that the result is positive. Stores the result in the first element.
   * 
   * @tparam Unsigned_integer_type A native unsigned integer type: unsigned int, bool, etc.
   * @param e First element.
   * @param a Second element.
   * @param m Third element.
   */
  template <typename Unsigned_integer_type, class = isUnsignedInteger<Unsigned_integer_type> >
  static void add_and_multiply_inplace_front(Unsigned_integer_type& e, Unsigned_integer_type a, Unsigned_integer_type m) {
    if constexpr (std::is_same_v<Unsigned_integer_type, bool>) {
      e = (e != a) && m;
    } else {
      e = add(e, a) ? get_value(m) : false;
    }
  }
  /**
   * @brief Adds the first element to the second one and multiplies the third one with it, that is
   * `((e + a) * m) % 2`, such that the result is positive. Stores the result in the third element.
   * 
   * @tparam Unsigned_integer_type A native unsigned integer type: unsigned int, bool, etc.
   * @param e First element.
   * @param a Second element.
   * @param m Third element.
   */
  template <typename Unsigned_integer_type, class = isUnsignedInteger<Unsigned_integer_type> >
  static void add_and_multiply_inplace_back(Unsigned_integer_type& e, Unsigned_integer_type a, Unsigned_integer_type m) {
    if constexpr (std::is_same_v<Unsigned_integer_type, bool>) {
      m = (e != a) && m;
    } else {
      m = add(e, a) ? get_value(m) : false;
    }
  }

  /**
   * @brief Returns true if the two given elements are equal in the field, false otherwise.
   * 
   * @tparam Unsigned_integer_type A native unsigned integer type: unsigned int, bool, etc.
   * @param e1 First element to compare.
   * @param e2 Second element to compare.
   * @return true If `e1 % 2 == e2 % 2`.
   * @return false Otherwise.
   */
  template <typename Unsigned_integer_type, class = isUnsignedInteger<Unsigned_integer_type> >
  static bool are_equal(Unsigned_integer_type e1, Unsigned_integer_type e2) {
    if constexpr (std::is_same_v<Unsigned_integer_type, bool>) {
      return e1 == e2;
    } else {
      return get_value(e1) == get_value(e2);
    }
  }

  /**
   * @brief Returns the inverse of the given element in the field.
   * 
   * @tparam Unsigned_integer_type A native unsigned integer type: unsigned int, bool, etc.
   * @param e Element to get the inverse from.
   * @return Inverse in the current field of `e % 2`.
   */
  template <typename Unsigned_integer_type, class = isUnsignedInteger<Unsigned_integer_type> >
  static Element get_inverse(Unsigned_integer_type e) {
    if constexpr (std::is_same_v<Unsigned_integer_type, bool>) {
      return e;
    } else {
      return get_value(e);
    }
  }
  /**
   * @brief For interface purposes with multi-fields. Returns the inverse together with the second argument.
   * 
   * @tparam Unsigned_integer_type  A native unsigned integer type: unsigned int, bool, etc.
   * @param e Element to get the inverse from.
   * @param productOfCharacteristics Some value.
   * @return Pair whose first element is the inverse of @p e and the second element is @p productOfCharacteristics.
   */
  template <typename Unsigned_integer_type, class = isUnsignedInteger<Unsigned_integer_type> >
  static std::pair<Element, Characteristic> get_partial_inverse(
      Unsigned_integer_type e, Characteristic productOfCharacteristics) {
    return {get_inverse(e), productOfCharacteristics};
  }

  /**
   * @brief Returns the additive identity of the field.
   * 
   * @return false.
   */
  static constexpr Element get_additive_identity() { return false; }
  /**
   * @brief Returns the multiplicative identity of the field.
   * 
   * @return true.
   */
  static constexpr Element get_multiplicative_identity() { return true; }
  /**
   * @brief For interface purposes with multi-fields. Returns the multiplicative identity of the field.
   * 
   * @param productOfCharacteristics Some value.
   * @return true.
   */
  static constexpr Element get_partial_multiplicative_identity(
      [[maybe_unused]] Characteristic productOfCharacteristics) {
    return true;
  }

  // static constexpr bool handles_only_z2() { return true; }
};

}  // namespace persistence_fields
}  // namespace Gudhi

#endif  // MATRIX_FIELD_Z2_OPERATORS_H_

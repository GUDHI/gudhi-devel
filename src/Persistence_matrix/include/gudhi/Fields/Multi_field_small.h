/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber, Clément Maria
 *
 *    Copyright (C) 2022-24 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file Multi_field_small.h
 * @author Hannah Schreiber, Clément Maria
 * @brief Contains the @ref Gudhi::persistence_fields::Multi_field_element_with_small_characteristics class.
 */

#ifndef MATRIX_FIELD_MULTI_SMALL_H_
#define MATRIX_FIELD_MULTI_SMALL_H_

#include <utility>
#include <vector>
#include <limits.h>
#include <numeric>

namespace Gudhi {
namespace persistence_fields {

/**
 * @class Multi_field_element_with_small_characteristics Multi_field_small.h gudhi/Fields/Multi_field_small.h
 * @ingroup persistence_fields
 *
 * @brief Class representing an element of a multi-field, such that the product of all characteristics fits into
 * the given @p Unsigned_integer_type template argument. The characteristics will corresponds to all prime numbers
 * in the interval given as other template arguments.
 *
 * @tparam minimum Interval closed lower bound.
 * @tparam maximum Interval closed upper bound.
 * @tparam Unsigned_integer_type A native unsigned integer type: unsigned int, long unsigned int, etc.
 * Will be used as the field element type.
 */
template <unsigned int minimum, unsigned int maximum, typename Unsigned_integer_type = unsigned int,
          class = std::enable_if_t<std::is_unsigned_v<Unsigned_integer_type> > >
class Multi_field_element_with_small_characteristics {
 public:
  using Element = Unsigned_integer_type; /**< Type for the elements in the field. */
  using Characteristic = Element;   /**< Type for the field characteristic. */
  template <class T>
  using isInteger = std::enable_if_t<std::is_integral_v<T> >;

  /**
   * @brief Default constructor. Sets the element to 0.
   */
  Multi_field_element_with_small_characteristics() : element_(0) {
    static_assert(maximum >= 2, "Characteristics have to be positive.");
    static_assert(minimum <= maximum, "The given interval is not valid.");
    static_assert(minimum != maximum || _is_prime(minimum), "The given interval does not contain a prime number.");
    static_assert(productOfAllCharacteristics_ != 1, "The given interval does not contain a prime number.");
  }
  /**
   * @brief Constructor setting the element to the given value.
   *
   * @param element Value of the element.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  Multi_field_element_with_small_characteristics(Integer_type element)
      : element_(_get_value(element)) {
    static_assert(maximum >= 2, "Characteristics has to be positive.");
    static_assert(minimum <= maximum, "The given interval is not valid.");
    static_assert(minimum != maximum || _is_prime(minimum), "The given interval does not contain a prime number.");
    static_assert(productOfAllCharacteristics_ != 1, "The given interval does not contain a prime number.");
  }
  /**
   * @brief Copy constructor.
   *
   * @param toCopy Element to copy.
   */
  Multi_field_element_with_small_characteristics(const Multi_field_element_with_small_characteristics& toCopy)
      : element_(toCopy.element_) {}
  /**
   * @brief Move constructor.
   *
   * @param toMove Element to move.
   */
  Multi_field_element_with_small_characteristics(Multi_field_element_with_small_characteristics&& toMove) noexcept
      : element_(std::exchange(toMove.element_, 0)) {}

  /**
   * @brief operator+=
   */
  friend void operator+=(Multi_field_element_with_small_characteristics& f1,
                         Multi_field_element_with_small_characteristics const& f2) {
    f1.element_ = _add(f1.element_, f2.element_);
  }
  /**
   * @brief operator+
   */
  friend Multi_field_element_with_small_characteristics operator+(
      Multi_field_element_with_small_characteristics f1, Multi_field_element_with_small_characteristics const& f2) {
    f1 += f2;
    return f1;
  }
  /**
   * @brief operator+=
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend void operator+=(Multi_field_element_with_small_characteristics& f, const Integer_type& v) {
    f.element_ = _add(f.element_, _get_value(v));
  }
  /**
   * @brief operator+
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Multi_field_element_with_small_characteristics operator+(Multi_field_element_with_small_characteristics f,
                                                                  const Integer_type& v) {
    f += v;
    return f;
  }
  /**
   * @brief operator+
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Integer_type operator+(const Integer_type& v, Multi_field_element_with_small_characteristics f) {
    f += v;
    return f.element_;
  }

  /**
   * @brief operator-=
   */
  friend void operator-=(Multi_field_element_with_small_characteristics& f1,
                         Multi_field_element_with_small_characteristics const& f2) {
    f1.element_ = _subtract(f1.element_, f2.element_);
  }
  /**
   * @brief operator-
   */
  friend Multi_field_element_with_small_characteristics operator-(
      Multi_field_element_with_small_characteristics f1, Multi_field_element_with_small_characteristics const& f2) {
    f1 -= f2;
    return f1;
  }
  /**
   * @brief operator-=
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend void operator-=(Multi_field_element_with_small_characteristics& f, const Integer_type& v) {
    f.element_ = _subtract(f.element_, _get_value(v));
  }
  /**
   * @brief operator-
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Multi_field_element_with_small_characteristics operator-(Multi_field_element_with_small_characteristics f,
                                                                  const Integer_type& v) {
    f -= v;
    return f;
  }
  /**
   * @brief operator-
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Integer_type operator-(const Integer_type& v, const Multi_field_element_with_small_characteristics& f) {
    return _subtract(_get_value(v), f.element_);
  }

  /**
   * @brief operator*=
   */
  friend void operator*=(Multi_field_element_with_small_characteristics& f1,
                         Multi_field_element_with_small_characteristics const& f2) {
    f1.element_ = _multiply(f1.element_, f2.element_);
  }
  /**
   * @brief operator*
   */
  friend Multi_field_element_with_small_characteristics operator*(
      Multi_field_element_with_small_characteristics f1, Multi_field_element_with_small_characteristics const& f2) {
    f1 *= f2;
    return f1;
  }
  /**
   * @brief operator*=
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend void operator*=(Multi_field_element_with_small_characteristics& f, const Integer_type& v) {
    f.element_ = _multiply(f.element_, _get_value(v));
  }
  /**
   * @brief operator*
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Multi_field_element_with_small_characteristics operator*(Multi_field_element_with_small_characteristics f,
                                                                  const Integer_type& v) {
    f *= v;
    return f;
  }
  /**
   * @brief operator*
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Integer_type operator*(const Integer_type& v, Multi_field_element_with_small_characteristics f) {
    f *= v;
    return f.element_;
  }

  /**
   * @brief operator==
   */
  friend bool operator==(const Multi_field_element_with_small_characteristics& f1,
                         const Multi_field_element_with_small_characteristics& f2) {
    return f1.element_ == f2.element_;
  }
  /**
   * @brief operator==
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator==(const Integer_type v, const Multi_field_element_with_small_characteristics& f) {
    return _get_value(v) == f.element_;
  }
  /**
   * @brief operator==
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator==(const Multi_field_element_with_small_characteristics& f, const Integer_type v) {
    return _get_value(v) == f.element_;
  }
  /**
   * @brief operator!=
   */
  friend bool operator!=(const Multi_field_element_with_small_characteristics& f1,
                         const Multi_field_element_with_small_characteristics& f2) {
    return !(f1 == f2);
  }
  /**
   * @brief operator!=
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator!=(const Integer_type v, const Multi_field_element_with_small_characteristics& f) {
    return !(v == f);
  }
  /**
   * @brief operator!=
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator!=(const Multi_field_element_with_small_characteristics& f, const Integer_type v) {
    return !(v == f);
  }

  /**
   * @brief Assign operator.
   */
  Multi_field_element_with_small_characteristics& operator=(Multi_field_element_with_small_characteristics other) {
    std::swap(element_, other.element_);
    return *this;
  }
  /**
   * @brief Assign operator.
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  Multi_field_element_with_small_characteristics& operator=(const Integer_type& value) {
    element_ = _get_value(value);
    return *this;
  }
  /**
   * @brief Swap operator.
   */
  friend void swap(Multi_field_element_with_small_characteristics& f1,
                   Multi_field_element_with_small_characteristics& f2) {
    std::swap(f1.element_, f2.element_);
  }

  /**
   * @brief Casts the element into an unsigned int.
   */
  operator unsigned int() const { return element_; }

  /**
   * @brief Returns the inverse of the element in the multi-field, see @cite boissonnat:hal-00922572.
   *
   * @return The inverse.
   */
  Multi_field_element_with_small_characteristics get_inverse() const {
    return get_partial_inverse(productOfAllCharacteristics_).first;
  }
  /**
   * @brief Returns the inverse of the element with respect to a sub-product of the characteristics in the multi-field,
   * see @cite boissonnat:hal-00922572.
   *
   * @param productOfCharacteristics Sub-product of the characteristics.
   * @return Pair of the inverse and the characteristic the inverse corresponds to.
   */
  std::pair<Multi_field_element_with_small_characteristics,Characteristic> get_partial_inverse(
      Characteristic productOfCharacteristics) const {
    Characteristic gcd = std::gcd(element_, productOfAllCharacteristics_);

    if (gcd == productOfCharacteristics)
      return {Multi_field_element_with_small_characteristics(), multiplicativeID_};  // partial inverse is 0

    Characteristic QT = productOfCharacteristics / gcd;

    const Element inv_qt = _get_inverse(element_, QT);

    auto res = get_partial_multiplicative_identity(QT);
    res *= inv_qt;

    return {res, QT};
  }

  /**
   * @brief Returns the additive identity of a field.
   *
   * @return The additive identity of a field.
   */
  static Multi_field_element_with_small_characteristics get_additive_identity() {
    return Multi_field_element_with_small_characteristics<minimum, maximum>();
  }
  /**
   * @brief Returns the multiplicative identity of a field.
   *
   * @return The multiplicative identity of a field.
   */
  static Multi_field_element_with_small_characteristics get_multiplicative_identity() {
    return Multi_field_element_with_small_characteristics<minimum, maximum>(multiplicativeID_);
  }
  /**
   * @brief Returns the partial multiplicative identity of the multi-field from the given product.
   * See @cite boissonnat:hal-00922572 for more details.
   *
   * @param productOfCharacteristics Product of the different characteristics to take into account in the multi-field.
   * @return The partial multiplicative identity of the multi-field.
   */
  static Multi_field_element_with_small_characteristics get_partial_multiplicative_identity(
      const Characteristic& productOfCharacteristics) {
    if (productOfCharacteristics == 0) {
      return Multi_field_element_with_small_characteristics<minimum, maximum>(multiplicativeID_);
    }
    Multi_field_element_with_small_characteristics<minimum, maximum> mult;
    for (Characteristic idx = 0; idx < primes_.size(); ++idx) {
      if ((productOfCharacteristics % primes_[idx]) == 0) {
        mult += partials_[idx];
      }
    }
    return mult;
  }
  /**
   * @brief Returns the product of all characteristics.
   *
   * @return The product of all characteristics.
   */
  static constexpr Characteristic get_characteristic() { return productOfAllCharacteristics_; }

  /**
   * @brief Returns the value of the element.
   *
   * @return Value of the element.
   */
  Element get_value() const { return element_; }

  // static constexpr bool handles_only_z2() { return false; }

 private:
  static constexpr bool _is_prime(const unsigned int p) {
    if (p <= 1) return false;
    if (p <= 3) return true;
    if (p % 2 == 0 || p % 3 == 0) return false;

    for (unsigned long i = 5; i * i <= p; i = i + 6)
      if (p % i == 0 || p % (i + 2) == 0) return false;

    return true;
  }
  static constexpr Element _multiply(Element a, Element b) {
    Element res = 0;
    Element temp_b = 0;

    if (b < a) std::swap(a, b);

    while (a != 0) {
      if (a & 1) {
        /* Add b to res, modulo m, without overflow */
        if (b >= productOfAllCharacteristics_ - res) res -= productOfAllCharacteristics_;
        res += b;
      }
      a >>= 1;

      /* Double b, modulo m */
      temp_b = b;
      if (b >= productOfAllCharacteristics_ - b) temp_b -= productOfAllCharacteristics_;
      b += temp_b;
    }
    return res;
  }
  static constexpr Element _add(Element element, Element v) {
    if (UINT_MAX - element < v) {
      // automatic unsigned integer overflow behaviour will make it work
      element += v;
      element -= productOfAllCharacteristics_;
      return element;
    }

    element += v;
    if (element >= productOfAllCharacteristics_) element -= productOfAllCharacteristics_;

    return element;
  }
  static constexpr Element _subtract(Element element, Element v) {
    if (element < v) {
      element += productOfAllCharacteristics_;
    }
    element -= v;

    return element;
  }
  static constexpr int _get_inverse(Element element, const Element mod) {
    // to solve: Ax + My = 1
    int M = mod;
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

    if (x < 0) x += mod;

    return x;
  }

  template <typename Integer_type, class = isInteger<Integer_type> >
  static constexpr Element _get_value(Integer_type e) {
    if constexpr (std::is_signed_v<Integer_type>){
      if (e < -static_cast<Integer_type>(productOfAllCharacteristics_)) e = e % productOfAllCharacteristics_;
      if (e < 0) return e += productOfAllCharacteristics_;
      return e < static_cast<Integer_type>(productOfAllCharacteristics_) ? e : e % productOfAllCharacteristics_;
    } else {
      return e < productOfAllCharacteristics_ ? e : e % productOfAllCharacteristics_;
    }
  }

  Element element_;
  static inline const std::vector<Characteristic> primes_ = []() {
    std::vector<Characteristic> res;
    for (Characteristic i = minimum; i <= maximum; ++i) {
      if (_is_prime(i)) {
        res.push_back(i);
      }
    }
    return res;
  }();
  static inline constexpr Characteristic productOfAllCharacteristics_ = []() {
    Characteristic res = 1;
    for (Characteristic i = minimum; i <= maximum; ++i) {
      if (_is_prime(i)) {
        res *= i;
      }
    }
    return res;
  }();
  static inline const std::vector<Characteristic> partials_ = []() {
    std::vector<Characteristic> res;

    if (productOfAllCharacteristics_ == 1) return res;

    for (Characteristic i = 0; i < primes_.size(); ++i) {
      Characteristic p = primes_[i];
      Characteristic base = productOfAllCharacteristics_ / p;
      Characteristic exp = p - 1;
      res.push_back(1);

      while (exp > 0) {
        // If exp is odd, multiply with result
        if (exp & 1) res.back() = _multiply(res.back(), base);
        // y must be even now
        exp = exp >> 1;
        base = _multiply(base, base);
      }
    }

    return res;
  }();
  // If I understood the paper well, multiplicativeID_ always equals to 1. But in Clement's code,
  // multiplicativeID_ is computed (see commented lambda function below). TODO: verify with Clement.
  static inline constexpr Element multiplicativeID_ = 1; /*= [](){
          unsigned int res = 0;
          for (unsigned int i = 0; i < partials_.size(); ++i){
                  res = (res + partials_[i]) % productOfAllCharacteristics_;
          }

          return res;
  }();*/
};

}  // namespace persistence_fields
}  // namespace Gudhi

#endif  // MATRIX_FIELD_MULTI_SMALL_H_

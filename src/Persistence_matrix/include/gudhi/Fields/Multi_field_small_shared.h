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
 * @file Multi_field_small_shared.h
 * @author Hannah Schreiber, Clément Maria
 * @brief Contains the @ref Gudhi::persistence_fields::Shared_multi_field_element_with_small_characteristics class.
 */

#ifndef MATRIX_FIELD_MULTI_SMALL_SHARED_H_
#define MATRIX_FIELD_MULTI_SMALL_SHARED_H_

#include <utility>
#include <vector>
#include <limits.h>
#include <stdexcept>
#include <numeric>

namespace Gudhi {
namespace persistence_fields {

/**
 * @class Shared_multi_field_element_with_small_characteristics Multi_field_small_shared.h \
 * gudhi/Fields/Multi_field_small_shared.h
 * @ingroup persistence_fields
 *
 * @brief Class representing an element of a multi-field, such that `productOfAllCharacteristics ^ 2` fits into
 * the given @p Unsigned_integer_type template argument. If each instantiation of the class can represent another
 * element, they all share the same characteristics. That is if the characteristics are set for one, they will be
 * set for all the others. The characteristics can be set before instantiating the elements with the static
 * @ref Shared_multi_field_element_with_small_characteristics::initialize method.
 *
 * @tparam Unsigned_integer_type A native unsigned integer type: unsigned int, long unsigned int, etc.
 * Will be used as the field element type.
 */
template <typename Unsigned_integer_type = unsigned int,
          class = std::enable_if_t<std::is_unsigned_v<Unsigned_integer_type> > >
class Shared_multi_field_element_with_small_characteristics {
 public:
  using Element = Unsigned_integer_type; /**< Type for the elements in the field. */
  using Characteristic = Element;   /**< Type for the field characteristic. */
  template <class T>
  using isInteger = std::enable_if_t<std::is_integral_v<T> >;

  /**
   * @brief Default constructor. Sets the element to 0.
   */
  Shared_multi_field_element_with_small_characteristics() : element_(0) {}
  /**
   * @brief Constructor setting the element to the given value.
   *
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if negative.
   * @param element Value of the element.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  Shared_multi_field_element_with_small_characteristics(Integer_type element) : element_(_get_value(element)) {}
  /**
   * @brief Copy constructor.
   *
   * @param toCopy Element to copy.
   */
  Shared_multi_field_element_with_small_characteristics(
      const Shared_multi_field_element_with_small_characteristics& toCopy)
      : element_(toCopy.element_) {}
  /**
   * @brief Move constructor.
   *
   * @param toMove Element to move.
   */
  Shared_multi_field_element_with_small_characteristics(
      Shared_multi_field_element_with_small_characteristics&& toMove) noexcept
      : element_(std::exchange(toMove.element_, 0)) {}

  /**
   * @brief Initialize the multi-field to the characteristics (primes) contained in the given interval.
   * Should be called first before constructing the field elements.
   * The characteristics must be small enough such that `productOfAllCharacteristics ^ 2` fits into an unsigned int.
   *
   * @param minimum Lowest value in the interval.
   * @param maximum Highest value in the interval.
   */
  static void initialize(unsigned int minimum, unsigned int maximum) {
    if (maximum < 2) throw std::invalid_argument("Characteristic must be strictly positive");
    if (minimum > maximum) throw std::invalid_argument("The given interval is not valid.");
    if (minimum == maximum && !_is_prime(minimum))
      throw std::invalid_argument("The given interval does not contain a prime number.");

    productOfAllCharacteristics_ = 1;
    primes_.clear();
    for (unsigned int i = minimum; i <= maximum; ++i) {
      if (_is_prime(i)) {
        primes_.push_back(i);
        productOfAllCharacteristics_ *= i;
      }
    }

    if (primes_.empty()) throw std::invalid_argument("The given interval does not contain a prime number.");

    partials_.resize(primes_.size());
    for (Characteristic i = 0; i < primes_.size(); ++i) {
      Characteristic p = primes_[i];
      Characteristic base = productOfAllCharacteristics_ / p;
      Characteristic exp = p - 1;
      partials_[i] = 1;

      while (exp > 0) {
        // If exp is odd, multiply with result
        if (exp & 1) partials_[i] = _multiply(partials_[i], base);
        // y must be even now
        exp = exp >> 1;  // y = y/2
        base = _multiply(base, base);
      }
    }

    // If I understood the paper well, multiplicativeID_ always equals to 1. But in Clement's code,
    // multiplicativeID_ is computed (see commented loop below). TODO: verify with Clement.
    //	for (unsigned int i = 0; i < partials_.size(); ++i){
    //		multiplicativeID_ = (multiplicativeID_ + partials_[i]) % productOfAllCharacteristics_;
    //	}
  }

  /**
   * @brief operator+=
   */
  friend void operator+=(Shared_multi_field_element_with_small_characteristics& f1,
                         Shared_multi_field_element_with_small_characteristics const& f2) {
    f1.element_ = _add(f1.element_, f2.element_);
  }
  /**
   * @brief operator+
   */
  friend Shared_multi_field_element_with_small_characteristics operator+(
      Shared_multi_field_element_with_small_characteristics f1,
      Shared_multi_field_element_with_small_characteristics const& f2) {
    f1 += f2;
    return f1;
  }
  /**
   * @brief operator+=
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend void operator+=(Shared_multi_field_element_with_small_characteristics& f, const Integer_type& v) {
    f.element_ = _add(f.element_, _get_value(v));
  }
  /**
   * @brief operator+
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Shared_multi_field_element_with_small_characteristics operator+(
      Shared_multi_field_element_with_small_characteristics f, const Integer_type& v) {
    f += v;
    return f;
  }
  /**
   * @brief operator+
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Integer_type operator+(const Integer_type& v, Shared_multi_field_element_with_small_characteristics f) {
    f += v;
    return f.element_;
  }

  /**
   * @brief operator-=
   */
  friend void operator-=(Shared_multi_field_element_with_small_characteristics& f1,
                         Shared_multi_field_element_with_small_characteristics const& f2) {
    f1.element_ = _subtract(f1.element_, f2.element_);
  }
  /**
   * @brief operator-
   */
  friend Shared_multi_field_element_with_small_characteristics operator-(
      Shared_multi_field_element_with_small_characteristics f1,
      Shared_multi_field_element_with_small_characteristics const& f2) {
    f1 -= f2;
    return f1;
  }
  /**
   * @brief operator-=
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend void operator-=(Shared_multi_field_element_with_small_characteristics& f, const Integer_type& v) {
    f.element_ = _subtract(f.element_, _get_value(v));
  }
  /**
   * @brief operator-
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Shared_multi_field_element_with_small_characteristics operator-(
      Shared_multi_field_element_with_small_characteristics f, const Integer_type& v) {
    f -= v;
    return f;
  }
  /**
   * @brief operator-
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Integer_type operator-(const Integer_type& v, const Shared_multi_field_element_with_small_characteristics& f) {
    return _subtract(_get_value(v), f.element_);
  }

  /**
   * @brief operator*=
   */
  friend void operator*=(Shared_multi_field_element_with_small_characteristics& f1,
                         Shared_multi_field_element_with_small_characteristics const& f2) {
    f1.element_ = _multiply(f1.element_, f2.element_);
  }
  /**
   * @brief operator*
   */
  friend Shared_multi_field_element_with_small_characteristics operator*(
      Shared_multi_field_element_with_small_characteristics f1,
      Shared_multi_field_element_with_small_characteristics const& f2) {
    f1 *= f2;
    return f1;
  }
  /**
   * @brief operator*=
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend void operator*=(Shared_multi_field_element_with_small_characteristics& f, const Integer_type& v) {
    f.element_ = _multiply(f.element_, _get_value(v));
  }
  /**
   * @brief operator*
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Shared_multi_field_element_with_small_characteristics operator*(
      Shared_multi_field_element_with_small_characteristics f, const Integer_type& v) {
    f *= v;
    return f;
  }
  /**
   * @brief operator*
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend Integer_type operator*(const Integer_type& v, Shared_multi_field_element_with_small_characteristics f) {
    f *= v;
    return f.element_;
  }

  /**
   * @brief operator==
   */
  friend bool operator==(const Shared_multi_field_element_with_small_characteristics& f1,
                         const Shared_multi_field_element_with_small_characteristics& f2) {
    return f1.element_ == f2.element_;
  }
  /**
   * @brief operator==
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator==(const Integer_type& v, const Shared_multi_field_element_with_small_characteristics& f) {
    return _get_value(v) == f.element_;
  }
  /**
   * @brief operator==
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator==(const Shared_multi_field_element_with_small_characteristics& f, const Integer_type& v) {
    return _get_value(v) == f.element_;
  }
  /**
   * @brief operator!=
   */
  friend bool operator!=(const Shared_multi_field_element_with_small_characteristics& f1,
                         const Shared_multi_field_element_with_small_characteristics& f2) {
    return !(f1 == f2);
  }
  /**
   * @brief operator!=
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator!=(const Integer_type v, const Shared_multi_field_element_with_small_characteristics& f) {
    return !(v == f);
  }
  /**
   * @brief operator!=
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  friend bool operator!=(const Shared_multi_field_element_with_small_characteristics& f, const Integer_type v) {
    return !(v == f);
  }

  /**
   * @brief Assign operator.
   */
  Shared_multi_field_element_with_small_characteristics& operator=(
      Shared_multi_field_element_with_small_characteristics other) {
    std::swap(element_, other.element_);
    return *this;
  }
  /**
   * @brief Assign operator.
   * 
   * @tparam Integer_type A native integer type. Should be able to contain the characteristic if signed.
   */
  template <typename Integer_type, class = isInteger<Integer_type> >
  Shared_multi_field_element_with_small_characteristics& operator=(const Integer_type& value) {
    element_ = _get_value(value);
    return *this;
  }
  /**
   * @brief Swap operator.
   */
  friend void swap(Shared_multi_field_element_with_small_characteristics& f1,
                   Shared_multi_field_element_with_small_characteristics& f2) {
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
  Shared_multi_field_element_with_small_characteristics get_inverse() const {
    return get_partial_inverse(productOfAllCharacteristics_).first;
  }
  /**
   * @brief Returns the inverse of the element with respect to a sub-product of the characteristics in the multi-field,
   * see @cite boissonnat:hal-00922572.
   *
   * @param productOfCharacteristics Sub-product of the characteristics.
   * @return Pair of the inverse and the characteristic the inverse corresponds to.
   */
  std::pair<Shared_multi_field_element_with_small_characteristics,Characteristic> get_partial_inverse(
      Characteristic productOfCharacteristics) const {
    Characteristic gcd = std::gcd(element_, productOfAllCharacteristics_);

    if (gcd == productOfCharacteristics)
      return {Shared_multi_field_element_with_small_characteristics(), multiplicativeID_};  // partial inverse is 0

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
  static Shared_multi_field_element_with_small_characteristics get_additive_identity() {
    return Shared_multi_field_element_with_small_characteristics();
  }
  /**
   * @brief Returns the multiplicative identity of a field.
   *
   * @return The multiplicative identity of a field.
   */
  static Shared_multi_field_element_with_small_characteristics get_multiplicative_identity() {
    return Shared_multi_field_element_with_small_characteristics(multiplicativeID_);
  }
  /**
   * @brief Returns the partial multiplicative identity of the multi-field from the given product.
   * See @cite boissonnat:hal-00922572 for more details.
   *
   * @param productOfCharacteristics Product of the different characteristics to take into account in the multi-field.
   * @return The partial multiplicative identity of the multi-field.
   */
  static Shared_multi_field_element_with_small_characteristics get_partial_multiplicative_identity(
      const Characteristic& productOfCharacteristics) {
    if (productOfCharacteristics == 0) {
      return Shared_multi_field_element_with_small_characteristics(multiplicativeID_);
    }
    Shared_multi_field_element_with_small_characteristics mult;
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
  static Characteristic get_characteristic() { return productOfAllCharacteristics_; }

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
  static Element _multiply(Element a, Element b) {
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
  static Element _add(Element element, Element v) {
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
  static Element _subtract(Element element, Element v) {
    if (element < v) {
      element += productOfAllCharacteristics_;
    }
    element -= v;

    return element;
  }
  static constexpr int _get_inverse(Element element, const Characteristic mod) {
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

  Element element_;                                          /**< Element. */
  static inline std::vector<Characteristic> primes_;         /**< All characteristics. */
  static inline Characteristic productOfAllCharacteristics_; /**< Product of all characteristics. */
  static inline std::vector<Characteristic> partials_;       /**< Partial products of the characteristics. */
  static inline constexpr Element multiplicativeID_ = 1;     /**< Multiplicative identity. */
};

}  // namespace persistence_fields
}  // namespace Gudhi

#endif  // MATRIX_FIELD_MULTI_SMALL_SHARED_H_

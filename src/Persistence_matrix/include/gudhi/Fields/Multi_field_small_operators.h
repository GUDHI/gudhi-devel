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
 * @file Multi_field_small_operators.h
 * @author Hannah Schreiber, Clément Maria
 * @brief Contains the @ref Gudhi::persistence_fields::Multi_field_operators_with_small_characteristics class.
 */

#ifndef MATRIX_FIELD_MULTI_SMALL_OPERATORS_H_
#define MATRIX_FIELD_MULTI_SMALL_OPERATORS_H_

#include <utility>
#include <vector>
#include <limits.h>
#include <stdexcept>
#include <numeric>

namespace Gudhi {
namespace persistence_fields {

/**
 * @class Multi_field_operators_with_small_characteristics Multi_field_small_operators.h \
 * gudhi/Fields/Multi_field_small_operators.h
 * @ingroup persistence_fields
 *
 * @brief Class defining operators for a multi-field with "consecutive" characteristic range, such that
 * `productOfAllCharacteristics ^ 2` fits into an unsigned int.
 */
class Multi_field_operators_with_small_characteristics 
{
 public:
  using Element = unsigned int;          /**< Type for the elements in the field. */
  using Characteristic = Element;   /**< Type for the field characteristic. */

  /**
   * @brief Default constructor, sets the product of all characteristics to 0.
   */
  Multi_field_operators_with_small_characteristics() : productOfAllCharacteristics_(0) /* , multiplicativeID_(1) */ 
  {}
  /**
   * @brief Constructor setting the characteristics to all prime numbers between the two given integers.
   * The product of all primes to the square has to fit into an unsigned int.
   * 
   * @param minCharacteristic Smallest value of a prime.
   * @param maxCharacteristic Highest value of a prime.
   */
  Multi_field_operators_with_small_characteristics(int minCharacteristic, int maxCharacteristic)
      : productOfAllCharacteristics_(0)  //, multiplicativeID_(1)
  {
    set_characteristic(minCharacteristic, maxCharacteristic);
  }
  /**
   * @brief Copy constructor.
   * 
   * @param toCopy Operators to copy.
   */
  Multi_field_operators_with_small_characteristics(const Multi_field_operators_with_small_characteristics& toCopy)
      : primes_(toCopy.primes_),
        productOfAllCharacteristics_(toCopy.productOfAllCharacteristics_),
        partials_(toCopy.partials_) /* ,
         multiplicativeID_(toCopy.multiplicativeID_) */
  {}
  /**
   * @brief Move constructor.
   * 
   * @param toMove Operators to move.
   */
  Multi_field_operators_with_small_characteristics(Multi_field_operators_with_small_characteristics&& toMove) noexcept
      : primes_(std::move(toMove.primes_)),
        productOfAllCharacteristics_(std::move(toMove.productOfAllCharacteristics_)),
        partials_(std::move(toMove.partials_)) /* ,
         multiplicativeID_(std::move(toMove.multiplicativeID_)) */
  {}

  /**
   * @brief Set the characteristics of the field, which are stored in a single value as a product of all of them.
   * The characteristics will be all prime numbers in the given interval.
   * The product of all primes to the square has to fit into an unsigned int.
   * 
   * @param minimum Smallest value of a prime.
   * @param maximum Highest value of a prime.
   */
  void set_characteristic(int minimum, int maximum) {
    if (maximum < 2) throw std::invalid_argument("Characteristic must be strictly positive");
    if (minimum > maximum) throw std::invalid_argument("The given interval is not valid.");
    if (minimum == maximum && !_is_prime(minimum))
      throw std::invalid_argument("The given interval does not contain a prime number.");

    productOfAllCharacteristics_ = 1;
    primes_.clear();
    for (unsigned int i = minimum; i <= static_cast<unsigned int>(maximum); ++i) {
      if (_is_prime(i)) {
        primes_.push_back(i);
        productOfAllCharacteristics_ *= i;
      }
    }

    if (primes_.empty()) throw std::invalid_argument("The given interval does not contain a prime number.");

    partials_.resize(primes_.size());
    for (unsigned int i = 0; i < primes_.size(); ++i) {
      unsigned int p = primes_[i];
      Characteristic base = productOfAllCharacteristics_ / p;
      unsigned int exp = p - 1;
      partials_[i] = 1;

      while (exp > 0) {
        // If exp is odd, multiply with result
        if (exp & 1) partials_[i] = _multiply(partials_[i], base, productOfAllCharacteristics_);
        // y must be even now
        exp = exp >> 1;  // y = y/2
        base = _multiply(base, base, productOfAllCharacteristics_);
      }
    }

    // If I understood the paper well, multiplicativeID_ always equals to 1. But in Clement's code,
    // multiplicativeID_ is computed (see commented loop below). TODO: verify with Clement.
    //	for (unsigned int i = 0; i < partials_.size(); ++i){
    //		multiplicativeID_ = (multiplicativeID_ + partials_[i]) % productOfAllCharacteristics_;
    //	}
  }
  /**
   * @brief Returns the current characteristics as the product of all of them.
   * 
   * @return The value of the current characteristic.
   */
  const Characteristic& get_characteristic() const { return productOfAllCharacteristics_; }

  /**
   * @brief Returns the value of an element in the field.
   * That is the positive value of the integer modulo the current characteristic.
   * 
   * @param e Integer to return the value from.
   * @return @p e modulo the current characteristic, such that the result is positive.
   */
  Element get_value(Element e) const {
    return e < productOfAllCharacteristics_ ? e : e % productOfAllCharacteristics_;
  }

  /**
   * @brief Returns the sum of two elements in the field.
   * 
   * @param e1 First element.
   * @param e2 Second element.
   * @return `(e1 + e2) % productOfAllCharacteristics`, such that the result is positive.
   */
  Element add(Element e1, Element e2) const {
    return _add(get_value(e1), get_value(e2), productOfAllCharacteristics_);
  }

  /**
   * @brief Stores in the first element the sum of two given elements in the field, that is
   * `(e1 + e2) % productOfAllCharacteristics`, such that the result is positive.
   * 
   * @param e1 First element.
   * @param e2 Second element.
   */
  void add_inplace(Element& e1, Element e2) const {
    e1 = _add(get_value(e1), get_value(e2), productOfAllCharacteristics_);
  }

  /**
   * @brief Returns the subtraction in the field of the first element by the second element.
   * 
   * @param e1 First element.
   * @param e2 Second element.
   * @return `(e1 - e2) % productOfAllCharacteristics`, such that the result is positive.
   */
  Element subtract(Element e1, Element e2) const {
    return _subtract(get_value(e1), get_value(e2), productOfAllCharacteristics_);
  }

  /**
   * @brief Stores in the first element the subtraction in the field of the first element by the second element,
   * that is `(e1 - e2) % productOfAllCharacteristics`, such that the result is positive.
   * 
   * @param e1 First element.
   * @param e2 Second element.
   */
  void subtract_inplace_front(Element& e1, Element e2) const {
    e1 = _subtract(get_value(e1), get_value(e2), productOfAllCharacteristics_);
  }
  /**
   * @brief Stores in the second element the subtraction in the field of the first element by the second element,
   * that is `(e1 - e2) % productOfAllCharacteristics`, such that the result is positive.
   * 
   * @param e1 First element.
   * @param e2 Second element.
   */
  void subtract_inplace_back(Element e1, Element& e2) const {
    e2 = _subtract(get_value(e1), get_value(e2), productOfAllCharacteristics_);
  }

  /**
   * @brief Returns the multiplication of two elements in the field.
   * 
   * @param e1 First element.
   * @param e2 Second element.
   * @return `(e1 * e2) % productOfAllCharacteristics`, such that the result is positive.
   */
  Element multiply(Element e1, Element e2) const {
    return _multiply(get_value(e1), get_value(e2), productOfAllCharacteristics_);
  }

  /**
   * @brief Stores in the first element the multiplication of two given elements in the field,
   * that is `(e1 * e2) % productOfAllCharacteristics`, such that the result is positive.
   * 
   * @param e1 First element.
   * @param e2 Second element.
   */
  void multiply_inplace(Element& e1, Element e2) const {
    e1 = _multiply(get_value(e1), get_value(e2), productOfAllCharacteristics_);
  }

  /**
   * @brief Multiplies the first element with the second one and adds the third one. Returns the result in the field.
   *
   * @warning Not overflow safe.
   * 
   * @param e First element.
   * @param m Second element.
   * @param a Third element.
   * @return `(e * m + a) % productOfAllCharacteristics`, such that the result is positive.
   */
  Element multiply_and_add(Element e, Element m, Element a) const { return get_value(e * m + a); }

  /**
   * @brief Multiplies the first element with the second one and adds the third one, that is
   * `(e * m + a) % productOfAllCharacteristics`, such that the result is positive.
   * Stores the result in the first element.
   *
   * @warning Not overflow safe.
   * 
   * @param e First element.
   * @param m Second element.
   * @param a Third element.
   */
  void multiply_and_add_inplace_front(Element& e, Element m, Element a) const {
    e = get_value(e * m + a);
  }
  /**
   * @brief Multiplies the first element with the second one and adds the third one, that is
   * `(e * m + a) % productOfAllCharacteristics`, such that the result is positive.
   * Stores the result in the third element.
   *
   * @warning Not overflow safe.
   * 
   * @param e First element.
   * @param m Second element.
   * @param a Third element.
   */
  void multiply_and_add_inplace_back(Element e, Element m, Element& a) const {
    a = get_value(e * m + a);
  }

  /**
   * @brief Adds the first element to the second one and multiplies the third one with it.
   * Returns the result in the field.
   *
   * @warning Not overflow safe.
   * 
   * @param e First element.
   * @param a Second element.
   * @param m Third element.
   * @return `((e + a) * m) % productOfAllCharacteristics`, such that the result is positive.
   */
  Element add_and_multiply(Element e, Element a, Element m) const { return get_value((e + a) * m); }

  /**
   * @brief Adds the first element to the second one and multiplies the third one with it, that is
   * `((e + a) * m) % productOfAllCharacteristics`, such that the result is positive.
   * Stores the result in the first element.
   *
   * @warning Not overflow safe.
   * 
   * @param e First element.
   * @param a Second element.
   * @param m Third element.
   */
  void add_and_multiply_inplace_front(Element& e, Element a, Element m) const {
    e = get_value((e + a) * m);
  }
  /**
   * @brief Adds the first element to the second one and multiplies the third one with it, that is
   * `((e + a) * m) % productOfAllCharacteristics`, such that the result is positive.
   * Stores the result in the third element.
   *
   * @warning Not overflow safe.
   * 
   * @param e First element.
   * @param a Second element.
   * @param m Third element.
   */
  void add_and_multiply_inplace_back(Element e, Element a, Element& m) const {
    m = get_value((e + a) * m);
  }

  /**
   * @brief Returns true if the two given elements are equal in the field, false otherwise.
   * 
   * @param e1 First element to compare.
   * @param e2 Second element to compare.
   * @return true If `e1 % productOfAllCharacteristics == e2 % productOfAllCharacteristics`.
   * @return false Otherwise.
   */
  bool are_equal(Element e1, Element e2) const { return get_value(e1) == get_value(e2); }

  /**
   * @brief Returns the inverse of the given element in the sense of @cite boissonnat:hal-00922572 with respect
   * to the product of all characteristics.
   * 
   * @param e Element to get the inverse from.
   * @return Inverse in the current field.
   */
  Element get_inverse(const Element& e) const {
    return get_partial_inverse(e, productOfAllCharacteristics_).first;
  }
  /**
   * @brief Returns the inverse of the given element in the multi-field corresponding to the given sub-product
   * of the product of all characteristics in the multi-field. See @cite boissonnat:hal-00922572 for more details.
   * 
   * @param e Element to get the inverse from.
   * @param productOfCharacteristics Product of the different characteristics to take into account in the multi-field.
   * @return Pair of the inverse of @p e and the characteristic the inverse is coming from.
   */
  std::pair<Element, Characteristic> get_partial_inverse(
      const Element& e, const Characteristic& productOfCharacteristics) const {
    Characteristic gcd = std::gcd(e, productOfAllCharacteristics_);

    if (gcd == productOfCharacteristics) return {0, get_multiplicative_identity()};  // partial inverse is 0

    Characteristic QT = productOfCharacteristics / gcd;

    const Characteristic inv_qt = _get_inverse(e, QT);

    auto res = get_partial_multiplicative_identity(QT);
    res = _multiply(res, inv_qt, productOfAllCharacteristics_);

    return {res, QT};
  }

  /**
   * @brief Returns the additive identity of a field.
   * 
   * @return The additive identity of a field.
   */
  static constexpr Element get_additive_identity() { return 0; }
  /**
   * @brief Returns the multiplicative identity of a field.
   * 
   * @return The multiplicative identity of a field.
   */
  static constexpr Element get_multiplicative_identity() { return 1; }
  // static Element get_multiplicative_identity(){ return multiplicativeID_; }
  /**
   * @brief Returns the partial multiplicative identity of the multi-field from the given product.
   * See @cite boissonnat:hal-00922572 for more details.
   * 
   * @param productOfCharacteristics Product of the different characteristics to take into account in the multi-field.
   * @return The partial multiplicative identity of the multi-field.
   */
  Element get_partial_multiplicative_identity(const Characteristic& productOfCharacteristics) const {
    if (productOfCharacteristics == 0) {
      return get_multiplicative_identity();
    }
    Element multIdentity = 0;
    for (unsigned int idx = 0; idx < primes_.size(); ++idx) {
      if ((productOfCharacteristics % primes_[idx]) == 0) {
        multIdentity = _add(multIdentity, partials_[idx], productOfAllCharacteristics_);
      }
    }
    return multIdentity;
  }

  // static constexpr bool handles_only_z2() { return false; }

  /**
   * @brief Assign operator.
   */
  Multi_field_operators_with_small_characteristics& operator=(Multi_field_operators_with_small_characteristics other) {
    primes_.swap(other.primes_);
    productOfAllCharacteristics_ = other.productOfAllCharacteristics_;
    partials_.swap(other.partials_);

    return *this;
  }
  /**
   * @brief Swap operator.
   */
  friend void swap(Multi_field_operators_with_small_characteristics& f1,
                   Multi_field_operators_with_small_characteristics& f2) {
    f1.primes_.swap(f2.primes_);
    std::swap(f1.productOfAllCharacteristics_, f2.productOfAllCharacteristics_);
    f1.partials_.swap(f2.partials_);
  }

 private:
  std::vector<unsigned int> primes_;                /**< All characteristics. */
  Characteristic productOfAllCharacteristics_; /**< Product of all characteristics. */
  std::vector<Characteristic> partials_;       /**< Partial products of the characteristics. */
  // static inline constexpr unsigned int multiplicativeID_ = 1;

  static Element _add(Element element, Element v, Characteristic characteristic);
  static Element _subtract(Element element, Element v, Characteristic characteristic);
  static Element _multiply(Element a, Element b, Characteristic characteristic);
  static constexpr long int _get_inverse(Element element, Characteristic mod);
  static constexpr bool _is_prime(const int p);
};

inline Multi_field_operators_with_small_characteristics::Element
Multi_field_operators_with_small_characteristics::_add(Element element, Element v,
                                                       Characteristic characteristic) {
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

inline Multi_field_operators_with_small_characteristics::Element
Multi_field_operators_with_small_characteristics::_subtract(Element element, Element v,
                                                             Characteristic characteristic) {
  if (element < v) {
    element += characteristic;
  }
  element -= v;

  return element;
}

inline Multi_field_operators_with_small_characteristics::Element
Multi_field_operators_with_small_characteristics::_multiply(Element a, Element b,
                                                            Characteristic characteristic) {
  Element res = 0;
  Element temp_b = 0;

  if (b < a) std::swap(a, b);

  while (a != 0) {
    if (a & 1) {
      /* Add b to res, modulo m, without overflow */
      if (b >= characteristic - res) res -= characteristic;
      res += b;
    }
    a >>= 1;

    /* Double b, modulo m */
    temp_b = b;
    if (b >= characteristic - b) temp_b -= characteristic;
    b += temp_b;
  }
  return res;
}

inline constexpr long int Multi_field_operators_with_small_characteristics::_get_inverse(Element element,
                                                                                         Characteristic mod) {
  // to solve: Ax + My = 1
  Element M = mod;
  Element A = element;
  long int y = 0, x = 1;
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

inline constexpr bool Multi_field_operators_with_small_characteristics::_is_prime(const int p) {
  if (p <= 1) return false;
  if (p <= 3) return true;
  if (p % 2 == 0 || p % 3 == 0) return false;

  for (long i = 5; i * i <= p; i = i + 6)
    if (p % i == 0 || p % (i + 2) == 0) return false;

  return true;
}

}  // namespace persistence_fields
}  // namespace Gudhi

#endif  // MATRIX_FIELD_MULTI_SMALL_OPERATORS_H_

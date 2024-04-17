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
 * @file Multi_field_operators.h
 * @author Hannah Schreiber, Clément Maria
 * @brief Contains the @ref Multi_field_operators class.
 */

#ifndef MATRIX_FIELD_MULTI_OPERATORS_H_
#define MATRIX_FIELD_MULTI_OPERATORS_H_

#include <utility>
#include <vector>
#include <gmpxx.h>
#include <stdexcept>

namespace Gudhi {
namespace persistence_fields {

/**
 * @class Multi_field_operators Multi_field_operators.h gudhi/Fields/Multi_field_operators.h
 * @ingroup persistence_fields
 *
 * @brief Class defining operators for a multi-field with "consecutive" charateristic range.
 */
class Multi_field_operators
{
 public:
  using element_type = mpz_class;             /**< Type for the elements in the field. */
  using characteristic_type = element_type;   /**< Type for the field characteristic. */

  /**
   * @brief Default constructor, sets the product of all characteristics to 0.
   */
  Multi_field_operators() : productOfAllCharacteristics_(0) /* , multiplicativeID_(1) */ 
  {}
  /**
   * @brief Constructor setting the characteristics to all prime numbers between the two given integers.
   * 
   * @param minCharacteristic Smallest value of a prime.
   * @param maxCharacteristic Heighest value of a prime.
   */
  Multi_field_operators(int minCharacteristic, int maxCharacteristic)
      : productOfAllCharacteristics_(0)  //, multiplicativeID_(1)
  {
    set_characteristic(minCharacteristic, maxCharacteristic);
  }
  /**
   * @brief Copy constructor.
   * 
   * @param toCopy Operators to copy.
   */
  Multi_field_operators(const Multi_field_operators& toCopy)
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
  Multi_field_operators(Multi_field_operators&& toMove) noexcept
      : primes_(std::move(toMove.primes_)),
        productOfAllCharacteristics_(std::move(toMove.productOfAllCharacteristics_)),
        partials_(std::move(toMove.partials_)) /* ,
         multiplicativeID_(std::move(toMove.multiplicativeID_)) */
  {}

  /**
   * @brief Set the characteristics of the field, which are stored in a single value as a product of all of them.
   * The characteristics will be all prime numbers in the given interval.
   * 
   * @param minimum Smallest value of a prime.
   * @param maximum Heighest value of a prime.
   */
  void set_characteristic(int minimum, int maximum) {
    if (maximum < 2) throw std::invalid_argument("Characteristic must be strictly positive");
    if (minimum > maximum) throw std::invalid_argument("The given interval is not valid.");
    if (minimum == maximum && !_is_prime(minimum))
      throw std::invalid_argument("The given interval does not contain a prime number.");

    unsigned int curr_prime = minimum;
    mpz_t tmp_prime;
    mpz_init_set_ui(tmp_prime, minimum);
    // test if min_prime is prime
    int is_prime = mpz_probab_prime_p(tmp_prime, 25);  // probabilistic primality test

    if (is_prime == 0) {  // min_prime is composite
      mpz_nextprime(tmp_prime, tmp_prime);
      curr_prime = mpz_get_ui(tmp_prime);
    }

    primes_.clear();
    while (curr_prime <= static_cast<unsigned int>(maximum)) {
      primes_.push_back(curr_prime);
      mpz_nextprime(tmp_prime, tmp_prime);
      curr_prime = mpz_get_ui(tmp_prime);
    }
    mpz_clear(tmp_prime);

    if (primes_.empty()) throw std::invalid_argument("The given interval does not contain a prime number.");

    productOfAllCharacteristics_ = 1;
    for (const unsigned int p : primes_) {
      productOfAllCharacteristics_ *= p;
    }

    partials_.resize(primes_.size());
    for (unsigned int i = 0; i < primes_.size(); ++i) {
      unsigned int p = primes_[i];
      partials_[i] = productOfAllCharacteristics_ / p;
      mpz_powm_ui(partials_[i].get_mpz_t(), partials_[i].get_mpz_t(), p - 1, productOfAllCharacteristics_.get_mpz_t());
    }

    // If I understood the paper well, multiplicativeID_ always equals to 1. But in Clement's code,
    // multiplicativeID_ is computed (see commented loop below). TODO: verify with Clement.
    // for (unsigned int i = 0; i < partials_.size(); ++i) {
    //   multiplicativeID_ = (multiplicativeID_ + partials_[i]) % productOfAllCharacteristics_;
    // }
  }
  /**
   * @brief Returns the current characteristics as the product of all of them.
   * 
   * @return The value of the current characteristic.
   */
  characteristic_type get_characteristic() const { return productOfAllCharacteristics_; }

  // TODO: verify if mpz_class is trivial to copy around or not.
  // If yes, remove the unnecessery references?
  // If not, do the operations in-place instead?

  /**
   * @brief Returns the value of an element in the field.
   * That is the positive value of the integer modulo the current characteristic.
   * 
   * @param e Integer to return the value from.
   * @return @p e modulo the current characteristic, such that the result is positive.
   */
  element_type get_value(element_type e) const {
    if (e >= productOfAllCharacteristics_ || e < -productOfAllCharacteristics_)
      mpz_mod(e.get_mpz_t(), e.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
    if (e < 0) mpz_add(e.get_mpz_t(), e.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
    return e;
  }

  /**
   * @brief Returns the sum of two elements in the field.
   * 
   * @param e1 First element.
   * @param e2 Second element.
   * @return `(e1 + e2) % productOfAllCharacteristics`, such that the result is positive.
   */
  element_type add(element_type e1, const element_type& e2) const {
    mpz_add(e1.get_mpz_t(), e1.get_mpz_t(), e2.get_mpz_t());
    return _get_value(e1);
  }

  /**
   * @brief Returns the substraction in the field of the first element by the second element.
   * 
   * @param e1 First element.
   * @param e2 Second element.
   * @return `(e1 - e2) % productOfAllCharacteristics`, such that the result is positive.
   */
  element_type substract(element_type e1, const element_type& e2) const {
    mpz_sub(e1.get_mpz_t(), e1.get_mpz_t(), e2.get_mpz_t());
    return _get_value(e1);
  }

  /**
   * @brief Returns the multiplication of two elements in the field.
   * 
   * @param e1 First element.
   * @param e2 Second element.
   * @return `(e1 * e2) % productOfAllCharacteristics`, such that the result is positive.
   */
  element_type multiply(element_type e1, const element_type& e2) const {
    mpz_mul(e1.get_mpz_t(), e1.get_mpz_t(), e2.get_mpz_t());
    return _get_value(e1);
  }

  /**
   * @brief Multiplies the first element with the second one and adds the third one. Returns the result in the field.
   * 
   * @param e First element.
   * @param m Second element.
   * @param a Third element.
   * @return `(e * m + a) % productOfAllCharacteristics`, such that the result is positive.
   */
  element_type multiply_and_add(element_type e, const element_type& m, const element_type& a) const {
    mpz_mul(e.get_mpz_t(), e.get_mpz_t(), m.get_mpz_t());
    mpz_add(e.get_mpz_t(), e.get_mpz_t(), a.get_mpz_t());
    return _get_value(e);
  }

  /**
   * @brief Adds the first element to the second one and multiplies the third one with it.
   * Returns the result in the field.
   * 
   * @param e First element.
   * @param a Second element.
   * @param m Third element.
   * @return `((e + a) * m) % productOfAllCharacteristics`, such that the result is positive.
   */
  element_type add_and_multiply(element_type e, const element_type& a, const element_type& m) const {
    mpz_add(e.get_mpz_t(), e.get_mpz_t(), a.get_mpz_t());
    mpz_mul(e.get_mpz_t(), e.get_mpz_t(), m.get_mpz_t());
    return _get_value(e);
  }

  /**
   * @brief Returns true if the two given elements are equal in the field, false otherwise.
   * 
   * @param e1 First element to compare.
   * @param e2 Second element to compare.
   * @return true If `e1 % productOfAllCharacteristics == e2 % productOfAllCharacteristics`.
   * @return false Otherwise.
   */
  bool are_equal(element_type e1, element_type e2) const { return _get_value(e1) == _get_value(e2); }

  /**
   * @brief Returns the inverse of the given element in the sense of @cite boissonnat:hal-00922572 with respect
   * to the product of all characteristics.
   * 
   * @param e Element to get the inverse from.
   * @return Inverse in the current field.
   */
  element_type get_inverse(const element_type& e) const {
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
  std::pair<element_type, characteristic_type> get_partial_inverse(
      const element_type& e, const characteristic_type& productOfCharacteristics) const {
    characteristic_type QR;
    mpz_gcd(QR.get_mpz_t(), e.get_mpz_t(), productOfCharacteristics.get_mpz_t());  // QR <- gcd(x,QS)

    if (QR == productOfCharacteristics) return {0, get_multiplicative_identity()};  // partial inverse is 0

    characteristic_type QT = productOfCharacteristics / QR;

    characteristic_type inv_qt;
    mpz_invert(inv_qt.get_mpz_t(), e.get_mpz_t(), QT.get_mpz_t());

    auto res = get_partial_multiplicative_identity(QT);
    mpz_mul(res.get_mpz_t(), res.get_mpz_t(), inv_qt.get_mpz_t());

    return {_get_value(res), QT};
  }

  /**
   * @brief Returns the additive identity of a field.
   * 
   * @return The additive identity of a field.
   */
  static element_type get_additive_identity() { return 0; }
  /**
   * @brief Returns the multiplicative identity of a field.
   * 
   * @return The multiplicative identity of a field.
   */
  static element_type get_multiplicative_identity() { return 1; }
  // static element_type get_multiplicative_identity(){ return multiplicativeID_; }
  /**
   * @brief Returns the partial multiplicative identity of the multi-field from the given product.
   * See @cite boissonnat:hal-00922572 for more details.
   * 
   * @param productOfCharacteristics Product of the different characteristics to take into account in the multi-field.
   * @return The partial multiplicative identity of the multi-field.
   */
  element_type get_partial_multiplicative_identity(const characteristic_type& productOfCharacteristics) const {
    if (productOfCharacteristics == 0) {
      return get_multiplicative_identity();
    }
    element_type multIdentity(0);
    for (unsigned int idx = 0; idx < primes_.size(); ++idx) {
      if ((productOfCharacteristics % primes_[idx]) == 0) {
        mpz_add(multIdentity.get_mpz_t(), multIdentity.get_mpz_t(), partials_[idx].get_mpz_t());
      }
    }
    return _get_value(multIdentity);
  }

  // static constexpr bool handles_only_z2() { return false; }

  /**
   * @brief Assign operator.
   */
  Multi_field_operators& operator=(Multi_field_operators other) {
    primes_.swap(other.primes_);
    productOfAllCharacteristics_ = other.productOfAllCharacteristics_;
    partials_.swap(other.partials_);

    return *this;
  }
  /**
   * @brief Swap operator.
   */
  friend void swap(Multi_field_operators& f1, Multi_field_operators& f2) {
    f1.primes_.swap(f2.primes_);
    std::swap(f1.productOfAllCharacteristics_, f2.productOfAllCharacteristics_);
    f1.partials_.swap(f2.partials_);
  }

 private:
  std::vector<unsigned int> primes_;                /**< All characteristics. */
  characteristic_type productOfAllCharacteristics_; /**< Product of all characteristics. */
  std::vector<characteristic_type> partials_;       /**< Partial products of the characteristics. */
  // static const element_type multiplicativeID_;

  static constexpr bool _is_prime(const int p) {
    if (p <= 1) return false;
    if (p <= 3) return true;
    if (p % 2 == 0 || p % 3 == 0) return false;

    for (long i = 5; i * i <= p; i = i + 6)
      if (p % i == 0 || p % (i + 2) == 0) return false;

    return true;
  }

  element_type& _get_value(element_type& e) const {
    if (e >= productOfAllCharacteristics_ || e < -productOfAllCharacteristics_)
      mpz_mod(e.get_mpz_t(), e.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
    if (e < 0) mpz_add(e.get_mpz_t(), e.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
    return e;
  }
};

}  // namespace persistence_fields
}  // namespace Gudhi

#endif  // MATRIX_FIELD_MULTI_OPERATORS_H_

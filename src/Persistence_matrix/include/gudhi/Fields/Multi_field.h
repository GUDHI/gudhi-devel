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
 * @file Multi_field.h
 * @author Hannah Schreiber, Clément Maria
 * @brief Contains the @ref Gudhi::persistence_fields::Multi_field_element class.
 */

#ifndef MATRIX_FIELD_MULTI_H_
#define MATRIX_FIELD_MULTI_H_

#include <utility>
#include <vector>
#include <gmpxx.h>
#include <stdexcept>

namespace Gudhi {
namespace persistence_fields {

/**
 * @class Multi_field_element Multi_field.h gudhi/Fields/Multi_field.h
 * @ingroup persistence_fields
 *
 * @brief Class representing an element of a multi-field.
 * The characteristics will corresponds to all prime numbers in the interval given as template.
 *
 * @tparam minimum Interval closed lower bound.
 * @tparam maximum Interval closed upper bound.
 */
template <unsigned int minimum, unsigned int maximum>
class Multi_field_element {
 public:
  using Element = mpz_class;           /**< Type for the elements in the field. */
  using Characteristic = Element; /**< Type for the field characteristic. */

  /**
   * @brief Default constructor. Sets the element to 0.
   */
  Multi_field_element();
  /**
   * @brief Constructor setting the element to the given value.
   *
   * @param element Value of the element.
   */
  Multi_field_element(const Element& element);
  /**
   * @brief Copy constructor.
   *
   * @param toCopy Element to copy.
   */
  Multi_field_element(const Multi_field_element& toCopy);
  /**
   * @brief Move constructor.
   *
   * @param toMove Element to move.
   */
  Multi_field_element(Multi_field_element&& toMove) noexcept;

  /**
   * @brief operator+=
   */
  friend void operator+=(Multi_field_element& f1, Multi_field_element const& f2) {
    f1.element_ += f2.element_;
    mpz_mod(f1.element_.get_mpz_t(), f1.element_.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
  }
  /**
   * @brief operator+
   */
  friend Multi_field_element operator+(Multi_field_element f1, Multi_field_element const& f2) {
    f1 += f2;
    return f1;
  }
  /**
   * @brief operator+=
   */
  friend void operator+=(Multi_field_element& f, const Element& v) {
    f.element_ += v;
    mpz_mod(f.element_.get_mpz_t(), f.element_.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
  }
  /**
   * @brief operator+
   */
  friend Multi_field_element operator+(Multi_field_element f, const Element& v) {
    f += v;
    return f;
  }
  /**
   * @brief operator+
   */
  friend Element operator+(Element v, Multi_field_element const& f) {
    v += f.element_;
    mpz_mod(v.get_mpz_t(), v.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
    return v;
  }

  /**
   * @brief operator-=
   */
  friend void operator-=(Multi_field_element& f1, Multi_field_element const& f2) {
    f1.element_ -= f2.element_;
    mpz_mod(f1.element_.get_mpz_t(), f1.element_.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
  }
  /**
   * @brief operator-
   */
  friend Multi_field_element operator-(Multi_field_element f1, Multi_field_element const& f2) {
    f1 -= f2;
    return f1;
  }
  /**
   * @brief operator-=
   */
  friend void operator-=(Multi_field_element& f, const Element& v) {
    f.element_ -= v;
    mpz_mod(f.element_.get_mpz_t(), f.element_.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
  }
  /**
   * @brief operator-
   */
  friend Multi_field_element operator-(Multi_field_element f, const Element& v) {
    f -= v;
    return f;
  }
  /**
   * @brief operator-
   */
  friend Element operator-(Element v, Multi_field_element const& f) {
    // Element e(v);
    if (v >= productOfAllCharacteristics_)
      mpz_mod(v.get_mpz_t(), v.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
    if (f.element_ > v) v += productOfAllCharacteristics_;
    v -= f.element_;
    return v;
  }

  /**
   * @brief operator*=
   */
  friend void operator*=(Multi_field_element& f1, Multi_field_element const& f2) {
    f1.element_ *= f2.element_;
    mpz_mod(f1.element_.get_mpz_t(), f1.element_.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
  }
  /**
   * @brief operator*
   */
  friend Multi_field_element operator*(Multi_field_element f1, Multi_field_element const& f2) {
    f1 *= f2;
    return f1;
  }
  /**
   * @brief operator*=
   */
  friend void operator*=(Multi_field_element& f, const Element& v) {
    f.element_ *= v;
    mpz_mod(f.element_.get_mpz_t(), f.element_.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
  }
  /**
   * @brief operator*
   */
  friend Multi_field_element operator*(Multi_field_element f, const Element& v) {
    f *= v;
    return f;
  }
  /**
   * @brief operator*
   */
  friend Element operator*(Element v, Multi_field_element const& f) {
    v *= f.element_;
    mpz_mod(v.get_mpz_t(), v.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
    return v;
  }

  /**
   * @brief operator==
   */
  friend bool operator==(const Multi_field_element& f1, const Multi_field_element& f2) {
    return f1.element_ == f2.element_;
  }
  /**
   * @brief operator==
   */
  friend bool operator==(const Element& v, const Multi_field_element& f) {
    if (v < productOfAllCharacteristics_) return v == f.element_;
    Element e(v);
    mpz_mod(e.get_mpz_t(), e.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
    return e == f.element_;
  }
  /**
   * @brief operator==
   */
  friend bool operator==(const Multi_field_element& f, const Element& v) {
    if (v < productOfAllCharacteristics_) return v == f.element_;
    Element e(v);
    mpz_mod(e.get_mpz_t(), e.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
    return e == f.element_;
  }
  /**
   * @brief operator!=
   */
  friend bool operator!=(const Multi_field_element& f1, const Multi_field_element& f2) { return !(f1 == f2); }
  /**
   * @brief operator!=
   */
  friend bool operator!=(const Element& v, const Multi_field_element& f) { return !(v == f); }
  /**
   * @brief operator!=
   */
  friend bool operator!=(const Multi_field_element& f, const Element& v) { return !(v == f); }

  /**
   * @brief Assign operator.
   */
  Multi_field_element& operator=(Multi_field_element other);
  /**
   * @brief Assign operator.
   */
  Multi_field_element& operator=(const Element& value);
  /**
   * @brief Swap operator.
   */
  friend void swap(Multi_field_element& f1, Multi_field_element& f2) { std::swap(f1.element_, f2.element_); }

  /**
   * @brief Casts the element into an unsigned int.
   */
  operator unsigned int() const;
  /**
   * @brief Casts the element into an mpz_class.
   */
  operator mpz_class() const;

  /**
   * @brief Returns the inverse of the element in the multi-field, see @cite boissonnat:hal-00922572.
   *
   * @return The inverse.
   */
  Multi_field_element get_inverse() const;
  /**
   * @brief Returns the inverse of the element with respect to a sub-product of the characteristics in the multi-field,
   * see @cite boissonnat:hal-00922572.
   *
   * @param productOfCharacteristics Sub-product of the characteristics.
   * @return Pair of the inverse and the characteristic the inverse corresponds to.
   */
  std::pair<Multi_field_element, Characteristic> get_partial_inverse(
      const Characteristic& productOfCharacteristics) const;

  /**
   * @brief Returns the additive identity of a field.
   *
   * @return The additive identity of a field.
   */
  static Multi_field_element get_additive_identity();
  /**
   * @brief Returns the multiplicative identity of a field.
   *
   * @return The multiplicative identity of a field.
   */
  static Multi_field_element get_multiplicative_identity();
  /**
   * @brief Returns the partial multiplicative identity of the multi-field from the given product.
   * See @cite boissonnat:hal-00922572 for more details.
   *
   * @param productOfCharacteristics Product of the different characteristics to take into account in the multi-field.
   * @return The partial multiplicative identity of the multi-field.
   */
  static Multi_field_element get_partial_multiplicative_identity(const Characteristic& productOfCharacteristics);
  /**
   * @brief Returns the product of all characteristics.
   *
   * @return The product of all characteristics.
   */
  static Characteristic get_characteristic();

  /**
   * @brief Returns the value of the element.
   *
   * @return Value of the element.
   */
  Element get_value() const;

  // static constexpr bool handles_only_z2() { return false; }

 private:
  Element element_;
  static inline const std::vector<unsigned int> primes_ = []() {
    std::vector<unsigned int> res;

    unsigned int curr_prime = minimum;
    mpz_t tmp_prime;
    mpz_init_set_ui(tmp_prime, minimum);
    // test if min_prime is prime
    int is_prime = mpz_probab_prime_p(tmp_prime, 25);  // probabilistic primality test

    if (is_prime == 0) {  // min_prime is composite
      mpz_nextprime(tmp_prime, tmp_prime);
      curr_prime = mpz_get_ui(tmp_prime);
    }

    while (curr_prime <= maximum) {
      res.push_back(curr_prime);
      mpz_nextprime(tmp_prime, tmp_prime);
      curr_prime = mpz_get_ui(tmp_prime);
    }
    mpz_clear(tmp_prime);

    return res;
  }();
  static inline const Characteristic productOfAllCharacteristics_ = []() {
    Characteristic res = 1;
    for (const auto p : primes_) {
      res *= p;
    }

    return res;
  }();
  static inline const std::vector<Characteristic> partials_ = []() {
    std::vector<Characteristic> res;

    if (productOfAllCharacteristics_ == 1) return res;

    for (unsigned int i = 0; i < primes_.size(); ++i) {
      unsigned int p = primes_[i];
      res.push_back(productOfAllCharacteristics_ / p);
      mpz_powm_ui(res.back().get_mpz_t(), res.back().get_mpz_t(), p - 1, productOfAllCharacteristics_.get_mpz_t());
    }

    return res;
  }();
  // If I understood the paper well, multiplicativeID_ always equals to 1. But in Clement's code,
  // multiplicativeID_ is computed (see commented lambda function below). TODO: verify with Clement.
  static inline const Element multiplicativeID_ = 1; /*[](){
           mpz_class res = 0;
           for (unsigned int i = 0; i < partials_.size(); ++i){
                   res = (res + partials_[i]) % productOfAllCharacteristics_;
           }

           return res;
   }();*/

  static constexpr bool _is_prime(const int p);
};

template <unsigned int minimum, unsigned int maximum>
inline Multi_field_element<minimum, maximum>::Multi_field_element() : element_(0) {
  static_assert(maximum >= 2, "Characteristics have to be positive.");
  static_assert(minimum <= maximum, "The given interval is not valid.");
  static_assert(minimum != maximum || _is_prime(minimum), "The given interval does not contain a prime number.");

  if (productOfAllCharacteristics_ == 1)
    throw std::runtime_error("The given interval does not contain a prime number.");
}

template <unsigned int minimum, unsigned int maximum>
inline Multi_field_element<minimum, maximum>::Multi_field_element(const Element& element) : element_(element) {
  static_assert(maximum >= 2, "Characteristics has to be positive.");
  static_assert(minimum <= maximum, "The given interval is not valid.");
  static_assert(minimum != maximum || _is_prime(minimum), "The given interval does not contain a prime number.");

  if (productOfAllCharacteristics_ == 1)
    throw std::runtime_error("The given interval does not contain a prime number.");

  mpz_mod(element_.get_mpz_t(), element_.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
}

template <unsigned int minimum, unsigned int maximum>
inline Multi_field_element<minimum, maximum>::Multi_field_element(const Multi_field_element<minimum, maximum>& toCopy)
    : element_(toCopy.element_) {}

template <unsigned int minimum, unsigned int maximum>
inline Multi_field_element<minimum, maximum>::Multi_field_element(
    Multi_field_element<minimum, maximum>&& toMove) noexcept
    : element_(std::move(toMove.element_)) {}

template <unsigned int minimum, unsigned int maximum>
inline Multi_field_element<minimum, maximum>& Multi_field_element<minimum, maximum>::operator=(
    Multi_field_element other) {
  std::swap(element_, other.element_);
  return *this;
}

template <unsigned int minimum, unsigned int maximum>
inline Multi_field_element<minimum, maximum>& Multi_field_element<minimum, maximum>::operator=(
    const Element& value) {
  mpz_mod(element_.get_mpz_t(), value.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
  return *this;
}

template <unsigned int minimum, unsigned int maximum>
inline Multi_field_element<minimum, maximum>::operator unsigned int() const {
  return element_.get_ui();
}

template <unsigned int minimum, unsigned int maximum>
inline Multi_field_element<minimum, maximum>::operator mpz_class() const {
  return element_;
}

template <unsigned int minimum, unsigned int maximum>
inline Multi_field_element<minimum, maximum> Multi_field_element<minimum, maximum>::get_inverse() const {
  return get_partial_inverse(productOfAllCharacteristics_).first;
}

template <unsigned int minimum, unsigned int maximum>
inline std::pair<Multi_field_element<minimum, maximum>,
                 typename Multi_field_element<minimum, maximum>::Characteristic>
Multi_field_element<minimum, maximum>::get_partial_inverse(const Characteristic& productOfCharacteristics) const {
  Characteristic QR;
  mpz_gcd(QR.get_mpz_t(), element_.get_mpz_t(), productOfCharacteristics.get_mpz_t());  // QR <- gcd(x,QS)

  if (QR == productOfCharacteristics) return {Multi_field_element(), multiplicativeID_};  // partial inverse is 0

  Characteristic QT = productOfCharacteristics / QR;

  Characteristic inv_qt;
  mpz_invert(inv_qt.get_mpz_t(), element_.get_mpz_t(), QT.get_mpz_t());

  auto res = get_partial_multiplicative_identity(QT);
  res *= inv_qt;

  return {res, QT};
}

template <unsigned int minimum, unsigned int maximum>
inline Multi_field_element<minimum, maximum> Multi_field_element<minimum, maximum>::get_additive_identity() {
  return Multi_field_element<minimum, maximum>();
}

template <unsigned int minimum, unsigned int maximum>
inline Multi_field_element<minimum, maximum> Multi_field_element<minimum, maximum>::get_multiplicative_identity() {
  return Multi_field_element<minimum, maximum>(multiplicativeID_);
}

template <unsigned int minimum, unsigned int maximum>
inline Multi_field_element<minimum, maximum> Multi_field_element<minimum, maximum>::get_partial_multiplicative_identity(
    const Characteristic& productOfCharacteristics) {
  if (productOfCharacteristics == 0) {
    return Multi_field_element<minimum, maximum>(multiplicativeID_);
  }
  Multi_field_element<minimum, maximum> mult;
  for (unsigned int idx = 0; idx < primes_.size(); ++idx) {
    if ((productOfCharacteristics % primes_[idx]) == 0) {
      mult += partials_[idx];
    }
  }
  return mult;
}

template <unsigned int minimum, unsigned int maximum>
inline typename Multi_field_element<minimum, maximum>::Characteristic
Multi_field_element<minimum, maximum>::get_characteristic() {
  return productOfAllCharacteristics_;
}

template <unsigned int minimum, unsigned int maximum>
inline typename Multi_field_element<minimum, maximum>::Element Multi_field_element<minimum, maximum>::get_value()
    const {
  return element_;
}

template <unsigned int minimum, unsigned int maximum>
inline constexpr bool Multi_field_element<minimum, maximum>::_is_prime(const int p) {
  if (p <= 1) return false;
  if (p <= 3) return true;
  if (p % 2 == 0 || p % 3 == 0) return false;

  for (long i = 5; i * i <= p; i = i + 6)
    if (p % i == 0 || p % (i + 2) == 0) return false;

  return true;
}

}  // namespace persistence_fields
}  // namespace Gudhi

#endif  // MATRIX_FIELD_MULTI_H_

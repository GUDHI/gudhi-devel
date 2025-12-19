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
 * @file Multi_field_shared.h
 * @author Hannah Schreiber, Clément Maria
 * @brief Contains the @ref Gudhi::persistence_fields::Shared_multi_field_element class.
 */

#ifndef MATRIX_FIELD_MULTI_SHARED_H_
#define MATRIX_FIELD_MULTI_SHARED_H_

#include <utility>
#include <vector>
#include <gmpxx.h>
#include <stdexcept>

#include <gudhi/Debug_utils.h>

namespace Gudhi {
namespace persistence_fields {

/**
 * @class Shared_multi_field_element Multi_field_shared.h gudhi/Fields/Multi_field_shared.h
 * @ingroup persistence_fields
 *
 * @brief Class representing an element of a multi-field. If each instantiation of the class can represent another
 * element, they all share the same characteristics. That is if the characteristics are set for one, they will be
 * set for all the others. The characteristics can be set before instantiating the elements with the static
 * @ref Shared_multi_field_element::initialize method.
 */
class Shared_multi_field_element
{
 public:
  using Element = mpz_class;      /**< Type for the elements in the field. */
  using Characteristic = Element; /**< Type for the field characteristic. */

  /**
   * @brief Default constructor. Sets the element to 0.
   */
  Shared_multi_field_element() : element_(0) {}

  /**
   * @brief Constructor setting the element to the given value.
   *
   * @param element Value of the element.
   */
  Shared_multi_field_element(Element element) : element_(std::move(element))
  {
    mpz_mod(element_.get_mpz_t(), element_.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
  }

  /**
   * @brief Initialize the multi-field to the characteristics (primes) contained in the given interval.
   * Should be called first before constructing the field elements.
   *
   * @param minimum Lowest value in the interval.
   * @param maximum Highest value in the interval.
   */
  static void initialize(unsigned int minimum, unsigned int maximum)
  {
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
    while (curr_prime <= maximum) {
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
    //	for (unsigned int i = 0; i < partials_.size(); ++i){
    //		multiplicativeID_ = (multiplicativeID_ + partials_[i]) % productOfAllCharacteristics_;
    //	}
  }

  /**
   * @brief operator+=
   */
  friend void operator+=(Shared_multi_field_element& f1, Shared_multi_field_element const& f2)
  {
    f1.element_ += f2.element_;
    mpz_mod(f1.element_.get_mpz_t(), f1.element_.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
  }

  /**
   * @brief operator+
   */
  friend Shared_multi_field_element operator+(Shared_multi_field_element f1, Shared_multi_field_element const& f2)
  {
    f1 += f2;
    return f1;
  }

  /**
   * @brief operator+=
   */
  friend void operator+=(Shared_multi_field_element& f, const Element& v)
  {
    f.element_ += v;
    mpz_mod(f.element_.get_mpz_t(), f.element_.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
  }

  /**
   * @brief operator+
   */
  friend Shared_multi_field_element operator+(Shared_multi_field_element f, const Element& v)
  {
    f += v;
    return f;
  }

  /**
   * @brief operator+
   */
  friend Element operator+(Element v, Shared_multi_field_element const& f)
  {
    v += f.element_;
    mpz_mod(v.get_mpz_t(), v.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
    return v;
  }

  /**
   * @brief operator-=
   */
  friend void operator-=(Shared_multi_field_element& f1, Shared_multi_field_element const& f2)
  {
    f1.element_ -= f2.element_;
    mpz_mod(f1.element_.get_mpz_t(), f1.element_.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
  }

  /**
   * @brief operator-
   */
  friend Shared_multi_field_element operator-(Shared_multi_field_element f1, Shared_multi_field_element const& f2)
  {
    f1 -= f2;
    return f1;
  }

  /**
   * @brief operator-=
   */
  friend void operator-=(Shared_multi_field_element& f, const Element& v)
  {
    f.element_ -= v;
    mpz_mod(f.element_.get_mpz_t(), f.element_.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
  }

  /**
   * @brief operator-
   */
  friend Shared_multi_field_element operator-(Shared_multi_field_element f, const Element& v)
  {
    f -= v;
    return f;
  }

  /**
   * @brief operator-
   */
  friend Element operator-(Element v, Shared_multi_field_element const& f)
  {
    if (v >= productOfAllCharacteristics_)
      mpz_mod(v.get_mpz_t(), v.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
    if (f.element_ > v) v += productOfAllCharacteristics_;
    v -= f.element_;
    return v;
  }

  /**
   * @brief operator*=
   */
  friend void operator*=(Shared_multi_field_element& f1, Shared_multi_field_element const& f2)
  {
    f1.element_ *= f2.element_;
    mpz_mod(f1.element_.get_mpz_t(), f1.element_.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
  }

  /**
   * @brief operator*
   */
  friend Shared_multi_field_element operator*(Shared_multi_field_element f1, Shared_multi_field_element const& f2)
  {
    f1 *= f2;
    return f1;
  }

  /**
   * @brief operator*=
   */
  friend void operator*=(Shared_multi_field_element& f, const Element& v)
  {
    f.element_ *= v;
    mpz_mod(f.element_.get_mpz_t(), f.element_.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
  }

  /**
   * @brief operator*
   */
  friend Shared_multi_field_element operator*(Shared_multi_field_element f, const Element& v)
  {
    f *= v;
    return f;
  }

  /**
   * @brief operator*
   */
  friend Element operator*(Element v, Shared_multi_field_element const& f)
  {
    v *= f.element_;
    mpz_mod(v.get_mpz_t(), v.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
    return v;
  }

  /**
   * @brief operator==
   */
  friend bool operator==(const Shared_multi_field_element& f1, const Shared_multi_field_element& f2)
  {
    return f1.element_ == f2.element_;
  }

  /**
   * @brief operator==
   */
  friend bool operator==(const Element& v, const Shared_multi_field_element& f)
  {
    if (v < productOfAllCharacteristics_) return v == f.element_;
    Element e(v);
    mpz_mod(e.get_mpz_t(), e.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
    return e == f.element_;
  }

  /**
   * @brief operator==
   */
  friend bool operator==(const Shared_multi_field_element& f, const Element& v)
  {
    if (v < productOfAllCharacteristics_) return v == f.element_;
    Element e(v);
    mpz_mod(e.get_mpz_t(), e.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
    return e == f.element_;
  }

  /**
   * @brief operator!=
   */
  friend bool operator!=(const Shared_multi_field_element& f1, const Shared_multi_field_element& f2)
  {
    return !(f1 == f2);
  }

  /**
   * @brief operator!=
   */
  friend bool operator!=(const Element& v, const Shared_multi_field_element& f) { return !(v == f); }

  /**
   * @brief operator!=
   */
  friend bool operator!=(const Shared_multi_field_element& f, const Element& v) { return !(v == f); }

  // /**
  //  * @brief Assign operator.
  //  */
  // Shared_multi_field_element& operator=(const Element& value)
  // {
  //   mpz_mod(element_.get_mpz_t(), value.get_mpz_t(), productOfAllCharacteristics_.get_mpz_t());
  //   return *this;
  // }

  /**
   * @brief Swap operator.
   */
  friend void swap(Shared_multi_field_element& f1, Shared_multi_field_element& f2) noexcept
  {
    std::swap(f1.element_, f2.element_);
  }

  /**
   * @brief Casts the element into an unsigned int.
   */
  operator unsigned int() const { return element_.get_ui(); }

  /**
   * @brief Casts the element into a mpz_class.
   */
  explicit operator mpz_class() const { return element_; }

  /**
   * @brief Returns the inverse of the element in the multi-field, see @cite boissonnat:hal-00922572.
   *
   * @return The inverse.
   */
  [[nodiscard]] Shared_multi_field_element get_inverse() const
  {
    return get_partial_inverse(productOfAllCharacteristics_).first;
  }

  /**
   * @brief Returns the inverse of the element with respect to a sub-product of the characteristics in the multi-field,
   * see @cite boissonnat:hal-00922572.
   *
   * @param productOfCharacteristics Sub-product of the characteristics.
   * @return Pair of the inverse and the characteristic the inverse corresponds to.
   */
  [[nodiscard]] std::pair<Shared_multi_field_element, Characteristic> get_partial_inverse(
      const Characteristic& productOfCharacteristics) const
  {
    GUDHI_CHECK(productOfCharacteristics >= 0 && productOfCharacteristics <= productOfAllCharacteristics_,
                "The given product is not the product of a subset of the current Multi-field characteristics.");

    Element QR;
    mpz_gcd(QR.get_mpz_t(), element_.get_mpz_t(), productOfCharacteristics.get_mpz_t());  // QR <- gcd(x,QS)

    if (QR == productOfCharacteristics)
      return {Shared_multi_field_element(), multiplicativeID_};  // partial inverse is 0

    Element QT = productOfCharacteristics / QR;

    Element inv_qt;
    mpz_invert(inv_qt.get_mpz_t(), element_.get_mpz_t(), QT.get_mpz_t());

    auto res = get_partial_multiplicative_identity(QT);
    res *= inv_qt;

    return {res, QT};
  }

  /**
   * @brief Returns the additive identity of a field.
   *
   * @return The additive identity of a field.
   */
  static Shared_multi_field_element get_additive_identity() { return {}; }

  /**
   * @brief Returns the multiplicative identity of a field.
   *
   * @return The multiplicative identity of a field.
   */
  static Shared_multi_field_element get_multiplicative_identity() { return {multiplicativeID_}; }

  /**
   * @brief Returns the partial multiplicative identity of the multi-field from the given product.
   * See @cite boissonnat:hal-00922572 for more details.
   *
   * @param productOfCharacteristics Product of the different characteristics to take into account in the multi-field.
   * @return The partial multiplicative identity of the multi-field.
   */
  static Shared_multi_field_element get_partial_multiplicative_identity(const Characteristic& productOfCharacteristics)
  {
    GUDHI_CHECK(productOfCharacteristics >= 0 && productOfCharacteristics <= productOfAllCharacteristics_,
                "The given product is not the product of a subset of the current Multi-field characteristics.");

    if (productOfCharacteristics == 0 || productOfCharacteristics == productOfAllCharacteristics_) {
      return get_multiplicative_identity();
    }
    Shared_multi_field_element mult;
    for (unsigned int idx = 0; idx < primes_.size(); ++idx) {
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
  [[nodiscard]] Element get_value() const { return element_; }

 private:
  Element element_;                                              /**< Element. */
  static inline std::vector<unsigned int> primes_;               /**< All characteristics. */
  static inline Characteristic productOfAllCharacteristics_ = 0; /**< Product of all characteristics. */
  static inline std::vector<Element> partials_;                  /**< Partial products of the characteristics. */
  static inline const Element multiplicativeID_ = 1;             /**< Multiplicative identity. */

  static constexpr bool _is_prime(const unsigned int p)
  {
    if (p <= 1) return false;
    if (p <= 3) return true;
    if (p % 2 == 0 || p % 3 == 0) return false;

    for (long i = 5; i * i <= p; i = i + 6)
      if (p % i == 0 || p % (i + 2) == 0) return false;

    return true;
  }
};

}  // namespace persistence_fields
}  // namespace Gudhi

#endif  // MATRIX_FIELD_MULTI_SHARED_H_

/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PERSISTENT_COHOMOLOGY_FIELD_ZP_H_
#define PERSISTENT_COHOMOLOGY_FIELD_ZP_H_

#include <utility>
#include <vector>
#include <stdexcept>

namespace Gudhi {

namespace persistent_cohomology {

/** \brief Structure representing the coefficient field \f$\mathbb{Z}/p\mathbb{Z}\f$
 *
 * \implements CoefficientField
 * \ingroup persistent_cohomology
 */
class Field_Zp {
 public:
  typedef int Element;

  Field_Zp()
      : Prime(0),
        inverse_() {
  }

  void init(int charac) {
    Prime = charac;

    // Check that the provided prime is less than the maximum allowed as int, calculation below, and 'plus_times_equal' function : 46337 ; i.e (max_prime-1)*max_prime <= INT_MAX
    if(Prime > 46337)
        throw std::invalid_argument("Maximum homology_coeff_field allowed value is 46337");

    // Check for primality
    if (Prime <= 1)
        throw std::invalid_argument("homology_coeff_field must be a prime number");

    inverse_.clear();
    inverse_.reserve(charac);
    inverse_.push_back(0);
    for (int i = 1; i < Prime; ++i) {
      int inv = 1;
      int mult = inv * i;
      while ( (mult % Prime) != 1) {
        ++inv;
        if(mult == Prime)
            throw std::invalid_argument("homology_coeff_field must be a prime number");
        mult = inv * i;
      }
      inverse_.push_back(inv);
    }
  }

  /** Set x <- x + w * y*/
  Element plus_times_equal(const Element& x, const Element& y, const Element& w) {
    assert(Prime > 0);  // division by zero + non negative values
    Element result = (x + w * y) % Prime;
    if (result < 0)
      result += Prime;
    return result;
  }

// operator= defined on Element

  /** Returns y * w */
  Element times(const Element& y, const Element& w) {
    return plus_times_equal(0, y, (Element)w);
  }

  Element plus_equal(const Element& x, const Element& y) {
    return plus_times_equal(x, y, (Element)1);
  }

  /** \brief Returns the additive idendity \f$0_{\Bbbk}\f$ of the field.*/
  Element additive_identity() const {
    return 0;
  }
  /** \brief Returns the multiplicative identity \f$1_{\Bbbk}\f$ of the field.*/
  Element multiplicative_identity(Element = 0) const {
    return 1;
  }
  /** Returns the inverse in the field. Modifies P. ??? */
  std::pair<Element, Element> inverse(Element x, Element P) {
    return std::pair<Element, Element>(inverse_[x], P);
  }  // <------ return the product of field characteristic for which x is invertible

  /** Returns -x * y.*/
  Element times_minus(Element x, Element y) {
    assert(Prime > 0);  // division by zero + non negative values
    Element out = (-x * y) % Prime;
    return (out < 0) ? out + Prime : out;
  }

  /** \brief Returns the characteristic \f$p\f$ of the field.*/
  int characteristic() const {
    return Prime;
  }

 private:
  int Prime;
  /** Property map Element -> Element, which associate to an element its inverse in the field.*/
  std::vector<Element> inverse_;
};

}  // namespace persistent_cohomology

}  // namespace Gudhi

#endif  // PERSISTENT_COHOMOLOGY_FIELD_ZP_H_

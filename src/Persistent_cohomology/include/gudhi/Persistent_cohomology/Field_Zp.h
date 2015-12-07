/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014  INRIA Sophia Antipolis-Méditerranée (France)
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SRC_PERSISTENT_COHOMOLOGY_INCLUDE_GUDHI_PERSISTENT_COHOMOLOGY_FIELD_ZP_H_
#define SRC_PERSISTENT_COHOMOLOGY_INCLUDE_GUDHI_PERSISTENT_COHOMOLOGY_FIELD_ZP_H_

#include <utility>
#include <vector>

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
        inverse_(),
        mult_id_all(1),
        add_id_all(0) {
  }

  void init(int charac) {
    assert(charac > 0);  // division by zero + non negative values
    Prime = charac;
    inverse_.clear();
    inverse_.reserve(charac);
    inverse_.push_back(0);
    for (int i = 1; i < Prime; ++i) {
      int inv = 1;
      while (((inv * i) % Prime) != 1)
        ++inv;
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
  const Element& additive_identity() const {
    return add_id_all;
  }
  /** \brief Returns the multiplicative identity \f$1_{\Bbbk}\f$ of the field.*/
  const Element& multiplicative_identity(Element = 0) const {
    return mult_id_all;
  }
  /** Returns the inverse in the field. Modifies P.*/
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
  const int& characteristic() const {
    return Prime;
  }

 private:
  int Prime;
  /** Property map Element -> Element, which associate to an element its inverse in the field.*/
  std::vector<Element> inverse_;
  const Element mult_id_all;
  const Element add_id_all;
};

}  // namespace persistent_cohomology

}  // namespace Gudhi

#endif  // SRC_PERSISTENT_COHOMOLOGY_INCLUDE_GUDHI_PERSISTENT_COHOMOLOGY_FIELD_ZP_H_

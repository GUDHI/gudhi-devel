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

  void init(uint8_t charac) {
    assert(charac != 0);  // division by zero
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
  void plus_times_equal(Element & x, Element y, Element w) {
    assert(Prime != 0);  // division by zero
    x = (x + w * y) % Prime;
  }

// operator= defined on Element

  /** Returns y * w */
  Element times(Element y, int w) {
    assert(Prime != 0);  // division by zero
    Element res = (y * w) % Prime;
    if (res < 0)
      return res + Prime;
    else
      return res;
  }

  void clear_coefficient(Element x) {
  }

  void plus_equal(Element & x, Element y) {
    assert(Prime != 0);  // division by zero
    x = ((x + y) % Prime);
  }

  /** \brief Returns the additive idendity \f$0_{\Bbbk}\f$ of the field.*/
  const Element& additive_identity() const {
    return add_id_all;
  }
  /** \brief Returns the multiplicative identity \f$1_{\Bbbk}\f$ of the field.*/
  const Element& multiplicative_identity(Element P = 0) const {
    return mult_id_all;
  }
  /** Returns the inverse in the field. Modifies P.*/
  std::pair<Element, Element> inverse(Element x, Element P) {
    return std::pair<Element, Element>(inverse_[x], P);
  }  // <------ return the product of field characteristic for which x is invertible

  /** Returns -x * y.*/
  Element times_minus(Element x, Element y) {
    assert(Prime != 0);  // division by zero
    Element out = (-x * y) % Prime;
    return (out < 0) ? out + Prime : out;
  }

  /** \brief Returns the characteristic \f$p\f$ of the field.*/
  const uint8_t& characteristic() const {
    return Prime;
  }

 private:
  uint8_t Prime;
  /** Property map Element -> Element, which associate to an element its inverse in the field.*/
  std::vector<Element> inverse_;
  const Element mult_id_all;
  const Element add_id_all;
};

}  // namespace persistent_cohomology

}  // namespace Gudhi

#endif  // SRC_PERSISTENT_COHOMOLOGY_INCLUDE_GUDHI_PERSISTENT_COHOMOLOGY_FIELD_ZP_H_

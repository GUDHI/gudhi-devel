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

#ifndef GUDHI_FIELD_ZP_H
#define GUDHI_FIELD_ZP_H

namespace Gudhi{

/** \brief Structure representing the coefficient field \f$\mathbb{Z}/p\mathbb{Z}\f$
  *
  * \implements CoefficientField
  * \ingroup persistent_cohomology
  */
class Field_Zp {
public:
typedef int Element;

Field_Zp()
: Prime(-1)
, inverse_() {}

void init(int charac ) {
  assert(charac <= 32768);
  Prime = charac;
  inverse_.clear();
  inverse_.reserve(charac);
  inverse_.push_back(0);
  for(int i=1 ; i<Prime ; ++i)
  { 
    int inv = 1; 
    while(((inv * i) % Prime) != 1) ++inv;
    inverse_.push_back(inv); 
  }
}

/** Set x <- x + w * y*/
void plus_times_equal ( Element & x, Element y, Element w ) 
{ x = (x + w * y) % Prime; }

// operator= defined on Element

/** Returns y * w */
Element times ( Element y, int w ) { 
  Element res = (y * w) % Prime;
  if(res < 0) return res+Prime;
  else return res; 
}

void clear_coefficient(Element x) {}

void plus_equal(Element & x, Element y) { x = ((x+y)%Prime); }

/** \brief Returns the additive idendity \f$0_{\Bbbk}\f$ of the field.*/
Element additive_identity () { return 0; }
/** \brief Returns the multiplicative identity \f$1_{\Bbbk}\f$ of the field.*/
Element multiplicative_identity ( Element P = 0) { return 1; }
/** Returns the inverse in the field. Modifies P.*/
std::pair<Element,Element> inverse ( Element x
                                   , Element P ) 
{ return std::pair<Element,Element>(inverse_[x],P); 
}  // <------ return the product of field characteristic for which x is invertible

/** Returns -x * y.*/
Element times_minus ( Element x, Element y ) 
{ 
  Element out = (-x * y) % Prime;
  return (out < 0) ? out + Prime : out; 
}


bool is_one ( Element x ) { return x == 1; }
bool is_zero ( Element x ) { return x == 0; }

//bool is_null()

/** \brief Returns the characteristic \f$p\f$ of the field.*/
Element characteristic() { return Prime; }

private:
  Element Prime;
/** Property map Element -> Element, which associate to an element its inverse in the field.*/
  std::vector< Element > inverse_;
};

}  // namespace GUDHI

#endif // GUDHI_FIELD_ZP_H

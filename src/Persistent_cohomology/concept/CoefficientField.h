 /*    This file is part of the Gudhi Library. The Gudhi library 
  *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
  *    library for computational topology.
  *
  *    Author(s):       Clément Maria
  *
  *    Copyright (C) 2014 Inria
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

/** \brief Concept describing the requirements for a class to represent 
  * a field of coefficients to compute persistent homology.
  */
struct CoefficientField {

/** \brief Type of element of the field.
  *
  * Must be Assignable. */
  typedef unspecified Element;

/** Default constructible. */
  CoefficientField();
  
  void init(Element charac);
  void init(Element charac_min, Element charac_max);

/** Return the characteristic of the field. */
  Element characteristic();
/** Return the element 1 of the field. */
  Element multiplicative_identity();
/** Return the element 0 of the field. */
  Element additive_identity();

/** Assign: x <- x + y */
  void plus_equal(Element x, Element y);

/** */
//... inverse()

  };
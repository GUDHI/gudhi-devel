/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
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
/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/** \brief Concept describing the requirements for a class to compute 
  * persistent homology. */
struct PersistentHomology {

/** \brief Type of filtered cell complex on which persistent homology
  * is computed.
  *
  * Must be a model of concept FilteredComplex.
  */
  typedef unspecified Filtered_complex;

/** \brief Type of coefficients to be used for computing persistent 
  * homology.
  *
  * Must be a model of concept CoefficientField.
  */
  typedef unspecified Coefficient_field;

 };

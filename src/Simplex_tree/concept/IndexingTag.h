/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/** \brief Concept describing an indexing scheme (see FilteredComplex) 
  * for applying 
  * continuous maps to a cell complex, and compute its persistent
  * homology.
  *
  * Must be `Gudhi::linear_indexing_tag`. 
  */
struct IndexingTag {};

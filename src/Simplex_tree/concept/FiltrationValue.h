/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/** \brief Value type for a filtration function on a cell complex.
  *
  * A <EM>filtration</EM> of a cell complex (see FilteredComplex) is
  * a function \f$f:\mathbf{K} \rightarrow \mathbb{R}\f$ satisfying \f$f(\tau)\leq 
  * f(\sigma)\f$ whenever \f$\tau \subseteq \sigma\f$. Ordering the simplices 
  * by increasing filtration values (breaking ties so as a simplex appears after 
  * its subsimplices of same filtration value) provides an indexing scheme 
  * (see IndexingTag).
  */
  struct FiltrationValue {
    /** \brief Operator < is a StrictWeakOrdering. */
    bool operator<(FiltrationValue f1, FiltrationValue f2);
  };

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

  /** \brief For multiparameter filtrations, this methods pushes a filtration value 
   * to the first moment, for operator< such that f < this. For instance, for a one critical filtration, with
   *  - this = (1,2)
   *  - x = (2,1)
   *
   * after calling this method, x should be equal to (2,2).
   * This function is called when using, e.g. `make_filtration_non_decreasing`, as the filtration of a simplex
   * has to be greater than the filtration of any of its faces.
   * */ 
  void push_to_least_common_upper_bound(const FiltrationValue f);

};

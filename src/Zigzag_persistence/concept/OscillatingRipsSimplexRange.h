/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CONCEPT_ZZ_OSCI_RIPS_RANGE_H_
#define CONCEPT_ZZ_OSCI_RIPS_RANGE_H_

/** @file OscillatingRipsSimplexRange.h
 * @brief Contains @ref Gudhi::zigzag_persistence::OscillatingRipsSimplexRange concept.
 */

#include <gudhi/Zigzag_persistence/oscillating_rips_iterators.h>
#include <gudhi/Simplex_tree.h>

namespace Gudhi {
namespace zigzag_persistence {

/**
 * @brief Class giving access to a range over the simplices in an oscillating Rips filtration 
 * in order of the filtration.
 *
 * A simplex has to be represented by a tuple of three elements:
 * the first is the simplex handle of the simplex in the given complex,
 * the second is the filtration value of the corresponding arrow,
 * the third is the direction of the arrow, i.e., indicates if the simplex is inserted or removed.
 */
class OscillatingRipsSimplexRange {
 public:
  /**
   * @brief Returns a range over the simplices in an oscillating Rips filtration 
   * in order of the filtration.
   * 
   * @param edgeStartIterator Begin iterator of the edge range.
   * @param edgeEndIterator End iterator of the edge range.
   * @param complex Structure storing the complex at each step.
   * @param maxDimension Maximal dimension of the expansion.
   * @return A range with begin() and end() methods.
   */
  static auto get_iterator_range(Oscillating_rips_edge_range<double>::Oscillating_rips_edge_iterator& edgeStartIterator,
                                 Oscillating_rips_edge_range<double>::Oscillating_rips_edge_iterator& edgeEndIterator,
                                 StableFilteredComplex& complex,
                                 int maxDimension);

  /**
   * @brief Returns a range over the simplices in an oscillating Rips filtration 
   * in order of the filtration.
   * 
   * @param edgeStartIterator Begin iterator of the edge range.
   * @param edgeEndIterator End iterator of the edge range.
   * @param complex Structure storing the complex at each step.
   * @param maxDimension Maximal dimension of the expansion.
   * @return A range with begin() and end() methods.
   */
  static auto get_iterator_range(std::vector<Zigzag_edge<double> >::iterator& edgeStartIterator,
                                 std::vector<Zigzag_edge<double> >::iterator& edgeEndIterator,
                                 StableFilteredComplex& complex,
                                 int maxDimension);
};

}  // namespace zigzag_persistence
}  // namespace Gudhi

#endif  // CONCEPT_ZZ_OSCI_RIPS_RANGE_H_

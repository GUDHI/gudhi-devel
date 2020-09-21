/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PERMUTAHEDRAL_REPRESENTATION_COMBINATION_ITERATOR_H_
#define PERMUTAHEDRAL_REPRESENTATION_COMBINATION_ITERATOR_H_

#include <vector>
#include <boost/range/iterator_range.hpp>

namespace Gudhi {

namespace coxeter_triangulation {

typedef unsigned uint;

/** \brief Class that allows the user to generate combinations of 
 *   k elements in a set of n elements.
 *  Based on the algorithm by Mifsud.
*/
class Combination_iterator : public boost::iterator_facade< Combination_iterator,
							    std::vector<uint> const,
							    boost::forward_traversal_tag> {
  typedef std::vector<uint> value_t;
  
protected:
  friend class boost::iterator_core_access;
  
  bool equal(Combination_iterator const& other) const {
    return (is_end_ && other.is_end_);
  }

  value_t const& dereference() const {
    return value_;
  }

  void increment() {
    if (value_[0] == n_ - k_) {
      is_end_ = true;
      return;
    }
    uint j = k_ - 1;
    if (value_[j] < n_ - 1) {
      value_[j]++;
      return;
    }
    for (; j > 0; --j)
      if (value_[j-1] < n_ - k_ + j-1) {
	value_[j-1]++;
	for (uint s = j; s < k_; s++)
	  value_[s] = value_[j-1] + s - (j-1);
	return;
      }
  }

public:
  
  Combination_iterator(const uint& n, const uint& k)
    :
    value_(k),
    is_end_(n == 0),
    n_(n),
    k_(k)
  {
    for (uint i = 0; i < k; ++i)
      value_[i] = i;
  }

  // Used for the creating an end iterator
  Combination_iterator() : is_end_(true), n_(0), k_(0) {}

  void reinitialize() {
    if (n_ > 0) {
      is_end_ = false;
      for (uint i = 0; i < n_; ++i)
      value_[i] = i;
    }
  }
  
protected:
  value_t value_; // the dereference value
  bool is_end_;   // is true when the current permutation is the final one 

  uint n_;
  uint k_;
};

} // namespace coxeter_triangulation

} // namespace Gudhi


#endif

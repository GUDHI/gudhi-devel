/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PERMUTAHEDRAL_REPRESENTATION_PERMUTATION_ITERATOR_H_
#define PERMUTAHEDRAL_REPRESENTATION_PERMUTATION_ITERATOR_H_

#include <vector>
#include <boost/range/iterator_range.hpp>

namespace Gudhi {

namespace coxeter_triangulation {

typedef unsigned uint;

/** \brief Class that allows the user to generate permutations.
 *   Based on the optimization of the Heap's algorithm by Sedgewick.
*/
class Permutation_iterator : public boost::iterator_facade< Permutation_iterator,
							    std::vector<uint> const,
							    boost::forward_traversal_tag> {
  typedef std::vector<uint> value_t;
  
protected:
  friend class boost::iterator_core_access;

  bool equal(Permutation_iterator const& other) const {
    return (is_end_ && other.is_end_);
  }

  value_t const& dereference() const {
    return value_;
  }

  void swap_two_indices(std::size_t i, std::size_t j) {
    uint t = value_[i];
    value_[i] = value_[j];
    value_[j] = t;
  }
  
  void elementary_increment() {
    uint j = 0;
    while (d_[j] == j+1) {
      d_[j] = 0;
      ++j;
    }
    if (j == n_ - 1) {
      is_end_ = true;
      return;
    }
    uint k = j+1;
    uint x = (k%2 ? d_[j] : 0);
    swap_two_indices(k, x);
    ++d_[j];
  }

  void elementary_increment_optim_3() {
    if (ct_ != 0) {
      --ct_;
      swap_two_indices(1 + (ct_%2), 0);
    }
    else {
      ct_ = 5;
      uint j = 2;
      while (d_[j] == j+1) {
	d_[j] = 0;
	++j;
      }
      if (j == n_ - 1) {
	is_end_ = true;
	return;
      }
      uint k = j+1;
      uint x = (k%2 ? d_[j] : 0);
      swap_two_indices(k, x);
      ++d_[j];
    }
  }
    
  void increment() {
    if (optim_3_)
      elementary_increment_optim_3();
    else
      elementary_increment();      
  }

public:
  
  Permutation_iterator(const uint& n)
    :
    value_(n),
    is_end_(n == 0),
    optim_3_(n >= 3),
    n_(n),
    d_(n),
    ct_(5)
  {
    for (uint i = 0; i < n; ++i) {
      value_[i] = i;
      d_[i] = 0;
    }
    if (n > 0)
      d_[n-1] = -1;
  }

  // Used for the creating an end iterator
  Permutation_iterator() : is_end_(true), n_(0) {}

  void reinitialize() {
    if (n_ > 0)
      is_end_ = false;
  }
  
protected:
  value_t value_; // the dereference value
  bool is_end_;   // is true when the current permutation is the final one 
  bool optim_3_;  // true if n>=3. for n >= 3, the algorithm is optimized

  uint n_;
  std::vector<uint> d_; // mix radix digits with radix [2 3 4 ... n-1 (sentinel=-1)]
  uint ct_;             // counter with values in {0,...,5} used in the n>=3 optimization.
};

} // namespace coxeter_triangulation

} // namespace Gudhi

#endif

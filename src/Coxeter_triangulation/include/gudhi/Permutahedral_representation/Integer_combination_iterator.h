/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PERMUTAHEDRAL_REPRESENTATION_INTEGER_COMBINATION_ITERATOR_H_
#define PERMUTAHEDRAL_REPRESENTATION_INTEGER_COMBINATION_ITERATOR_H_

#include <vector>
#include <boost/range/iterator_range.hpp>

namespace Gudhi {

namespace coxeter_triangulation {

typedef unsigned uint;

/** \brief Class that allows the user to generate combinations of 
 *   k elements in a set of n elements.
 *  Based on the algorithm by Mifsud.
*/
class Integer_combination_iterator : public boost::iterator_facade< Integer_combination_iterator,
                                                                    std::vector<uint> const,
                                                                    boost::forward_traversal_tag> {
  using value_t = std::vector<uint>;
  
 private:
  friend class boost::iterator_core_access;
  
  bool equal(Integer_combination_iterator const& other) const {
    return (is_end_ && other.is_end_);
  }

  value_t const& dereference() const {
    return value_;
  }

  void increment() {
    uint j1 = 0;
    uint s = 0;
    while (value_[j1] == 0 && j1 < k_)
      j1++;
    uint j2 = j1+1;
    while (value_[j2] == bounds_[j2]) {
      if (bounds_[j2] != 0) {
	s += value_[j1];
	value_[j1] = 0;
	j1 = j2;
      }
      j2++;
    }
    if (j2 >= k_) {
      is_end_ = true;
      return;
    }
    s += value_[j1] - 1;
    value_[j1] = 0;
    value_[j2]++;
    uint i = 0;
    while (s >= bounds_[i]) {
      value_[i] = bounds_[i];
      s -= bounds_[i];
      i++;
    }
    value_[i++] = s;
  }

public:
  template <class Bound_range>
  Integer_combination_iterator(const uint& n, const uint& k, const Bound_range& bounds)
    :
    value_(k+2),
    is_end_(n == 0 || k == 0),
    n_(n),
    k_(k)
  {
    bounds_.reserve(k+2);
    uint sum_radices = 0;
    for (auto b: bounds) {
      bounds_.push_back(b);
      sum_radices += b;      
    }
    bounds_.push_back(2);
    bounds_.push_back(1);
    if (n > sum_radices) {
      is_end_ = true;
      return;
    }
    uint i = 0;
    uint s = n;
    while (s >= bounds_[i]) {
      value_[i] = bounds_[i];
      s -= bounds_[i];
      i++;
    }
    value_[i++] = s;

    while (i < k_)
      value_[i++] = 0;    
    value_[k] = 1;
    value_[k+1] = 0;
  }

  // Used for the creating an end iterator
  Integer_combination_iterator() : is_end_(true), n_(0), k_(0) {}

 private:
  value_t value_; // the dereference value
  bool is_end_;   // is true when the current integer combination is the final one 

  uint n_;
  uint k_;
  std::vector<uint> bounds_;
};

} // namespace coxeter_triangulation

} // namespace Gudhi

#endif

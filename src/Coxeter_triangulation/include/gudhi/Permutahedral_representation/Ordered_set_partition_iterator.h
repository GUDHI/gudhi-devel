/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PERMUTAHEDRAL_REPRESENTATION_ORDERED_SET_PARTITION_ITERATOR_H_
#define PERMUTAHEDRAL_REPRESENTATION_ORDERED_SET_PARTITION_ITERATOR_H_

#include <vector>
#include <limits>

#include <gudhi/Permutahedral_representation/Permutation_iterator.h>
#include <gudhi/Permutahedral_representation/Set_partition_iterator.h>

#include <boost/range/iterator_range.hpp>

namespace Gudhi {

namespace coxeter_triangulation {

typedef unsigned uint;

/** \brief Class that represents an ordered set partition of a set {0,...,n-1} in k parts as
 *         a pair of an unordered set partition given in lexicographic order and
 *         a permutation of the parts.
 */
struct Ordered_set_partition {
  Set_partition_iterator s_it_;
  Permutation_iterator p_it_;

  // Ordered_set_partition(const Set_partition_iterator& s_it, const Permutation_iterator& p_it)
  //   : s_it_(s_it), p_it_(p_it) {}

  const std::vector<uint> operator[](const uint& i) const { return (*s_it_)[(*p_it_)[i]]; }

  std::size_t size() const { return s_it_->size(); }
};

/** \brief Class that allows the user to generate set partitions of a set {0,...,n-1} in k parts.
 *
 */
class Ordered_set_partition_iterator
    : public boost::iterator_facade<Ordered_set_partition_iterator, Ordered_set_partition const,
                                    boost::forward_traversal_tag> {
  using value_t = Ordered_set_partition;

 private:
  friend class boost::iterator_core_access;

  bool equal(Ordered_set_partition_iterator const& other) const { return (is_end_ && other.is_end_); }

  value_t const& dereference() const { return value_; }

  void increment() {
    if (++value_.p_it_ == p_end_) {
      if (++value_.s_it_ == s_end_) {
        is_end_ = true;
        return;
      } else
        value_.p_it_.reinitialize();
    }
  }

 public:
  Ordered_set_partition_iterator(const uint& n, const uint& k)
      : value_({Set_partition_iterator(n, k), Permutation_iterator(k)}), is_end_(n == 0) {}

  // Used for the creating an end iterator
  Ordered_set_partition_iterator() : is_end_(true) {}

  void reinitialize() {
    is_end_ = false;
    value_.p_it_.reinitialize();
    value_.s_it_.reinitialize();
  }

 private:
  Set_partition_iterator s_end_;  // Set partition iterator and the corresponding end iterator
  Permutation_iterator p_end_;    // Permutation iterator and the corresponding end iterator
  value_t value_;                 // the dereference value
  bool is_end_;                   // is true when the current permutation is the final one
};

}  // namespace coxeter_triangulation

}  // namespace Gudhi

#endif

/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria, Hannah Schreiber
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef GUDHI_SIMPLE_UNION_FIND_H_
#define GUDHI_SIMPLE_UNION_FIND_H_

#include <cstddef>  // std::size_t
#include <vector>

#include <boost/pending/disjoint_sets.hpp>

namespace Gudhi {

namespace persistent_cohomology {

template <typename T>
struct Simple_union_find {
 public:
  Simple_union_find(int num_sets = 0)
      : ds_rank_(num_sets), ds_parent_(num_sets), dsets_(ds_rank_.data(), ds_parent_.data()) {}

  void resize(std::size_t n) {
    ds_rank_.resize(n);
    ds_parent_.resize(n);
    dsets_ = boost::disjoint_sets<int*, T*>(ds_rank_.data(), ds_parent_.data());
  }

  void make_set(T key) { dsets_.make_set(key); }

  T find_set(T key) { return dsets_.find_set(key); }

  void merge_sets(T k1, T k2) { dsets_.link(k1, k2); }

  T get_parent(T key) const { return ds_parent_[key]; }

 private:
  std::vector<int> ds_rank_;
  std::vector<T> ds_parent_;
  boost::disjoint_sets<int*, T*> dsets_;
};

}  // namespace persistent_cohomology

}  // namespace Gudhi

#endif  // GUDHI_SIMPLE_UNION_FIND_H_

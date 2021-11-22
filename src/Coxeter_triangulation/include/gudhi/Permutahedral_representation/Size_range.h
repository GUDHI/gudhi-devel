/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PERMUTAHEDRAL_REPRESENTATION_SIZE_RANGE_H_
#define PERMUTAHEDRAL_REPRESENTATION_SIZE_RANGE_H_

#include <cstdlib>  // for std::size_t

#include <boost/range/iterator_range.hpp>

namespace Gudhi {

namespace coxeter_triangulation {

/** \brief Auxillary iterator class for sizes of parts in an ordered set partition.
 */
template <class T_it>
class Size_iterator
    : public boost::iterator_facade<Size_iterator<T_it>, std::size_t const, boost::forward_traversal_tag> {
  friend class boost::iterator_core_access;

 private:
  bool equal(Size_iterator const& other) const { return (is_end_ && other.is_end_); }

  std::size_t const& dereference() const { return value_; }

  void increment() {
    if (++t_it_ == t_end_) {
      is_end_ = true;
      return;
    }
    value_ = t_it_->size() - 1;
  }

 public:
  Size_iterator(const T_it& t_begin, const T_it& t_end) : t_it_(t_begin), t_end_(t_end), is_end_(t_begin == t_end) {
    if (!is_end_) value_ = t_it_->size() - 1;
  }

 private:
  T_it t_it_, t_end_;
  bool is_end_;
  std::size_t value_;
};

template <class T>
class Size_range {
  const T& t_;

 public:
  typedef Size_iterator<typename T::const_iterator> iterator;

  Size_range(const T& t) : t_(t) {}

  std::size_t operator[](std::size_t i) const { return t_[i].size() - 1; }

  iterator begin() const { return iterator(t_.begin(), t_.end()); }

  iterator end() const { return iterator(t_.end(), t_.end()); }
};

}  // namespace coxeter_triangulation

}  // namespace Gudhi

#endif

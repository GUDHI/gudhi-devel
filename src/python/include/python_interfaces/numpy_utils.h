/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef INCLUDE_NUMPY_UTILS_PYTHON_H_
#define INCLUDE_NUMPY_UTILS_PYTHON_H_

#include <cstddef>
#include <type_traits>
#include <vector>

#include <boost/iterator/iterator_facade.hpp>
#include <nanobind/ndarray.h>

template <class T, typename... Shape>
inline auto _wrap_as_numpy_array(std::vector<T> *tensor, Shape... shapes)
{
  return nanobind::ndarray<nanobind::numpy, T>(
      tensor->data(), {static_cast<std::size_t>(shapes)...}, nanobind::capsule(tensor, [](void *p) noexcept {
        delete reinterpret_cast<std::vector<T> *>(p);
      }));
}

// tensor has to be declared with 'new []'
template <class T, typename... Shape>
inline auto _wrap_as_numpy_array(T *tensor, Shape... shapes)
{
  return nanobind::ndarray<nanobind::numpy, T>(
      tensor, {static_cast<std::size_t>(shapes)...}, nanobind::capsule(tensor, [](void *p) noexcept {
        delete[] reinterpret_cast<T *>(p);
      }));
}

template <typename T, class = std::enable_if<std::is_arithmetic_v<T> > >
class Numpy_span
{
 public:
  using value_type = const T;
  using const_iterator = value_type *;
  using iterator = const_iterator;
  using difference_type = std::ptrdiff_t;
  using size_type = std::size_t;

  Numpy_span(const nanobind::ndarray<const T, nanobind::ndim<1> > &array)
      : begin_(array.data()), end_(begin_ + array.shape(0)) {};

  Numpy_span(value_type *begin, value_type *end) : begin_(begin), end_(end) {};

  iterator begin() const noexcept { return begin_; }

  iterator end() const noexcept { return end_; }

  size_type size() const { return end_ - begin_; }

  bool empty() const { return end_ == begin_; }

 private:
  value_type *begin_;
  value_type *end_;
};

template <typename T, class = std::enable_if<std::is_arithmetic_v<T> > >
class Numpy_2d_span
{
 public:
  using value_type = const T;
  using difference_type = std::ptrdiff_t;
  using size_type = std::size_t;

  class iterator : public boost::iterator_facade<iterator, value_type *, boost::forward_traversal_tag, value_type *>
  {
   public:
    iterator(T *curr, T *end, std::size_t stride) : curr_(curr), end_(end), stride_(stride)
    {
      if (curr_ > end_) curr_ = end_;
    }

    iterator(T *end) : curr_(end), end_(end), stride_(0) {}

   private:
    friend class boost::iterator_core_access;

    bool equal(iterator const &other) const { return curr_ == other.curr_ && end_ == other.end_; }

    T const *dereference() const { return curr_; }

    void increment()
    {
      curr_ += stride_;
      if (curr_ > end_) curr_ = end_;
    }

    value_type *curr_;
    value_type *end_;
    const std::size_t stride_;
  };
  using const_iterator = iterator;

  Numpy_2d_span(const nanobind::ndarray<const T, nanobind::ndim<2> > &array)
      : begin_(array.data()), end_(begin_ + (array.shape(0) * array.shape(1))), stride_(array.shape(1)) {};

  Numpy_2d_span(value_type *begin, value_type *end, std::size_t stride) : begin_(begin), end_(end), stride_(stride) {};

  iterator begin() const noexcept { return iterator(begin_, end_, stride_); }

  iterator end() const noexcept { return iterator(end_); }

  size_type size() const { return (end_ - begin_) / stride_; }

  bool empty() const { return end_ == begin_; }

 private:
  value_type *begin_;
  value_type *end_;
  const std::size_t stride_;
};

#endif  // INCLUDE_NUMPY_UTILS_PYTHON_H_

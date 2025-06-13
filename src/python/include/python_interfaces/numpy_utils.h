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
  using Array = nanobind::ndarray<const T, nanobind::ndim<2> >;
  using value_type = const T;
  using difference_type = std::ptrdiff_t;
  using size_type = std::size_t;
  using const_reference = value_type *;
  using const_iterator = value_type *;

  using iterator = const_iterator;

  Numpy_2d_span(const Array &array) : array_view_(array.view()) {};

  iterator begin() const noexcept { return &array_view_(0, 0); }

  iterator end() const noexcept { return begin() + array_view_.shape(0) * array_view_.shape(1); }

  size_type size() const { return array_view_.shape(0); }

  bool empty() const { return array_view_.size() == 0; }

  const_reference operator[](size_type pos) const { return &array_view_(pos, 0); }

 private:
  using View = decltype(std::declval<Array>().view());

  View array_view_;
};

template <typename T, class = std::enable_if<std::is_arithmetic_v<T> > >
class Numpy_array_element_iterator
    : public boost::iterator_facade<Numpy_array_element_iterator<T>, const T, boost::forward_traversal_tag, const T &>
{
 public:
  using value_type = const T;
  using difference_type = std::ptrdiff_t;
  using size_type = std::size_t;
  using const_reference = value_type &;

  Numpy_array_element_iterator(value_type *curr, value_type *end, std::size_t stride)
      : curr_(curr), end_(end), stride_(stride)
  {
    if (curr_ > end_) curr_ = end_;
  }

  Numpy_array_element_iterator(value_type *end) : curr_(end), end_(end), stride_(0) {}

 private:
  friend class boost::iterator_core_access;

  bool equal(Numpy_array_element_iterator const &other) const { return curr_ == other.curr_ && end_ == other.end_; }

  const_reference dereference() const { return *curr_; }

  void increment()
  {
    curr_ += stride_;
    if (curr_ > end_) curr_ = end_;
  }

  value_type *curr_;
  value_type *end_;
  size_type stride_;
};

#endif  // INCLUDE_NUMPY_UTILS_PYTHON_H_

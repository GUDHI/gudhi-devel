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

#include <nanobind/ndarray.h>

#include <boost/iterator/iterator_facade.hpp>
#include <boost/range/iterator_range_core.hpp>

#include <gudhi/Debug_utils.h>

template <class T, typename... Shape>
inline auto _wrap_as_numpy_array(std::vector<T> &&tensor, Shape... shapes)
{
  std::vector<T> *tensor_ptr = new std::vector<T>(std::move(tensor));
  return nanobind::ndarray<nanobind::numpy, T>(
      tensor_ptr->data(), {static_cast<std::size_t>(shapes)...}, nanobind::capsule(tensor_ptr, [](void *p) noexcept {
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

template <class T, std::size_t I>
inline auto _wrap_as_numpy_array(std::vector<std::array<T, I> > &&tensor)
{
  std::vector<std::array<T, I> > *tensor_ptr = new std::vector<std::array<T, I> >(std::move(tensor));
  return nanobind::ndarray<nanobind::numpy, T>(
      tensor_ptr->data(), {tensor_ptr->size(), I}, nanobind::capsule(tensor_ptr, [](void *p) noexcept {
        delete reinterpret_cast<std::vector<std::array<T, I> > *>(p);
      }));
}

template <typename T, class = std::enable_if<std::is_arithmetic_v<T> > >
class Numpy_span
{
 public:
  using value_type = T;
  using const_pointer = value_type const *;
  using const_iterator = const_pointer;
  using iterator = const_iterator;
  using difference_type = std::ptrdiff_t;
  using size_type = std::size_t;

  Numpy_span(const nanobind::ndarray<const T, nanobind::ndim<1>, nanobind::any_contig> &array)
      : begin_(array.data()), end_(begin_ + array.shape(0)) {};

  Numpy_span(const_pointer begin, const_pointer end) : begin_(begin), end_(end) {};

  const_iterator begin() const noexcept { return begin_; }

  const_iterator end() const noexcept { return end_; }

  size_type size() const { return end_ - begin_; }

  bool empty() const { return end_ == begin_; }

 private:
  const_pointer begin_;
  const_pointer end_;
};

template <typename T, class = std::enable_if<std::is_arithmetic_v<T> > >
class Numpy_array_element_iterator
    : public boost::iterator_facade<Numpy_array_element_iterator<T>, const T, boost::random_access_traversal_tag>
{
 private:
  using Base = boost::iterator_facade<Numpy_array_element_iterator<T>, const T, boost::random_access_traversal_tag>;

 public:
  using const_pointer = T const *;
  using const_reference = typename Base::reference;        // T const&
  using difference_type = typename Base::difference_type;  // std::ptrdiff_t
  using size_type = std::size_t;

  Numpy_array_element_iterator(const_pointer curr, const_pointer end, std::size_t stride)
      : curr_(curr), end_(end), stride_(stride)
  {
    if (curr_ > end_) curr_ = end_;
  }

  Numpy_array_element_iterator(const_pointer end, std::size_t stride) : curr_(end), end_(end), stride_(stride) {}

  // necessary for boost::iterator_range to be able to use operator[].
  // Seg fails otherwise for some reasons
  const_reference operator[](difference_type n) const { return *(curr_ + n * stride_); }

 private:
  friend class boost::iterator_core_access;

  bool equal(Numpy_array_element_iterator const &other) const
  {
    GUDHI_CHECK(end_ == other.end_, "The two iterators are not set on the same range and are therefore not comparable");
    return curr_ == other.curr_;
  }

  const_reference dereference() const { return *curr_; }

  void increment()
  {
    curr_ += stride_;
    GUDHI_CHECK(curr_ <= end_, "Incrementation went out of bound.");
  }

  void decrement() { curr_ -= stride_; }

  void advance(difference_type n)
  {
    curr_ += n * stride_;
    GUDHI_CHECK(curr_ <= end_, "Incrementation went out of bound.");
  }

  difference_type distance_to(const Numpy_array_element_iterator &other) const
  {
    if (stride_ == 0) return 0;  // in case of empty ranges
    return (other.curr_ - curr_) / static_cast<difference_type>(stride_);
  }

  const_pointer curr_;
  const_pointer end_;
  size_type stride_;
};

template <typename T, class View, class = std::enable_if<std::is_arithmetic_v<T> > >
boost::iterator_range<Numpy_array_element_iterator<T> > make_element_range(const T *start,
                                                                           const View &view,
                                                                           bool nonTransposed = true)
{
  auto end1 = start + view.shape(nonTransposed) * view.stride(nonTransposed);
  return boost::iterator_range<Numpy_array_element_iterator<T> >(
      Numpy_array_element_iterator(start, end1, view.stride(nonTransposed)),
      Numpy_array_element_iterator(end1, view.stride(nonTransposed)));
}

template <typename T, class = std::enable_if<std::is_arithmetic_v<T> > >
class Numpy_2d_span
{
 public:
  using Array = nanobind::ndarray<const T, nanobind::ndim<2> >;
  using value_type = T;
  using const_pointer = T const *;
  using difference_type = std::ptrdiff_t;
  using size_type = std::size_t;
  using const_iterator = Numpy_array_element_iterator<T>;
  using iterator = const_iterator;
  using const_reference = boost::iterator_range<Numpy_array_element_iterator<T> >;

  Numpy_2d_span(const Array &array) : array_view_(array.view()) {};

  const_iterator begin() const noexcept { return iterator(get_start_ptr(), get_end_ptr(), array_view_.stride(0)); }

  const_iterator end() const noexcept { return iterator(get_end_ptr(), array_view_.stride(0)); }

  size_type size() const { return array_view_.shape(0); }

  bool empty() const { return array_view_.size() == 0; }

  const_reference operator[](size_type pos) const { return make_element_range(&array_view_(pos, 0), array_view_); }

 private:
  using View = decltype(std::declval<Array>().view());

  View array_view_;

  const_pointer get_start_ptr() const { return &array_view_(0, 0); }

  const_pointer get_end_ptr() const { return get_start_ptr() + array_view_.shape(0) * array_view_.shape(1); }
};

#endif  // INCLUDE_NUMPY_UTILS_PYTHON_H_

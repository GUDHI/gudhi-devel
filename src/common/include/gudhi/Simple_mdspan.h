/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber, David Loiseaux
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file simple_mdspan.h
 * @author Hannah Schreiber, David Loiseaux
 * @brief Contains the @ref Gudhi::Simple_mdspan class.
 */

#ifndef GUDHI_SIMPLE_MDSPAN_H_
#define GUDHI_SIMPLE_MDSPAN_H_

#include <cstddef>      // std::size_t
#include <type_traits>  // std::remove_cv_t
#include <limits>
#include <initializer_list>
#include <vector>

#include <gudhi/Debug_utils.h>

namespace Gudhi {

template <typename T>
class Simple_mdspan
{
 public:
  using index_type = std::size_t;
  using extents_type = std::vector<index_type>;
  using element_type = T;
  using value_type = std::remove_cv_t<T>;
  using size_type = extents_type::size_type;
  using rank_type = std::size_t;
  using data_handle_type = T*;
  using reference = T&;

  constexpr Simple_mdspan() : ptr_(nullptr) {}

  constexpr Simple_mdspan(const Simple_mdspan& rhs) = default;
  constexpr Simple_mdspan(Simple_mdspan&& rhs) = default;

  template <class... IndexTypes>
  constexpr explicit Simple_mdspan(data_handle_type ptr, IndexTypes... exts)
      : Simple_mdspan(ptr, {static_cast<index_type>(exts)...})
  {}

  template <class IndexRange = std::initializer_list<index_type> >
  constexpr Simple_mdspan(data_handle_type ptr, const IndexRange& exts) : ptr_(ptr), exts_(exts.begin(), exts.end())
  {
    GUDHI_CHECK(ptr != nullptr || exts_.empty() || exts_[0] == 0, "Given pointer is not properly initialized.");
    if (!exts_.empty()) _initialize_strides();
  }

  constexpr Simple_mdspan& operator=(const Simple_mdspan& rhs) = default;
  constexpr Simple_mdspan& operator=(Simple_mdspan&& rhs) = default;

  //version with [] not possible before C++23
  template <class... IndexTypes>
  constexpr reference operator()(IndexTypes... indices) const
  {
    return operator[]({static_cast<index_type>(indices)...});
  }

  template <class IndexRange = std::initializer_list<index_type> >
  constexpr reference operator[](const IndexRange& indices) const
  {
    GUDHI_CHECK(indices.size() == exts_.size(), "Wrong number of parameters.");

    data_handle_type data = ptr_;
    auto it = indices.begin();
    GUDHI_CHECK_code(unsigned int i = 0);
    for (auto stride : ext_shifts_) {
      GUDHI_CHECK_code(GUDHI_CHECK(*it < exts_[i], "Out of bound index."));
      data += (stride * (*it));
      ++it;
      GUDHI_CHECK_code(++i);
    }

    return *data;
  }

  //replaces mapping() from the original mdspan
  template <class... IndexTypes>
  constexpr index_type get_index(IndexTypes... indices) const
  {
    return get_index({static_cast<index_type>(indices)...});
  }

  template <class IndexRange = std::initializer_list<index_type> >
  constexpr index_type get_index(const IndexRange& indices) const
  {
    data_handle_type data = &operator[](indices);
    return data - ptr_;
  }

  constexpr rank_type rank() noexcept { return exts_.size(); }

  constexpr rank_type rank_dynamic() noexcept { return exts_.size(); }

  static constexpr std::size_t static_extent(rank_type r) noexcept { return std::numeric_limits<std::size_t>::max(); }

  constexpr index_type extent(rank_type r) const
  {
    GUDHI_CHECK(r < exts_.size(), "Out of bound index.");
    return exts_[r];
  }

  constexpr size_type size() const noexcept
  {
    if (exts_.empty()) return 0;
    return ext_shifts_[0] * exts_[0];
  }

  constexpr bool empty() const noexcept { return exts_.empty(); }

  constexpr index_type stride(rank_type r) const
  {
    GUDHI_CHECK(r < ext_shifts_.size(), "Stride out of bound.");
    return ext_shifts_[r];
  }

  constexpr const extents_type& extents() const noexcept { return exts_; }

  constexpr const data_handle_type& data_handle() const noexcept { return ptr_; }

  // if is_unique() is true for all possible instantiations of this class
  static constexpr bool is_always_unique() { return true; }

  // if is_exhaustive() is true for all possible instantiations of this class
  static constexpr bool is_always_exhaustive() { return true; }

  // if is_strided() is true for all possible instantiations of this class
  static constexpr bool is_always_strided() { return true; }

  // unicity of the mapping (i,j,k,...) -> real index
  constexpr bool is_unique() const { return true; }

  // if all real indices have a preimage in form (i,j,k,...)
  constexpr bool is_exhaustive() const { return true; }

  // if distance in memory is constant between two values in same rank
  constexpr bool is_strided() const { return true; }

  friend constexpr void swap(Simple_mdspan& x, Simple_mdspan& y) noexcept
  {
    std::swap(x.ptr_, y.ptr_);
    x.exts_.swap(y.exts_);
    x.ext_shifts_.swap(y.ext_shifts_);
  }

  //as not everything is computed at compile time as for mdspan, update is usually faster than reconstructing everytime.
  void update_extent(rank_type r, index_type new_value){
    GUDHI_CHECK(r < exts_.size(), "Index out of bound.");
    exts_[r] = new_value;
    _update_strides(r);
  }
  //for update_extent to make sense, as resizing the vector can move it in the memory
  void update_data(data_handle_type ptr){
    GUDHI_CHECK(ptr != nullptr, "Null pointer not valid input.");
    ptr_ = ptr;
  }

 private:
  data_handle_type ptr_;
  extents_type exts_;
  extents_type ext_shifts_;

  void _initialize_strides()
  {
    ext_shifts_.resize(exts_.size());
    ext_shifts_.back() = 1;
    for (auto i = exts_.size() - 1; i > 0; --i) {
      ext_shifts_[i - 1] = ext_shifts_[i] * exts_[i];
    }
  }

  void _update_strides(rank_type start)
  {
    for (auto i = start; i > 0; --i) {
      ext_shifts_[i - 1] = ext_shifts_[i] * exts_[i];
    }
  }
};

}  // namespace Gudhi

#endif  // GUDHI_SIMPLE_MDSPAN_H_
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
 * @private
 * @file simple_mdspan.h
 * @author Hannah Schreiber, David Loiseaux
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

/**
 * @private
 * @brief Reproduces the behaviour of C++23 `std::layout_right` class.
 */
class layout_right
{
 public:
  class mapping
  {
   public:
    using index_type = std::size_t;
    using extents_type = std::vector<index_type>;
    using size_type = typename extents_type::size_type;
    using rank_type = std::size_t;
    using layout_type = layout_right;

    // constructors
    mapping() noexcept = default;
    mapping(const mapping&) noexcept = default;

    mapping(const extents_type& exts) noexcept : exts_(exts)
    {
      if (!exts_.empty()) _initialize_strides();
    }

    mapping& operator=(const mapping&) noexcept = default;

    // observers
    constexpr const extents_type& extents() const noexcept { return exts_; }

    index_type required_span_size() const noexcept
    {
      if (exts_.empty()) return 0;
      return ext_shifts_[0] * exts_[0];
    }

    template <class... Indices>
    constexpr index_type operator()(Indices... indices) const
    {
      return operator()({static_cast<index_type>(indices)...});
    }

    template <class IndexRange = std::initializer_list<index_type> >
    constexpr index_type operator()(const IndexRange& indices) const
    {
      GUDHI_CHECK(indices.size() == exts_.size(), "Wrong number of parameters.");

      index_type newIndex = 0;
      auto it = indices.begin();
      GUDHI_CHECK_code(unsigned int i = 0);
      for (auto stride : ext_shifts_) {
        GUDHI_CHECK_code(GUDHI_CHECK(*it < exts_[i], "Out of bound index."));
        newIndex += (stride * (*it));
        ++it;
        GUDHI_CHECK_code(++i);
      }

      return newIndex;
    }

    static constexpr bool is_always_unique() noexcept { return true; }

    static constexpr bool is_always_exhaustive() noexcept { return true; }

    static constexpr bool is_always_strided() noexcept { return true; }

    static constexpr bool is_unique() noexcept { return true; }

    static constexpr bool is_exhaustive() noexcept { return true; }

    static constexpr bool is_strided() noexcept { return true; }

    index_type stride(rank_type r) const
    {
      GUDHI_CHECK(r < ext_shifts_.size(), "Stride out of bound.");
      return ext_shifts_[r];
    }

    friend bool operator==(const mapping& m1, const mapping& m2) noexcept { return m1.exts_ == m2.exts_; }

    friend void swap(mapping& m1, mapping& m2) noexcept
    {
      m1.exts_.swap(m2.exts_);
      m1.ext_shifts_.swap(m2.ext_shifts_);
    }

    // as not everything is computed at compile time as for mdspan, update is usually faster than reconstructing
    // everytime.
    void update_extent(rank_type r, index_type new_value)
    {
      GUDHI_CHECK(r < exts_.size(), "Index out of bound.");
      exts_[r] = new_value;
      _update_strides(r);
    }

   private:
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
};

/**
 * @private
 * @brief Simplified version of C++23 `std::mdspan` class that compiles with C++17.
 *
 * Main differences:
 * - extends are all dynamic,
 * - there is no Extend class, everything is managed by the mapper class instead,
 * - there is no AccessorPolicy template: the container pointed by the stored pointer is assumed to be vector-like,
 * i.e., continuous and, e.g., you can do `ptr_ + 2` to access the third element,
 * - `object[i,j,k,...]` is replaced by either `object(i,j,k,...)` or `object[{i,j,k,...}]`, as C++17 does not
 * allow more than one argument for `operator[]`,
 * - two additional methods: `update_extent` and `update_data` to avoid recalculating the helpers in the mapping class
 * at each size modification of the underlying container. In the original `std::mdspan` most of the work is done
 * at compile time, so it is usually fine to reconstruct a view everytime needed. That is not the case here.
 */
template <typename T, class LayoutPolicy = layout_right>
class Simple_mdspan
{
 public:
  using layout_type = LayoutPolicy;
  using mapping_type = typename LayoutPolicy::mapping;
  using extents_type = typename mapping_type::extents_type;
  using element_type = T;
  using value_type = std::remove_cv_t<T>;
  using index_type = typename mapping_type::index_type;
  using size_type = typename mapping_type::size_type;
  using rank_type = typename mapping_type::rank_type;
  using data_handle_type = T*;
  using reference = T&;

  Simple_mdspan() : ptr_(nullptr) {}

  Simple_mdspan(const Simple_mdspan& rhs) = default;
  Simple_mdspan(Simple_mdspan&& rhs) = default;

  template <class... IndexTypes>
  explicit Simple_mdspan(data_handle_type ptr, IndexTypes... exts)
      : Simple_mdspan(ptr, {static_cast<index_type>(exts)...})
  {}

  template <class IndexRange = std::initializer_list<index_type> >
  Simple_mdspan(data_handle_type ptr, const IndexRange& exts) : ptr_(ptr), map_(extents_type(exts.begin(), exts.end()))
  {
    GUDHI_CHECK(ptr != nullptr || empty() || *(exts.begin()) == 0, "Given pointer is not properly initialized.");
  }

  Simple_mdspan(data_handle_type ptr, const mapping_type& m) : ptr_(ptr), map_(m) {}

  Simple_mdspan& operator=(const Simple_mdspan& rhs) = default;
  Simple_mdspan& operator=(Simple_mdspan&& rhs) = default;

  // version with [] not possible before C++23
  template <class... IndexTypes>
  constexpr reference operator()(IndexTypes... indices) const
  {
    return operator[]({static_cast<index_type>(indices)...});
  }

  template <class IndexRange = std::initializer_list<index_type> >
  reference operator[](const IndexRange& indices) const
  {
    return *(ptr_ + map_(indices));
  }

  constexpr rank_type rank() noexcept { return map_.extents().size(); }

  constexpr rank_type rank_dynamic() noexcept { return map_.extents().size(); }

  static constexpr std::size_t static_extent(rank_type r) noexcept { return std::numeric_limits<std::size_t>::max(); }

  constexpr index_type extent(rank_type r) const
  {
    GUDHI_CHECK(r < map_.extents().size(), "Out of bound index.");
    return map_.extents()[r];
  }

  constexpr size_type size() const noexcept { return map_.required_span_size(); }

  constexpr bool empty() const noexcept { return map_.required_span_size() == 0; }

  constexpr index_type stride(rank_type r) const { return map_.stride(r); }

  constexpr const extents_type& extents() const noexcept { return map_.extents(); }

  constexpr const data_handle_type& data_handle() const noexcept { return ptr_; }

  constexpr const mapping_type& mapping() const noexcept { return map_; }

  // if is_unique() is true for all possible instantiations of this class
  static constexpr bool is_always_unique() { return mapping_type::is_always_unique(); }

  // if is_exhaustive() is true for all possible instantiations of this class
  static constexpr bool is_always_exhaustive() { return mapping_type::is_always_exhaustive(); }

  // if is_strided() is true for all possible instantiations of this class
  static constexpr bool is_always_strided() { return mapping_type::is_always_strided(); }

  // unicity of the mapping (i,j,k,...) -> real index
  constexpr bool is_unique() const { return map_.is_unique(); }

  // if all real indices have a preimage in form (i,j,k,...)
  constexpr bool is_exhaustive() const { return map_.is_exhaustive(); }

  // if distance in memory is constant between two values in same rank
  constexpr bool is_strided() const { return map_.is_strided(); }

  friend constexpr void swap(Simple_mdspan& x, Simple_mdspan& y) noexcept
  {
    std::swap(x.ptr_, y.ptr_);
    swap(x.map_, y.map_);
  }

  // as not everything is computed at compile time as for mdspan, update is usually faster than reconstructing
  // everytime.
  void update_extent(rank_type r, index_type new_value) { map_.update_extent(r, new_value); }

  // for update_extent to make sense, as resizing the vector can move it in the memory
  void update_data(data_handle_type ptr)
  {
    GUDHI_CHECK(ptr != nullptr, "Null pointer not valid input.");
    ptr_ = ptr;
  }

 private:
  data_handle_type ptr_;
  mapping_type map_;
};

}  // namespace Gudhi

#endif  // GUDHI_SIMPLE_MDSPAN_H_
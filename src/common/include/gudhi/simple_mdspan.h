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
#include <stdexcept>
#include <type_traits>  // std::remove_cv_t, std::make_unsigned_t, std::integral_constant
#include <limits>
#include <initializer_list>
#include <utility>
#include <array>

#include <gudhi/Debug_utils.h>

namespace Gudhi {

inline constexpr std::size_t dynamic_extent = std::numeric_limits<std::size_t>::max();

template <class IndexType, std::size_t... Extents>
class extents;

namespace detail {

  template <std::size_t v>
  struct is_dynamic : std::integral_constant<std::size_t, 0> {};

  template <>
  struct is_dynamic<dynamic_extent> : std::integral_constant<std::size_t, 1> {};

  template <std::size_t I, class T>
  struct dynamic_count;

  template <std::size_t I, std::size_t first, std::size_t... Tail>
  struct dynamic_count<I, std::integer_sequence<std::size_t, first, Tail...> >
      : std::integral_constant<std::size_t,
                              (is_dynamic<first>::value +
                                dynamic_count<I - 1, std::integer_sequence<std::size_t, Tail...>>::value)> {
  };

  template <std::size_t first, std::size_t... Tail>
  struct dynamic_count<0, std::integer_sequence<std::size_t, first, Tail...> >
      : std::integral_constant<std::size_t, is_dynamic<first>::value> {
  };

  template <std::size_t first>
  struct dynamic_count<0, std::integer_sequence<std::size_t, first> >
      : std::integral_constant<std::size_t, is_dynamic<first>::value> {};

  template <std::size_t I, class T>
  struct extent_value;

  template <std::size_t I, std::size_t first, std::size_t... Tail>
  struct extent_value<I, std::integer_sequence<std::size_t, first, Tail...>>
      : extent_value<I - 1, std::integer_sequence<std::size_t, Tail...>> {};

  template <std::size_t first, std::size_t... Tail>
  struct extent_value<0, std::integer_sequence<std::size_t, first, Tail...>>
      : std::integral_constant<std::size_t, first> {};

  template <class T, T I, T N, T... integers>
  struct dynamic_value_sequence {
    using type = typename dynamic_value_sequence<T, I + 1, N, integers..., dynamic_extent>::type;
  };

  template <class T, T N, T... integers>
  struct dynamic_value_sequence<T, N, N, integers...> {
    using type = std::integer_sequence<T, integers...>;
  };

  template <class IndexType, std::size_t... Pack>
  constexpr auto dynamic_value_extents(std::integer_sequence<std::size_t, Pack...>)
  {
    return extents<IndexType, Pack...>();
  };

  template <class IndexType, std::size_t Rank>
  constexpr auto dynamic_value_extents_value =
      dynamic_value_extents<IndexType>((typename Gudhi::detail::dynamic_value_sequence<std::size_t, 0, Rank>::type{}));

}  // namespace detail

template <class IndexType, std::size_t Rank>
using dextents = decltype(detail::dynamic_value_extents_value<IndexType, Rank>);

/**
 * @private
 * @brief Reproduces the behaviour of C++23 `std::extents` class.
 */
template <class IndexType, std::size_t... Extents>
class extents
{
 public:
  using index_type = IndexType;
  using size_type = std::make_unsigned_t<index_type>;
  using rank_type = std::size_t;

  // observers of the multidimensional index space
  static constexpr rank_type rank() noexcept { return sizeof...(Extents); }

  static constexpr rank_type rank_dynamic() noexcept
  {
    return detail::dynamic_count<rank() - 1, std::integer_sequence<std::size_t, Extents...> >::value;
  }

  static constexpr std::size_t static_extent(rank_type r) noexcept
  {
    std::array<std::size_t, sizeof...(Extents)> exts{Extents...};
    return exts[r];
  }

  constexpr index_type extent(rank_type r) const noexcept
  {
    if (dynamic_extent_shifts_[r] < 0) return static_extent(r);
    return dynamic_extents_[dynamic_extent_shifts_[r]];
  }

  void update_dynamic_extent(rank_type r, index_type i){
    if (dynamic_extent_shifts_[r] < 0) throw std::invalid_argument("Given rank is not dynamic.");
    dynamic_extents_[dynamic_extent_shifts_[r]] = i;
  }

  // constructors
  constexpr extents() noexcept : dynamic_extents_(), dynamic_extent_shifts_(_init_shifts()) {}

  template <class OtherIndexType, std::size_t... OtherExtents>
  constexpr explicit extents(const extents<OtherIndexType, OtherExtents...>& other) noexcept
      : dynamic_extents_(), dynamic_extent_shifts_(_init_shifts())
  {
    for (rank_type r = 0; r < rank(); ++r) {
      if (dynamic_extent_shifts_[r] >= 0) dynamic_extents_[dynamic_extent_shifts_[r]] = other.extent(r);
    }
  }

  template <class... OtherIndexTypes>
  constexpr explicit extents(OtherIndexTypes... extents) noexcept
      : dynamic_extents_{static_cast<IndexType>(extents)...}, dynamic_extent_shifts_(_init_shifts())
  {}

  template <class OtherIndexType, std::size_t N>
  constexpr explicit extents(const std::array<OtherIndexType, N>& other) noexcept
      : dynamic_extents_{other}, dynamic_extent_shifts_(_init_shifts())
  {}

  // comparison operators
  template <class OtherIndexType, std::size_t... OtherExtents>
  friend constexpr bool operator==(const extents& e1, const extents<OtherIndexType, OtherExtents...>& e2) noexcept
  {
    if (e1.rank() != e2.rank()) return false;
    for (rank_type r = 0; r < rank(); ++r) {
      if (e1.extent(r) != e2.extent(r)) return false;
    }
    return true;
  }

  friend void swap(extents& e1, extents& e2) noexcept
  {
    e1.dynamic_extents_.swap(e2.dynamic_extents_);
    e1.dynamic_extent_shifts_.swap(e2.dynamic_extent_shifts_);
  }

  friend std::ostream &operator<<(std::ostream &stream, const extents &e)
  {
    stream << "[ " << sizeof...(Extents) << " ] ";
    ((stream << Extents << ' '), ...);

    return stream;
  }

 private:
  std::array<index_type, rank_dynamic()> dynamic_extents_;
  std::array<int, rank()> dynamic_extent_shifts_;

  static constexpr std::array<int, rank()> _init_shifts()
  {
    std::array<std::size_t, sizeof...(Extents)> exts{Extents...};
    std::array<int, rank()> res = {};
    std::size_t index = 0;
    for (rank_type i = 0; i < rank(); ++i) {
      if (exts[i] == dynamic_extent) {
        res[i] = index;
        ++index;
      } else {
        res[i] = -1;
      }
    }
    return res;
  }
};

// Does not seem to work with C++17(?) because the use of 'dextents' is not explicit enough:
// "trailing return type ‘Gudhi::dextents<long unsigned int, sizeof... (Integrals)>’ of deduction guide is not a
// specialization of ‘Gudhi::extents<IndexType, Extents>’"
// Or does someone knows a workaround...?
// template<class... Integrals>
// explicit extents(Integrals...) -> dextents<std::size_t, sizeof...(Integrals)>;

/**
 * @private
 * @brief Reproduces the behaviour of C++23 `std::layout_right` class.
 */
class layout_right
{
 public:
  template<class Extents>
  class mapping
  {
   public:
    using extents_type = Extents;
    using index_type = typename extents_type::index_type;
    using size_type = typename extents_type::size_type;
    using rank_type = typename extents_type::rank_type;
    using layout_type = layout_right;

    // constructors
    mapping() noexcept = default;
    mapping(const mapping&) noexcept = default;

    mapping(const extents_type& exts) noexcept : exts_(exts)
    {
      if constexpr (extents_type::rank() != 0) _initialize_strides();
    }

    mapping& operator=(const mapping&) noexcept = default;

    // observers
    constexpr const extents_type& extents() const noexcept { return exts_; }

    index_type required_span_size() const noexcept
    {
      if constexpr (extents_type::rank() == 0) return 0;
      else return ext_shifts_[0] * exts_.extent(0);
    }

    template <class... Indices>
    constexpr index_type operator()(Indices... indices) const
    {
      return operator()({static_cast<index_type>(indices)...});
    }

    template <class IndexRange = std::initializer_list<index_type> >
    constexpr index_type operator()(const IndexRange& indices) const
    {
      GUDHI_CHECK(indices.size() == extents_type::rank(), "Wrong number of parameters.");

      index_type newIndex = 0;
      auto it = indices.begin();
      GUDHI_CHECK_code(unsigned int i = 0);
      for (auto stride : ext_shifts_) {
        GUDHI_CHECK_code(GUDHI_CHECK(*it < exts_.extent(i), "Out of bound index."));
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
      swap(m1.exts_, m2.exts_);
      m1.ext_shifts_.swap(m2.ext_shifts_);
    }

    // update can be faster than reconstructing everytime if only relatively small r's are updated.
    void update_extent(rank_type r, index_type new_value)
    {
      GUDHI_CHECK(r < extents_type::rank(), "Index out of bound.");
      exts_.update_dynamic_extent(r, new_value);
      _update_strides(r);
    }

   private:
    extents_type exts_;
    std::array<index_type,extents_type::rank()> ext_shifts_;

    constexpr void _initialize_strides()
    {
      ext_shifts_[extents_type::rank() - 1] = 1;
      for (auto i = extents_type::rank() - 1; i > 0; --i) {
        ext_shifts_[i - 1] = ext_shifts_[i] * exts_.extent(i);
      }
    }

    constexpr void _update_strides(rank_type start)
    {
      for (auto i = start; i > 0; --i) {
        ext_shifts_[i - 1] = ext_shifts_[i] * exts_.extent(i);
      }
    }
  };
};

/**
 * @private
 * @brief Simplified version of C++23 `std::mdspan` class that compiles with C++17.
 *
 * Main differences:
 * - there is no AccessorPolicy template: the container pointed by the stored pointer is assumed to be vector-like,
 * i.e., continuous and, e.g., you can do `ptr_ + 2` to access the third element,
 * - does not implement any "submdspan" methods (C++26),
 * - `object[i,j,k,...]` is replaced by either `object(i,j,k,...)` or `object[{i,j,k,...}]`, as C++17 does not
 * allow more than one argument for `operator[]`,
 * - two additional methods: `update_extent` and `update_data` to avoid recalculating the helpers in the mapping class
 * at each size modification of the underlying container, when the update is trivial (i.e. when only rank 0 is
 * modified, which happens often in our use case). 
 */
template <typename T, class Extents, class LayoutPolicy = layout_right>
class Simple_mdspan
{
 public:
  using layout_type = LayoutPolicy;
  using mapping_type = typename LayoutPolicy::template mapping<Extents>;
  using extents_type = Extents;
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
      : ptr_(ptr), map_(extents_type(exts...))
  {
    GUDHI_CHECK(ptr != nullptr || empty() || Extents::rank() == 0, "Given pointer is not properly initialized.");
  }

  template <class OtherIndexType, size_t N>
  constexpr explicit Simple_mdspan(data_handle_type ptr, const std::array<OtherIndexType, N>& exts)
      : ptr_(ptr), map_(extents_type(exts))
  {
    GUDHI_CHECK(ptr != nullptr || empty() || Extents::rank() == 0, "Given pointer is not properly initialized.");
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

  constexpr rank_type rank() noexcept { return map_.extents().rank(); }

  constexpr rank_type rank_dynamic() noexcept { return map_.extents().rank_dynamic(); }

  static constexpr std::size_t static_extent(rank_type r) noexcept { return std::numeric_limits<std::size_t>::max(); }

  constexpr index_type extent(rank_type r) const
  {
    GUDHI_CHECK(r < map_.extents().rank(), "Out of bound index.");
    return map_.extents().extent(r);
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

template <class CArray>
Simple_mdspan(CArray&)
    -> Simple_mdspan<std::remove_all_extents_t<CArray>, Gudhi::extents<std::size_t, std::extent_v<CArray, 0>>>;

template <class Pointer>
Simple_mdspan(Pointer&&)
    -> Simple_mdspan<std::remove_pointer_t<std::remove_reference_t<Pointer>>, Gudhi::extents<std::size_t>>;

template <class ElementType, class... Integrals>
explicit Simple_mdspan(ElementType*, Integrals...)
    -> Simple_mdspan<ElementType, Gudhi::dextents<std::size_t, sizeof...(Integrals)>>;

template <class ElementType, class OtherIndexType, std::size_t N>
Simple_mdspan(ElementType*, const std::array<OtherIndexType, N>&)
    -> Simple_mdspan<ElementType, Gudhi::dextents<std::size_t, N>>;

template <class ElementType, class IndexType, std::size_t... ExtentsPack>
Simple_mdspan(ElementType*, const Gudhi::extents<IndexType, ExtentsPack...>&)
    -> Simple_mdspan<ElementType, Gudhi::extents<IndexType, ExtentsPack...>>;

template <class ElementType, class MappingType>
Simple_mdspan(ElementType*, const MappingType&)
    -> Simple_mdspan<ElementType, typename MappingType::extents_type, typename MappingType::layout_type>;

}  // namespace Gudhi

#endif  // GUDHI_SIMPLE_MDSPAN_H_
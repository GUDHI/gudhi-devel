/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022-24 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file unordered_set_column.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_matrix::Unordered_set_column class.
 * Also defines the std::hash method for @ref Gudhi::persistence_matrix::Unordered_set_column.
 */

#ifndef PM_UNORDERED_SET_COLUMN_H
#define PM_UNORDERED_SET_COLUMN_H

#include <stdexcept>
#include <type_traits>
#include <algorithm>  // std::lexicographical_compare
#include <utility>    //std::swap, std::move & std::exchange
#include <set>
#include <vector>

#include <boost/iterator/indirect_iterator.hpp>
#if BOOST_VERSION >= 108100
#include <boost/unordered/unordered_flat_set.hpp>
#else
#include <unordered_set>
#endif

#include <gudhi/Persistence_matrix/allocators/entry_constructors.h>

namespace Gudhi {
namespace persistence_matrix {

// For unordered_set container. Outside of Unordered_set_column because of a msvc bug who can't compile properly
// unordered_flat_set if the hash method is nested.
template <class Entry>
struct EntryPointerHash {
  size_t operator()(const Entry* c) const { return std::hash<Entry>()(*c); }
};

template <class Entry>
struct EntryPointerEq {
  bool operator()(const Entry* c1, const Entry* c2) const { return *c1 == *c2; }
};

/**
 * @class Unordered_set_column unordered_set_column.h gudhi/Persistence_matrix/columns/unordered_set_column.h
 * @ingroup persistence_matrix
 *
 * @brief Column class following the @ref PersistenceMatrixColumn concept.
 *
 * Column based on an unordered set structure. The entries are not ordered, but only non-zero values
 * are stored uniquely in the underlying container. When adding an entry range into it, the given entry range
 * also does not need to be ordered (contrary to most other column types).
 *
 * @tparam Master_matrix An instantiation of @ref Matrix from which all types and options are deduced.
 */
template <class Master_matrix>
class Unordered_set_column : public Master_matrix::Row_access_option,
                             public Master_matrix::Column_dimension_option,
                             public Master_matrix::Chain_column_option
{
 public:
  using Master = Master_matrix;
  using Index = typename Master_matrix::Index;
  using ID_index = typename Master_matrix::ID_index;
  using Dimension = typename Master_matrix::Dimension;
  using Field_element = typename Master_matrix::Element;
  using Entry = typename Master_matrix::Matrix_entry;
  using Column_settings = typename Master_matrix::Column_settings;

 private:
  using Field_operators = typename Master_matrix::Field_operators;
  using Entry_constructor = typename Master_matrix::Entry_constructor;

  struct EntryPointerComp {
    bool operator()(const Entry* c1, const Entry* c2) const { return *c1 < *c2; }
  };

#if BOOST_VERSION >= 108100
  using Column_support = boost::unordered_flat_set<Entry*, EntryPointerHash<Entry>, EntryPointerEq<Entry>>;
#else
  using Column_support = std::unordered_set<Entry*, EntryPointerHash<Entry>, EntryPointerEq<Entry>>;
#endif

 public:
  using iterator = boost::indirect_iterator<typename Column_support::iterator>;
  using const_iterator = boost::indirect_iterator<typename Column_support::const_iterator>;
  using Content_range = std::vector<Entry>;

  Unordered_set_column(Column_settings* colSettings = nullptr);
  template <class Container = typename Master_matrix::Boundary>
  Unordered_set_column(const Container& nonZeroRowIndices, Column_settings* colSettings);
  template <class Container = typename Master_matrix::Boundary, class Row_container>
  Unordered_set_column(Index columnIndex,
                       const Container& nonZeroRowIndices,
                       Row_container* rowContainer,
                       Column_settings* colSettings);
  template <class Container = typename Master_matrix::Boundary,
            class = std::enable_if_t<!std::is_arithmetic_v<Container>>>
  Unordered_set_column(const Container& nonZeroRowIndices, Dimension dimension, Column_settings* colSettings);
  template <class Container = typename Master_matrix::Boundary,
            class Row_container,
            class = std::enable_if_t<!std::is_arithmetic_v<Container>>>
  Unordered_set_column(Index columnIndex,
                       const Container& nonZeroRowIndices,
                       Dimension dimension,
                       Row_container* rowContainer,
                       Column_settings* colSettings);
  Unordered_set_column(ID_index idx, Dimension dimension, Column_settings* colSettings);
  Unordered_set_column(ID_index idx,
                       Field_element e,
                       Dimension dimension,
                       Column_settings* colSettings);
  template <class Row_container>
  Unordered_set_column(Index columnIndex,
                       ID_index idx,
                       Dimension dimension,
                       Row_container* rowContainer,
                       Column_settings* colSettings);
  template <class Row_container>
  Unordered_set_column(Index columnIndex,
                       ID_index idx,
                       Field_element e,
                       Dimension dimension,
                       Row_container* rowContainer,
                       Column_settings* colSettings);
  Unordered_set_column(const Unordered_set_column& column, Column_settings* colSettings = nullptr);
  template <class Row_container>
  Unordered_set_column(const Unordered_set_column& column,
                       Index columnIndex,
                       Row_container* rowContainer,
                       Column_settings* colSettings = nullptr);
  Unordered_set_column(Unordered_set_column&& column) noexcept;
  ~Unordered_set_column();

  std::vector<Field_element> get_content(int columnLength = -1) const;
  bool is_non_zero(ID_index rowIndex) const;
  [[nodiscard]] bool is_empty() const;
  [[nodiscard]] std::size_t size() const;

  template <class Row_index_map>
  void reorder(const Row_index_map& valueMap,
               [[maybe_unused]] Index columnIndex = Master_matrix::template get_null_value<Index>());
  void clear();
  void clear(ID_index rowIndex);

  ID_index get_pivot() const;
  Field_element get_pivot_value() const;

  iterator begin() noexcept;
  const_iterator begin() const noexcept;
  iterator end() noexcept;
  const_iterator end() const noexcept;

  Content_range get_non_zero_content_range() const;

  template <class Entry_range>
  Unordered_set_column& operator+=(const Entry_range& column);
  Unordered_set_column& operator+=(Unordered_set_column& column);

  Unordered_set_column& operator*=(const Field_element& v);

  // this = v * this + column
  template <class Entry_range>
  Unordered_set_column& multiply_target_and_add(const Field_element& val, const Entry_range& column);
  Unordered_set_column& multiply_target_and_add(const Field_element& val, Unordered_set_column& column);
  // this = this + column * v
  template <class Entry_range>
  Unordered_set_column& multiply_source_and_add(const Entry_range& column, const Field_element& val);
  Unordered_set_column& multiply_source_and_add(Unordered_set_column& column, const Field_element& val);

  void push_back(const Entry& entry);

  friend bool operator==(const Unordered_set_column& c1, const Unordered_set_column& c2)
  {
    if (&c1 == &c2) return true;
    if (c1.column_.size() != c2.column_.size()) return false;

    for (Entry* entry : c1.column_) {
      auto it = c2.column_.find(entry);
      if (it == c2.column_.end()) return false;
      if (Master_matrix::get_element(**it) != Master_matrix::get_element(*entry)) return false;
    }
    return true;
  }

  friend bool operator<(const Unordered_set_column& c1, const Unordered_set_column& c2)
  {
    if (&c1 == &c2) return false;

    auto comp = [](const Entry* n1, const Entry* n2) -> bool {
      Index r1 = Master_matrix::get_row_index(*n1);
      Index r2 = Master_matrix::get_row_index(*n2);
      Field_element e1 = Master_matrix::get_element(*n1);
      Field_element e2 = Master_matrix::get_element(*n2);

      if (r1 != r2) return r1 < r2;
      if (e1 != e2) return e1 < e2;

      return false;
    };

    std::set<Entry*, decltype(comp)> entries1(comp), entries2(comp);
    entries1.insert(c1.column_.begin(), c1.column_.end());
    entries2.insert(c2.column_.begin(), c2.column_.end());

    return std::lexicographical_compare(entries1.begin(), entries1.end(), entries2.begin(), entries2.end(), comp);
  }

  // Disabled with row access.
  Unordered_set_column& operator=(const Unordered_set_column& other);
  Unordered_set_column& operator=(Unordered_set_column&& other) noexcept;

  friend void swap(Unordered_set_column& col1, Unordered_set_column& col2) noexcept
  {
    swap(static_cast<typename Master_matrix::Row_access_option&>(col1),
         static_cast<typename Master_matrix::Row_access_option&>(col2));
    swap(static_cast<typename Master_matrix::Column_dimension_option&>(col1),
         static_cast<typename Master_matrix::Column_dimension_option&>(col2));
    swap(static_cast<typename Master_matrix::Chain_column_option&>(col1),
         static_cast<typename Master_matrix::Chain_column_option&>(col2));
    col1.column_.swap(col2.column_);
    std::swap(col1.operators_, col2.operators_);
    std::swap(col1.entryPool_, col2.entryPool_);
  }

 private:
  using RA_opt = typename Master_matrix::Row_access_option;
  using Dim_opt = typename Master_matrix::Column_dimension_option;
  using Chain_opt = typename Master_matrix::Chain_column_option;

  Column_support column_;
  Field_operators const* operators_;
  Entry_constructor* entryPool_;

  void _delete_entry(typename Column_support::iterator& it);
  Entry* _insert_entry(ID_index rowIndex, const Field_element& value);
  template <class Entry_range>
  bool _add(const Entry_range& column);
  template <class Entry_range>
  bool _multiply_target_and_add(const Field_element& val, const Entry_range& column);
  template <class Entry_range>
  bool _multiply_source_and_add(const Entry_range& column, const Field_element& val);
  template <class Entry_range, typename F1, typename F2>
  bool _generic_add(const Entry_range& source, F1&& process_source, F2&& update_target);
};

template <class Master_matrix>
inline Unordered_set_column<Master_matrix>::Unordered_set_column(Column_settings* colSettings)
    : RA_opt(),
      Dim_opt(),
      Chain_opt(),
      operators_(Master_matrix::get_operator_ptr(colSettings)),
      entryPool_(colSettings == nullptr ? nullptr : &(colSettings->entryConstructor))
{}

template <class Master_matrix>
template <class Container>
inline Unordered_set_column<Master_matrix>::Unordered_set_column(const Container& nonZeroRowIndices,
                                                                 Column_settings* colSettings)
    : Unordered_set_column(nonZeroRowIndices,
                           nonZeroRowIndices.size() == 0 ? 0 : nonZeroRowIndices.size() - 1,
                           colSettings)
{
  static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type,
                "Constructor not available for chain columns, please specify the dimension of the chain.");
}

template <class Master_matrix>
template <class Container, class Row_container>
inline Unordered_set_column<Master_matrix>::Unordered_set_column(Index columnIndex,
                                                                 const Container& nonZeroRowIndices,
                                                                 Row_container* rowContainer,
                                                                 Column_settings* colSettings)
    : Unordered_set_column(columnIndex,
                           nonZeroRowIndices,
                           nonZeroRowIndices.size() == 0 ? 0 : nonZeroRowIndices.size() - 1,
                           rowContainer,
                           colSettings)
{
  static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type,
                "Constructor not available for chain columns, please specify the dimension of the chain.");
}

template <class Master_matrix>
template <class Container, class>
inline Unordered_set_column<Master_matrix>::Unordered_set_column(const Container& nonZeroRowIndices,
                                                                 Dimension dimension,
                                                                 Column_settings* colSettings)
    : RA_opt(),
      Dim_opt(dimension),
      Chain_opt(nonZeroRowIndices.begin() == nonZeroRowIndices.end()
                    ? Master_matrix::template get_null_value<ID_index>()
                    : Master_matrix::get_row_index(*std::prev(nonZeroRowIndices.end()))),
      column_(nonZeroRowIndices.size()),
      operators_(Master_matrix::get_operator_ptr(colSettings)),
      entryPool_(&(colSettings->entryConstructor))
{
  for (const auto& id : nonZeroRowIndices) {
    _insert_entry(Master_matrix::get_row_index(id),
                  Master_matrix::get_coefficient_value(Master_matrix::get_element(id), operators_));
  }
}

template <class Master_matrix>
template <class Container, class Row_container, class>
inline Unordered_set_column<Master_matrix>::Unordered_set_column(Index columnIndex,
                                                                 const Container& nonZeroRowIndices,
                                                                 Dimension dimension,
                                                                 Row_container* rowContainer,
                                                                 Column_settings* colSettings)
    : RA_opt(columnIndex, rowContainer),
      Dim_opt(dimension),
      Chain_opt(nonZeroRowIndices.begin() == nonZeroRowIndices.end()
                    ? Master_matrix::template get_null_value<ID_index>()
                    : Master_matrix::get_row_index(*std::prev(nonZeroRowIndices.end()))),
      column_(nonZeroRowIndices.size()),
      operators_(Master_matrix::get_operator_ptr(colSettings)),
      entryPool_(&(colSettings->entryConstructor))
{
  for (const auto& id : nonZeroRowIndices) {
    _insert_entry(Master_matrix::get_row_index(id),
                  Master_matrix::get_coefficient_value(Master_matrix::get_element(id), operators_));
  }
}

template <class Master_matrix>
inline Unordered_set_column<Master_matrix>::Unordered_set_column(ID_index idx,
                                                                 Dimension dimension,
                                                                 Column_settings* colSettings)
    : RA_opt(),
      Dim_opt(dimension),
      Chain_opt(idx),
      column_(1),
      operators_(nullptr),
      entryPool_(&(colSettings->entryConstructor))
{
  static_assert(Master_matrix::Option_list::is_z2,
                "Constructor not available for Zp != Z2. Please specify the coefficient.");
  _insert_entry(idx, 1);
}

template <class Master_matrix>
inline Unordered_set_column<Master_matrix>::Unordered_set_column(ID_index idx,
                                                                 Field_element e,
                                                                 Dimension dimension,
                                                                 Column_settings* colSettings)
    : RA_opt(),
      Dim_opt(dimension),
      Chain_opt(idx),
      column_(1),
      operators_(&(colSettings->operators)),
      entryPool_(&(colSettings->entryConstructor))
{
  static_assert(!Master_matrix::Option_list::is_z2,
                "Constructor not available for Zp == Z2. Please do not specify any coefficient.");
  _insert_entry(idx, operators_->get_value(e));
}

template <class Master_matrix>
template <class Row_container>
inline Unordered_set_column<Master_matrix>::Unordered_set_column(Index columnIndex,
                                                                 ID_index idx,
                                                                 Dimension dimension,
                                                                 Row_container* rowContainer,
                                                                 Column_settings* colSettings)
    : RA_opt(columnIndex, rowContainer),
      Dim_opt(dimension),
      Chain_opt(idx),
      column_(1),
      operators_(nullptr),
      entryPool_(&(colSettings->entryConstructor))
{
  static_assert(Master_matrix::Option_list::is_z2,
                "Constructor not available for Zp != Z2. Please specify the coefficient.");
  _insert_entry(idx, 1);
}

template <class Master_matrix>
template <class Row_container>
inline Unordered_set_column<Master_matrix>::Unordered_set_column(Index columnIndex,
                                                                 ID_index idx,
                                                                 Field_element e,
                                                                 Dimension dimension,
                                                                 Row_container* rowContainer,
                                                                 Column_settings* colSettings)
    : RA_opt(columnIndex, rowContainer),
      Dim_opt(dimension),
      Chain_opt(idx),
      column_(1),
      operators_(&(colSettings->operators)),
      entryPool_(&(colSettings->entryConstructor))
{
  static_assert(!Master_matrix::Option_list::is_z2,
                "Constructor not available for Zp == Z2. Please do not specify any coefficient.");
  _insert_entry(idx, operators_->get_value(e));
}

template <class Master_matrix>
inline Unordered_set_column<Master_matrix>::Unordered_set_column(const Unordered_set_column& column,
                                                                 Column_settings* colSettings)
    : RA_opt(),
      Dim_opt(static_cast<const Dim_opt&>(column)),
      Chain_opt(static_cast<const Chain_opt&>(column)),
      column_(column.column_.bucket_count()),
      operators_(colSettings == nullptr ? column.operators_ : Master_matrix::get_operator_ptr(colSettings)),
      entryPool_(colSettings == nullptr ? column.entryPool_ : &(colSettings->entryConstructor))
{
  static_assert(!Master_matrix::Option_list::has_row_access,
                "Simple copy constructor not available when row access option enabled. Please specify the new column "
                "index and the row container.");

  for (const Entry* entry : column.column_) {
    _insert_entry(entry->get_row_index(), entry->get_element());
  }
}

template <class Master_matrix>
template <class Row_container>
inline Unordered_set_column<Master_matrix>::Unordered_set_column(const Unordered_set_column& column,
                                                                 Index columnIndex,
                                                                 Row_container* rowContainer,
                                                                 Column_settings* colSettings)
    : RA_opt(columnIndex, rowContainer),
      Dim_opt(static_cast<const Dim_opt&>(column)),
      Chain_opt(static_cast<const Chain_opt&>(column)),
      column_(column.column_.bucket_count()),
      operators_(colSettings == nullptr ? column.operators_ : Master_matrix::get_operator_ptr(colSettings)),
      entryPool_(colSettings == nullptr ? column.entryPool_ : &(colSettings->entryConstructor))
{
  for (const Entry* entry : column.column_) {
    _insert_entry(entry->get_row_index(), entry->get_element());
  }
}

template <class Master_matrix>
inline Unordered_set_column<Master_matrix>::Unordered_set_column(Unordered_set_column&& column) noexcept
    : RA_opt(std::move(static_cast<RA_opt&>(column))),
      Dim_opt(std::move(static_cast<Dim_opt&>(column))),
      Chain_opt(std::move(static_cast<Chain_opt&>(column))),
      column_(std::move(column.column_)),
      operators_(std::exchange(column.operators_, nullptr)),
      entryPool_(std::exchange(column.entryPool_, nullptr))
{}

template <class Master_matrix>
inline Unordered_set_column<Master_matrix>::~Unordered_set_column()
{
  for (auto* entry : column_) {
    if constexpr (Master_matrix::Option_list::has_row_access) RA_opt::unlink(entry);
    entryPool_->destroy(entry);
  }
}

template <class Master_matrix>
inline std::vector<typename Unordered_set_column<Master_matrix>::Field_element>
Unordered_set_column<Master_matrix>::get_content(int columnLength) const
{
  if (columnLength < 0 && column_.size() > 0)
    columnLength = (*std::max_element(column_.begin(), column_.end(), EntryPointerComp()))->get_row_index() + 1;
  else if (columnLength < 0)
    return std::vector<Field_element>();

  std::vector<Field_element> container(columnLength, 0);
  for (auto it = column_.begin(); it != column_.end(); ++it) {
    if ((*it)->get_row_index() < static_cast<ID_index>(columnLength)) {
      container[(*it)->get_row_index()] = Master_matrix::get_element(**it);
    }
  }
  return container;
}

template <class Master_matrix>
inline bool Unordered_set_column<Master_matrix>::is_non_zero(ID_index rowIndex) const
{
  Entry entry(rowIndex);
  return column_.find(&entry) != column_.end();
}

template <class Master_matrix>
inline bool Unordered_set_column<Master_matrix>::is_empty() const
{
  return column_.empty();
}

template <class Master_matrix>
inline std::size_t Unordered_set_column<Master_matrix>::size() const
{
  return column_.size();
}

template <class Master_matrix>
template <class Row_index_map>
inline void Unordered_set_column<Master_matrix>::reorder(const Row_index_map& valueMap,
                                                         [[maybe_unused]] Index columnIndex)
{
  static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type,
                "Method not available for chain columns.");

  Column_support newSet;

  for (Entry* entry : column_) {
    if constexpr (Master_matrix::Option_list::has_row_access) {
      RA_opt::unlink(entry);
      if (columnIndex != Master_matrix::template get_null_value<Index>()) entry->set_column_index(columnIndex);
    }
    entry->set_row_index(valueMap.at(entry->get_row_index()));
    newSet.insert(entry);
    if constexpr (Master_matrix::Option_list::has_row_access &&
                  Master_matrix::Option_list::has_intrusive_rows)  // intrusive list
      RA_opt::insert_entry(entry->get_row_index(), entry);
  }

  // when row is a set, all entries have to be deleted first, to avoid colliding when inserting
  if constexpr (Master_matrix::Option_list::has_row_access && !Master_matrix::Option_list::has_intrusive_rows) {  // set
    for (Entry* entry : newSet) {
      RA_opt::insert_entry(entry->get_row_index(), entry);
    }
  }

  column_.swap(newSet);
}

template <class Master_matrix>
inline void Unordered_set_column<Master_matrix>::clear()
{
  static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type,
                "Method not available for chain columns as a base element should not be empty.");

  for (auto* entry : column_) {
    if constexpr (Master_matrix::Option_list::has_row_access) RA_opt::unlink(entry);
    entryPool_->destroy(entry);
  }

  column_.clear();
}

template <class Master_matrix>
inline void Unordered_set_column<Master_matrix>::clear(ID_index rowIndex)
{
  static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type,
                "Method not available for chain columns.");

  auto entry = entryPool_->construct(rowIndex);
  auto it = column_.find(entry);
  if (it != column_.end()) {
    _delete_entry(it);
  }
  entryPool_->destroy(entry);
}

template <class Master_matrix>
inline typename Unordered_set_column<Master_matrix>::ID_index Unordered_set_column<Master_matrix>::get_pivot() const
{
  static_assert(Master_matrix::isNonBasic,
                "Method not available for base columns.");  // could technically be, but is the notion useful then?

  if constexpr (Master_matrix::Option_list::is_of_boundary_type) {
    if (column_.empty()) return Master_matrix::template get_null_value<ID_index>();
    // linear search could be avoided with storing the pivot. But even then, some modifications of the column requires
    // the max, so not clear how much it is worth it.
    return (*std::max_element(column_.begin(), column_.end(), EntryPointerComp()))->get_row_index();
  } else {
    return Chain_opt::_get_pivot();
  }
}

template <class Master_matrix>
inline typename Unordered_set_column<Master_matrix>::Field_element
Unordered_set_column<Master_matrix>::get_pivot_value() const
{
  static_assert(Master_matrix::isNonBasic,
                "Method not available for base columns.");  // could technically be, but is the notion useful then?

  if constexpr (Master_matrix::Option_list::is_z2) {
    return 1;
  } else {
    if constexpr (Master_matrix::Option_list::is_of_boundary_type) {
      if (column_.empty()) return 0;
      return (*std::max_element(column_.begin(), column_.end(), EntryPointerComp()))->get_element();
    } else {
      if (Chain_opt::_get_pivot() == Master_matrix::template get_null_value<ID_index>()) return Field_element();
      for (const Entry* entry : column_) {
        if (entry->get_row_index() == Chain_opt::_get_pivot()) return entry->get_element();
      }
      return Field_element();  // should never happen if chain column is used properly
    }
  }
}

template <class Master_matrix>
inline typename Unordered_set_column<Master_matrix>::iterator Unordered_set_column<Master_matrix>::begin() noexcept
{
  return column_.begin();
}

template <class Master_matrix>
inline typename Unordered_set_column<Master_matrix>::const_iterator Unordered_set_column<Master_matrix>::begin()
    const noexcept
{
  return column_.begin();
}

template <class Master_matrix>
inline typename Unordered_set_column<Master_matrix>::iterator Unordered_set_column<Master_matrix>::end() noexcept
{
  return column_.end();
}

template <class Master_matrix>
inline typename Unordered_set_column<Master_matrix>::const_iterator Unordered_set_column<Master_matrix>::end()
    const noexcept
{
  return column_.end();
}

template <class Master_matrix>
inline typename Unordered_set_column<Master_matrix>::Content_range
Unordered_set_column<Master_matrix>::get_non_zero_content_range() const
{
  Content_range res(column_.size());
  std::size_t i = 0;
  for (const auto& entry : column_) res[i++] = *entry;
  std::sort(res.begin(), res.end());
  return res;
}

template <class Master_matrix>
template <class Entry_range>
inline Unordered_set_column<Master_matrix>& Unordered_set_column<Master_matrix>::operator+=(const Entry_range& column)
{
  static_assert((!Master_matrix::isNonBasic || std::is_same_v<Entry_range, Unordered_set_column>),
                "For boundary columns, the range has to be a column of same type to help ensure the validity of the "
                "base element.");  // could be removed, if we give the responsibility to the user.
  static_assert((!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type),
                "For chain columns, the given column cannot be constant.");

  _add(column);

  return *this;
}

template <class Master_matrix>
inline Unordered_set_column<Master_matrix>& Unordered_set_column<Master_matrix>::operator+=(
    Unordered_set_column& column)
{
  if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type) {
    // assumes that the addition never zeros out this column.
    if (_add(column)) {
      Chain_opt::_swap_pivots(column);
      Dim_opt::_swap_dimension(column);
    }
  } else {
    _add(column);
  }

  return *this;
}

template <class Master_matrix>
inline Unordered_set_column<Master_matrix>& Unordered_set_column<Master_matrix>::operator*=(const Field_element& v)
{
  Field_element val = Master_matrix::get_coefficient_value(v, operators_);

  if (val == Field_operators::get_additive_identity()) {
    if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type) {
      throw std::invalid_argument("A chain column should not be multiplied by 0.");
    } else {
      clear();
    }
    return *this;
  }

  if (val == Field_operators::get_multiplicative_identity()) return *this;

  // multiply_inplace needs a non-const reference to element, so even if Z2 never reaches here, it won't compile
  // without the constexpr, as we are not storing a dummy value just for this purpose.
  if constexpr (!Master_matrix::Option_list::is_z2) {
    for (Entry* entry : column_) {
      operators_->multiply_inplace(entry->get_element(), val);
      if constexpr (Master_matrix::Option_list::has_row_access) RA_opt::update_entry(*entry);
    }
  }

  return *this;
}

template <class Master_matrix>
template <class Entry_range>
inline Unordered_set_column<Master_matrix>& Unordered_set_column<Master_matrix>::multiply_target_and_add(
    const Field_element& val,
    const Entry_range& column)
{
  static_assert((!Master_matrix::isNonBasic || std::is_same_v<Entry_range, Unordered_set_column>),
                "For boundary columns, the range has to be a column of same type to help ensure the validity of the "
                "base element.");  // could be removed, if we give the responsibility to the user.
  static_assert((!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type),
                "For chain columns, the given column cannot be constant.");

  _multiply_target_and_add(Master_matrix::get_coefficient_value(val, operators_), column);

  return *this;
}

template <class Master_matrix>
inline Unordered_set_column<Master_matrix>& Unordered_set_column<Master_matrix>::multiply_target_and_add(
    const Field_element& val,
    Unordered_set_column& column)
{
  if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type) {
    // assumes that the addition never zeros out this column.
    if (_multiply_target_and_add(Master_matrix::get_coefficient_value(val, operators_), column)) {
      Chain_opt::_swap_pivots(column);
      Dim_opt::_swap_dimension(column);
    }
  } else {
    _multiply_target_and_add(Master_matrix::get_coefficient_value(val, operators_), column);
  }

  return *this;
}

template <class Master_matrix>
template <class Entry_range>
inline Unordered_set_column<Master_matrix>& Unordered_set_column<Master_matrix>::multiply_source_and_add(
    const Entry_range& column,
    const Field_element& val)
{
  static_assert((!Master_matrix::isNonBasic || std::is_same_v<Entry_range, Unordered_set_column>),
                "For boundary columns, the range has to be a column of same type to help ensure the validity of the "
                "base element.");  // could be removed, if we give the responsibility to the user.
  static_assert((!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type),
                "For chain columns, the given column cannot be constant.");

  _multiply_source_and_add(column, Master_matrix::get_coefficient_value(val, operators_));

  return *this;
}

template <class Master_matrix>
inline Unordered_set_column<Master_matrix>& Unordered_set_column<Master_matrix>::multiply_source_and_add(
    Unordered_set_column& column,
    const Field_element& val)
{
  if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type) {
    // assumes that the addition never zeros out this column.
    if (_multiply_source_and_add(column, Master_matrix::get_coefficient_value(val, operators_))) {
      Chain_opt::_swap_pivots(column);
      Dim_opt::_swap_dimension(column);
    }
  } else {
    _multiply_source_and_add(column, Master_matrix::get_coefficient_value(val, operators_));
  }

  return *this;
}

template <class Master_matrix>
inline void Unordered_set_column<Master_matrix>::push_back(const Entry& entry)
{
  static_assert(Master_matrix::Option_list::is_of_boundary_type, "`push_back` is not available for Chain matrices.");

  GUDHI_CHECK(entry.get_row_index() > get_pivot(), "The new row index has to be higher than the current pivot.");

  _insert_entry(entry.get_row_index(), entry.get_element());
}

template <class Master_matrix>
inline Unordered_set_column<Master_matrix>& Unordered_set_column<Master_matrix>::operator=(
    const Unordered_set_column& other)
{
  static_assert(!Master_matrix::Option_list::has_row_access, "= assignment not enabled with row access option.");

  // to avoid destroying the column when building from it-self in the for loop below...
  if (this == &other) return *this;

  Dim_opt::operator=(other);
  Chain_opt::operator=(other);

  for (auto* entry : column_) {
    if constexpr (Master_matrix::Option_list::has_row_access) RA_opt::unlink(entry);
    entryPool_->destroy(entry);
  }
  column_.clear();

  operators_ = other.operators_;
  entryPool_ = other.entryPool_;

  for (const Entry* entry : other.column_) {
    _insert_entry(entry->get_row_index(), entry->get_element());
  }

  return *this;
}

template <class Master_matrix>
inline Unordered_set_column<Master_matrix>& Unordered_set_column<Master_matrix>::operator=(
    Unordered_set_column&& other) noexcept
{
  static_assert(!Master_matrix::Option_list::has_row_access, "= assignment not enabled with row access option.");

  // to avoid destroying the column before building from it-self...
  if (&column_ == &(other.column_)) return *this;

  Dim_opt::operator=(std::move(other));
  Chain_opt::operator=(std::move(other));

  for (auto* entry : column_) {
    if (entry != nullptr) entryPool_->destroy(entry);
  }

  column_ = std::move(other.column_);
  operators_ = std::exchange(other.operators_, nullptr);
  entryPool_ = std::exchange(other.entryPool_, nullptr);

  return *this;
}

template <class Master_matrix>
inline void Unordered_set_column<Master_matrix>::_delete_entry(typename Column_support::iterator& it)
{
  if constexpr (Master_matrix::Option_list::has_row_access) RA_opt::unlink(*it);
  entryPool_->destroy(*it);
  auto tmp = it++;
  // it = column_.erase(it);
  column_.erase(tmp);
}

template <class Master_matrix>
inline typename Unordered_set_column<Master_matrix>::Entry* Unordered_set_column<Master_matrix>::_insert_entry(
    ID_index rowIndex,
    const Field_element& value)
{
  Entry* newEntry;
  if constexpr (Master_matrix::Option_list::has_row_access) {
    newEntry = entryPool_->construct(RA_opt::get_column_index(), rowIndex);
  } else {
    newEntry = entryPool_->construct(rowIndex);
  }
  newEntry->set_element(value);
  column_.insert(newEntry);
  if constexpr (Master_matrix::Option_list::has_row_access) RA_opt::insert_entry(rowIndex, newEntry);
  return newEntry;
}

template <class Master_matrix>
template <class Entry_range>
inline bool Unordered_set_column<Master_matrix>::_add(const Entry_range& column)
{
  if (column.begin() == column.end()) return false;
  if (column_.empty()) {  // chain should never enter here.
    for (const Entry& entry : column) {
      _insert_entry(entry.get_row_index(), entry.get_element());
    }
    return true;
  }

  return _generic_add(
      column,
      [&](const Entry& oldEntry, Entry* newEntry) {
        newEntry->set_element(oldEntry.get_element());
      },
      [&](Entry* targetEntry, const Entry& sourceEntry) {
        if constexpr (!Master_matrix::Option_list::is_z2)
          operators_->add_inplace(targetEntry->get_element(), sourceEntry.get_element());
      });
}

template <class Master_matrix>
template <class Entry_range>
inline bool Unordered_set_column<Master_matrix>::_multiply_target_and_add(const Field_element& val,
                                                                          const Entry_range& column)
{
  // because the column is unordered, I don't see a way to do both operations in one go
  // without guarantees on the entry range...
  operator*=(val);
  return _add(column) || val == Field_operators::get_additive_identity();
}

template <class Master_matrix>
template <class Entry_range>
inline bool Unordered_set_column<Master_matrix>::_multiply_source_and_add(const Entry_range& column,
                                                                          const Field_element& val)
{
  if (val == Field_operators::get_additive_identity() || column.begin() == column.end()) {
    return false;
  }

  if (val == Field_operators::get_multiplicative_identity()) {
    return _add(column);
  }

  // multiply_inplace needs a non-const reference to element, so even if Z2 never reaches here, it won't compile
  // without the constexpr, as we are not storing a dummy value just for this purpose.
  if constexpr (!Master_matrix::Option_list::is_z2) {
    return _generic_add(
        column,
        [&](const Entry& oldEntry, Entry* newEntry) {
          newEntry->set_element(oldEntry.get_element());
          operators_->multiply_inplace(newEntry->get_element(), val);
        },
        [&](Entry* targetEntry, const Entry& sourceEntry) {
          operators_->multiply_and_add_inplace_back(sourceEntry.get_element(), val, targetEntry->get_element());
        });
  } else {
    return false;  // we should never arrive here, just to suppress the warning
  }
}

template <class Master_matrix>
template <class Entry_range, typename F1, typename F2>
inline bool Unordered_set_column<Master_matrix>::_generic_add(const Entry_range& source,
                                                              [[maybe_unused]] F1&& process_source,
                                                              [[maybe_unused]] F2&& update_target)
{
  bool pivotIsZeroed = false;

  for (const Entry& entry : source) {
    Entry* newEntry;
    if constexpr (Master_matrix::Option_list::has_row_access) {
      newEntry = entryPool_->construct(RA_opt::get_column_index(), entry.get_row_index());
    } else {
      newEntry = entryPool_->construct(entry.get_row_index());
    }
    auto res = column_.insert(newEntry);
    if (res.second) {
      std::forward<F1>(process_source)(entry, newEntry);
      if constexpr (Master_matrix::Option_list::has_row_access) RA_opt::insert_entry(entry.get_row_index(), newEntry);
    } else {
      entryPool_->destroy(newEntry);
      if constexpr (Master_matrix::Option_list::is_z2) {
        if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type) {
          if (entry.get_row_index() == Chain_opt::_get_pivot()) pivotIsZeroed = true;
        }
        _delete_entry(res.first);
      } else {
        std::forward<F2>(update_target)(*res.first, entry);
        if ((*res.first)->get_element() == Field_operators::get_additive_identity()) {
          if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type) {
            if ((*res.first)->get_row_index() == Chain_opt::_get_pivot()) pivotIsZeroed = true;
          }
          _delete_entry(res.first);
        } else {
          if constexpr (Master_matrix::Option_list::has_row_access) RA_opt::update_entry(**res.first);
        }
      }
    }
  }

  return pivotIsZeroed;
}

}  // namespace persistence_matrix
}  // namespace Gudhi

/**
 * @ingroup persistence_matrix
 *
 * @brief Hash method for @ref Gudhi::persistence_matrix::Unordered_set_column.
 *
 * @tparam Master_matrix Template parameter of @ref Gudhi::persistence_matrix::Unordered_set_column.
 * @tparam Entry_constructor Template parameter of @ref Gudhi::persistence_matrix::Unordered_set_column.
 */
template <class Master_matrix>
struct std::hash<Gudhi::persistence_matrix::Unordered_set_column<Master_matrix>> {
  std::size_t operator()(const Gudhi::persistence_matrix::Unordered_set_column<Master_matrix>& column) const
  {
    // can't use Gudhi::persistence_matrix::hash_column because unordered
    std::size_t seed = 0;
    for (const auto& entry : column) {
      seed ^= std::hash<unsigned int>()(entry.get_row_index() * static_cast<unsigned int>(entry.get_element()));
    }
    return seed;
  }
};

#endif  // PM_UNORDERED_SET_COLUMN_H

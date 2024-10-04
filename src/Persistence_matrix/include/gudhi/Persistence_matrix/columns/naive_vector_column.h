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
 * @file naive_vector_column.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_matrix::Naive_vector_column class.
 * Also defines the std::hash method for @ref Gudhi::persistence_matrix::Naive_vector_column.
 */

#ifndef PM_NAIVE_VECTOR_COLUMN_H
#define PM_NAIVE_VECTOR_COLUMN_H

#include <vector>
#include <stdexcept>
#include <type_traits>
#include <algorithm>  //binary_search
#include <utility>    //std::swap, std::move & std::exchange

#include <boost/iterator/indirect_iterator.hpp>

#include <gudhi/Persistence_matrix/allocators/entry_constructors.h>
#include <gudhi/Persistence_matrix/columns/column_utilities.h>

namespace Gudhi {
namespace persistence_matrix {

/**
 * @class Naive_vector_column naive_vector_column.h gudhi/Persistence_matrix/columns/naive_vector_column.h
 * @ingroup persistence_matrix
 *
 * @brief Column class following the @ref PersistenceMatrixColumn concept.
 *
 * Column based on a vector structure. The entries are always ordered by row index and only non-zero values
 * are stored uniquely in the underlying container.
 *
 * @tparam Master_matrix An instantiation of @ref Matrix from which all types and options are deduced.
 * @tparam Entry_constructor Factory of @ref Entry classes.
 */
template <class Master_matrix>
class Naive_vector_column : public Master_matrix::Row_access_option,
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
  using Column_support = std::vector<Entry*>;
  using Entry_constructor = typename Master_matrix::Entry_constructor;

 public:
  using iterator = boost::indirect_iterator<typename Column_support::iterator>;
  using const_iterator = boost::indirect_iterator<typename Column_support::const_iterator>;
  using reverse_iterator = boost::indirect_iterator<typename Column_support::reverse_iterator>;
  using const_reverse_iterator = boost::indirect_iterator<typename Column_support::const_reverse_iterator>;

  Naive_vector_column(Column_settings* colSettings = nullptr);
  template <class Container = typename Master_matrix::Boundary>
  Naive_vector_column(const Container& nonZeroRowIndices, Column_settings* colSettings);
  template <class Container = typename Master_matrix::Boundary, class Row_container>
  Naive_vector_column(Index columnIndex,
                      const Container& nonZeroRowIndices,
                      Row_container* rowContainer,
                      Column_settings* colSettings);
  template <class Container = typename Master_matrix::Boundary>
  Naive_vector_column(const Container& nonZeroChainRowIndices, Dimension dimension, Column_settings* colSettings);
  template <class Container = typename Master_matrix::Boundary, class Row_container>
  Naive_vector_column(Index columnIndex,
                      const Container& nonZeroChainRowIndices,
                      Dimension dimension,
                      Row_container* rowContainer,
                      Column_settings* colSettings);
  Naive_vector_column(const Naive_vector_column& column, Column_settings* colSettings = nullptr);
  template <class Row_container>
  Naive_vector_column(const Naive_vector_column& column,
                      Index columnIndex,
                      Row_container* rowContainer,
                      Column_settings* colSettings = nullptr);
  Naive_vector_column(Naive_vector_column&& column) noexcept;
  ~Naive_vector_column();

  std::vector<Field_element> get_content(int columnLength = -1) const;
  bool is_non_zero(ID_index rowIndex) const;
  bool is_empty() const;
  std::size_t size() const;

  template <class Row_index_map>
  void reorder(const Row_index_map& valueMap, [[maybe_unused]] Index columnIndex = -1);
  void clear();
  void clear(ID_index rowIndex);

  ID_index get_pivot() const;
  Field_element get_pivot_value() const;

  iterator begin() noexcept;
  const_iterator begin() const noexcept;
  iterator end() noexcept;
  const_iterator end() const noexcept;
  reverse_iterator rbegin() noexcept;
  const_reverse_iterator rbegin() const noexcept;
  reverse_iterator rend() noexcept;
  const_reverse_iterator rend() const noexcept;

  template <class Entry_range>
  Naive_vector_column& operator+=(const Entry_range& column);
  Naive_vector_column& operator+=(Naive_vector_column& column);

  Naive_vector_column& operator*=(unsigned int v);

  // this = v * this + column
  template <class Entry_range>
  Naive_vector_column& multiply_target_and_add(const Field_element& val, const Entry_range& column);
  Naive_vector_column& multiply_target_and_add(const Field_element& val, Naive_vector_column& column);
  // this = this + column * v
  template <class Entry_range>
  Naive_vector_column& multiply_source_and_add(const Entry_range& column, const Field_element& val);
  Naive_vector_column& multiply_source_and_add(Naive_vector_column& column, const Field_element& val);

  void push_back(const Entry& entry);

  friend bool operator==(const Naive_vector_column& c1, const Naive_vector_column& c2) {
    if (&c1 == &c2) return true;
    if (c1.column_.size() != c2.column_.size()) return false;

    auto it1 = c1.column_.begin();
    auto it2 = c2.column_.begin();
    while (it1 != c1.column_.end() && it2 != c2.column_.end()) {
      if constexpr (Master_matrix::Option_list::is_z2) {
        if ((*it1)->get_row_index() != (*it2)->get_row_index()) return false;
      } else {
        if ((*it1)->get_row_index() != (*it2)->get_row_index() || (*it1)->get_element() != (*it2)->get_element())
          return false;
      }
      ++it1;
      ++it2;
    }
    return true;
  }
  friend bool operator<(const Naive_vector_column& c1, const Naive_vector_column& c2) {
    if (&c1 == &c2) return false;

    auto it1 = c1.column_.begin();
    auto it2 = c2.column_.begin();
    while (it1 != c1.column_.end() && it2 != c2.column_.end()) {
      if ((*it1)->get_row_index() != (*it2)->get_row_index()) return (*it1)->get_row_index() < (*it2)->get_row_index();
      if constexpr (!Master_matrix::Option_list::is_z2) {
        if ((*it1)->get_element() != (*it2)->get_element()) return (*it1)->get_element() < (*it2)->get_element();
      }
      ++it1;
      ++it2;
    }
    return it2 != c2.column_.end();
  }

  // Disabled with row access.
  Naive_vector_column& operator=(const Naive_vector_column& other);

  friend void swap(Naive_vector_column& col1, Naive_vector_column& col2) {
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
  Field_operators* operators_;
  Entry_constructor* entryPool_;

  template <class Column, class Entry_iterator, typename F1, typename F2, typename F3, typename F4>
  friend void _generic_merge_entry_to_column(Column& targetColumn,
                                             Entry_iterator& itSource,
                                             typename Column::Column_support::iterator& itTarget,
                                             F1&& process_target,
                                             F2&& process_source,
                                             F3&& update_target1,
                                             F4&& update_target2,
                                             bool& pivotIsZeroed);
  template <class Column, class Entry_range, typename F1, typename F2, typename F3, typename F4, typename F5>
  friend bool _generic_add_to_column(const Entry_range& source,
                                     Column& targetColumn,
                                     F1&& process_target,
                                     F2&& process_source,
                                     F3&& update_target1,
                                     F4&& update_target2,
                                     F5&& finish_target);

  void _delete_entry(Entry* entry);
  void _delete_entry(typename Column_support::iterator& it);
  Entry* _insert_entry(const Field_element& value, ID_index rowIndex, Column_support& column);
  void _insert_entry(ID_index rowIndex, Column_support& column);
  void _update_entry(const Field_element& value, ID_index rowIndex, Index position);
  void _update_entry(ID_index rowIndex, Index position);
  template <class Entry_range>
  bool _add(const Entry_range& column);
  template <class Entry_range>
  bool _multiply_target_and_add(const Field_element& val, const Entry_range& column);
  template <class Entry_range>
  bool _multiply_source_and_add(const Entry_range& column, const Field_element& val);
};

template <class Master_matrix>
inline Naive_vector_column<Master_matrix>::Naive_vector_column(Column_settings* colSettings)
    : RA_opt(),
      Dim_opt(),
      Chain_opt(),
      operators_(nullptr),
      entryPool_(colSettings == nullptr ? nullptr : &(colSettings->entryConstructor))
{
  if (operators_ == nullptr && entryPool_ == nullptr) return; // to allow default constructor which gives a dummy column
  if constexpr (!Master_matrix::Option_list::is_z2) {
    operators_ = &(colSettings->operators);
  }
}

template <class Master_matrix>
template <class Container>
inline Naive_vector_column<Master_matrix>::Naive_vector_column(const Container& nonZeroRowIndices,
                                                               Column_settings* colSettings)
    : RA_opt(),
      Dim_opt(nonZeroRowIndices.size() == 0 ? 0 : nonZeroRowIndices.size() - 1),
      Chain_opt(),
      column_(nonZeroRowIndices.size(), nullptr),
      operators_(nullptr),
      entryPool_(&(colSettings->entryConstructor))
{
  static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type,
                "Constructor not available for chain columns, please specify the dimension of the chain.");

  Index i = 0;
  if constexpr (Master_matrix::Option_list::is_z2) {
    for (ID_index id : nonZeroRowIndices) {
      _update_entry(id, i++);
    }
  } else {
    operators_ = &(colSettings->operators);
    for (const auto& p : nonZeroRowIndices) {
      _update_entry(operators_->get_value(p.second), p.first, i++);
    }
  }
}

template <class Master_matrix>
template <class Container, class Row_container>
inline Naive_vector_column<Master_matrix>::Naive_vector_column(Index columnIndex,
                                                               const Container& nonZeroRowIndices,
                                                               Row_container* rowContainer,
                                                               Column_settings* colSettings)
    : RA_opt(columnIndex, rowContainer),
      Dim_opt(nonZeroRowIndices.size() == 0 ? 0 : nonZeroRowIndices.size() - 1),
      Chain_opt([&] {
        if constexpr (Master_matrix::Option_list::is_z2) {
          return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : *std::prev(nonZeroRowIndices.end());
        } else {
          return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : std::prev(nonZeroRowIndices.end())->first;
        }
      }()),
      column_(nonZeroRowIndices.size(), nullptr),
      operators_(nullptr),
      entryPool_(&(colSettings->entryConstructor))
{
  static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type,
                "Constructor not available for chain columns, please specify the dimension of the chain.");

  Index i = 0;
  if constexpr (Master_matrix::Option_list::is_z2) {
    for (ID_index id : nonZeroRowIndices) {
      _update_entry(id, i++);
    }
  } else {
    operators_ = &(colSettings->operators);
    for (const auto& p : nonZeroRowIndices) {
      _update_entry(operators_->get_value(p.second), p.first, i++);
    }
  }
}

template <class Master_matrix>
template <class Container>
inline Naive_vector_column<Master_matrix>::Naive_vector_column(const Container& nonZeroRowIndices,
                                                               Dimension dimension,
                                                               Column_settings* colSettings)
    : RA_opt(),
      Dim_opt(dimension),
      Chain_opt([&] {
        if constexpr (Master_matrix::Option_list::is_z2) {
          return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : *std::prev(nonZeroRowIndices.end());
        } else {
          return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : std::prev(nonZeroRowIndices.end())->first;
        }
      }()),
      column_(nonZeroRowIndices.size(), nullptr),
      operators_(nullptr),
      entryPool_(&(colSettings->entryConstructor))
{
  Index i = 0;
  if constexpr (Master_matrix::Option_list::is_z2) {
    for (ID_index id : nonZeroRowIndices) {
      _update_entry(id, i++);
    }
  } else {
    operators_ = &(colSettings->operators);
    for (const auto& p : nonZeroRowIndices) {
      _update_entry(operators_->get_value(p.second), p.first, i++);
    }
  }
}

template <class Master_matrix>
template <class Container, class Row_container>
inline Naive_vector_column<Master_matrix>::Naive_vector_column(Index columnIndex,
                                                               const Container& nonZeroRowIndices,
                                                               Dimension dimension,
                                                               Row_container* rowContainer,
                                                               Column_settings* colSettings)
    : RA_opt(columnIndex, rowContainer),
      Dim_opt(dimension),
      Chain_opt([&] {
        if constexpr (Master_matrix::Option_list::is_z2) {
          return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : *std::prev(nonZeroRowIndices.end());
        } else {
          return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : std::prev(nonZeroRowIndices.end())->first;
        }
      }()),
      column_(nonZeroRowIndices.size(), nullptr),
      operators_(nullptr),
      entryPool_(&(colSettings->entryConstructor))
{
  Index i = 0;
  if constexpr (Master_matrix::Option_list::is_z2) {
    for (ID_index id : nonZeroRowIndices) {
      _update_entry(id, i++);
    }
  } else {
    operators_ = &(colSettings->operators);
    for (const auto& p : nonZeroRowIndices) {
      _update_entry(operators_->get_value(p.second), p.first, i++);
    }
  }
}

template <class Master_matrix>
inline Naive_vector_column<Master_matrix>::Naive_vector_column(const Naive_vector_column& column,
                                                               Column_settings* colSettings)
    : RA_opt(),
      Dim_opt(static_cast<const Dim_opt&>(column)),
      Chain_opt(static_cast<const Chain_opt&>(column)),
      column_(column.column_.size(), nullptr),
      operators_(colSettings == nullptr ? column.operators_ : nullptr),
      entryPool_(colSettings == nullptr ? column.entryPool_ : &(colSettings->entryConstructor))
{
  static_assert(!Master_matrix::Option_list::has_row_access,
                "Simple copy constructor not available when row access option enabled. Please specify the new column "
                "index and the row container.");

  if constexpr (!Master_matrix::Option_list::is_z2) {
    if (colSettings != nullptr) operators_ = &(colSettings->operators);
  }

  Index i = 0;
  for (const Entry* entry : column.column_) {
    if constexpr (Master_matrix::Option_list::is_z2) {
      _update_entry(entry->get_row_index(), i++);
    } else {
      _update_entry(entry->get_element(), entry->get_row_index(), i++);
    }
  }
}

template <class Master_matrix>
template <class Row_container>
inline Naive_vector_column<Master_matrix>::Naive_vector_column(const Naive_vector_column& column,
                                                               Index columnIndex,
                                                               Row_container* rowContainer,
                                                               Column_settings* colSettings)
    : RA_opt(columnIndex, rowContainer),
      Dim_opt(static_cast<const Dim_opt&>(column)),
      Chain_opt(static_cast<const Chain_opt&>(column)),
      column_(column.column_.size(), nullptr),
      operators_(colSettings == nullptr ? column.operators_ : nullptr),
      entryPool_(colSettings == nullptr ? column.entryPool_ : &(colSettings->entryConstructor))
{
  if constexpr (!Master_matrix::Option_list::is_z2) {
    if (colSettings != nullptr) operators_ = &(colSettings->operators);
  }

  Index i = 0;
  for (const Entry* entry : column.column_) {
    if constexpr (Master_matrix::Option_list::is_z2) {
      _update_entry(entry->get_row_index(), i++);
    } else {
      _update_entry(entry->get_element(), entry->get_row_index(), i++);
    }
  }
}

template <class Master_matrix>
inline Naive_vector_column<Master_matrix>::Naive_vector_column(Naive_vector_column&& column) noexcept
    : RA_opt(std::move(static_cast<RA_opt&>(column))),
      Dim_opt(std::move(static_cast<Dim_opt&>(column))),
      Chain_opt(std::move(static_cast<Chain_opt&>(column))),
      column_(std::move(column.column_)),
      operators_(std::exchange(column.operators_, nullptr)),
      entryPool_(std::exchange(column.entryPool_, nullptr))
{}

template <class Master_matrix>
inline Naive_vector_column<Master_matrix>::~Naive_vector_column()
{
  for (auto* entry : column_) {
    _delete_entry(entry);
  }
}

template <class Master_matrix>
inline std::vector<typename Naive_vector_column<Master_matrix>::Field_element>
Naive_vector_column<Master_matrix>::get_content(int columnLength) const
{
  if (columnLength < 0 && column_.size() > 0)
    columnLength = column_.back()->get_row_index() + 1;
  else if (columnLength < 0)
    return std::vector<Field_element>();

  std::vector<Field_element> container(columnLength, 0);
  for (auto it = column_.begin(); it != column_.end() && (*it)->get_row_index() < static_cast<ID_index>(columnLength);
       ++it) {
    if constexpr (Master_matrix::Option_list::is_z2) {
      container[(*it)->get_row_index()] = 1;
    } else {
      container[(*it)->get_row_index()] = (*it)->get_element();
    }
  }
  return container;
}

template <class Master_matrix>
inline bool Naive_vector_column<Master_matrix>::is_non_zero(ID_index rowIndex) const
{
  Entry entry(rowIndex);
  return std::binary_search(column_.begin(), column_.end(), &entry,
                            [](const Entry* a, const Entry* b) { return a->get_row_index() < b->get_row_index(); });
}

template <class Master_matrix>
inline bool Naive_vector_column<Master_matrix>::is_empty() const
{
  return column_.empty();
}

template <class Master_matrix>
inline std::size_t Naive_vector_column<Master_matrix>::size() const
{
  return column_.size();
}

template <class Master_matrix>
template <class Row_index_map>
inline void Naive_vector_column<Master_matrix>::reorder(const Row_index_map& valueMap,
                                                        [[maybe_unused]] Index columnIndex)
{
  static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type,
                "Method not available for chain columns.");

  for (Entry* entry : column_) {
    if constexpr (Master_matrix::Option_list::has_row_access) {
      RA_opt::unlink(entry);
      if (columnIndex != static_cast<Index>(-1)) entry->set_column_index(columnIndex);
    }
    entry->set_row_index(valueMap.at(entry->get_row_index()));
    if constexpr (Master_matrix::Option_list::has_intrusive_rows && Master_matrix::Option_list::has_row_access)
      RA_opt::insert_entry(entry->get_row_index(), entry);
  }

  // all entries have to be deleted first, to avoid problem with insertion when row is a set
  if constexpr (!Master_matrix::Option_list::has_intrusive_rows && Master_matrix::Option_list::has_row_access) {
    for (Entry* entry : column_) {
      RA_opt::insert_entry(entry->get_row_index(), entry);
    }
  }

  std::sort(column_.begin(), column_.end(), [](const Entry* c1, const Entry* c2) { return *c1 < *c2; });
}

template <class Master_matrix>
inline void Naive_vector_column<Master_matrix>::clear()
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
inline void Naive_vector_column<Master_matrix>::clear(ID_index rowIndex)
{
  static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type,
                "Method not available for chain columns.");

  auto it = column_.begin();
  while (it != column_.end() && (*it)->get_row_index() != rowIndex) ++it;
  if (it != column_.end()) {
    _delete_entry(*it);
    column_.erase(it);
  }
}

template <class Master_matrix>
inline typename Naive_vector_column<Master_matrix>::ID_index Naive_vector_column<Master_matrix>::get_pivot() const
{
  static_assert(Master_matrix::isNonBasic,
                "Method not available for base columns.");  // could technically be, but is the notion useful then?

  if constexpr (Master_matrix::Option_list::is_of_boundary_type) {
    return column_.empty() ? -1 : column_.back()->get_row_index();
  } else {
    return Chain_opt::get_pivot();
  }
}

template <class Master_matrix>
inline typename Naive_vector_column<Master_matrix>::Field_element Naive_vector_column<Master_matrix>::get_pivot_value()
    const
{
  static_assert(Master_matrix::isNonBasic,
                "Method not available for base columns.");  // could technically be, but is the notion useful then?

  if constexpr (Master_matrix::Option_list::is_z2) {
    return 1;
  } else {
    if constexpr (Master_matrix::Option_list::is_of_boundary_type) {
      return column_.empty() ? Field_element() : column_.back()->get_element();
    } else {
      if (Chain_opt::get_pivot() == static_cast<ID_index>(-1)) return Field_element();
      for (const Entry* entry : column_) {
        if (entry->get_row_index() == Chain_opt::get_pivot()) return entry->get_element();
      }
      return Field_element();  // should never happen if chain column is used properly
    }
  }
}

template <class Master_matrix>
inline typename Naive_vector_column<Master_matrix>::iterator Naive_vector_column<Master_matrix>::begin() noexcept
{
  return column_.begin();
}

template <class Master_matrix>
inline typename Naive_vector_column<Master_matrix>::const_iterator Naive_vector_column<Master_matrix>::begin()
    const noexcept
{
  return column_.begin();
}

template <class Master_matrix>
inline typename Naive_vector_column<Master_matrix>::iterator Naive_vector_column<Master_matrix>::end() noexcept
{
  return column_.end();
}

template <class Master_matrix>
inline typename Naive_vector_column<Master_matrix>::const_iterator Naive_vector_column<Master_matrix>::end()
    const noexcept
{
  return column_.end();
}

template <class Master_matrix>
inline typename Naive_vector_column<Master_matrix>::reverse_iterator
Naive_vector_column<Master_matrix>::rbegin() noexcept
{
  return column_.rbegin();
}

template <class Master_matrix>
inline typename Naive_vector_column<Master_matrix>::const_reverse_iterator Naive_vector_column<Master_matrix>::rbegin()
    const noexcept
{
  return column_.rbegin();
}

template <class Master_matrix>
inline typename Naive_vector_column<Master_matrix>::reverse_iterator
Naive_vector_column<Master_matrix>::rend() noexcept
{
  return column_.rend();
}

template <class Master_matrix>
inline typename Naive_vector_column<Master_matrix>::const_reverse_iterator Naive_vector_column<Master_matrix>::rend()
    const noexcept
{
  return column_.rend();
}

template <class Master_matrix>
template <class Entry_range>
inline Naive_vector_column<Master_matrix>& Naive_vector_column<Master_matrix>::operator+=(const Entry_range& column)
{
  static_assert((!Master_matrix::isNonBasic || std::is_same_v<Entry_range, Naive_vector_column>),
                "For boundary columns, the range has to be a column of same type to help ensure the validity of the "
                "base element.");  // could be removed, if we give the responsibility to the user.
  static_assert((!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type),
                "For chain columns, the given column cannot be constant.");

  _add(column);

  return *this;
}

template <class Master_matrix>
inline Naive_vector_column<Master_matrix>& Naive_vector_column<Master_matrix>::operator+=(Naive_vector_column& column)
{
  if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type) {
    // assumes that the addition never zeros out this column.
    if (_add(column)) {
      Chain_opt::swap_pivots(column);
      Dim_opt::swap_dimension(column);
    }
  } else {
    _add(column);
  }

  return *this;
}

template <class Master_matrix>
inline Naive_vector_column<Master_matrix>& Naive_vector_column<Master_matrix>::operator*=(unsigned int v)
{
  if constexpr (Master_matrix::Option_list::is_z2) {
    if (v % 2 == 0) {
      if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type) {
        throw std::invalid_argument("A chain column should not be multiplied by 0.");
      } else {
        clear();
      }
    }
  } else {
    Field_element val = operators_->get_value(v);

    if (val == Field_operators::get_additive_identity()) {
      if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type) {
        throw std::invalid_argument("A chain column should not be multiplied by 0.");
      } else {
        clear();
      }
      return *this;
    }

    if (val == Field_operators::get_multiplicative_identity()) return *this;

    for (Entry* entry : column_) {
      operators_->multiply_inplace(entry->get_element(), val);
      if constexpr (Master_matrix::Option_list::has_row_access) RA_opt::update_entry(*entry);
    }
  }

  return *this;
}

template <class Master_matrix>
template <class Entry_range>
inline Naive_vector_column<Master_matrix>& Naive_vector_column<Master_matrix>::multiply_target_and_add(
    const Field_element& val, const Entry_range& column)
{
  static_assert((!Master_matrix::isNonBasic || std::is_same_v<Entry_range, Naive_vector_column>),
                "For boundary columns, the range has to be a column of same type to help ensure the validity of the "
                "base element.");  // could be removed, if we give the responsibility to the user.
  static_assert((!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type),
                "For chain columns, the given column cannot be constant.");

  if constexpr (Master_matrix::Option_list::is_z2) {
    if (val) {
      _add(column);
    } else {
      clear();
      _add(column);
    }
  } else {
    _multiply_target_and_add(val, column);
  }

  return *this;
}

template <class Master_matrix>
inline Naive_vector_column<Master_matrix>& Naive_vector_column<Master_matrix>::multiply_target_and_add(
    const Field_element& val, Naive_vector_column& column)
{
  if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type) {
    // assumes that the addition never zeros out this column.
    if constexpr (Master_matrix::Option_list::is_z2) {
      if (val) {
        if (_add(column)) {
          Chain_opt::swap_pivots(column);
          Dim_opt::swap_dimension(column);
        }
      } else {
        throw std::invalid_argument("A chain column should not be multiplied by 0.");
      }
    } else {
      if (_multiply_target_and_add(val, column)) {
        Chain_opt::swap_pivots(column);
        Dim_opt::swap_dimension(column);
      }
    }
  } else {
    if constexpr (Master_matrix::Option_list::is_z2) {
      if (val) {
        _add(column);
      } else {
        clear();
        _add(column);
      }
    } else {
      _multiply_target_and_add(val, column);
    }
  }

  return *this;
}

template <class Master_matrix>
template <class Entry_range>
inline Naive_vector_column<Master_matrix>& Naive_vector_column<Master_matrix>::multiply_source_and_add(
    const Entry_range& column, const Field_element& val)
{
  static_assert((!Master_matrix::isNonBasic || std::is_same_v<Entry_range, Naive_vector_column>),
                "For boundary columns, the range has to be a column of same type to help ensure the validity of the "
                "base element.");  // could be removed, if we give the responsibility to the user.
  static_assert((!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type),
                "For chain columns, the given column cannot be constant.");

  if constexpr (Master_matrix::Option_list::is_z2) {
    if (val) {
      _add(column);
    }
  } else {
    _multiply_source_and_add(column, val);
  }

  return *this;
}

template <class Master_matrix>
inline Naive_vector_column<Master_matrix>& Naive_vector_column<Master_matrix>::multiply_source_and_add(
    Naive_vector_column& column, const Field_element& val)
{
  if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type) {
    // assumes that the addition never zeros out this column.
    if constexpr (Master_matrix::Option_list::is_z2) {
      if (val) {
        if (_add(column)) {
          Chain_opt::swap_pivots(column);
          Dim_opt::swap_dimension(column);
        }
      }
    } else {
      if (_multiply_source_and_add(column, val)) {
        Chain_opt::swap_pivots(column);
        Dim_opt::swap_dimension(column);
      }
    }
  } else {
    if constexpr (Master_matrix::Option_list::is_z2) {
      if (val) {
        _add(column);
      }
    } else {
      _multiply_source_and_add(column, val);
    }
  }

  return *this;
}

template <class Master_matrix>
inline void Naive_vector_column<Master_matrix>::push_back(const Entry& entry)
{
  static_assert(Master_matrix::Option_list::is_of_boundary_type, "`push_back` is not available for Chain matrices.");

  if constexpr (Master_matrix::Option_list::is_z2) {
    _insert_entry(entry.get_row_index(), column_);
  } else {
    _insert_entry(entry.get_element(), entry.get_row_index(), column_);
  }
}

template <class Master_matrix>
inline Naive_vector_column<Master_matrix>& Naive_vector_column<Master_matrix>::operator=(
    const Naive_vector_column& other)
{
  static_assert(!Master_matrix::Option_list::has_row_access, "= assignment not enabled with row access option.");

  Dim_opt::operator=(other);
  Chain_opt::operator=(other);

  auto tmpPool = entryPool_;
  entryPool_ = other.entryPool_;

  while (column_.size() > other.column_.size()) {
    if (column_.back() == nullptr) {
      if constexpr (Master_matrix::Option_list::has_row_access) RA_opt::unlink(column_.back());
      tmpPool->destroy(column_.back());
    }
    column_.pop_back();
  }

  column_.resize(other.column_.size(), nullptr);
  Index i = 0;
  for (const Entry* entry : other.column_) {
    if (column_[i] != nullptr) {
      if constexpr (Master_matrix::Option_list::has_row_access) RA_opt::unlink(column_[i]);
      tmpPool->destroy(column_[i]);
    }
    if constexpr (Master_matrix::Option_list::is_z2) {
      _update_entry(entry->get_row_index(), i++);
    } else {
      _update_entry(entry->get_element(), entry->get_row_index(), i++);
    }
  }

  operators_ = other.operators_;

  return *this;
}

template <class Master_matrix>
inline void Naive_vector_column<Master_matrix>::_delete_entry(Entry* entry)
{
  if constexpr (Master_matrix::Option_list::has_row_access) RA_opt::unlink(entry);
  entryPool_->destroy(entry);
}

template <class Master_matrix>
inline void Naive_vector_column<Master_matrix>::_delete_entry(typename Column_support::iterator& it)
{
  _delete_entry(*it);
  ++it;
}

template <class Master_matrix>
inline typename Naive_vector_column<Master_matrix>::Entry* Naive_vector_column<Master_matrix>::_insert_entry(
    const Field_element& value, ID_index rowIndex, Column_support& column)
{
  if constexpr (Master_matrix::Option_list::has_row_access) {
    Entry* newEntry = entryPool_->construct(RA_opt::columnIndex_, rowIndex);
    newEntry->set_element(value);
    column.push_back(newEntry);
    RA_opt::insert_entry(rowIndex, newEntry);
    return newEntry;
  } else {
    Entry* newEntry = entryPool_->construct(rowIndex);
    column.push_back(newEntry);
    newEntry->set_element(value);
    return newEntry;
  }
}

template <class Master_matrix>
inline void Naive_vector_column<Master_matrix>::_insert_entry(ID_index rowIndex, Column_support& column)
{
  if constexpr (Master_matrix::Option_list::has_row_access) {
    Entry* newEntry = entryPool_->construct(RA_opt::columnIndex_, rowIndex);
    column.push_back(newEntry);
    RA_opt::insert_entry(rowIndex, newEntry);
  } else {
    column.push_back(entryPool_->construct(rowIndex));
  }
}

template <class Master_matrix>
inline void Naive_vector_column<Master_matrix>::_update_entry(const Field_element& value,
                                                             ID_index rowIndex,
                                                             Index position)
{
  if constexpr (Master_matrix::Option_list::has_row_access) {
    Entry* newEntry = entryPool_->construct(RA_opt::columnIndex_, rowIndex);
    newEntry->set_element(value);
    column_[position] = newEntry;
    RA_opt::insert_entry(rowIndex, newEntry);
  } else {
    column_[position] = entryPool_->construct(rowIndex);
    column_[position]->set_element(value);
  }
}

template <class Master_matrix>
inline void Naive_vector_column<Master_matrix>::_update_entry(ID_index rowIndex, Index position)
{
  if constexpr (Master_matrix::Option_list::has_row_access) {
    Entry* newEntry = entryPool_->construct(RA_opt::columnIndex_, rowIndex);
    column_[position] = newEntry;
    RA_opt::insert_entry(rowIndex, newEntry);
  } else {
    column_[position] = entryPool_->construct(rowIndex);
  }
}

template <class Master_matrix>
template <class Entry_range>
inline bool Naive_vector_column<Master_matrix>::_add(const Entry_range& column)
{
  if (column.begin() == column.end()) return false;
  if (column_.empty()) {  // chain should never enter here.
    column_.resize(column.size());
    Index i = 0;
    for (const Entry& entry : column) {
      if constexpr (Master_matrix::Option_list::is_z2) {
        _update_entry(entry.get_row_index(), i++);
      } else {
        _update_entry(entry.get_element(), entry.get_row_index(), i++);
      }
    }
    return true;
  }

  Column_support newColumn;
  newColumn.reserve(column_.size() + column.size());  // safe upper bound

  auto pivotIsZeroed = _generic_add_to_column(
      column, 
      *this, 
      [&](Entry* entryTarget) { newColumn.push_back(entryTarget); },
      [&](typename Entry_range::const_iterator& itSource,
          [[maybe_unused]] const typename Column_support::iterator& itTarget) {
        if constexpr (Master_matrix::Option_list::is_z2) {
          _insert_entry(itSource->get_row_index(), newColumn);
        } else {
          _insert_entry(itSource->get_element(), itSource->get_row_index(), newColumn);
        }
      },
      [&](Field_element& targetElement, typename Entry_range::const_iterator& itSource) {
        if constexpr (!Master_matrix::Option_list::is_z2)
          operators_->add_inplace(targetElement, itSource->get_element());
      },
      [&](Entry* entryTarget) { newColumn.push_back(entryTarget); },
      [&](typename Column_support::iterator& itTarget) {
        while (itTarget != column_.end()) {
          newColumn.push_back(*itTarget);
          itTarget++;
        }
      });

  column_.swap(newColumn);

  return pivotIsZeroed;
}

template <class Master_matrix>
template <class Entry_range>
inline bool Naive_vector_column<Master_matrix>::_multiply_target_and_add(const Field_element& val,
                                                                         const Entry_range& column)
{
  if (val == 0u) {
    if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type) {
      throw std::invalid_argument("A chain column should not be multiplied by 0.");
      // this would not only mess up the base, but also the pivots stored.
    } else {
      clear();
    }
  }
  if (column_.empty()) {  // chain should never enter here.
    column_.resize(column.size());
    Index i = 0;
    for (const Entry& entry : column) {
      if constexpr (Master_matrix::Option_list::is_z2) {
        _update_entry(entry.get_row_index(), i++);
      } else {
        _update_entry(entry.get_element(), entry.get_row_index(), i++);
      }
    }
    return true;
  }

  Column_support newColumn;
  newColumn.reserve(column_.size() + column.size());  // safe upper bound

  auto pivotIsZeroed = _generic_add_to_column(
      column,
      *this,
      [&](Entry* entryTarget) {
        operators_->multiply_inplace(entryTarget->get_element(), val);
        if constexpr (Master_matrix::Option_list::has_row_access) RA_opt::update_entry(*entryTarget);
        newColumn.push_back(entryTarget);
      },
      [&](typename Entry_range::const_iterator& itSource, const typename Column_support::iterator& itTarget) {
        _insert_entry(itSource->get_element(), itSource->get_row_index(), newColumn);
      },
      [&](Field_element& targetElement, typename Entry_range::const_iterator& itSource) {
        operators_->multiply_and_add_inplace_front(targetElement, val, itSource->get_element());
      },
      [&](Entry* entryTarget) { newColumn.push_back(entryTarget); },
      [&](typename Column_support::iterator& itTarget) {
        while (itTarget != column_.end()) {
          operators_->multiply_inplace((*itTarget)->get_element(), val);
          if constexpr (Master_matrix::Option_list::has_row_access) RA_opt::update_entry(**itTarget);
          newColumn.push_back(*itTarget);
          itTarget++;
        }
      });

  column_.swap(newColumn);

  return pivotIsZeroed;
}

template <class Master_matrix>
template <class Entry_range>
inline bool Naive_vector_column<Master_matrix>::_multiply_source_and_add(const Entry_range& column,
                                                                         const Field_element& val)
{
  if (val == 0u || column.begin() == column.end()) {
    return false;
  }

  Column_support newColumn;
  newColumn.reserve(column_.size() + column.size());  // safe upper bound

  auto pivotIsZeroed = _generic_add_to_column(
      column,
      *this,
      [&](Entry* entryTarget) { newColumn.push_back(entryTarget); },
      [&](typename Entry_range::const_iterator& itSource, const typename Column_support::iterator& itTarget) {
        Entry* newEntry = _insert_entry(itSource->get_element(), itSource->get_row_index(), newColumn);
        operators_->multiply_inplace(newEntry->get_element(), val);
        if constexpr (Master_matrix::Option_list::has_row_access) RA_opt::update_entry(*newEntry);
      },
      [&](Field_element& targetElement, typename Entry_range::const_iterator& itSource) {
        operators_->multiply_and_add_inplace_back(itSource->get_element(), val, targetElement);
      },
      [&](Entry* entryTarget) { newColumn.push_back(entryTarget); },
      [&](typename Column_support::iterator& itTarget) {
        while (itTarget != column_.end()) {
          newColumn.push_back(*itTarget);
          itTarget++;
        }
      });

  column_.swap(newColumn);

  return pivotIsZeroed;
}

}  // namespace persistence_matrix
}  // namespace Gudhi

/**
 * @ingroup persistence_matrix
 *
 * @brief Hash method for @ref Gudhi::persistence_matrix::Naive_vector_column.
 *
 * @tparam Master_matrix Template parameter of @ref Gudhi::persistence_matrix::Naive_vector_column.
 * @tparam Entry_constructor Template parameter of @ref Gudhi::persistence_matrix::Naive_vector_column.
 */
template <class Master_matrix>
struct std::hash<Gudhi::persistence_matrix::Naive_vector_column<Master_matrix> > {
  std::size_t operator()(const Gudhi::persistence_matrix::Naive_vector_column<Master_matrix>& column) const {
    return Gudhi::persistence_matrix::hash_column(column);
  }
};

#endif  // PM_NAIVE_VECTOR_COLUMN_H

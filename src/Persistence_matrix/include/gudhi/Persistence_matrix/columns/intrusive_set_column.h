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
 * @file intrusive_set_column.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_matrix::Intrusive_set_column class.
 * Also defines the std::hash method for @ref Gudhi::persistence_matrix::Intrusive_set_column.
 */

#ifndef PM_INTRUSIVE_SET_COLUMN_H
#define PM_INTRUSIVE_SET_COLUMN_H

#include <vector>
#include <stdexcept>
#include <type_traits>
#include <utility>  //std::swap, std::move & std::exchange

#include <boost/intrusive/set.hpp>

#include <gudhi/Debug_utils.h>
#include <gudhi/Persistence_matrix/allocators/entry_constructors.h>
#include <gudhi/Persistence_matrix/columns/column_utilities.h>

namespace Gudhi {
namespace persistence_matrix {

/**
 * @class Intrusive_set_column intrusive_set_column.h gudhi/Persistence_matrix/columns/intrusive_set_column.h
 * @ingroup persistence_matrix
 *
 * @brief Column class following the @ref PersistenceMatrixColumn concept.
 *
 * Column based on a intrusive set structure. The entries are ordered by row index and only non-zero values
 * are stored uniquely in the underlying container.
 *
 * @tparam Master_matrix An instantiation of @ref Matrix from which all types and options are deduced.
 */
template <class Master_matrix>
class Intrusive_set_column : public Master_matrix::Row_access_option,
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
  using Column_support =
      boost::intrusive::set<Entry,
                            boost::intrusive::constant_time_size<false>,
                            boost::intrusive::base_hook<typename Master_matrix::Base_hook_matrix_set_column> >;
  using Entry_constructor = typename Master_matrix::Entry_constructor;

 public:
  using iterator = typename Column_support::iterator;
  using const_iterator = typename Column_support::const_iterator;
  using reverse_iterator = typename Column_support::reverse_iterator;
  using const_reverse_iterator = typename Column_support::const_reverse_iterator;
  using Content_range = const Column_support&;

  Intrusive_set_column(Column_settings* colSettings = nullptr);
  template <class Container = typename Master_matrix::Boundary>
  Intrusive_set_column(const Container& nonZeroRowIndices, Column_settings* colSettings);
  template <class Container = typename Master_matrix::Boundary, class Row_container>
  Intrusive_set_column(Index columnIndex,
                       const Container& nonZeroRowIndices,
                       Row_container* rowContainer,
                       Column_settings* colSettings);
  template <class Container = typename Master_matrix::Boundary,
            class = std::enable_if_t<!std::is_arithmetic_v<Container> > >
  Intrusive_set_column(const Container& nonZeroRowIndices, Dimension dimension, Column_settings* colSettings);
  template <class Container = typename Master_matrix::Boundary,
            class Row_container,
            class = std::enable_if_t<!std::is_arithmetic_v<Container> > >
  Intrusive_set_column(Index columnIndex,
                       const Container& nonZeroRowIndices,
                       Dimension dimension,
                       Row_container* rowContainer,
                       Column_settings* colSettings);
  Intrusive_set_column(ID_index idx, Dimension dimension, Column_settings* colSettings);
  Intrusive_set_column(ID_index idx,
                       Field_element e,
                       Dimension dimension,
                       Column_settings* colSettings);
  template <class Row_container>
  Intrusive_set_column(Index columnIndex,
                       ID_index idx,
                       Dimension dimension,
                       Row_container* rowContainer,
                       Column_settings* colSettings);
  template <class Row_container>
  Intrusive_set_column(Index columnIndex,
                       ID_index idx,
                       Field_element e,
                       Dimension dimension,
                       Row_container* rowContainer,
                       Column_settings* colSettings);
  Intrusive_set_column(const Intrusive_set_column& column, Column_settings* colSettings = nullptr);
  template <class Row_container>
  Intrusive_set_column(const Intrusive_set_column& column,
                       Index columnIndex,
                       Row_container* rowContainer,
                       Column_settings* colSettings = nullptr);
  Intrusive_set_column(Intrusive_set_column&& column) noexcept;
  ~Intrusive_set_column();

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
  reverse_iterator rbegin() noexcept;
  const_reverse_iterator rbegin() const noexcept;
  reverse_iterator rend() noexcept;
  const_reverse_iterator rend() const noexcept;

  Content_range get_non_zero_content_range() const;

  template <class Entry_range>
  Intrusive_set_column& operator+=(const Entry_range& column);
  Intrusive_set_column& operator+=(Intrusive_set_column& column);

  Intrusive_set_column& operator*=(const Field_element& v);

  // this = v * this + column
  template <class Entry_range>
  Intrusive_set_column& multiply_target_and_add(const Field_element& val, const Entry_range& column);
  Intrusive_set_column& multiply_target_and_add(const Field_element& val, Intrusive_set_column& column);
  // this = this + column * v
  template <class Entry_range>
  Intrusive_set_column& multiply_source_and_add(const Entry_range& column, const Field_element& val);
  Intrusive_set_column& multiply_source_and_add(Intrusive_set_column& column, const Field_element& val);

  void push_back(const Entry& entry);

  friend bool operator==(const Intrusive_set_column& c1, const Intrusive_set_column& c2)
  {
    if (&c1 == &c2) return true;

    return std::equal(c1.column_.begin(),
                      c1.column_.end(),
                      c2.column_.begin(),
                      c2.column_.end(),
                      [](const Entry& e1, const Entry& e2) {
                        return e1.get_row_index() == e2.get_row_index() && e1.get_element() == e2.get_element();
                      });
  }

  friend bool operator<(const Intrusive_set_column& c1, const Intrusive_set_column& c2)
  {
    if (&c1 == &c2) return false;

    return std::lexicographical_compare(c1.column_.begin(),
                                        c1.column_.end(),
                                        c2.column_.begin(),
                                        c2.column_.end(),
                                        [](const Entry& e1, const Entry& e2) {
                                          if (e1.get_row_index() != e2.get_row_index())
                                            return e1.get_row_index() < e2.get_row_index();
                                          if (e1.get_element() != e2.get_element())
                                            return e1.get_element() < e2.get_element();
                                          return false;
                                        });
  }

  // Disabled with row access.
  Intrusive_set_column& operator=(const Intrusive_set_column& other);
  Intrusive_set_column& operator=(Intrusive_set_column&& other) noexcept;

  friend void swap(Intrusive_set_column& col1, Intrusive_set_column& col2) noexcept
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

  // Cloner object function for boost intrusive container
  struct New_cloner {
    New_cloner(Entry_constructor* entryPool) : entryPool_(entryPool) {};

    Entry* operator()(const Entry& clone_this) { return entryPool_->construct(clone_this); }

    Entry_constructor* entryPool_;
  };

  // The disposer object function for boost intrusive container
  struct Delete_disposer {
    Delete_disposer() = default;
    Delete_disposer(Intrusive_set_column* col) : col_(col) {};

    void operator()(Entry* delete_this)
    {
      if constexpr (Master_matrix::Option_list::has_row_access) col_->unlink(delete_this);
      col_->entryPool_->destroy(delete_this);
    }

    Intrusive_set_column* col_;
  };

  Column_support column_;
  Field_operators const* operators_;
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
  template <class Column, class Entry_range>
  friend bool _add_to_column(const Entry_range& source, Column& targetColumn);
  template <class Column, class Entry_range>
  friend bool _multiply_target_and_add_to_column(const typename Column::Field_element& val,
                                                 const Entry_range& source,
                                                 Column& targetColumn);
  template <class Column, class Entry_range>
  friend bool _multiply_source_and_add_to_column(const typename Column::Field_element& val,
                                                 const Entry_range& source,
                                                 Column& targetColumn);

  void _delete_entry(iterator& it);
  Entry* _insert_entry(const iterator& position, ID_index rowIndex, const Field_element& value);
  template <class Entry_range>
  bool _add(const Entry_range& column);
  template <class Entry_range>
  bool _multiply_target_and_add(const Field_element& val, const Entry_range& column);
  template <class Entry_range>
  bool _multiply_source_and_add(const Entry_range& column, const Field_element& val);
};

template <class Master_matrix>
inline Intrusive_set_column<Master_matrix>::Intrusive_set_column(Column_settings* colSettings)
    : RA_opt(),
      Dim_opt(),
      Chain_opt(),
      operators_(Master_matrix::get_operator_ptr(colSettings)),
      entryPool_(colSettings == nullptr ? nullptr : &(colSettings->entryConstructor))
{}

template <class Master_matrix>
template <class Container>
inline Intrusive_set_column<Master_matrix>::Intrusive_set_column(const Container& nonZeroRowIndices,
                                                                 Column_settings* colSettings)
    : Intrusive_set_column(nonZeroRowIndices,
                           nonZeroRowIndices.size() == 0 ? 0 : nonZeroRowIndices.size() - 1,
                           colSettings)
{
  static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type,
                "Constructor not available for chain columns, please specify the dimension of the chain.");
}

template <class Master_matrix>
template <class Container, class Row_container>
inline Intrusive_set_column<Master_matrix>::Intrusive_set_column(Index columnIndex,
                                                                 const Container& nonZeroRowIndices,
                                                                 Row_container* rowContainer,
                                                                 Column_settings* colSettings)
    : Intrusive_set_column(columnIndex,
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
inline Intrusive_set_column<Master_matrix>::Intrusive_set_column(const Container& nonZeroRowIndices,
                                                                 Dimension dimension,
                                                                 Column_settings* colSettings)
    : RA_opt(),
      Dim_opt(dimension),
      Chain_opt(nonZeroRowIndices.begin() == nonZeroRowIndices.end()
                    ? Master_matrix::template get_null_value<ID_index>()
                    : Master_matrix::get_row_index(*std::prev(nonZeroRowIndices.end()))),
      operators_(Master_matrix::get_operator_ptr(colSettings)),
      entryPool_(&(colSettings->entryConstructor))
{
  for (const auto& id : nonZeroRowIndices) {
    _insert_entry(column_.end(),
                  Master_matrix::get_row_index(id),
                  Master_matrix::get_coefficient_value(Master_matrix::get_element(id), operators_));
  }
}

template <class Master_matrix>
template <class Container, class Row_container, class>
inline Intrusive_set_column<Master_matrix>::Intrusive_set_column(Index columnIndex,
                                                                 const Container& nonZeroRowIndices,
                                                                 Dimension dimension,
                                                                 Row_container* rowContainer,
                                                                 Column_settings* colSettings)
    : RA_opt(columnIndex, rowContainer),
      Dim_opt(dimension),
      Chain_opt(nonZeroRowIndices.begin() == nonZeroRowIndices.end()
                    ? Master_matrix::template get_null_value<ID_index>()
                    : Master_matrix::get_row_index(*std::prev(nonZeroRowIndices.end()))),
      operators_(Master_matrix::get_operator_ptr(colSettings)),
      entryPool_(&(colSettings->entryConstructor))
{
  for (const auto& id : nonZeroRowIndices) {
    _insert_entry(column_.end(),
                  Master_matrix::get_row_index(id),
                  Master_matrix::get_coefficient_value(Master_matrix::get_element(id), operators_));
  }
}

template <class Master_matrix>
inline Intrusive_set_column<Master_matrix>::Intrusive_set_column(ID_index idx,
                                                                 Dimension dimension,
                                                                 Column_settings* colSettings)
    : RA_opt(), Dim_opt(dimension), Chain_opt(idx), operators_(nullptr), entryPool_(&(colSettings->entryConstructor))
{
  static_assert(Master_matrix::Option_list::is_z2,
                "Constructor not available for Zp != Z2. Please specify the coefficient.");
  _insert_entry(column_.end(), idx, 1);
}

template <class Master_matrix>
inline Intrusive_set_column<Master_matrix>::Intrusive_set_column(ID_index idx,
                                                                 Field_element e,
                                                                 Dimension dimension,
                                                                 Column_settings* colSettings)
    : RA_opt(),
      Dim_opt(dimension),
      Chain_opt(idx),
      operators_(&(colSettings->operators)),
      entryPool_(&(colSettings->entryConstructor))
{
  static_assert(!Master_matrix::Option_list::is_z2,
                "Constructor not available for Zp == Z2. Please do not specify any coefficient.");
  _insert_entry(column_.end(), idx, operators_->get_value(e));
}

template <class Master_matrix>
template <class Row_container>
inline Intrusive_set_column<Master_matrix>::Intrusive_set_column(Index columnIndex,
                                                                 ID_index idx,
                                                                 Dimension dimension,
                                                                 Row_container* rowContainer,
                                                                 Column_settings* colSettings)
    : RA_opt(columnIndex, rowContainer),
      Dim_opt(dimension),
      Chain_opt(idx),
      operators_(nullptr),
      entryPool_(&(colSettings->entryConstructor))
{
  static_assert(Master_matrix::Option_list::is_z2,
                "Constructor not available for Zp != Z2. Please specify the coefficient.");
  _insert_entry(column_.end(), idx, 1);
}

template <class Master_matrix>
template <class Row_container>
inline Intrusive_set_column<Master_matrix>::Intrusive_set_column(Index columnIndex,
                                                                 ID_index idx,
                                                                 Field_element e,
                                                                 Dimension dimension,
                                                                 Row_container* rowContainer,
                                                                 Column_settings* colSettings)
    : RA_opt(columnIndex, rowContainer),
      Dim_opt(dimension),
      Chain_opt(idx),
      operators_(&(colSettings->operators)),
      entryPool_(&(colSettings->entryConstructor))
{
  static_assert(!Master_matrix::Option_list::is_z2,
                "Constructor not available for Zp == Z2. Please do not specify any coefficient.");
  _insert_entry(column_.end(), idx, operators_->get_value(e));
}

template <class Master_matrix>
inline Intrusive_set_column<Master_matrix>::Intrusive_set_column(const Intrusive_set_column& column,
                                                                 Column_settings* colSettings)
    : RA_opt(),
      Dim_opt(static_cast<const Dim_opt&>(column)),
      Chain_opt(static_cast<const Chain_opt&>(column)),
      operators_(colSettings == nullptr ? column.operators_ : Master_matrix::get_operator_ptr(colSettings)),
      entryPool_(colSettings == nullptr ? column.entryPool_ : &(colSettings->entryConstructor))
{
  static_assert(!Master_matrix::Option_list::has_row_access,
                "Simple copy constructor not available when row access option enabled. Please specify the new column "
                "index and the row container.");

  column_.clone_from(column.column_, New_cloner(entryPool_), Delete_disposer(this));
}

template <class Master_matrix>
template <class Row_container>
inline Intrusive_set_column<Master_matrix>::Intrusive_set_column(const Intrusive_set_column& column,
                                                                 Index columnIndex,
                                                                 Row_container* rowContainer,
                                                                 Column_settings* colSettings)
    : RA_opt(columnIndex, rowContainer),
      Dim_opt(static_cast<const Dim_opt&>(column)),
      Chain_opt(static_cast<const Chain_opt&>(column)),
      operators_(colSettings == nullptr ? column.operators_ : Master_matrix::get_operator_ptr(colSettings)),
      entryPool_(colSettings == nullptr ? column.entryPool_ : &(colSettings->entryConstructor))
{
  for (const Entry& entry : column.column_) {
    _insert_entry(column_.end(), entry.get_row_index(), entry.get_element());
  }
}

template <class Master_matrix>
inline Intrusive_set_column<Master_matrix>::Intrusive_set_column(Intrusive_set_column&& column) noexcept
    : RA_opt(std::move(static_cast<RA_opt&>(column))),
      Dim_opt(std::move(static_cast<Dim_opt&>(column))),
      Chain_opt(std::move(static_cast<Chain_opt&>(column))),
      column_(std::move(column.column_)),
      operators_(std::exchange(column.operators_, nullptr)),
      entryPool_(std::exchange(column.entryPool_, nullptr))
{}

template <class Master_matrix>
inline Intrusive_set_column<Master_matrix>::~Intrusive_set_column()
{
  column_.clear_and_dispose(Delete_disposer(this));
}

template <class Master_matrix>
inline std::vector<typename Intrusive_set_column<Master_matrix>::Field_element>
Intrusive_set_column<Master_matrix>::get_content(int columnLength) const
{
  if (columnLength < 0 && column_.size() > 0)
    columnLength = column_.rbegin()->get_row_index() + 1;
  else if (columnLength < 0)
    return std::vector<Field_element>();

  std::vector<Field_element> container(columnLength);
  for (auto it = column_.begin(); it != column_.end() && it->get_row_index() < static_cast<ID_index>(columnLength);
       ++it) {
    container[it->get_row_index()] = Master_matrix::get_element(*it);
  }
  return container;
}

template <class Master_matrix>
inline bool Intrusive_set_column<Master_matrix>::is_non_zero(ID_index rowIndex) const
{
  return column_.find(Entry(rowIndex)) != column_.end();
}

template <class Master_matrix>
inline bool Intrusive_set_column<Master_matrix>::is_empty() const
{
  return column_.empty();
}

template <class Master_matrix>
inline std::size_t Intrusive_set_column<Master_matrix>::size() const
{
  return column_.size();
}

template <class Master_matrix>
template <class Row_index_map>
inline void Intrusive_set_column<Master_matrix>::reorder(const Row_index_map& valueMap,
                                                         [[maybe_unused]] Index columnIndex)
{
  static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type,
                "Method not available for chain columns.");

  Column_support newSet;

  if constexpr (Master_matrix::Option_list::has_row_access) {
    for (auto it = column_.begin(); it != column_.end();) {
      Entry* newEntry = entryPool_->construct(
          columnIndex == Master_matrix::template get_null_value<Index>() ? RA_opt::get_column_index() : columnIndex,
          valueMap.at(it->get_row_index()));
      newEntry->set_element(it->get_element());
      newSet.insert(newSet.end(), *newEntry);
      _delete_entry(it);                                             // increases it
      if constexpr (Master_matrix::Option_list::has_intrusive_rows)  // intrusive list
        RA_opt::insert_entry(newEntry->get_row_index(), newEntry);
    }

    // when row is a set, all entries have to be deleted first, to avoid colliding when inserting
    if constexpr (!Master_matrix::Option_list::has_intrusive_rows) {  // set
      for (Entry& entry : newSet) {
        RA_opt::insert_entry(entry.get_row_index(), &entry);
      }
    }
  } else {
    for (auto it = column_.begin(); it != column_.end();) {
      Entry* newEntry = entryPool_->construct(valueMap.at(it->get_row_index()));
      newEntry->set_element(it->get_element());
      newSet.insert(newSet.end(), *newEntry);
      _delete_entry(it);  // increases it
    }
  }

  column_.swap(newSet);
}

template <class Master_matrix>
inline void Intrusive_set_column<Master_matrix>::clear()
{
  static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type,
                "Method not available for chain columns as a base element should not be empty.");

  column_.clear_and_dispose(Delete_disposer(this));
}

template <class Master_matrix>
inline void Intrusive_set_column<Master_matrix>::clear(ID_index rowIndex)
{
  static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type,
                "Method not available for chain columns.");

  auto it = column_.find(Entry(rowIndex));
  if (it != column_.end()) {
    _delete_entry(it);
  }
}

template <class Master_matrix>
inline typename Intrusive_set_column<Master_matrix>::ID_index Intrusive_set_column<Master_matrix>::get_pivot() const
{
  static_assert(Master_matrix::isNonBasic,
                "Method not available for base columns.");  // could technically be, but is the notion useful then?

  if constexpr (Master_matrix::Option_list::is_of_boundary_type) {
    if (column_.empty()) return Master_matrix::template get_null_value<ID_index>();
    return column_.rbegin()->get_row_index();
  } else {
    return Chain_opt::_get_pivot();
  }
}

template <class Master_matrix>
inline typename Intrusive_set_column<Master_matrix>::Field_element
Intrusive_set_column<Master_matrix>::get_pivot_value() const
{
  static_assert(Master_matrix::isNonBasic,
                "Method not available for base columns.");  // could technically be, but is the notion useful then?

  if constexpr (Master_matrix::Option_list::is_z2) {
    return 1;
  } else {
    if constexpr (Master_matrix::Option_list::is_of_boundary_type) {
      if (column_.empty()) return 0;
      return column_.rbegin()->get_element();
    } else {
      if (Chain_opt::_get_pivot() == Master_matrix::template get_null_value<ID_index>()) return 0;
      auto it = column_.find(Entry(Chain_opt::_get_pivot()));
      GUDHI_CHECK(it != column_.end(),
                  "Intrusive_set_column::get_pivot_value - Pivot not found only if the column was misused.");
      return it->get_element();
    }
  }
}

template <class Master_matrix>
inline typename Intrusive_set_column<Master_matrix>::iterator Intrusive_set_column<Master_matrix>::begin() noexcept
{
  return column_.begin();
}

template <class Master_matrix>
inline typename Intrusive_set_column<Master_matrix>::const_iterator Intrusive_set_column<Master_matrix>::begin()
    const noexcept
{
  return column_.begin();
}

template <class Master_matrix>
inline typename Intrusive_set_column<Master_matrix>::iterator Intrusive_set_column<Master_matrix>::end() noexcept
{
  return column_.end();
}

template <class Master_matrix>
inline typename Intrusive_set_column<Master_matrix>::const_iterator Intrusive_set_column<Master_matrix>::end()
    const noexcept
{
  return column_.end();
}

template <class Master_matrix>
inline typename Intrusive_set_column<Master_matrix>::reverse_iterator
Intrusive_set_column<Master_matrix>::rbegin() noexcept
{
  return column_.rbegin();
}

template <class Master_matrix>
inline typename Intrusive_set_column<Master_matrix>::const_reverse_iterator
Intrusive_set_column<Master_matrix>::rbegin() const noexcept
{
  return column_.rbegin();
}

template <class Master_matrix>
inline typename Intrusive_set_column<Master_matrix>::reverse_iterator
Intrusive_set_column<Master_matrix>::rend() noexcept
{
  return column_.rend();
}

template <class Master_matrix>
inline typename Intrusive_set_column<Master_matrix>::const_reverse_iterator Intrusive_set_column<Master_matrix>::rend()
    const noexcept
{
  return column_.rend();
}

template <class Master_matrix>
inline typename Intrusive_set_column<Master_matrix>::Content_range
Intrusive_set_column<Master_matrix>::get_non_zero_content_range() const
{
  return column_;
}

template <class Master_matrix>
template <class Entry_range>
inline Intrusive_set_column<Master_matrix>& Intrusive_set_column<Master_matrix>::operator+=(const Entry_range& column)
{
  static_assert((!Master_matrix::isNonBasic || std::is_same_v<Entry_range, Intrusive_set_column>),
                "For boundary columns, the range has to be a column of same type to help ensure the validity of the "
                "base element.");  // could be removed, if we give the responsibility to the user.
  static_assert((!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type),
                "For chain columns, the given column cannot be constant.");

  _add(column);

  return *this;
}

template <class Master_matrix>
inline Intrusive_set_column<Master_matrix>& Intrusive_set_column<Master_matrix>::operator+=(
    Intrusive_set_column& column)
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
inline Intrusive_set_column<Master_matrix>& Intrusive_set_column<Master_matrix>::operator*=(const Field_element& v)
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
    for (Entry& entry : column_) {
      operators_->multiply_inplace(entry.get_element(), val);
      if constexpr (Master_matrix::Option_list::has_row_access) RA_opt::update_entry(entry);
    }
  }

  return *this;
}

template <class Master_matrix>
template <class Entry_range>
inline Intrusive_set_column<Master_matrix>& Intrusive_set_column<Master_matrix>::multiply_target_and_add(
    const Field_element& val,
    const Entry_range& column)
{
  static_assert((!Master_matrix::isNonBasic || std::is_same_v<Entry_range, Intrusive_set_column>),
                "For boundary columns, the range has to be a column of same type to help ensure the validity of the "
                "base element.");  // could be removed, if we give the responsibility to the user.
  static_assert((!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type),
                "For chain columns, the given column cannot be constant.");

  _multiply_target_and_add(Master_matrix::get_coefficient_value(val, operators_), column);

  return *this;
}

template <class Master_matrix>
inline Intrusive_set_column<Master_matrix>& Intrusive_set_column<Master_matrix>::multiply_target_and_add(
    const Field_element& val,
    Intrusive_set_column& column)
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
inline Intrusive_set_column<Master_matrix>& Intrusive_set_column<Master_matrix>::multiply_source_and_add(
    const Entry_range& column,
    const Field_element& val)
{
  static_assert((!Master_matrix::isNonBasic || std::is_same_v<Entry_range, Intrusive_set_column>),
                "For boundary columns, the range has to be a column of same type to help ensure the validity of the "
                "base element.");  // could be removed, if we give the responsibility to the user.
  static_assert((!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type),
                "For chain columns, the given column cannot be constant.");

  _multiply_source_and_add(column, Master_matrix::get_coefficient_value(val, operators_));

  return *this;
}

template <class Master_matrix>
inline Intrusive_set_column<Master_matrix>& Intrusive_set_column<Master_matrix>::multiply_source_and_add(
    Intrusive_set_column& column,
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
inline void Intrusive_set_column<Master_matrix>::push_back(const Entry& entry)
{
  static_assert(Master_matrix::Option_list::is_of_boundary_type, "`push_back` is not available for Chain matrices.");

  GUDHI_CHECK(entry.get_row_index() > get_pivot(), "The new row index has to be higher than the current pivot.");

  _insert_entry(column_.end(), entry.get_row_index(), entry.get_element());
}

template <class Master_matrix>
inline Intrusive_set_column<Master_matrix>& Intrusive_set_column<Master_matrix>::operator=(
    const Intrusive_set_column& other)
{
  static_assert(!Master_matrix::Option_list::has_row_access, "= assignment not enabled with row access option.");

  // otherwise the column will be destroyed before copying itself...
  if (this == &other) return *this;

  Dim_opt::operator=(other);
  Chain_opt::operator=(other);

  // order is important
  column_.clear_and_dispose(Delete_disposer(this));
  operators_ = other.operators_;
  entryPool_ = other.entryPool_;
  column_.clone_from(other.column_, New_cloner(entryPool_), Delete_disposer(this));

  return *this;
}

template <class Master_matrix>
inline Intrusive_set_column<Master_matrix>& Intrusive_set_column<Master_matrix>::operator=(
    Intrusive_set_column&& other) noexcept
{
  static_assert(!Master_matrix::Option_list::has_row_access, "= assignment not enabled with row access option.");

  // to avoid destroying the column before building from it-self...
  if (&column_ == &(other.column_)) return *this;

  Dim_opt::operator=(std::move(other));
  Chain_opt::operator=(std::move(other));

  column_.clear_and_dispose(Delete_disposer(this));

  operators_ = std::exchange(other.operators_, nullptr);
  entryPool_ = std::exchange(other.entryPool_, nullptr);
  column_ = std::move(other.column_);

  return *this;
}

template <class Master_matrix>
inline void Intrusive_set_column<Master_matrix>::_delete_entry(iterator& it)
{
  it = column_.erase_and_dispose(it, Delete_disposer(this));
}

template <class Master_matrix>
inline typename Intrusive_set_column<Master_matrix>::Entry* Intrusive_set_column<Master_matrix>::_insert_entry(
    const iterator& position,
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
  column_.insert(position, *newEntry);
  if constexpr (Master_matrix::Option_list::has_row_access) RA_opt::insert_entry(rowIndex, newEntry);
  return newEntry;
}

template <class Master_matrix>
template <class Entry_range>
inline bool Intrusive_set_column<Master_matrix>::_add(const Entry_range& column)
{
  return _add_to_column(column, *this);
}

template <class Master_matrix>
template <class Entry_range>
inline bool Intrusive_set_column<Master_matrix>::_multiply_target_and_add(const Field_element& val,
                                                                          const Entry_range& column)
{
  return _multiply_target_and_add_to_column(val, column, *this);
}

template <class Master_matrix>
template <class Entry_range>
inline bool Intrusive_set_column<Master_matrix>::_multiply_source_and_add(const Entry_range& column,
                                                                          const Field_element& val)
{
  return _multiply_source_and_add_to_column(val, column, *this);
}

}  // namespace persistence_matrix
}  // namespace Gudhi

/**
 * @ingroup persistence_matrix
 *
 * @brief Hash method for @ref Gudhi::persistence_matrix::Intrusive_set_column.
 *
 * @tparam Master_matrix Template parameter of @ref Gudhi::persistence_matrix::Intrusive_set_column.
 * @tparam Entry_constructor Template parameter of @ref Gudhi::persistence_matrix::Intrusive_set_column.
 */
template <class Master_matrix>
struct std::hash<Gudhi::persistence_matrix::Intrusive_set_column<Master_matrix> > {
  std::size_t operator()(const Gudhi::persistence_matrix::Intrusive_set_column<Master_matrix>& column) const
  {
    return Gudhi::persistence_matrix::hash_column(column);
  }
};

#endif  // PM_INTRUSIVE_SET_COLUMN_H

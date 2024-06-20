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
 * @file vector_column.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Vector_column class.
 * Also defines the std::hash method for @ref Vector_column.
 */

#ifndef PM_VECTOR_COLUMN_H
#define PM_VECTOR_COLUMN_H

#include <cstddef>
#include <vector>
#include <stdexcept>
#include <type_traits>
#include <algorithm>      //binary_search
#include <unordered_set>
#include <utility>        //std::swap, std::move & std::exchange

#include <boost/iterator/indirect_iterator.hpp>
#include <gudhi/Persistence_matrix/columns/column_utilities.h>

namespace Gudhi {
namespace persistence_matrix {

/**
 * @class Vector_column vector_column.h gudhi/Persistence_matrix/columns/vector_column.h
 * @ingroup persistence_matrix
 *
 * @brief Column class following the @ref PersistenceMatrixColumn concept.
 *
 * Column based on a vector structure. The cells are always ordered by row index, but cells are removed by 
 * @ref PersistenceMatrixColumn::clear(PersistenceMatrixOptions::index_type rowIndex) "clear(index)" in a lazy way,
 * so erased values can still be in the underlying container.
 * On the other hand, two cells will never have the same row index.
 * 
 * @tparam Master_matrix An instanciation of @ref Matrix from which all types and options are deduced.
 * @tparam Cell_constructor Factory of @ref Cell classes.
 */
template <class Master_matrix>
class Vector_column : public Master_matrix::Row_access_option,
                      public Master_matrix::Column_dimension_option,
                      public Master_matrix::Chain_column_option 
{
 public:
  using Master = Master_matrix;
  using index = typename Master_matrix::index;
  using id_index = typename Master_matrix::id_index;
  using dimension_type = typename Master_matrix::dimension_type;
  using Field_element_type = typename Master_matrix::element_type;
  using Cell = typename Master_matrix::Cell_type;
  using Column_settings = typename Master_matrix::Column_settings;

 private:
  using Field_operators = typename Master_matrix::Field_operators;
  using Column_type = std::vector<Cell*>;
  using Cell_constructor = typename Master_matrix::Cell_constructor;

 public:
  using iterator = boost::indirect_iterator<typename Column_type::iterator>;
  using const_iterator = boost::indirect_iterator<typename Column_type::const_iterator>;
  using reverse_iterator = boost::indirect_iterator<typename Column_type::reverse_iterator>;
  using const_reverse_iterator = boost::indirect_iterator<typename Column_type::const_reverse_iterator>;

  Vector_column(Column_settings* colSettings = nullptr);
  template <class Container_type = typename Master_matrix::boundary_type>
  Vector_column(const Container_type& nonZeroRowIndices, Column_settings* colSettings);
  template <class Container_type = typename Master_matrix::boundary_type, class Row_container_type>
  Vector_column(index columnIndex, 
                const Container_type& nonZeroRowIndices, 
                Row_container_type* rowContainer,
                Column_settings* colSettings);
  template <class Container_type = typename Master_matrix::boundary_type>
  Vector_column(const Container_type& nonZeroChainRowIndices, 
                dimension_type dimension, 
                Column_settings* colSettings);
  template <class Container_type = typename Master_matrix::boundary_type, class Row_container_type>
  Vector_column(index columnIndex, 
                const Container_type& nonZeroChainRowIndices, 
                dimension_type dimension,
                Row_container_type* rowContainer, 
                Column_settings* colSettings);
  Vector_column(const Vector_column& column, 
                Column_settings* colSettings = nullptr);
  template <class Row_container_type>
  Vector_column(const Vector_column& column, 
                index columnIndex, 
                Row_container_type* rowContainer,
                Column_settings* colSettings = nullptr);
  Vector_column(Vector_column&& column) noexcept;
  ~Vector_column();

  std::vector<Field_element_type> get_content(int columnLength = -1) const;
  bool is_non_zero(id_index rowIndex) const;
  bool is_empty() const;
  std::size_t size() const;

  template <class Map_type>
  void reorder(const Map_type& valueMap, [[maybe_unused]] index columnIndex = -1);
  void clear();
  // do not clear a cell to 0 if the cell was already 0, otherwise size/is_empty will be wrong.
  void clear(id_index rowIndex);

  id_index get_pivot();
  Field_element_type get_pivot_value();

  iterator begin() noexcept;
  const_iterator begin() const noexcept;
  iterator end() noexcept;
  const_iterator end() const noexcept;
  reverse_iterator rbegin() noexcept;
  const_reverse_iterator rbegin() const noexcept;
  reverse_iterator rend() noexcept;
  const_reverse_iterator rend() const noexcept;

  template <class Cell_range>
  Vector_column& operator+=(const Cell_range& column);
  Vector_column& operator+=(Vector_column& column);

  Vector_column& operator*=(unsigned int v);

  // this = v * this + column
  template <class Cell_range>
  Vector_column& multiply_target_and_add(const Field_element_type& val, const Cell_range& column);
  Vector_column& multiply_target_and_add(const Field_element_type& val, Vector_column& column);
  // this = this + column * v
  template <class Cell_range>
  Vector_column& multiply_source_and_add(const Cell_range& column, const Field_element_type& val);
  Vector_column& multiply_source_and_add(Vector_column& column, const Field_element_type& val);

  std::size_t compute_hash_value();

  friend bool operator==(const Vector_column& c1, const Vector_column& c2) {
    if (&c1 == &c2) return true;
    if (c1.erasedValues_.empty() && c2.erasedValues_.empty() && c1.column_.size() != c2.column_.size()) return false;

    auto it1 = c1.column_.begin();
    auto it2 = c2.column_.begin();
    while (it1 != c1.column_.end() && it2 != c2.column_.end()) {
      if constexpr (!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type) {
        while (it1 != c1.column_.end() && c1.erasedValues_.find((*it1)->get_row_index()) != c1.erasedValues_.end())
          ++it1;
        while (it2 != c2.column_.end() && c2.erasedValues_.find((*it2)->get_row_index()) != c2.erasedValues_.end())
          ++it2;
        if (it1 == c1.column_.end() || it2 == c2.column_.end()) break;
      }
      if constexpr (Master_matrix::Option_list::is_z2) {
        if ((*it1)->get_row_index() != (*it2)->get_row_index()) return false;
      } else {
        if ((*it1)->get_row_index() != (*it2)->get_row_index() || (*it1)->get_element() != (*it2)->get_element())
          return false;
      }
      ++it1;
      ++it2;
    }

    if constexpr (!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type) {
      while (it1 != c1.column_.end() && c1.erasedValues_.find((*it1)->get_row_index()) != c1.erasedValues_.end()) ++it1;
      while (it2 != c2.column_.end() && c2.erasedValues_.find((*it2)->get_row_index()) != c2.erasedValues_.end()) ++it2;
      return it2 == c2.column_.end() && it1 == c1.column_.end();
    } else {
      return true;
    }
  }
  friend bool operator<(const Vector_column& c1, const Vector_column& c2) {
    if (&c1 == &c2) return false;

    auto it1 = c1.column_.begin();
    auto it2 = c2.column_.begin();
    while (it1 != c1.column_.end() && it2 != c2.column_.end()) {
      if constexpr (!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type) {
        while (it1 != c1.column_.end() && c1.erasedValues_.find((*it1)->get_row_index()) != c1.erasedValues_.end())
          ++it1;
        while (it2 != c2.column_.end() && c2.erasedValues_.find((*it2)->get_row_index()) != c2.erasedValues_.end())
          ++it2;
        if (it1 == c1.column_.end() || it2 == c2.column_.end()) break;
      }

      if ((*it1)->get_row_index() != (*it2)->get_row_index()) return (*it1)->get_row_index() < (*it2)->get_row_index();
      if constexpr (!Master_matrix::Option_list::is_z2) {
        if ((*it1)->get_element() != (*it2)->get_element()) return (*it1)->get_element() < (*it2)->get_element();
      }
      ++it1;
      ++it2;
    }
    if constexpr (!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type) {
      while (it1 != c1.column_.end() && c1.erasedValues_.find((*it1)->get_row_index()) != c1.erasedValues_.end()) ++it1;
      while (it2 != c2.column_.end() && c2.erasedValues_.find((*it2)->get_row_index()) != c2.erasedValues_.end()) ++it2;
    }
    return it2 != c2.column_.end();
  }

  // Disabled with row access.
  Vector_column& operator=(const Vector_column& other);

  friend void swap(Vector_column& col1, Vector_column& col2) {
    swap(static_cast<typename Master_matrix::Row_access_option&>(col1),
         static_cast<typename Master_matrix::Row_access_option&>(col2));
    swap(static_cast<typename Master_matrix::Column_dimension_option&>(col1),
         static_cast<typename Master_matrix::Column_dimension_option&>(col2));
    swap(static_cast<typename Master_matrix::Chain_column_option&>(col1),
         static_cast<typename Master_matrix::Chain_column_option&>(col2));
    col1.column_.swap(col2.column_);
    col1.erasedValues_.swap(col2.erasedValues_);
    std::swap(col1.operators_, col2.operators_);
    std::swap(col1.cellPool_, col2.cellPool_);
  }

 private:
  using ra_opt = typename Master_matrix::Row_access_option;
  using dim_opt = typename Master_matrix::Column_dimension_option;
  using chain_opt = typename Master_matrix::Chain_column_option;

  Column_type column_;
  std::unordered_set<id_index> erasedValues_;  // TODO: test other containers? Useless when clear(index) is never
                                               // called, how much is it worth it?
  Field_operators* operators_;
  Cell_constructor* cellPool_;

  template <class Column_type, class Cell_iterator, typename F1, typename F2, typename F3, typename F4>
  friend void _generic_merge_cell_to_column(Column_type& targetColumn,
                                            Cell_iterator& itSource,
                                            typename Column_type::Column_type::iterator& itTarget,
                                            F1&& process_target,
                                            F2&& process_source,
                                            F3&& update_target1,
                                            F4&& update_target2,
                                            bool& pivotIsZeroed);

  void _delete_cell(Cell* cell);
  void _delete_cell(typename Column_type::iterator& it);
  Cell* _insert_cell(const Field_element_type& value, id_index rowIndex, Column_type& column);
  void _insert_cell(id_index rowIndex, Column_type& column);
  void _update_cell(const Field_element_type& value, id_index rowIndex, index position);
  void _update_cell(id_index rowIndex, index position);
  template <class Cell_range>
  bool _add(const Cell_range& column);
  template <class Cell_range>
  bool _multiply_target_and_add(const Field_element_type& val, const Cell_range& column);
  template <class Cell_range>
  bool _multiply_source_and_add(const Cell_range& column, const Field_element_type& val);
  template <class Cell_range, typename F1, typename F2, typename F3, typename F4>
  bool _generic_add(const Cell_range& source,
                    F1&& process_target,
                    F2&& process_source,
                    F3&& update_target1,
                    F4&& update_target2);
};

template <class Master_matrix>
inline Vector_column<Master_matrix>::Vector_column(Column_settings* colSettings)
    : ra_opt(), dim_opt(), chain_opt(), operators_(nullptr), cellPool_(colSettings == nullptr ? nullptr : &(colSettings->cellConstructor))
{
  if (operators_ == nullptr && cellPool_ == nullptr) return;  //to allow default constructor which gives a dummy column
  if constexpr (!Master_matrix::Option_list::is_z2){
    operators_ = &(colSettings->operators);
  }
}

template <class Master_matrix>
template <class Container_type>
inline Vector_column<Master_matrix>::Vector_column(const Container_type& nonZeroRowIndices,
                                                                     Column_settings* colSettings)
    : ra_opt(),
      dim_opt(nonZeroRowIndices.size() == 0 ? 0 : nonZeroRowIndices.size() - 1),
      chain_opt(),
      column_(nonZeroRowIndices.size(), nullptr),
      operators_(nullptr),
      cellPool_(&(colSettings->cellConstructor))
{
  static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type,
                "Constructor not available for chain columns, please specify the dimension of the chain.");

  if constexpr (!Master_matrix::Option_list::is_z2){
    operators_ = &(colSettings->operators);
  }

  index i = 0;
  if constexpr (Master_matrix::Option_list::is_z2) {
    for (id_index id : nonZeroRowIndices) {
      _update_cell(id, i++);
    }
  } else {
    for (const auto& p : nonZeroRowIndices) {
      _update_cell(operators_->get_value(p.second), p.first, i++);
    }
  }
}

template <class Master_matrix>
template <class Container_type, class Row_container_type>
inline Vector_column<Master_matrix>::Vector_column(index columnIndex,
                                                                     const Container_type& nonZeroRowIndices,
                                                                     Row_container_type* rowContainer,
                                                                     Column_settings* colSettings)
    : ra_opt(columnIndex, rowContainer),
      dim_opt(nonZeroRowIndices.size() == 0 ? 0 : nonZeroRowIndices.size() - 1),
      chain_opt([&] {
        if constexpr (Master_matrix::Option_list::is_z2) {
          return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : *std::prev(nonZeroRowIndices.end());
        } else {
          return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : std::prev(nonZeroRowIndices.end())->first;
        }
      }()),
      column_(nonZeroRowIndices.size(), nullptr),
      operators_(nullptr),
      cellPool_(&(colSettings->cellConstructor))
{
  static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type,
                "Constructor not available for chain columns, please specify the dimension of the chain.");

  if constexpr (!Master_matrix::Option_list::is_z2){
    operators_ = &(colSettings->operators);
  }

  index i = 0;
  if constexpr (Master_matrix::Option_list::is_z2) {
    for (id_index id : nonZeroRowIndices) {
      _update_cell(id, i++);
    }
  } else {
    for (const auto& p : nonZeroRowIndices) {
      _update_cell(operators_->get_value(p.second), p.first, i++);
    }
  }
}

template <class Master_matrix>
template <class Container_type>
inline Vector_column<Master_matrix>::Vector_column(const Container_type& nonZeroRowIndices,
                                                                     dimension_type dimension,
                                                                     Column_settings* colSettings)
    : ra_opt(),
      dim_opt(dimension),
      chain_opt([&] {
        if constexpr (Master_matrix::Option_list::is_z2) {
          return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : *std::prev(nonZeroRowIndices.end());
        } else {
          return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : std::prev(nonZeroRowIndices.end())->first;
        }
      }()),
      column_(nonZeroRowIndices.size(), nullptr),
      operators_(nullptr),
      cellPool_(&(colSettings->cellConstructor))
{
  if constexpr (!Master_matrix::Option_list::is_z2){
    operators_ = &(colSettings->operators);
  }

  index i = 0;
  if constexpr (Master_matrix::Option_list::is_z2) {
    for (id_index id : nonZeroRowIndices) {
      _update_cell(id, i++);
    }
  } else {
    for (const auto& p : nonZeroRowIndices) {
      _update_cell(operators_->get_value(p.second), p.first, i++);
    }
  }
}

template <class Master_matrix>
template <class Container_type, class Row_container_type>
inline Vector_column<Master_matrix>::Vector_column(
    index columnIndex, 
    const Container_type& nonZeroRowIndices, 
    dimension_type dimension,
    Row_container_type* rowContainer, 
    Column_settings* colSettings)
    : ra_opt(columnIndex, rowContainer),
      dim_opt(dimension),
      chain_opt([&] {
        if constexpr (Master_matrix::Option_list::is_z2) {
          return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : *std::prev(nonZeroRowIndices.end());
        } else {
          return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : std::prev(nonZeroRowIndices.end())->first;
        }
      }()),
      column_(nonZeroRowIndices.size(), nullptr),
      operators_(nullptr),
      cellPool_(&(colSettings->cellConstructor))
{
  if constexpr (!Master_matrix::Option_list::is_z2){
    operators_ = &(colSettings->operators);
  }

  index i = 0;
  if constexpr (Master_matrix::Option_list::is_z2) {
    for (id_index id : nonZeroRowIndices) {
      _update_cell(id, i++);
    }
  } else {
    for (const auto& p : nonZeroRowIndices) {
      _update_cell(operators_->get_value(p.second), p.first, i++);
    }
  }
}

template <class Master_matrix>
inline Vector_column<Master_matrix>::Vector_column(const Vector_column& column,
                                                                     Column_settings* colSettings)
    : ra_opt(),
      dim_opt(static_cast<const dim_opt&>(column)),
      chain_opt(static_cast<const chain_opt&>(column)),
      column_(column.column_.size(), nullptr),
      erasedValues_(column.erasedValues_),
      operators_(colSettings == nullptr ? column.operators_ : nullptr),
      cellPool_(colSettings == nullptr ? column.cellPool_ : &(colSettings->cellConstructor))
{
  static_assert(!Master_matrix::Option_list::has_row_access,
                "Simple copy constructor not available when row access option enabled. Please specify the new column "
                "index and the row container.");

  if constexpr (!Master_matrix::Option_list::is_z2){
    if (colSettings != nullptr) operators_ = &(colSettings->operators);
  }

  index i = 0;
  for (const Cell* cell : column.column_) {
    if constexpr (Master_matrix::Option_list::is_z2) {
      _update_cell(cell->get_row_index(), i++);
    } else {
      _update_cell(cell->get_element(), cell->get_row_index(), i++);
    }
  }
}

template <class Master_matrix>
template <class Row_container_type>
inline Vector_column<Master_matrix>::Vector_column(const Vector_column& column, index columnIndex,
                                                                     Row_container_type* rowContainer,
                                                                     Column_settings* colSettings)
    : ra_opt(columnIndex, rowContainer),
      dim_opt(static_cast<const dim_opt&>(column)),
      chain_opt(static_cast<const chain_opt&>(column)),
      column_(column.column_.size(), nullptr),
      erasedValues_(column.erasedValues_),
      operators_(colSettings == nullptr ? column.operators_ : nullptr),
      cellPool_(colSettings == nullptr ? column.cellPool_ : &(colSettings->cellConstructor))
{
  if constexpr (!Master_matrix::Option_list::is_z2){
    if (colSettings != nullptr) operators_ = &(colSettings->operators);
  }

  index i = 0;
  for (const Cell* cell : column.column_) {
    if constexpr (Master_matrix::Option_list::is_z2) {
      _update_cell(cell->get_row_index(), i++);
    } else {
      _update_cell(cell->get_element(), cell->get_row_index(), i++);
    }
  }
}

template <class Master_matrix>
inline Vector_column<Master_matrix>::Vector_column(Vector_column&& column) noexcept
    : ra_opt(std::move(static_cast<ra_opt&>(column))),
      dim_opt(std::move(static_cast<dim_opt&>(column))),
      chain_opt(std::move(static_cast<chain_opt&>(column))),
      column_(std::move(column.column_)),
      erasedValues_(std::move(column.erasedValues_)),
      operators_(std::exchange(column.operators_, nullptr)),
      cellPool_(std::exchange(column.cellPool_, nullptr)) 
{}

template <class Master_matrix>
inline Vector_column<Master_matrix>::~Vector_column() 
{
  for (auto* cell : column_) {
    _delete_cell(cell);
  }
}

template <class Master_matrix>
inline std::vector<typename Vector_column<Master_matrix>::Field_element_type>
Vector_column<Master_matrix>::get_content(int columnLength) const 
{
  if (columnLength < 0 && column_.size() > 0)
    columnLength = column_.back()->get_row_index() + 1;
  else if (columnLength < 0)
    return std::vector<Field_element_type>();

  std::vector<Field_element_type> container(columnLength, 0);
  for (auto it = column_.begin(); it != column_.end() && (*it)->get_row_index() < static_cast<id_index>(columnLength);
       ++it) {
    if constexpr (!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type) {
      if (erasedValues_.find((*it)->get_row_index()) != erasedValues_.end()) continue;
    }
    if constexpr (Master_matrix::Option_list::is_z2) {
      container[(*it)->get_row_index()] = 1;
    } else {
      container[(*it)->get_row_index()] = (*it)->get_element();
    }
  }
  return container;
}

template <class Master_matrix>
inline bool Vector_column<Master_matrix>::is_non_zero(id_index rowIndex) const 
{
  if constexpr (!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type)
    if (erasedValues_.find(rowIndex) != erasedValues_.end()) return false;

  Cell cell(rowIndex);
  return std::binary_search(column_.begin(), column_.end(), &cell,
                            [](const Cell* a, const Cell* b) { return a->get_row_index() < b->get_row_index(); });
}

template <class Master_matrix>
inline bool Vector_column<Master_matrix>::is_empty() const 
{
  if constexpr (!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type) {
    return column_.size() == erasedValues_.size();  // assumes that erasedValues is always a subset of column_, which is
                                                    // wrong if someone cleared an non exitsing value...
  } else {
    return column_.empty();
  }
}

template <class Master_matrix>
inline std::size_t Vector_column<Master_matrix>::size() const 
{
  if constexpr (!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type) {
    return column_.size() - erasedValues_.size();  // assumes that erasedValues is always a subset of column_, which is
                                                   // wrong if someone cleared an non exitsing value...
  } else {
    return column_.size();
  }
}

template <class Master_matrix>
template <class Map_type>
inline void Vector_column<Master_matrix>::reorder(const Map_type& valueMap,
                                                                    [[maybe_unused]] index columnIndex) 
{
  static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type,
                "Method not available for chain columns.");

  if (erasedValues_.empty()) {  // to avoid useless push_backs.
    for (Cell* cell : column_) {
      if constexpr (Master_matrix::Option_list::has_row_access) {
        ra_opt::unlink(cell);
        if (columnIndex != static_cast<index>(-1)) cell->set_column_index(columnIndex);
      }
      cell->set_row_index(valueMap.at(cell->get_row_index()));
      if constexpr (Master_matrix::Option_list::has_intrusive_rows && Master_matrix::Option_list::has_row_access)
        ra_opt::insert_cell(cell->get_row_index(), cell);
    }

    // all cells have to be deleted first, to avoid problem with insertion when row is a set
    if constexpr (!Master_matrix::Option_list::has_intrusive_rows && Master_matrix::Option_list::has_row_access) {
      for (Cell* cell : column_) {
        ra_opt::insert_cell(cell->get_row_index(), cell);
      }
    }

    std::sort(column_.begin(), column_.end(), [](const Cell* c1, const Cell* c2) { return *c1 < *c2; });
  } else {
    Column_type newColumn;
    for (Cell* cell : column_) {
      if (erasedValues_.find(cell->get_row_index()) == erasedValues_.end()) {
        if constexpr (Master_matrix::Option_list::has_row_access) {
          ra_opt::unlink(cell);
          if (columnIndex != static_cast<index>(-1)) cell->set_column_index(columnIndex);
        }
        cell->set_row_index(valueMap.at(cell->get_row_index()));
        newColumn.push_back(cell);
        if constexpr (Master_matrix::Option_list::has_intrusive_rows && Master_matrix::Option_list::has_row_access)
          ra_opt::insert_cell(cell->get_row_index(), cell);
      } else {
        _delete_cell(cell);
      }
    }
    // all cells have to be deleted first, to avoid problem with insertion when row is a set
    if constexpr (!Master_matrix::Option_list::has_intrusive_rows && Master_matrix::Option_list::has_row_access) {
      for (Cell* cell : column_) {
        ra_opt::insert_cell(cell->get_row_index(), cell);
      }
    }
    std::sort(newColumn.begin(), newColumn.end(), [](const Cell* c1, const Cell* c2) { return *c1 < *c2; });
    erasedValues_.clear();
    column_.swap(newColumn);
  }
}

template <class Master_matrix>
inline void Vector_column<Master_matrix>::clear() 
{
  static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type,
                "Method not available for chain columns as a base element should not be empty.");

  for (auto* cell : column_) {
    if constexpr (Master_matrix::Option_list::has_row_access) ra_opt::unlink(cell);
    cellPool_->destroy(cell);
  }

  column_.clear();
  erasedValues_.clear();
}

template <class Master_matrix>
inline void Vector_column<Master_matrix>::clear(id_index rowIndex) 
{
  static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type,
                "Method not available for chain columns.");

  erasedValues_.insert(rowIndex);
}

template <class Master_matrix>
inline typename Vector_column<Master_matrix>::id_index
Vector_column<Master_matrix>::get_pivot() 
{
  static_assert(Master_matrix::isNonBasic,
                "Method not available for base columns.");  // could technically be, but is the notion usefull then?

  if constexpr (Master_matrix::Option_list::is_of_boundary_type) {
    if (column_.empty()) return -1;
    if (erasedValues_.empty()) return column_.back()->get_row_index();

    auto it = erasedValues_.find(column_.back()->get_row_index());
    while (!column_.empty() && it != erasedValues_.end()) {
      erasedValues_.erase(it);
      _delete_cell(column_.back());
      column_.pop_back();
      if (!column_.empty()) it = erasedValues_.find(column_.back()->get_row_index());
    }

    if (column_.empty()) return -1;
    return column_.back()->get_row_index();
  } else {
    return chain_opt::get_pivot();
  }
}

template <class Master_matrix>
inline typename Vector_column<Master_matrix>::Field_element_type
Vector_column<Master_matrix>::get_pivot_value() 
{
  static_assert(Master_matrix::isNonBasic,
                "Method not available for base columns.");  // could technically be, but is the notion usefull then?

  if constexpr (Master_matrix::Option_list::is_z2) {
    return 1;
  } else {
    if constexpr (Master_matrix::Option_list::is_of_boundary_type) {
      if (column_.empty()) return 0;
      if (erasedValues_.empty()) return column_.back()->get_element();

      auto it = erasedValues_.find(column_.back()->get_row_index());
      while (!column_.empty() && it != erasedValues_.end()) {
        erasedValues_.erase(it);
        _delete_cell(column_.back());
        column_.pop_back();
        if (!column_.empty()) it = erasedValues_.find(column_.back()->get_row_index());
      }

      if (column_.empty()) return 0;
      return column_.back()->get_element();
    } else {
      if (chain_opt::get_pivot() == static_cast<id_index>(-1)) return Field_element_type();
      for (const Cell* cell : column_) {
        if (cell->get_row_index() == chain_opt::get_pivot()) return cell->get_element();
      }
      return Field_element_type();  // should never happen if chain column is used properly
    }
  }
}

template <class Master_matrix>
inline typename Vector_column<Master_matrix>::iterator
Vector_column<Master_matrix>::begin() noexcept 
{
  return column_.begin();
}

template <class Master_matrix>
inline typename Vector_column<Master_matrix>::const_iterator
Vector_column<Master_matrix>::begin() const noexcept 
{
  return column_.begin();
}

template <class Master_matrix>
inline typename Vector_column<Master_matrix>::iterator
Vector_column<Master_matrix>::end() noexcept 
{
  return column_.end();
}

template <class Master_matrix>
inline typename Vector_column<Master_matrix>::const_iterator
Vector_column<Master_matrix>::end() const noexcept 
{
  return column_.end();
}

template <class Master_matrix>
inline typename Vector_column<Master_matrix>::reverse_iterator
Vector_column<Master_matrix>::rbegin() noexcept 
{
  return column_.rbegin();
}

template <class Master_matrix>
inline typename Vector_column<Master_matrix>::const_reverse_iterator
Vector_column<Master_matrix>::rbegin() const noexcept 
{
  return column_.rbegin();
}

template <class Master_matrix>
inline typename Vector_column<Master_matrix>::reverse_iterator
Vector_column<Master_matrix>::rend() noexcept 
{
  return column_.rend();
}

template <class Master_matrix>
inline typename Vector_column<Master_matrix>::const_reverse_iterator
Vector_column<Master_matrix>::rend() const noexcept 
{
  return column_.rend();
}

template <class Master_matrix>
template <class Cell_range>
inline Vector_column<Master_matrix>& Vector_column<Master_matrix>::operator+=(
    const Cell_range& column) 
{
  static_assert((!Master_matrix::isNonBasic || std::is_same_v<Cell_range, Vector_column>),
                "For boundary columns, the range has to be a column of same type to help ensure the validity of the "
                "base element.");  // could be removed, if we give the responsability to the user.
  static_assert((!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type),
                "For chain columns, the given column cannot be constant.");

  _add(column);

  return *this;
}

template <class Master_matrix>
inline Vector_column<Master_matrix>& Vector_column<Master_matrix>::operator+=(
    Vector_column& column) 
{
  if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type) {
    // assumes that the addition never zeros out this column.
    if (_add(column)) {
      chain_opt::swap_pivots(column);
      dim_opt::swap_dimension(column);
    }
  } else {
    _add(column);
  }

  return *this;
}

template <class Master_matrix>
inline Vector_column<Master_matrix>& Vector_column<Master_matrix>::operator*=(
    unsigned int v) 
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
    Field_element_type val = operators_->get_value(v);

    if (val == Field_operators::get_additive_identity()) {
      if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type) {
        throw std::invalid_argument("A chain column should not be multiplied by 0.");
      } else {
        clear();
      }
      return *this;
    }

    if (val == Field_operators::get_multiplicative_identity()) return *this;

    for (Cell* cell : column_) {
      operators_->multiply_inplace(cell->get_element(), val);
      if constexpr (Master_matrix::Option_list::has_row_access) ra_opt::update_cell(*cell);
    }
  }

  return *this;
}

template <class Master_matrix>
template <class Cell_range>
inline Vector_column<Master_matrix>& Vector_column<Master_matrix>::multiply_target_and_add(
    const Field_element_type& val, const Cell_range& column) 
{
  static_assert((!Master_matrix::isNonBasic || std::is_same_v<Cell_range, Vector_column>),
                "For boundary columns, the range has to be a column of same type to help ensure the validity of the "
                "base element.");  // could be removed, if we give the responsability to the user.
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
inline Vector_column<Master_matrix>& Vector_column<Master_matrix>::multiply_target_and_add(
    const Field_element_type& val, Vector_column& column) 
{
  if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type) {
    // assumes that the addition never zeros out this column.
    if constexpr (Master_matrix::Option_list::is_z2) {
      if (val) {
        if (_add(column)) {
          chain_opt::swap_pivots(column);
          dim_opt::swap_dimension(column);
        }
      } else {
        throw std::invalid_argument("A chain column should not be multiplied by 0.");
      }
    } else {
      if (_multiply_target_and_add(val, column)) {
        chain_opt::swap_pivots(column);
        dim_opt::swap_dimension(column);
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
template <class Cell_range>
inline Vector_column<Master_matrix>& Vector_column<Master_matrix>::multiply_source_and_add(
    const Cell_range& column, const Field_element_type& val) 
{
  static_assert((!Master_matrix::isNonBasic || std::is_same_v<Cell_range, Vector_column>),
                "For boundary columns, the range has to be a column of same type to help ensure the validity of the "
                "base element.");  // could be removed, if we give the responsability to the user.
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
inline Vector_column<Master_matrix>& Vector_column<Master_matrix>::multiply_source_and_add(
    Vector_column& column, const Field_element_type& val) 
{
  if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type) {
    // assumes that the addition never zeros out this column.
    if constexpr (Master_matrix::Option_list::is_z2) {
      if (val) {
        if (_add(column)) {
          chain_opt::swap_pivots(column);
          dim_opt::swap_dimension(column);
        }
      }
    } else {
      if (_multiply_source_and_add(column, val)) {
        chain_opt::swap_pivots(column);
        dim_opt::swap_dimension(column);
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
inline Vector_column<Master_matrix>& Vector_column<Master_matrix>::operator=(
    const Vector_column& other) 
{
  static_assert(!Master_matrix::Option_list::has_row_access, "= assignement not enabled with row access option.");

  dim_opt::operator=(other);
  chain_opt::operator=(other);

  auto tmpPool = cellPool_;
  cellPool_ = other.cellPool_;

  while (column_.size() > other.column_.size()) {
    if (column_.back() != nullptr) {
      if constexpr (Master_matrix::Option_list::has_row_access) ra_opt::unlink(column_.back());
      tmpPool->destroy(column_.back());
    }
    column_.pop_back();
  }

  column_.resize(other.column_.size(), nullptr);
  index i = 0;
  for (const Cell* cell : other.column_) {
    if (column_[i] != nullptr) {
      if constexpr (Master_matrix::Option_list::has_row_access) ra_opt::unlink(column_[i]);
      tmpPool->destroy(column_[i]);
    }
    if constexpr (Master_matrix::Option_list::is_z2) {
      _update_cell(cell->get_row_index(), i++);
    } else {
      _update_cell(cell->get_element(), cell->get_row_index(), i++);
    }
  }
  erasedValues_ = other.erasedValues_;
  operators_ = other.operators_;

  return *this;
}

template <class Master_matrix>
inline std::size_t Vector_column<Master_matrix>::compute_hash_value()
{
  std::size_t seed = 0;
  for (Cell* cell : column_) {
    if (erasedValues_.find(cell->get_row_index()) == erasedValues_.end()){
      seed ^= std::hash<unsigned int>()(cell->get_row_index() * static_cast<unsigned int>(cell->get_element())) +
              0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
  }
  return seed;
}

template <class Master_matrix>
inline void Vector_column<Master_matrix>::_delete_cell(Cell* cell) 
{
  if constexpr (Master_matrix::Option_list::has_row_access) ra_opt::unlink(cell);
  cellPool_->destroy(cell);
}

template <class Master_matrix>
inline void Vector_column<Master_matrix>::_delete_cell(typename Column_type::iterator& it)
{
  _delete_cell(*it);
  ++it;
}

template <class Master_matrix>
inline typename Vector_column<Master_matrix>::Cell* Vector_column<Master_matrix>::_insert_cell(
    const Field_element_type& value, id_index rowIndex, Column_type& column)
{
  if constexpr (Master_matrix::Option_list::has_row_access) {
    Cell* newCell = cellPool_->construct(ra_opt::columnIndex_, rowIndex);
    newCell->set_element(value);
    column.push_back(newCell);
    ra_opt::insert_cell(rowIndex, newCell);
    return newCell;
  } else {
    Cell* newCell = cellPool_->construct(rowIndex);
    newCell->set_element(value);
    column.push_back(newCell);
    return newCell;
  }
}

template <class Master_matrix>
inline void Vector_column<Master_matrix>::_insert_cell(id_index rowIndex, Column_type& column) 
{
  if constexpr (Master_matrix::Option_list::has_row_access) {
    Cell* newCell = cellPool_->construct(ra_opt::columnIndex_, rowIndex);
    column.push_back(newCell);
    ra_opt::insert_cell(rowIndex, newCell);
  } else {
    Cell* newCell = cellPool_->construct(rowIndex);
    column.push_back(newCell);
  }
}

template <class Master_matrix>
inline void Vector_column<Master_matrix>::_update_cell(const Field_element_type& value,
                                                                         id_index rowIndex, index position) 
{
  if constexpr (Master_matrix::Option_list::has_row_access) {
    Cell* newCell = cellPool_->construct(ra_opt::columnIndex_, rowIndex);
    newCell->set_element(value);
    column_[position] = newCell;
    ra_opt::insert_cell(rowIndex, newCell);
  } else {
    column_[position] = cellPool_->construct(rowIndex);
    column_[position]->set_element(value);
  }
}

template <class Master_matrix>
inline void Vector_column<Master_matrix>::_update_cell(id_index rowIndex, index position) 
{
  if constexpr (Master_matrix::Option_list::has_row_access) {
    Cell* newCell = cellPool_->construct(ra_opt::columnIndex_, rowIndex);
    column_[position] = newCell;
    ra_opt::insert_cell(rowIndex, newCell);
  } else {
    column_[position] = cellPool_->construct(rowIndex);
  }
}

template <class Master_matrix>
template <class Cell_range>
inline bool Vector_column<Master_matrix>::_add(const Cell_range& column) 
{
  if (column.begin() == column.end()) return false;
  if (column_.empty()) {  // chain should never enter here.
    column_.resize(column.size());
    index i = 0;
    for (const Cell& cell : column) {
      if constexpr (Master_matrix::Option_list::is_z2) {
        _update_cell(cell.get_row_index(), i++);
      } else {
        _update_cell(cell.get_element(), cell.get_row_index(), i++);
      }
    }
    return true;
  }

  Column_type newColumn;
  newColumn.reserve(column_.size() + column.size());  // safe upper bound

  auto pivotIsZeroed = _generic_add(
      column,
      [&](Cell* cellTarget) { newColumn.push_back(cellTarget); },
      [&](typename Cell_range::const_iterator& itSource,
          [[maybe_unused]] const typename Column_type::iterator& itTarget) {
        if constexpr (Master_matrix::Option_list::is_z2) {
          _insert_cell(itSource->get_row_index(), newColumn);
        } else {
          _insert_cell(itSource->get_element(), itSource->get_row_index(), newColumn);
        }
      },
      [&](Field_element_type& targetElement, typename Cell_range::const_iterator& itSource) {
        if constexpr (!Master_matrix::Option_list::is_z2)
          operators_->add_inplace(targetElement, itSource->get_element());
      },
      [&](Cell* cellTarget) { newColumn.push_back(cellTarget); });

  column_.swap(newColumn);

  return pivotIsZeroed;
}

template <class Master_matrix>
template <class Cell_range>
inline bool Vector_column<Master_matrix>::_multiply_target_and_add(const Field_element_type& val,
                                                                              const Cell_range& column) 
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
    index i = 0;
    for (const Cell& cell : column) {
      if constexpr (Master_matrix::Option_list::is_z2) {
        _update_cell(cell.get_row_index(), i++);
      } else {
        _update_cell(cell.get_element(), cell.get_row_index(), i++);
      }
    }
    if constexpr (std::is_same_v<Cell_range, Vector_column<Master_matrix> >) erasedValues_ = column.erasedValues_;
    return true;
  }

  Column_type newColumn;
  newColumn.reserve(column_.size() + column.size());  // safe upper bound

  auto pivotIsZeroed = _generic_add(
      column,
      [&](Cell* cellTarget) {
        operators_->multiply_inplace(cellTarget->get_element(), val);
        if constexpr (Master_matrix::Option_list::has_row_access) ra_opt::update_cell(*cellTarget);
        newColumn.push_back(cellTarget);
      },
      [&](typename Cell_range::const_iterator& itSource, const typename Column_type::iterator& itTarget) {
        _insert_cell(itSource->get_element(), itSource->get_row_index(), newColumn);
      },
      [&](Field_element_type& targetElement, typename Cell_range::const_iterator& itSource) {
        operators_->multiply_and_add_inplace_front(targetElement, val, itSource->get_element());
      },
      [&](Cell* cellTarget) { newColumn.push_back(cellTarget); });

  column_.swap(newColumn);

  return pivotIsZeroed;
}

template <class Master_matrix>
template <class Cell_range>
inline bool Vector_column<Master_matrix>::_multiply_source_and_add(const Cell_range& column,
                                                                              const Field_element_type& val) 
{
  if (val == 0u || column.begin() == column.end()) {
    return false;
  }

  Column_type newColumn;
  newColumn.reserve(column_.size() + column.size());  // safe upper bound

  auto pivotIsZeroed = _generic_add(
      column,
      [&](Cell* cellTarget) { newColumn.push_back(cellTarget); },
      [&](typename Cell_range::const_iterator& itSource, const typename Column_type::iterator& itTarget) {
        Cell* newCell = _insert_cell(itSource->get_element(), itSource->get_row_index(), newColumn);
        operators_->multiply_inplace(newCell->get_element(), val);
        if constexpr (Master_matrix::Option_list::has_row_access) ra_opt::update_cell(*newCell);
      },
      [&](Field_element_type& targetElement, typename Cell_range::const_iterator& itSource) {
        operators_->multiply_and_add_inplace_back(itSource->get_element(), val, targetElement);
      },
      [&](Cell* cellTarget) { newColumn.push_back(cellTarget); });

  column_.swap(newColumn);

  return pivotIsZeroed;
}

template <class Master_matrix>
template <class Cell_range, typename F1, typename F2, typename F3, typename F4>
inline bool Vector_column<Master_matrix>::_generic_add(const Cell_range& column,
                                                       F1&& process_target,
                                                       F2&& process_source,
                                                       F3&& update_target1,
                                                       F4&& update_target2)
{
  auto updateTargetIterator = [&](typename Column_type::iterator& itTarget){
    if constexpr (!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type){
      while (itTarget != column_.end() && erasedValues_.find((*itTarget)->get_row_index()) != erasedValues_.end()) {
        _delete_cell(*itTarget);
        ++itTarget;
      }
    }
  };
  auto updateSourceIterator = [&](typename Cell_range::const_iterator& itSource){
    if constexpr (std::is_same_v<Cell_range, Vector_column<Master_matrix> > &&
                  (!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type)) {
      while (itSource != column.end() &&
             column.erasedValues_.find(itSource->get_row_index()) != column.erasedValues_.end())
        ++itSource;
    }
  };

  bool pivotIsZeroed = false;

  auto itTarget = column_.begin();
  auto itSource = column.begin();
  while (itTarget != column_.end() && itSource != column.end()) {
    updateTargetIterator(itTarget);
    updateSourceIterator(itSource);
    if (itTarget == column_.end() || itSource == column.end()) break;

    _generic_merge_cell_to_column(*this,
                                   itSource, itTarget,
                                   process_target, process_source, update_target1, update_target2,
                                   pivotIsZeroed);
  }

  while (itSource != column.end()) {
    updateSourceIterator(itSource);
    if (itSource == column.end()) break;

    process_source(itSource, column_.end());
    ++itSource;
  }

  while (itTarget != column_.end()) {
    updateTargetIterator(itTarget);
    if (itTarget == column_.end()) break;

    process_target(*itTarget);
    ++itTarget;
  }

  erasedValues_.clear();

  return pivotIsZeroed;
}

}  // namespace persistence_matrix
}  // namespace Gudhi

/**
 * @ingroup persistence_matrix
 *
 * @brief Hash method for @ref Gudhi::persistence_matrix::Vector_column.
 * 
 * @tparam Master_matrix Template parameter of @ref Gudhi::persistence_matrix::Vector_column.
 * @tparam Cell_constructor Template parameter of @ref Gudhi::persistence_matrix::Vector_column.
 */
template <class Master_matrix>
struct std::hash<Gudhi::persistence_matrix::Vector_column<Master_matrix> > 
{
  std::size_t operator()(const Gudhi::persistence_matrix::Vector_column<Master_matrix>& column) const {
    return column.compute_hash_value();
  }
};

#endif  // PM_VECTOR_COLUMN_H

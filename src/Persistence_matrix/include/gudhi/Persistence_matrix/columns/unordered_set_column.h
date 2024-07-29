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

#include <vector>
#include <stdexcept>
#include <type_traits>
#include <set>
#include <utility>  //std::swap, std::move & std::exchange

#include <boost/iterator/indirect_iterator.hpp>
#if BOOST_VERSION >= 108100
#include <boost/unordered/unordered_flat_set.hpp>
#else
#include <unordered_set>
#endif

#include <gudhi/Persistence_matrix/allocators/cell_constructors.h>

namespace Gudhi {
namespace persistence_matrix {

// For unordered_set container. Outside of Unordered_set_column because of a msvc bug who can't compile properly
// unordered_flat_set if the hash method is nested.
template <class Cell>
struct CellPointerHash
{
  size_t operator()(const Cell* c) const { return std::hash<Cell>()(*c); }
};
template <class Cell>
struct CellPointerEq
{
  bool operator()(const Cell* c1, const Cell* c2) const { return *c1 == *c2; }
};

/**
 * @class Unordered_set_column unordered_set_column.h gudhi/Persistence_matrix/columns/unordered_set_column.h
 * @ingroup persistence_matrix
 *
 * @brief Column class following the @ref PersistenceMatrixColumn concept.
 *
 * Column based on an unordered set structure. The cells are not ordered, but only non-zero values
 * are stored uniquely in the underlying container. When adding a cell range into it, the given cell range
 * also does not need to be ordered (contrary to most other column types).
 *
 * @tparam Master_matrix An instantiation of @ref Matrix from which all types and options are deduced.
 * @tparam Cell_constructor Factory of @ref Cell classes.
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
  using Cell = typename Master_matrix::Matrix_cell;
  using Column_settings = typename Master_matrix::Column_settings;

 private:
  using Field_operators = typename Master_matrix::Field_operators;
  using Cell_constructor = typename Master_matrix::Cell_constructor;

  struct CellPointerComp {
    bool operator()(const Cell* c1, const Cell* c2) const { return *c1 < *c2; }
  };

#if BOOST_VERSION >= 108100
  using Column_support = boost::unordered_flat_set<Cell*, CellPointerHash<Cell>, CellPointerEq<Cell>>;
#else
  using Column_support = std::unordered_set<Cell*, CellPointerHash<Cell>, CellPointerEq<Cell>>;
#endif

 public:
  using iterator = boost::indirect_iterator<typename Column_support::iterator>;
  using const_iterator = boost::indirect_iterator<typename Column_support::const_iterator>;

  Unordered_set_column(Column_settings* colSettings = nullptr);
  template <class Container = typename Master_matrix::Boundary>
  Unordered_set_column(const Container& nonZeroRowIndices, Column_settings* colSettings);
  template <class Container = typename Master_matrix::Boundary, class Row_container>
  Unordered_set_column(Index columnIndex,
                       const Container& nonZeroRowIndices,
                       Row_container* rowContainer,
                       Column_settings* colSettings);
  template <class Container = typename Master_matrix::Boundary>
  Unordered_set_column(const Container& nonZeroChainRowIndices, Dimension dimension, Column_settings* colSettings);
  template <class Container = typename Master_matrix::Boundary, class Row_container>
  Unordered_set_column(Index columnIndex,
                       const Container& nonZeroChainRowIndices,
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

  template <class Cell_range>
  Unordered_set_column& operator+=(const Cell_range& column);
  Unordered_set_column& operator+=(Unordered_set_column& column);

  Unordered_set_column& operator*=(unsigned int v);

  // this = v * this + column
  template <class Cell_range>
  Unordered_set_column& multiply_target_and_add(const Field_element& val, const Cell_range& column);
  Unordered_set_column& multiply_target_and_add(const Field_element& val, Unordered_set_column& column);
  // this = this + column * v
  template <class Cell_range>
  Unordered_set_column& multiply_source_and_add(const Cell_range& column, const Field_element& val);
  Unordered_set_column& multiply_source_and_add(Unordered_set_column& column, const Field_element& val);

  friend bool operator==(const Unordered_set_column& c1, const Unordered_set_column& c2) {
    if (&c1 == &c2) return true;
    if (c1.column_.size() != c2.column_.size()) return false;

    for (Cell* cell : c1.column_) {
      auto it = c2.column_.find(cell);
      if (it == c2.column_.end()) return false;
      if constexpr (!Master_matrix::Option_list::is_z2)
        if ((*it)->get_element() != cell->get_element()) return false;
    }
    return true;
  }
  friend bool operator<(const Unordered_set_column& c1, const Unordered_set_column& c2) {
    if (&c1 == &c2) return false;

    using ID_index = Unordered_set_column<Master_matrix>::ID_index;
    using Cell_rep =
        typename std::conditional<Master_matrix::Option_list::is_z2,
                                  ID_index,
                                  std::pair<ID_index, unsigned int>
                                 >::type;

    auto it1 = c1.column_.begin();
    auto it2 = c2.column_.begin();
    std::set<Cell_rep> cells1, cells2;
    while (it1 != c1.column_.end() && it2 != c2.column_.end()) {
      if constexpr (Master_matrix::Option_list::is_z2) {
        cells1.insert((*it1)->get_row_index());
        cells2.insert((*it2)->get_row_index());
      } else {
        cells1.emplace((*it1)->get_row_index(), (*it1)->get_element());
        cells2.emplace((*it2)->get_row_index(), (*it2)->get_element());
      }
      ++it1;
      ++it2;
    }
    while (it1 != c1.column_.end()) {
      if constexpr (Master_matrix::Option_list::is_z2) {
        cells1.insert((*it1)->get_row_index());
      } else {
        cells1.emplace((*it1)->get_row_index(), (*it1)->get_element());
      }
      ++it1;
    }
    while (it2 != c2.column_.end()) {
      if constexpr (Master_matrix::Option_list::is_z2) {
        cells2.insert((*it2)->get_row_index());
      } else {
        cells2.emplace((*it2)->get_row_index(), (*it2)->get_element());
      }
      ++it2;
    }
    return cells1 < cells2;
  }

  // Disabled with row access.
  Unordered_set_column& operator=(const Unordered_set_column& other);

  friend void swap(Unordered_set_column& col1, Unordered_set_column& col2) {
    swap(static_cast<typename Master_matrix::Row_access_option&>(col1),
         static_cast<typename Master_matrix::Row_access_option&>(col2));
    swap(static_cast<typename Master_matrix::Column_dimension_option&>(col1),
         static_cast<typename Master_matrix::Column_dimension_option&>(col2));
    swap(static_cast<typename Master_matrix::Chain_column_option&>(col1),
         static_cast<typename Master_matrix::Chain_column_option&>(col2));
    col1.column_.swap(col2.column_);
    std::swap(col1.operators_, col2.operators_);
    std::swap(col1.cellPool_, col2.cellPool_);
  }

 private:
  using RA_opt = typename Master_matrix::Row_access_option;
  using Dim_opt = typename Master_matrix::Column_dimension_option;
  using Chain_opt = typename Master_matrix::Chain_column_option;

  Column_support column_;
  Field_operators* operators_;
  Cell_constructor* cellPool_;

  void _delete_cell(typename Column_support::iterator& it);
  Cell* _insert_cell(const Field_element& value, ID_index rowIndex);
  void _insert_cell(ID_index rowIndex);
  template <class Cell_range>
  bool _add(const Cell_range& column);
  template <class Cell_range>
  bool _multiply_target_and_add(const Field_element& val, const Cell_range& column);
  template <class Cell_range>
  bool _multiply_source_and_add(const Cell_range& column, const Field_element& val);
  template <class Cell_range, typename F1, typename F2>
  bool _generic_add(const Cell_range& source, F1&& process_source, F2&& update_target);
};

template <class Master_matrix>
inline Unordered_set_column<Master_matrix>::Unordered_set_column(Column_settings* colSettings)
    : RA_opt(),
      Dim_opt(),
      Chain_opt(),
      operators_(nullptr),
      cellPool_(colSettings == nullptr ? nullptr : &(colSettings->cellConstructor))
{
  if (operators_ == nullptr && cellPool_ == nullptr) return;  // to allow default constructor which gives a dummy column
  if constexpr (!Master_matrix::Option_list::is_z2) {
    operators_ = &(colSettings->operators);
  }
}

template <class Master_matrix>
template <class Container>
inline Unordered_set_column<Master_matrix>::Unordered_set_column(const Container& nonZeroRowIndices,
                                                                 Column_settings* colSettings)
    : RA_opt(),
      Dim_opt(nonZeroRowIndices.size() == 0 ? 0 : nonZeroRowIndices.size() - 1),
      Chain_opt(),
      column_(nonZeroRowIndices.size()),
      operators_(nullptr),
      cellPool_(&(colSettings->cellConstructor))
{
  static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type,
                "Constructor not available for chain columns, please specify the dimension of the chain.");

  if constexpr (Master_matrix::Option_list::is_z2) {
    for (ID_index id : nonZeroRowIndices) {
      _insert_cell(id);
    }
  } else {
    operators_ = &(colSettings->operators);
    for (const auto& p : nonZeroRowIndices) {
      _insert_cell(operators_->get_value(p.second), p.first);
    }
  }
}

template <class Master_matrix>
template <class Container, class Row_container>
inline Unordered_set_column<Master_matrix>::Unordered_set_column(Index columnIndex,
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
      column_(nonZeroRowIndices.size()),
      operators_(nullptr),
      cellPool_(&(colSettings->cellConstructor))
{
  static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type,
                "Constructor not available for chain columns, please specify the dimension of the chain.");

  if constexpr (Master_matrix::Option_list::is_z2) {
    for (ID_index id : nonZeroRowIndices) {
      _insert_cell(id);
    }
  } else {
    operators_ = &(colSettings->operators);
    for (const auto& p : nonZeroRowIndices) {
      _insert_cell(operators_->get_value(p.second), p.first);
    }
  }
}

template <class Master_matrix>
template <class Container>
inline Unordered_set_column<Master_matrix>::Unordered_set_column(const Container& nonZeroRowIndices,
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
      column_(nonZeroRowIndices.size()),
      operators_(nullptr),
      cellPool_(&(colSettings->cellConstructor))
{
  if constexpr (Master_matrix::Option_list::is_z2) {
    for (ID_index id : nonZeroRowIndices) {
      _insert_cell(id);
    }
  } else {
    operators_ = &(colSettings->operators);
    for (const auto& p : nonZeroRowIndices) {
      _insert_cell(operators_->get_value(p.second), p.first);
    }
  }
}

template <class Master_matrix>
template <class Container, class Row_container>
inline Unordered_set_column<Master_matrix>::Unordered_set_column(Index columnIndex,
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
      column_(nonZeroRowIndices.size()),
      operators_(nullptr),
      cellPool_(&(colSettings->cellConstructor))
{
  if constexpr (Master_matrix::Option_list::is_z2) {
    for (ID_index id : nonZeroRowIndices) {
      _insert_cell(id);
    }
  } else {
    operators_ = &(colSettings->operators);
    for (const auto& p : nonZeroRowIndices) {
      _insert_cell(operators_->get_value(p.second), p.first);
    }
  }
}

template <class Master_matrix>
inline Unordered_set_column<Master_matrix>::Unordered_set_column(const Unordered_set_column& column,
                                                                 Column_settings* colSettings)
    : RA_opt(),
      Dim_opt(static_cast<const Dim_opt&>(column)),
      Chain_opt(static_cast<const Chain_opt&>(column)),
      column_(column.column_.bucket_count()),
      operators_(colSettings == nullptr ? column.operators_ : nullptr),
      cellPool_(colSettings == nullptr ? column.cellPool_ : &(colSettings->cellConstructor))
{
  static_assert(!Master_matrix::Option_list::has_row_access,
                "Simple copy constructor not available when row access option enabled. Please specify the new column "
                "index and the row container.");

  if constexpr (!Master_matrix::Option_list::is_z2) {
    if (colSettings != nullptr) operators_ = &(colSettings->operators);
  }

  for (const Cell* cell : column.column_) {
    if constexpr (Master_matrix::Option_list::is_z2) {
      _insert_cell(cell->get_row_index());
    } else {
      _insert_cell(cell->get_element(), cell->get_row_index());
    }
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
      operators_(colSettings == nullptr ? column.operators_ : nullptr),
      cellPool_(colSettings == nullptr ? column.cellPool_ : &(colSettings->cellConstructor))
{
  if constexpr (!Master_matrix::Option_list::is_z2) {
    if (colSettings != nullptr) operators_ = &(colSettings->operators);
  }

  for (const Cell* cell : column.column_) {
    if constexpr (Master_matrix::Option_list::is_z2) {
      _insert_cell(cell->get_row_index());
    } else {
      _insert_cell(cell->get_element(), cell->get_row_index());
    }
  }
}

template <class Master_matrix>
inline Unordered_set_column<Master_matrix>::Unordered_set_column(Unordered_set_column&& column) noexcept
    : RA_opt(std::move(static_cast<RA_opt&>(column))),
      Dim_opt(std::move(static_cast<Dim_opt&>(column))),
      Chain_opt(std::move(static_cast<Chain_opt&>(column))),
      column_(std::move(column.column_)),
      operators_(std::exchange(column.operators_, nullptr)),
      cellPool_(std::exchange(column.cellPool_, nullptr))
{}

template <class Master_matrix>
inline Unordered_set_column<Master_matrix>::~Unordered_set_column()
{
  for (auto* cell : column_) {
    if constexpr (Master_matrix::Option_list::has_row_access) RA_opt::unlink(cell);
    cellPool_->destroy(cell);
  }
}

template <class Master_matrix>
inline std::vector<typename Unordered_set_column<Master_matrix>::Field_element>
Unordered_set_column<Master_matrix>::get_content(int columnLength) const
{
  if (columnLength < 0 && column_.size() > 0)
    columnLength = (*std::max_element(column_.begin(), column_.end(), CellPointerComp()))->get_row_index() + 1;
  else if (columnLength < 0)
    return std::vector<Field_element>();

  std::vector<Field_element> container(columnLength, 0);
  for (auto it = column_.begin(); it != column_.end(); ++it) {
    if ((*it)->get_row_index() < static_cast<ID_index>(columnLength)) {
      if constexpr (Master_matrix::Option_list::is_z2) {
        container[(*it)->get_row_index()] = 1;
      } else {
        container[(*it)->get_row_index()] = (*it)->get_element();
      }
    }
  }
  return container;
}

template <class Master_matrix>
inline bool Unordered_set_column<Master_matrix>::is_non_zero(ID_index rowIndex) const
{
  Cell cell(rowIndex);
  return column_.find(&cell) != column_.end();
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

  for (Cell* cell : column_) {
    if constexpr (Master_matrix::Option_list::has_row_access) {
      RA_opt::unlink(cell);
      if (columnIndex != static_cast<Index>(-1)) cell->set_column_index(columnIndex);
    }
    cell->set_row_index(valueMap.at(cell->get_row_index()));
    newSet.insert(cell);
    if constexpr (Master_matrix::Option_list::has_row_access &&
                  Master_matrix::Option_list::has_intrusive_rows)  // intrusive list
      RA_opt::insert_cell(cell->get_row_index(), cell);
  }

  // when row is a set, all cells have to be deleted first, to avoid colliding when inserting
  if constexpr (Master_matrix::Option_list::has_row_access && !Master_matrix::Option_list::has_intrusive_rows) {  // set
    for (Cell* cell : newSet) {
      RA_opt::insert_cell(cell->get_row_index(), cell);
    }
  }

  column_.swap(newSet);
}

template <class Master_matrix>
inline void Unordered_set_column<Master_matrix>::clear()
{
  static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type,
                "Method not available for chain columns as a base element should not be empty.");

  for (auto* cell : column_) {
    if constexpr (Master_matrix::Option_list::has_row_access) RA_opt::unlink(cell);
    cellPool_->destroy(cell);
  }

  column_.clear();
}

template <class Master_matrix>
inline void Unordered_set_column<Master_matrix>::clear(ID_index rowIndex)
{
  static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type,
                "Method not available for chain columns.");

  auto cell = cellPool_->construct(rowIndex);
  auto it = column_.find(cell);
  if (it != column_.end()) {
    _delete_cell(it);
  }
  cellPool_->destroy(cell);
}

template <class Master_matrix>
inline typename Unordered_set_column<Master_matrix>::ID_index Unordered_set_column<Master_matrix>::get_pivot() const
{
  static_assert(Master_matrix::isNonBasic,
                "Method not available for base columns.");  // could technically be, but is the notion useful then?

  if constexpr (Master_matrix::Option_list::is_of_boundary_type) {
    if (column_.empty()) return -1;
    // linear search could be avoided with storing the pivot. But even then, some modifications of the column requires
    // the max, so not clear how much it is worth it.
    return (*std::max_element(column_.begin(), column_.end(), CellPointerComp()))->get_row_index();
  } else {
    return Chain_opt::get_pivot();
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
      return (*std::max_element(column_.begin(), column_.end(), CellPointerComp()))->get_element();
    } else {
      if (Chain_opt::get_pivot() == static_cast<ID_index>(-1)) return Field_element();
      for (const Cell* cell : column_) {
        if (cell->get_row_index() == Chain_opt::get_pivot()) return cell->get_element();
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
template <class Cell_range>
inline Unordered_set_column<Master_matrix>& Unordered_set_column<Master_matrix>::operator+=(const Cell_range& column)
{
  static_assert((!Master_matrix::isNonBasic || std::is_same_v<Cell_range, Unordered_set_column>),
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
      Chain_opt::swap_pivots(column);
      Dim_opt::swap_dimension(column);
    }
  } else {
    _add(column);
  }

  return *this;
}

template <class Master_matrix>
inline Unordered_set_column<Master_matrix>& Unordered_set_column<Master_matrix>::operator*=(unsigned int v)
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

    for (Cell* cell : column_) {
      operators_->multiply_inplace(cell->get_element(), val);
      if constexpr (Master_matrix::Option_list::has_row_access) RA_opt::update_cell(*cell);
    }
  }

  return *this;
}

template <class Master_matrix>
template <class Cell_range>
inline Unordered_set_column<Master_matrix>& Unordered_set_column<Master_matrix>::multiply_target_and_add(
    const Field_element& val, const Cell_range& column)
{
  static_assert((!Master_matrix::isNonBasic || std::is_same_v<Cell_range, Unordered_set_column>),
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
inline Unordered_set_column<Master_matrix>& Unordered_set_column<Master_matrix>::multiply_target_and_add(
    const Field_element& val, Unordered_set_column& column)
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
template <class Cell_range>
inline Unordered_set_column<Master_matrix>& Unordered_set_column<Master_matrix>::multiply_source_and_add(
    const Cell_range& column, const Field_element& val)
{
  static_assert((!Master_matrix::isNonBasic || std::is_same_v<Cell_range, Unordered_set_column>),
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
inline Unordered_set_column<Master_matrix>& Unordered_set_column<Master_matrix>::multiply_source_and_add(
    Unordered_set_column& column, const Field_element& val)
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
inline Unordered_set_column<Master_matrix>& Unordered_set_column<Master_matrix>::operator=(
    const Unordered_set_column& other)
{
  static_assert(!Master_matrix::Option_list::has_row_access, "= assignment not enabled with row access option.");

  Dim_opt::operator=(other);
  Chain_opt::operator=(other);

  for (auto* cell : column_) {
    if constexpr (Master_matrix::Option_list::has_row_access) RA_opt::unlink(cell);
    cellPool_->destroy(cell);
  }
  column_.clear();

  operators_ = other.operators_;
  cellPool_ = other.cellPool_;

  for (const Cell* cell : other.column_) {
    if constexpr (Master_matrix::Option_list::is_z2) {
      _insert_cell(cell->get_row_index());
    } else {
      _insert_cell(cell->get_element(), cell->get_row_index());
    }
  }

  return *this;
}

template <class Master_matrix>
inline void Unordered_set_column<Master_matrix>::_delete_cell(typename Column_support::iterator& it)
{
  if constexpr (Master_matrix::Option_list::has_row_access) RA_opt::unlink(*it);
  cellPool_->destroy(*it);
  auto tmp = it++;
  // it = column_.erase(it);
  column_.erase(tmp);
}

template <class Master_matrix>
inline typename Unordered_set_column<Master_matrix>::Cell* Unordered_set_column<Master_matrix>::_insert_cell(
    const Field_element& value, ID_index rowIndex)
{
  if constexpr (Master_matrix::Option_list::has_row_access) {
    Cell* newCell = cellPool_->construct(RA_opt::columnIndex_, rowIndex);
    newCell->set_element(value);
    column_.insert(newCell);
    RA_opt::insert_cell(rowIndex, newCell);
    return newCell;
  } else {
    Cell* newCell = cellPool_->construct(rowIndex);
    newCell->set_element(value);
    column_.insert(newCell);
    return newCell;
  }
}

template <class Master_matrix>
inline void Unordered_set_column<Master_matrix>::_insert_cell(ID_index rowIndex)
{
  if constexpr (Master_matrix::Option_list::has_row_access) {
    Cell* newCell = cellPool_->construct(RA_opt::columnIndex_, rowIndex);
    column_.insert(newCell);
    RA_opt::insert_cell(rowIndex, newCell);
  } else {
    Cell* newCell = cellPool_->construct(rowIndex);
    column_.insert(newCell);
  }
}

template <class Master_matrix>
template <class Cell_range>
inline bool Unordered_set_column<Master_matrix>::_add(const Cell_range& column)
{
  return _generic_add(
      column,
      [&](const Cell& oldCell, Cell* newCell) {
        if constexpr (!Master_matrix::Option_list::is_z2) newCell->set_element(oldCell.get_element());
      },
      [&](Cell* targetCell, const Cell& sourceCell) {
        if constexpr (!Master_matrix::Option_list::is_z2)
          operators_->add_inplace(targetCell->get_element(), sourceCell.get_element());
      });
}

template <class Master_matrix>
template <class Cell_range>
inline bool Unordered_set_column<Master_matrix>::_multiply_target_and_add(const Field_element& val,
                                                                          const Cell_range& column)
{
  if (val == 0u) {
    if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type) {
      throw std::invalid_argument("A chain column should not be multiplied by 0.");
      // this would not only mess up the base, but also the pivots stored.
    } else {
      clear();
      for (const Cell& v : column) {
        _insert_cell(v.get_element(), v.get_row_index());
      }
      return true;
    }
  }

  // because the column is unordered, I don't see a way to do both operations in one go
  // without guarantees on the cell range...
  operator*=(val);
  return _add(column);
}

template <class Master_matrix>
template <class Cell_range>
inline bool Unordered_set_column<Master_matrix>::_multiply_source_and_add(const Cell_range& column,
                                                                          const Field_element& val)
{
  if (val == 0u) {
    return false;
  }

  return _generic_add(
      column,
      [&](const Cell& oldCell, Cell* newCell) {
        newCell->set_element(oldCell.get_element());
        operators_->multiply_inplace(newCell->get_element(), val);
      },
      [&](Cell* targetCell, const Cell& sourceCell) {
        operators_->multiply_and_add_inplace_back(sourceCell.get_element(), val, targetCell->get_element());
      });
}

template <class Master_matrix>
template <class Cell_range, typename F1, typename F2>
inline bool Unordered_set_column<Master_matrix>::_generic_add(const Cell_range& source,
                                                              F1&& process_source,
                                                              F2&& update_target)
{
  bool pivotIsZeroed = false;

  for (const Cell& cell : source) {
    Cell* newCell;
    if constexpr (Master_matrix::Option_list::has_row_access) {
      newCell = cellPool_->construct(RA_opt::columnIndex_, cell.get_row_index());
    } else {
      newCell = cellPool_->construct(cell.get_row_index());
    }
    auto res = column_.insert(newCell);
    if (res.second) {
      process_source(cell, newCell);
      if constexpr (Master_matrix::Option_list::has_row_access) RA_opt::insert_cell(cell.get_row_index(), newCell);
    } else {
      cellPool_->destroy(newCell);
      if constexpr (Master_matrix::Option_list::is_z2) {
        if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type) {
          if (cell.get_row_index() == Chain_opt::get_pivot()) pivotIsZeroed = true;
        }
        _delete_cell(res.first);
      } else {
        update_target(*res.first, cell);
        if ((*res.first)->get_element() == Field_operators::get_additive_identity()) {
          if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type) {
            if ((*res.first)->get_row_index() == Chain_opt::get_pivot()) pivotIsZeroed = true;
          }
          _delete_cell(res.first);
        } else {
          if constexpr (Master_matrix::Option_list::has_row_access) RA_opt::update_cell(**res.first);
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
 * @tparam Cell_constructor Template parameter of @ref Gudhi::persistence_matrix::Unordered_set_column.
 */
template <class Master_matrix>
struct std::hash<Gudhi::persistence_matrix::Unordered_set_column<Master_matrix>> {
  std::size_t operator()(const Gudhi::persistence_matrix::Unordered_set_column<Master_matrix>& column) const {
    // can't use Gudhi::persistence_matrix::hash_column because unordered
    std::size_t seed = 0;
    for (const auto& cell : column) {
      seed ^= std::hash<unsigned int>()(cell.get_row_index() * static_cast<unsigned int>(cell.get_element()));
    }
    return seed;
  }
};

#endif  // PM_UNORDERED_SET_COLUMN_H

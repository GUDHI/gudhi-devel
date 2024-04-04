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
 * @file heap_column.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Heap_column class. Also defines the std::hash method for @ref Heap_column.
 */

#ifndef PM_HEAP_COLUMN_H
#define PM_HEAP_COLUMN_H

#include <vector>
#include <stdexcept>
#include <type_traits>
#include <algorithm>  //binary_search
#include <utility>    //std::swap, std::move & std::exchange

#include <boost/iterator/indirect_iterator.hpp>

#include <gudhi/Persistence_matrix/columns/cell_constructors.h>

namespace Gudhi {
namespace persistence_matrix {

/**
 * @brief Column class following the @ref PersistenceMatrixColumn concept. Not compatible with row access.
 *
 * Column based on a heap structure. The heap is represented as a vector sorted as a heap. The top of the heap is
 * the cell with the biggest row index. The sum of two columns is lazy: the content of the source is simply inserted
 * into the heap of the target. Therefore the underlying vector can contain several cells with the same row index.
 * The real value of a cell at a row index corresponds to the sum in the coeffcient field of all values with same
 * row index.
 * 
 * @tparam Master_matrix An instanciation of @ref Matrix from which all types and options are deduced.
 * @tparam Cell_constructor Factory of @ref Cell classes.
 */
template <class Master_matrix, class Cell_constructor = New_cell_constructor<typename Master_matrix::Cell_type> >
class Heap_column : public Master_matrix::Column_dimension_option, public Master_matrix::Chain_column_option 
{
 public:
  using Master = Master_matrix;
  using Field_operators = typename Master_matrix::Field_operators;
  using Field_element_type = typename Master_matrix::element_type;
  using index = typename Master_matrix::index;
  using id_index = typename Master_matrix::id_index;
  using dimension_type = typename Master_matrix::dimension_type;
  using Cell = typename Master_matrix::Cell_type;
  using Column_type = std::vector<Cell*>;
  using iterator = boost::indirect_iterator<typename Column_type::iterator>;
  using const_iterator = boost::indirect_iterator<typename Column_type::const_iterator>;
  using reverse_iterator = boost::indirect_iterator<typename Column_type::reverse_iterator>;
  using const_reverse_iterator = boost::indirect_iterator<typename Column_type::const_reverse_iterator>;

  Heap_column(Field_operators* operators = nullptr, Cell_constructor* cellConstructor = nullptr);
  template <class Container_type = typename Master_matrix::boundary_type>
  Heap_column(const Container_type& nonZeroRowIndices, Field_operators* operators, Cell_constructor* cellConstructor);
  template <class Container_type = typename Master_matrix::boundary_type>
  Heap_column(const Container_type& nonZeroChainRowIndices, 
              dimension_type dimension, 
              Field_operators* operators,
              Cell_constructor* cellConstructor);
  Heap_column(const Heap_column& column, Field_operators* operators = nullptr,
              Cell_constructor* cellConstructor = nullptr);
  Heap_column(Heap_column&& column) noexcept;
  ~Heap_column();

  // just for the sake of the interface
  // row containers and column index are ignored as row access is not implemented for heap columns
  template <class Container_type = typename Master_matrix::boundary_type, class Row_container_type>
  Heap_column(index columnIndex, 
              const Container_type& nonZeroRowIndices, 
              Row_container_type* rowContainer,
              Field_operators* operators, 
              Cell_constructor* cellConstructor);
  template <class Container_type = typename Master_matrix::boundary_type, class Row_container_type>
  Heap_column(index columnIndex, 
              const Container_type& nonZeroChainRowIndices, 
              dimension_type dimension,
              Row_container_type* rowContainer, 
              Field_operators* operators, 
              Cell_constructor* cellConstructor);
  template <class Row_container_type>
  Heap_column(const Heap_column& column, 
              index columnIndex, 
              Row_container_type* rowContainer,
              Field_operators* operators = nullptr, 
              Cell_constructor* cellConstructor = nullptr);

  std::vector<Field_element_type> get_content(int columnLength = -1) const;
  bool is_non_zero(id_index rowIndex) const;
  bool is_empty();
  std::size_t size() const;

  template <class Map_type>
  void reorder(const Map_type& valueMap, [[maybe_unused]] index columnIndex = -1);
  void clear();
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
  Heap_column& operator+=(const Cell_range& column);
  Heap_column& operator+=(Heap_column& column);

  Heap_column& operator*=(unsigned int v);

  // this = v * this + column
  template <class Cell_range>
  Heap_column& multiply_and_add(const Field_element_type& val, const Cell_range& column);
  Heap_column& multiply_and_add(const Field_element_type& val, Heap_column& column);
  // this = this + column * v
  template <class Cell_range>
  Heap_column& multiply_and_add(const Cell_range& column, const Field_element_type& val);
  Heap_column& multiply_and_add(Heap_column& column, const Field_element_type& val);

  friend bool operator==(const Heap_column& c1, const Heap_column& c2) {
    if (&c1 == &c2) return true;
    return c1.get_content() == c2.get_content();
  }
  friend bool operator<(const Heap_column& c1, const Heap_column& c2) {
    if (&c1 == &c2) return false;

    auto cont1 = c1.get_content();
    auto cont2 = c2.get_content();

    auto it1 = cont1.begin();
    auto it2 = cont2.begin();
    while (it1 != cont1.end() && it2 != cont2.end()) {
      if (*it1 == 0u && *it2 != 0u) return false;
      if (*it1 != 0u && *it2 == 0u) return true;
      if (*it1 != *it2) return *it1 < *it2;
      ++it1;
      ++it2;
    }
    return it2 != cont2.end();
  }

  // void set_operators(Field_operators* operators){ operators_ = operators; }

  // Disabled with row access.
  Heap_column& operator=(const Heap_column& other);

  friend void swap(Heap_column& col1, Heap_column& col2) {
    swap(static_cast<typename Master_matrix::Column_dimension_option&>(col1),
         static_cast<typename Master_matrix::Column_dimension_option&>(col2));
    swap(static_cast<typename Master_matrix::Chain_column_option&>(col1),
         static_cast<typename Master_matrix::Chain_column_option&>(col2));
    col1.column_.swap(col2.column_);
    std::swap(col1.insertsSinceLastPrune_, col2.insertsSinceLastPrune_);
    std::swap(col1.operators_, col2.operators_);
    std::swap(col1.cellPool_, col2.cellPool_);
  }

 private:
  using dim_opt = typename Master_matrix::Column_dimension_option;
  using chain_opt = typename Master_matrix::Chain_column_option;

  struct {
    bool operator()(const Cell* c1, const Cell* c2) const { return *c1 < *c2; }
  } cellPointerComp_;

  Column_type column_;
  unsigned int insertsSinceLastPrune_;
  Field_operators* operators_;
  Cell_constructor* cellPool_;

  void _prune();
  Cell* _pop_pivot();
  template <class Cell_range>
  bool _add(const Cell_range& column);
  template <class Cell_range>
  bool _multiply_and_add(const Field_element_type& val, const Cell_range& column);
  template <class Cell_range>
  bool _multiply_and_add(const Cell_range& column, const Field_element_type& val);

  void _verifyCellConstructor() {
    if (cellPool_ == nullptr) {
      if constexpr (std::is_same_v<Cell_constructor, New_cell_constructor<typename Master_matrix::Cell_type> >) {
        cellPool_ = &Master_matrix::defaultCellConstructor;
      } else {
        throw std::invalid_argument("Cell constructor pointer cannot be null.");
      }
    }
  }
};

template <class Master_matrix, class Cell_constructor>
inline Heap_column<Master_matrix, Cell_constructor>::Heap_column(Field_operators* operators,
                                                                 Cell_constructor* cellConstructor)
    : dim_opt(), chain_opt(), insertsSinceLastPrune_(0), operators_(operators) 
{
  if (operators_ == nullptr && cellPool_ == nullptr) return;  //to allow default constructor which gives a dummy column
  _verifyCellConstructor();
}

template <class Master_matrix, class Cell_constructor>
template <class Container_type>
inline Heap_column<Master_matrix, Cell_constructor>::Heap_column(const Container_type& nonZeroRowIndices,
                                                                 Field_operators* operators,
                                                                 Cell_constructor* cellConstructor)
    : dim_opt(nonZeroRowIndices.size() == 0 ? 0 : nonZeroRowIndices.size() - 1),
      chain_opt(),
      column_(nonZeroRowIndices.size(), nullptr),
      insertsSinceLastPrune_(0),
      operators_(operators),
      cellPool_(cellConstructor) 
{
  static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type,
                "Constructor not available for chain columns, please specify the dimension of the chain.");

  _verifyCellConstructor();

  index i = 0;
  if constexpr (Master_matrix::Option_list::is_z2) {
    for (id_index id : nonZeroRowIndices) {
      column_[i++] = cellPool_->construct(id);
    }
  } else {
    for (const auto& p : nonZeroRowIndices) {
      column_[i] = cellPool_->construct(p.first);
      column_[i++]->set_element(operators_->get_value(p.second));
    }
  }
  std::make_heap(column_.begin(), column_.end(), cellPointerComp_);
}

template <class Master_matrix, class Cell_constructor>
template <class Container_type>
inline Heap_column<Master_matrix, Cell_constructor>::Heap_column(const Container_type& nonZeroRowIndices,
                                                                 dimension_type dimension, 
                                                                 Field_operators* operators,
                                                                 Cell_constructor* cellConstructor)
    : dim_opt(dimension),
      chain_opt([&] {
        if constexpr (Master_matrix::Option_list::is_z2) {
          return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : *std::prev(nonZeroRowIndices.end());
        } else {
          return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : std::prev(nonZeroRowIndices.end())->first;
        }
      }()),
      column_(nonZeroRowIndices.size(), nullptr),
      insertsSinceLastPrune_(0),
      operators_(operators),
      cellPool_(cellConstructor) 
{
  _verifyCellConstructor();

  index i = 0;
  if constexpr (Master_matrix::Option_list::is_z2) {
    for (id_index id : nonZeroRowIndices) {
      column_[i++] = cellPool_->construct(id);
    }
  } else {
    for (const auto& p : nonZeroRowIndices) {
      column_[i] = cellPool_->construct(p.first);
      column_[i++]->set_element(operators_->get_value(p.second));
    }
  }
  std::make_heap(column_.begin(), column_.end(), cellPointerComp_);
}

template <class Master_matrix, class Cell_constructor>
inline Heap_column<Master_matrix, Cell_constructor>::Heap_column(const Heap_column& column, 
                                                                 Field_operators* operators,
                                                                 Cell_constructor* cellConstructor)
    : dim_opt(static_cast<const dim_opt&>(column)),
      chain_opt(static_cast<const chain_opt&>(column)),
      column_(column.column_.size(), nullptr),
      insertsSinceLastPrune_(0),
      operators_(operators == nullptr ? column.operators_ : operators),
      cellPool_(cellConstructor == nullptr ? column.cellPool_ : cellConstructor) 
{
  static_assert(!Master_matrix::Option_list::has_row_access,
                "Simple copy constructor not available when row access option enabled. Please specify the new column "
                "index and the row container.");

  index i = 0;
  for (const Cell* cell : column.column_) {
    if constexpr (Master_matrix::Option_list::is_z2) {
      column_[i++] = cellPool_->construct(cell->get_row_index());
    } else {
      column_[i] = cellPool_->construct(cell->get_row_index());
      column_[i++]->set_element(cell->get_element());
    }
  }
  // column.column_ already ordered as a heap, so no need of make_heap.
}

template <class Master_matrix, class Cell_constructor>
inline Heap_column<Master_matrix, Cell_constructor>::Heap_column(Heap_column&& column) noexcept
    : dim_opt(std::move(static_cast<dim_opt&>(column))),
      chain_opt(std::move(static_cast<chain_opt&>(column))),
      column_(std::move(column.column_)),
      insertsSinceLastPrune_(std::exchange(column.insertsSinceLastPrune_, 0)),
      operators_(std::exchange(column.operators_, nullptr)),
      cellPool_(std::exchange(column.cellPool_, nullptr)) 
{}

template <class Master_matrix, class Cell_constructor>
template <class Container_type, class Row_container_type>
inline Heap_column<Master_matrix, Cell_constructor>::Heap_column(index columnIndex,
                                                                 const Container_type& nonZeroRowIndices,
                                                                 Row_container_type* rowContainer,
                                                                 Field_operators* operators,
                                                                 Cell_constructor* cellConstructor)
    : dim_opt(nonZeroRowIndices.size() == 0 ? 0 : nonZeroRowIndices.size() - 1),
      chain_opt([&] {
        if constexpr (Master_matrix::Option_list::is_z2) {
          return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : *std::prev(nonZeroRowIndices.end());
        } else {
          return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : std::prev(nonZeroRowIndices.end())->first;
        }
      }()),
      column_(nonZeroRowIndices.size(), nullptr),
      insertsSinceLastPrune_(0),
      operators_(operators),
      cellPool_(cellConstructor) 
{
  static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type,
                "Constructor not available for chain columns, please specify the dimension of the chain.");

  _verifyCellConstructor();

  index i = 0;
  if constexpr (Master_matrix::Option_list::is_z2) {
    for (id_index id : nonZeroRowIndices) {
      column_[i++] = cellPool_->construct(id);
    }
  } else {
    for (const auto& p : nonZeroRowIndices) {
      column_[i] = cellPool_->construct(p.first);
      column_[i++]->set_element(operators_->get_value(p.second));
    }
  }
  std::make_heap(column_.begin(), column_.end(), cellPointerComp_);
}

template <class Master_matrix, class Cell_constructor>
template <class Container_type, class Row_container_type>
inline Heap_column<Master_matrix, Cell_constructor>::Heap_column(
    index columnIndex, 
    const Container_type& nonZeroRowIndices, 
    dimension_type dimension,
    Row_container_type* rowContainer, 
    Field_operators* operators, 
    Cell_constructor* cellConstructor)
    : dim_opt(dimension),
      chain_opt([&] {
        if constexpr (Master_matrix::Option_list::is_z2) {
          return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : *std::prev(nonZeroRowIndices.end());
        } else {
          return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : std::prev(nonZeroRowIndices.end())->first;
        }
      }()),
      column_(nonZeroRowIndices.size(), nullptr),
      insertsSinceLastPrune_(0),
      operators_(operators),
      cellPool_(cellConstructor) 
{
  _verifyCellConstructor();

  index i = 0;
  if constexpr (Master_matrix::Option_list::is_z2) {
    for (id_index id : nonZeroRowIndices) {
      column_[i++] = cellPool_->construct(id);
    }
  } else {
    for (const auto& p : nonZeroRowIndices) {
      column_[i] = cellPool_->construct(p.first);
      column_[i++]->set_element(operators_->get_value(p.second));
    }
  }
  std::make_heap(column_.begin(), column_.end(), cellPointerComp_);
}

template <class Master_matrix, class Cell_constructor>
template <class Row_container_type>
inline Heap_column<Master_matrix, Cell_constructor>::Heap_column(const Heap_column& column, index columnIndex,
                                                                 Row_container_type* rowContainer,
                                                                 Field_operators* operators,
                                                                 Cell_constructor* cellConstructor)
    : dim_opt(static_cast<const dim_opt&>(column)),
      chain_opt(static_cast<const chain_opt&>(column)),
      column_(column.column_.size(), nullptr),
      insertsSinceLastPrune_(0),
      operators_(operators == nullptr ? column.operators_ : operators),
      cellPool_(cellConstructor == nullptr ? column.cellPool_ : cellConstructor) 
{
  index i = 0;
  for (const Cell* cell : column.column_) {
    if constexpr (Master_matrix::Option_list::is_z2) {
      column_[i++] = cellPool_->construct(cell->get_row_index());
    } else {
      column_[i] = cellPool_->construct(cell->get_row_index());
      column_[i++]->set_element(cell->get_element());
    }
  }
  // column.column_ already ordered as a heap, so no need of make_heap.
}

template <class Master_matrix, class Cell_constructor>
inline Heap_column<Master_matrix, Cell_constructor>::~Heap_column() 
{
  for (auto* cell : column_) {
    cellPool_->destroy(cell);
  }
}

template <class Master_matrix, class Cell_constructor>
inline std::vector<typename Heap_column<Master_matrix, Cell_constructor>::Field_element_type>
Heap_column<Master_matrix, Cell_constructor>::get_content(int columnLength) const 
{
  bool pivotLength = (columnLength < 0);
  if (columnLength < 0 && column_.size() > 0)
    columnLength = column_.front()->get_row_index() + 1;
  else if (columnLength < 0)
    return std::vector<Field_element_type>();

  std::vector<Field_element_type> container(columnLength, 0);
  for (auto it = column_.begin(); it != column_.end(); ++it) {
    if ((*it)->get_row_index() < static_cast<id_index>(columnLength)) {
      if constexpr (Master_matrix::Option_list::is_z2) {
        container[(*it)->get_row_index()] = !container[(*it)->get_row_index()];
      } else {
        container[(*it)->get_row_index()] = operators_->add(container[(*it)->get_row_index()], (*it)->get_element());
      }
    }
  }

  if (pivotLength) {
    while (!container.empty() && container.back() == 0u) container.pop_back();
  }

  return container;
}

template <class Master_matrix, class Cell_constructor>
inline bool Heap_column<Master_matrix, Cell_constructor>::is_non_zero(id_index rowIndex) const 
{
  if constexpr (Master_matrix::Option_list::is_z2) {
    bool c = false;
    for (const Cell* cell : column_) {
      if (cell->get_row_index() == rowIndex) c = !c;
    }
    return c;
  } else {
    Field_element_type c(0);
    for (const Cell* cell : column_) {
      if (cell->get_row_index() == rowIndex) c = operators_->add(c, cell->get_element());
    }
    return c != Field_operators::get_additive_identity();
  }
}

template <class Master_matrix, class Cell_constructor>
inline bool Heap_column<Master_matrix, Cell_constructor>::is_empty() 
{
  Cell* pivot = _pop_pivot();
  if (pivot != nullptr) {
    column_.push_back(pivot);
    std::push_heap(column_.begin(), column_.end(), cellPointerComp_);
    return false;
  }
  return true;
}

template <class Master_matrix, class Cell_constructor>
inline std::size_t Heap_column<Master_matrix, Cell_constructor>::size() const 
{
  return column_.size();
}

template <class Master_matrix, class Cell_constructor>
template <class Map_type>
inline void Heap_column<Master_matrix, Cell_constructor>::reorder(const Map_type& valueMap,
                                                                  [[maybe_unused]] index columnIndex) 
{
  static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type,
                "Method not available for chain columns.");

  Column_type tempCol;
  Cell* pivot = _pop_pivot();
  while (pivot != nullptr) {
    pivot->set_row_index(valueMap.at(pivot->get_row_index()));
    tempCol.push_back(pivot);
    pivot = _pop_pivot();
  }
  column_.swap(tempCol);
  std::make_heap(column_.begin(), column_.end(), cellPointerComp_);

  insertsSinceLastPrune_ = 0;
}

template <class Master_matrix, class Cell_constructor>
inline void Heap_column<Master_matrix, Cell_constructor>::clear() 
{
  static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type,
                "Method not available for chain columns as a base element should not be empty.");

  for (auto* cell : column_) {
    cellPool_->destroy(cell);
  }

  column_.clear();
  insertsSinceLastPrune_ = 0;
}

template <class Master_matrix, class Cell_constructor>
inline void Heap_column<Master_matrix, Cell_constructor>::clear(id_index rowIndex) 
{
  static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type,
                "Method not available for chain columns.");

  column_.erase(std::remove_if(column_.begin(), column_.end(),
                               [rowIndex](const Cell* c) { return c->get_row_index() == rowIndex; }),
                column_.end());
  std::make_heap(column_.begin(), column_.end(), cellPointerComp_);
}

template <class Master_matrix, class Cell_constructor>
inline typename Heap_column<Master_matrix, Cell_constructor>::id_index
Heap_column<Master_matrix, Cell_constructor>::get_pivot() 
{
  static_assert(Master_matrix::isNonBasic,
                "Method not available for base columns.");  // could technically be, but is the notion usefull then?

  if constexpr (Master_matrix::Option_list::is_of_boundary_type) {
    Cell* pivot = _pop_pivot();
    if (pivot != nullptr) {
      column_.push_back(pivot);
      std::push_heap(column_.begin(), column_.end(), cellPointerComp_);
      return pivot->get_row_index();
    }
    return -1;
  } else {
    return chain_opt::get_pivot();
  }
}

template <class Master_matrix, class Cell_constructor>
inline typename Heap_column<Master_matrix, Cell_constructor>::Field_element_type
Heap_column<Master_matrix, Cell_constructor>::get_pivot_value() 
{
  static_assert(Master_matrix::isNonBasic,
                "Method not available for base columns.");  // could technically be, but is the notion usefull then?

  if constexpr (Master_matrix::Option_list::is_z2) {
    return 1;
  } else {
    if constexpr (Master_matrix::Option_list::is_of_boundary_type) {
      Cell* pivot = _pop_pivot();
      if (pivot != nullptr) {
        column_.push_back(pivot);
        std::push_heap(column_.begin(), column_.end(), cellPointerComp_);
        return pivot->get_element();
      }
      return 0;
    } else {
      Field_element_type sum(0);
      if (chain_opt::get_pivot() == -1) return sum;
      for (const Cell* cell : column_) {
        if (cell->get_row_index() == chain_opt::get_pivot()) sum = operators_->add(sum, cell->get_element());
      }
      return sum;  // should not be 0 if properly used.
    }
  }
}

template <class Master_matrix, class Cell_constructor>
inline typename Heap_column<Master_matrix, Cell_constructor>::iterator
Heap_column<Master_matrix, Cell_constructor>::begin() noexcept 
{
  return column_.begin();
}

template <class Master_matrix, class Cell_constructor>
inline typename Heap_column<Master_matrix, Cell_constructor>::const_iterator
Heap_column<Master_matrix, Cell_constructor>::begin() const noexcept 
{
  return column_.begin();
}

template <class Master_matrix, class Cell_constructor>
inline typename Heap_column<Master_matrix, Cell_constructor>::iterator
Heap_column<Master_matrix, Cell_constructor>::end() noexcept 
{
  return column_.end();
}

template <class Master_matrix, class Cell_constructor>
inline typename Heap_column<Master_matrix, Cell_constructor>::const_iterator
Heap_column<Master_matrix, Cell_constructor>::end() const noexcept 
{
  return column_.end();
}

template <class Master_matrix, class Cell_constructor>
inline typename Heap_column<Master_matrix, Cell_constructor>::reverse_iterator
Heap_column<Master_matrix, Cell_constructor>::rbegin() noexcept 
{
  return column_.rbegin();
}

template <class Master_matrix, class Cell_constructor>
inline typename Heap_column<Master_matrix, Cell_constructor>::const_reverse_iterator
Heap_column<Master_matrix, Cell_constructor>::rbegin() const noexcept 
{
  return column_.rbegin();
}

template <class Master_matrix, class Cell_constructor>
inline typename Heap_column<Master_matrix, Cell_constructor>::reverse_iterator
Heap_column<Master_matrix, Cell_constructor>::rend() noexcept 
{
  return column_.rend();
}

template <class Master_matrix, class Cell_constructor>
inline typename Heap_column<Master_matrix, Cell_constructor>::const_reverse_iterator
Heap_column<Master_matrix, Cell_constructor>::rend() const noexcept 
{
  return column_.rend();
}

template <class Master_matrix, class Cell_constructor>
template <class Cell_range>
inline Heap_column<Master_matrix, Cell_constructor>& Heap_column<Master_matrix, Cell_constructor>::operator+=(
    const Cell_range& column) 
{
  static_assert((!Master_matrix::isNonBasic || std::is_same_v<Cell_range, Heap_column>),
                "For boundary columns, the range has to be a column of same type to help ensure the validity of the "
                "base element.");  // could be removed, if we give the responsability to the user.
  static_assert((!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type),
                "For chain columns, the given column cannot be constant.");

  _add(column);

  return *this;
}

template <class Master_matrix, class Cell_constructor>
inline Heap_column<Master_matrix, Cell_constructor>& Heap_column<Master_matrix, Cell_constructor>::operator+=(
    Heap_column& column) 
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

template <class Master_matrix, class Cell_constructor>
inline Heap_column<Master_matrix, Cell_constructor>& Heap_column<Master_matrix, Cell_constructor>::operator*=(
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
      cell->get_element() = operators_->multiply(cell->get_element(), val);
    }
  }

  return *this;
}

template <class Master_matrix, class Cell_constructor>
template <class Cell_range>
inline Heap_column<Master_matrix, Cell_constructor>& Heap_column<Master_matrix, Cell_constructor>::multiply_and_add(
    const Field_element_type& val, const Cell_range& column) 
{
  static_assert((!Master_matrix::isNonBasic || std::is_same_v<Cell_range, Heap_column>),
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
    _multiply_and_add(val, column);
  }

  return *this;
}

template <class Master_matrix, class Cell_constructor>
inline Heap_column<Master_matrix, Cell_constructor>& Heap_column<Master_matrix, Cell_constructor>::multiply_and_add(
    const Field_element_type& val, Heap_column& column) 
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
      if (_multiply_and_add(val, column)) {
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
      _multiply_and_add(val, column);
    }
  }

  return *this;
}

template <class Master_matrix, class Cell_constructor>
template <class Cell_range>
inline Heap_column<Master_matrix, Cell_constructor>& Heap_column<Master_matrix, Cell_constructor>::multiply_and_add(
    const Cell_range& column, const Field_element_type& val) 
{
  static_assert((!Master_matrix::isNonBasic || std::is_same_v<Cell_range, Heap_column>),
                "For boundary columns, the range has to be a column of same type to help ensure the validity of the "
                "base element.");  // could be removed, if we give the responsability to the user.
  static_assert((!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type),
                "For chain columns, the given column cannot be constant.");

  if constexpr (Master_matrix::Option_list::is_z2) {
    if (val) {
      _add(column);
    }
  } else {
    _multiply_and_add(column, val);
  }

  return *this;
}

template <class Master_matrix, class Cell_constructor>
inline Heap_column<Master_matrix, Cell_constructor>& Heap_column<Master_matrix, Cell_constructor>::multiply_and_add(
    Heap_column& column, const Field_element_type& val) 
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
      if (_multiply_and_add(column, val)) {
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
      _multiply_and_add(column, val);
    }
  }

  return *this;
}

template <class Master_matrix, class Cell_constructor>
inline Heap_column<Master_matrix, Cell_constructor>& Heap_column<Master_matrix, Cell_constructor>::operator=(
    const Heap_column& other) 
{
  static_assert(!Master_matrix::Option_list::has_row_access, "= assignement not enabled with row access option.");

  dim_opt::operator=(other);
  chain_opt::operator=(other);

  while (column_.size() > other.column_.size()) {
    if (column_.back() != nullptr) cellPool_->destroy(column_.back());
    column_.pop_back();
  }

  column_.resize(other.column_.size(), nullptr);
  index i = 0;
  for (const Cell* cell : other.column_) {
    if (column_[i] != nullptr) {
      cellPool_->destroy(column_[i]);
    }
    if constexpr (Master_matrix::Option_list::is_z2) {
      column_[i++] = other.cellPool_->construct(cell->get_row_index());
    } else {
      column_[i] = other.cellPool_->construct(cell->get_row_index());
      column_[i++]->set_element(cell->get_element());
    }
  }
  insertsSinceLastPrune_ = other.insertsSinceLastPrune_;
  operators_ = other.operators_;
  cellPool_ = other.cellPool_;

  return *this;
}

template <class Master_matrix, class Cell_constructor>
inline void Heap_column<Master_matrix, Cell_constructor>::_prune() 
{
  if (insertsSinceLastPrune_ == 0) return;

  Column_type tempCol;
  Cell* pivot = _pop_pivot();
  while (pivot != nullptr) {
    tempCol.push_back(pivot);
    pivot = _pop_pivot();
  }
  column_.swap(tempCol);
  std::make_heap(column_.begin(), column_.end(), cellPointerComp_);

  insertsSinceLastPrune_ = 0;
}

template <class Master_matrix, class Cell_constructor>
inline typename Heap_column<Master_matrix, Cell_constructor>::Cell*
Heap_column<Master_matrix, Cell_constructor>::_pop_pivot() 
{
  if (column_.empty()) {
    return nullptr;
  }

  Cell* pivot = column_.front();
  std::pop_heap(column_.begin(), column_.end(), cellPointerComp_);
  column_.pop_back();
  if constexpr (Master_matrix::Option_list::is_z2) {
    while (!column_.empty() && column_.front()->get_row_index() == pivot->get_row_index()) {
      std::pop_heap(column_.begin(), column_.end(), cellPointerComp_);
      column_.pop_back();

      if (column_.empty()) {
        return nullptr;
      }
      pivot = column_.front();
      std::pop_heap(column_.begin(), column_.end(), cellPointerComp_);
      column_.pop_back();
    }
  } else {
    while (!column_.empty() && column_.front()->get_row_index() == pivot->get_row_index()) {
      pivot->get_element() = operators_->add(pivot->get_element(), column_.front()->get_element());
      std::pop_heap(column_.begin(), column_.end(), cellPointerComp_);
      column_.pop_back();
    }

    if (pivot->get_element() == Field_operators::get_additive_identity()) return _pop_pivot();
  }

  return pivot;
}

template <class Master_matrix, class Cell_constructor>
template <class Cell_range>
inline bool Heap_column<Master_matrix, Cell_constructor>::_add(const Cell_range& column) 
{
  if (column.begin() == column.end()) return false;
  if (column_.empty()) {  // chain should never enter here.
    column_.resize(column.size());
    index i = 0;
    for (const Cell& cell : column) {
      if constexpr (Master_matrix::Option_list::is_z2) {
        column_[i++] = cellPool_->construct(cell.get_row_index());
      } else {
        column_[i] = cellPool_->construct(cell.get_row_index());
        column_[i++]->set_element(cell.get_element());
      }
    }
    insertsSinceLastPrune_ = column_.size();
    return true;
  }

  Field_element_type pivotVal(1);

  if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type)
    pivotVal = get_pivot_value();

  for (const Cell& cell : column) {
    ++insertsSinceLastPrune_;
    if constexpr (Master_matrix::Option_list::is_z2) {
      if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type) {
        if (cell.get_row_index() == chain_opt::get_pivot()) pivotVal = !pivotVal;
      }
      column_.push_back(cellPool_->construct(cell.get_row_index()));
    } else {
      if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type) {
        if (cell.get_row_index() == chain_opt::get_pivot()) pivotVal = operators_->add(pivotVal, cell.get_element());
      }
      column_.push_back(cellPool_->construct(cell.get_row_index()));
      column_.back()->set_element(cell.get_element());
    }
    std::push_heap(column_.begin(), column_.end(), cellPointerComp_);
  }

  if (2 * insertsSinceLastPrune_ > column_.size()) _prune();

  return pivotVal == Field_operators::get_additive_identity();
}

template <class Master_matrix, class Cell_constructor>
template <class Cell_range>
inline bool Heap_column<Master_matrix, Cell_constructor>::_multiply_and_add(const Field_element_type& val,
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
      column_[i] = cellPool_->construct(cell.get_row_index());
      column_[i++]->set_element(cell.get_element());
    }
    insertsSinceLastPrune_ = column_.size();
    return true;
  }

  Field_element_type pivotVal(0);

  for (Cell* cell : column_) {
    cell->get_element() = operators_->multiply(cell->get_element(), val);
    if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type) {
      if (cell->get_row_index() == chain_opt::get_pivot()) pivotVal = operators_->add(pivotVal, cell->get_element());
    }
  }

  for (const Cell& cell : column) {
    ++insertsSinceLastPrune_;
    if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type) {
      if (cell.get_row_index() == chain_opt::get_pivot()) pivotVal = operators_->add(pivotVal, cell.get_element());
    }
    column_.push_back(cellPool_->construct(cell.get_row_index()));
    column_.back()->set_element(cell.get_element());
    std::push_heap(column_.begin(), column_.end(), cellPointerComp_);
  }

  if (2 * insertsSinceLastPrune_ > column_.size()) _prune();

  return pivotVal == Field_operators::get_additive_identity();
}

template <class Master_matrix, class Cell_constructor>
template <class Cell_range>
inline bool Heap_column<Master_matrix, Cell_constructor>::_multiply_and_add(const Cell_range& column,
                                                                            const Field_element_type& val) 
{
  if (val == 0u || column.begin() == column.end()) {
    return false;
  }
  if (column_.empty()) {  // chain should never enter here.
    column_.resize(column.size());
    index i = 0;
    for (const Cell& cell : column) {
      column_[i] = cellPool_->construct(cell.get_row_index());
      column_[i++]->set_element(cell.get_element());
    }
    insertsSinceLastPrune_ = column_.size();
    return true;
  }

  Field_element_type pivotVal(1);

  if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type)
    pivotVal = get_pivot_value();

  for (const Cell& cell : column) {
    ++insertsSinceLastPrune_;
    if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type) {
      if (cell.get_row_index() == chain_opt::get_pivot())
        pivotVal = operators_->multiply_and_add(cell.get_element(), val, pivotVal);
    }
    column_.push_back(cellPool_->construct(cell.get_row_index()));
    column_.back()->set_element(operators_->multiply(cell.get_element(), val));
    std::push_heap(column_.begin(), column_.end(), cellPointerComp_);
  }

  if (2 * insertsSinceLastPrune_ > column_.size()) _prune();

  return pivotVal == 0u;
}

}  // namespace persistence_matrix
}  // namespace Gudhi

/**
 * @brief Hash method for @ref Gudhi::persistence_matrix::Heap_column.
 * 
 * @tparam Master_matrix Template parameter of @ref Gudhi::persistence_matrix::Heap_column.
 * @tparam Cell_constructor Template parameter of @ref Gudhi::persistence_matrix::Heap_column.
 */
template <class Master_matrix, class Cell_constructor>
struct std::hash<Gudhi::persistence_matrix::Heap_column<Master_matrix, Cell_constructor> > 
{
  size_t operator()(const Gudhi::persistence_matrix::Heap_column<Master_matrix, Cell_constructor>& column) const {
    std::size_t seed = 0;
    unsigned int i = 0;
    for (bool val : column.get_content()) {
      seed ^= std::hash<unsigned int>()(i++ * val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
  }
};

#endif  // PM_HEAP_COLUMN_H

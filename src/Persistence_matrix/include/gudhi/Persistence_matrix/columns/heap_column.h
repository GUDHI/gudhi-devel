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
 * @brief Contains the @ref Gudhi::persistence_matrix::Heap_column class. Also defines the std::hash method
 * for @ref Gudhi::persistence_matrix::Heap_column.
 */

#ifndef PM_HEAP_COLUMN_H
#define PM_HEAP_COLUMN_H

#include <vector>
#include <stdexcept>
#include <type_traits>
#include <algorithm>  //binary_search
#include <utility>    //std::swap, std::move & std::exchange

#include <boost/iterator/indirect_iterator.hpp>

#include <gudhi/Persistence_matrix/allocators/entry_constructors.h>

namespace Gudhi {
namespace persistence_matrix {

/**
 * @class Heap_column heap_column.h gudhi/Persistence_matrix/columns/heap_column.h
 * @ingroup persistence_matrix
 *
 * @brief Column class following the @ref PersistenceMatrixColumn concept. Not compatible with row access.
 *
 * Column based on a heap structure. The heap is represented as a vector sorted as a heap. The top of the heap is
 * the entry with the biggest row index. The sum of two columns is lazy: the content of the source is simply inserted
 * into the heap of the target. Therefore the underlying vector can contain several entries with the same row index.
 * The real value of an entry at a row index corresponds to the sum in the coefficient field of all values with same
 * row index. Additionally, the given entry range added into the heap does not need to be somehow ordered.
 *
 * @tparam Master_matrix An instantiation of @ref Matrix from which all types and options are deduced.
 * @tparam Entry_constructor Factory of @ref Entry classes.
 */
template <class Master_matrix>
class Heap_column : public Master_matrix::Column_dimension_option, public Master_matrix::Chain_column_option
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

  Heap_column(Column_settings* colSettings = nullptr);
  template <class Container = typename Master_matrix::Boundary>
  Heap_column(const Container& nonZeroRowIndices, Column_settings* colSettings);
  template <class Container = typename Master_matrix::Boundary>
  Heap_column(const Container& nonZeroChainRowIndices, Dimension dimension, Column_settings* colSettings);
  Heap_column(const Heap_column& column, Column_settings* colSettings = nullptr);
  Heap_column(Heap_column&& column) noexcept;
  ~Heap_column();

  // just for the sake of the interface
  // row containers and column index are ignored as row access is not implemented for heap columns
  template <class Container = typename Master_matrix::Boundary, class Row_container>
  Heap_column(Index columnIndex,
              const Container& nonZeroRowIndices,
              Row_container* rowContainer,
              Column_settings* colSettings);
  template <class Container = typename Master_matrix::Boundary, class Row_container>
  Heap_column(Index columnIndex,
              const Container& nonZeroChainRowIndices,
              Dimension dimension,
              Row_container* rowContainer,
              Column_settings* colSettings);
  template <class Row_container>
  Heap_column(const Heap_column& column,
              Index columnIndex,
              Row_container* rowContainer,
              Column_settings* colSettings = nullptr);

  std::vector<Field_element> get_content(int columnLength = -1) const;
  bool is_non_zero(ID_index rowIndex) const;
  bool is_empty();
  std::size_t size() const;

  template <class Row_index_map>
  void reorder(const Row_index_map& valueMap, [[maybe_unused]] Index columnIndex = -1);
  void clear();
  void clear(ID_index rowIndex);

  ID_index get_pivot();
  Field_element get_pivot_value();

  iterator begin() noexcept;
  const_iterator begin() const noexcept;
  iterator end() noexcept;
  const_iterator end() const noexcept;
  reverse_iterator rbegin() noexcept;
  const_reverse_iterator rbegin() const noexcept;
  reverse_iterator rend() noexcept;
  const_reverse_iterator rend() const noexcept;

  template <class Entry_range>
  Heap_column& operator+=(const Entry_range& column);
  Heap_column& operator+=(Heap_column& column);

  Heap_column& operator*=(unsigned int v);

  // this = v * this + column
  template <class Entry_range>
  Heap_column& multiply_target_and_add(const Field_element& val, const Entry_range& column);
  Heap_column& multiply_target_and_add(const Field_element& val, Heap_column& column);
  // this = this + column * v
  template <class Entry_range>
  Heap_column& multiply_source_and_add(const Entry_range& column, const Field_element& val);
  Heap_column& multiply_source_and_add(Heap_column& column, const Field_element& val);

  void push_back(const Entry& entry);

  std::size_t compute_hash_value();

  friend bool operator==(const Heap_column& c1, const Heap_column& c2) {
    if (&c1 == &c2) return true;

    Heap_column cc1(c1), cc2(c2);
    Entry* p1 = cc1._pop_pivot();
    Entry* p2 = cc2._pop_pivot();
    while (p1 != nullptr && p2 != nullptr) {
      if (p1->get_row_index() != p2->get_row_index()) {
        c1.entryPool_->destroy(p1);
        c2.entryPool_->destroy(p2);
        return false;
      }
      if constexpr (!Master_matrix::Option_list::is_z2) {
        if (p1->get_element() != p2->get_element()) {
          c1.entryPool_->destroy(p1);
          c2.entryPool_->destroy(p2);
          return false;
        }
      }
      c1.entryPool_->destroy(p1);
      c2.entryPool_->destroy(p2);
      p1 = cc1._pop_pivot();
      p2 = cc2._pop_pivot();
    }

    if (p1 == nullptr && p2 == nullptr) return true;
    if (p1 != nullptr) {
      c1.entryPool_->destroy(p1);
      return false;
    }
    c2.entryPool_->destroy(p2);
    return false;
  }
  friend bool operator<(const Heap_column& c1, const Heap_column& c2) {
    if (&c1 == &c2) return false;

    // lexicographical order but starting from last value and not first
    Heap_column cc1(c1), cc2(c2);
    Entry* p1 = cc1._pop_pivot();
    Entry* p2 = cc2._pop_pivot();
    while (p1 != nullptr && p2 != nullptr) {
      if (p1->get_row_index() > p2->get_row_index()) {
        c1.entryPool_->destroy(p1);
        c2.entryPool_->destroy(p2);
        return false;
      }
      if (p1->get_row_index() < p2->get_row_index()) {
        c1.entryPool_->destroy(p1);
        c2.entryPool_->destroy(p2);
        return true;
      }
      if constexpr (!Master_matrix::Option_list::is_z2) {
        if (p1->get_element() > p2->get_element()) {
          c1.entryPool_->destroy(p1);
          c2.entryPool_->destroy(p2);
          return false;
        }
        if (p1->get_element() < p2->get_element()) {
          c1.entryPool_->destroy(p1);
          c2.entryPool_->destroy(p2);
          return true;
        }
      }
      c1.entryPool_->destroy(p1);
      c2.entryPool_->destroy(p2);
      p1 = cc1._pop_pivot();
      p2 = cc2._pop_pivot();
    }

    if (p2 == nullptr) {
      c1.entryPool_->destroy(p1);
      return false;
    }
    c2.entryPool_->destroy(p2);
    return true;
  }

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
    std::swap(col1.entryPool_, col2.entryPool_);
  }

 private:
  using Dim_opt = typename Master_matrix::Column_dimension_option;
  using Chain_opt = typename Master_matrix::Chain_column_option;

  struct EntryPointerComp {
    bool operator()(const Entry* c1, const Entry* c2) const { return *c1 < *c2; }
  } entryPointerComp_;

  Column_support column_;
  unsigned int insertsSinceLastPrune_;
  Field_operators* operators_;
  Entry_constructor* entryPool_;

  void _prune();
  Entry* _pop_pivot();
  template <class Entry_range>
  bool _add(const Entry_range& column);
  template <class Entry_range>
  bool _multiply_target_and_add(const Field_element& val, const Entry_range& column);
  template <class Entry_range>
  bool _multiply_source_and_add(const Entry_range& column, const Field_element& val);
};

template <class Master_matrix>
inline Heap_column<Master_matrix>::Heap_column(Column_settings* colSettings)
    : Dim_opt(),
      Chain_opt(),
      insertsSinceLastPrune_(0),
      operators_(nullptr),
      entryPool_(colSettings == nullptr ? nullptr : &(colSettings->entryConstructor))
{
  if (colSettings == nullptr) return;  // to allow default constructor which gives a dummy column

  if constexpr (!Master_matrix::Option_list::is_z2) {
    operators_ = &(colSettings->operators);
  }
}

template <class Master_matrix>
template <class Container>
inline Heap_column<Master_matrix>::Heap_column(const Container& nonZeroRowIndices, Column_settings* colSettings)
    : Dim_opt(nonZeroRowIndices.size() == 0 ? 0 : nonZeroRowIndices.size() - 1),
      Chain_opt(),
      column_(nonZeroRowIndices.size(), nullptr),
      insertsSinceLastPrune_(0),
      operators_(nullptr),
      entryPool_(&(colSettings->entryConstructor))
{
  static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type,
                "Constructor not available for chain columns, please specify the dimension of the chain.");

  if constexpr (!Master_matrix::Option_list::is_z2) {
    operators_ = &(colSettings->operators);
  }

  Index i = 0;
  if constexpr (Master_matrix::Option_list::is_z2) {
    for (ID_index id : nonZeroRowIndices) {
      column_[i++] = entryPool_->construct(id);
    }
  } else {
    for (const auto& p : nonZeroRowIndices) {
      column_[i] = entryPool_->construct(p.first);
      column_[i++]->set_element(operators_->get_value(p.second));
    }
  }
  std::make_heap(column_.begin(), column_.end(), entryPointerComp_);
}

template <class Master_matrix>
template <class Container>
inline Heap_column<Master_matrix>::Heap_column(const Container& nonZeroRowIndices,
                                               Dimension dimension,
                                               Column_settings* colSettings)
    : Dim_opt(dimension),
      Chain_opt([&] {
        if constexpr (Master_matrix::Option_list::is_z2) {
          return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : *std::prev(nonZeroRowIndices.end());
        } else {
          return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : std::prev(nonZeroRowIndices.end())->first;
        }
      }()),
      column_(nonZeroRowIndices.size(), nullptr),
      insertsSinceLastPrune_(0),
      operators_(nullptr),
      entryPool_(&(colSettings->entryConstructor))
{
  Index i = 0;
  if constexpr (Master_matrix::Option_list::is_z2) {
    for (ID_index id : nonZeroRowIndices) {
      column_[i++] = entryPool_->construct(id);
    }
  } else {
    operators_ = &(colSettings->operators);
    for (const auto& p : nonZeroRowIndices) {
      column_[i] = entryPool_->construct(p.first);
      column_[i++]->set_element(operators_->get_value(p.second));
    }
  }
  std::make_heap(column_.begin(), column_.end(), entryPointerComp_);
}

template <class Master_matrix>
inline Heap_column<Master_matrix>::Heap_column(const Heap_column& column, Column_settings* colSettings)
    : Dim_opt(static_cast<const Dim_opt&>(column)),
      Chain_opt(static_cast<const Chain_opt&>(column)),
      column_(column.column_.size(), nullptr),
      insertsSinceLastPrune_(0),
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
      column_[i++] = entryPool_->construct(entry->get_row_index());
    } else {
      column_[i] = entryPool_->construct(entry->get_row_index());
      column_[i++]->set_element(entry->get_element());
    }
  }
  // column.column_ already ordered as a heap, so no need of make_heap.
}

template <class Master_matrix>
inline Heap_column<Master_matrix>::Heap_column(Heap_column&& column) noexcept
    : Dim_opt(std::move(static_cast<Dim_opt&>(column))),
      Chain_opt(std::move(static_cast<Chain_opt&>(column))),
      column_(std::move(column.column_)),
      insertsSinceLastPrune_(std::exchange(column.insertsSinceLastPrune_, 0)),
      operators_(std::exchange(column.operators_, nullptr)),
      entryPool_(std::exchange(column.entryPool_, nullptr))
{}

template <class Master_matrix>
template <class Container, class Row_container>
inline Heap_column<Master_matrix>::Heap_column(Index columnIndex,
                                               const Container& nonZeroRowIndices,
                                               Row_container* rowContainer,
                                               Column_settings* colSettings)
    : Dim_opt(nonZeroRowIndices.size() == 0 ? 0 : nonZeroRowIndices.size() - 1),
      Chain_opt([&] {
        if constexpr (Master_matrix::Option_list::is_z2) {
          return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : *std::prev(nonZeroRowIndices.end());
        } else {
          return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : std::prev(nonZeroRowIndices.end())->first;
        }
      }()),
      column_(nonZeroRowIndices.size(), nullptr),
      insertsSinceLastPrune_(0),
      operators_(nullptr),
      entryPool_(&(colSettings->entryConstructor))
{
  static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type,
                "Constructor not available for chain columns, please specify the dimension of the chain.");

  Index i = 0;
  if constexpr (Master_matrix::Option_list::is_z2) {
    for (ID_index id : nonZeroRowIndices) {
      column_[i++] = entryPool_->construct(id);
    }
  } else {
    operators_ = &(colSettings->operators);
    for (const auto& p : nonZeroRowIndices) {
      column_[i] = entryPool_->construct(p.first);
      column_[i++]->set_element(operators_->get_value(p.second));
    }
  }
  std::make_heap(column_.begin(), column_.end(), entryPointerComp_);
}

template <class Master_matrix>
template <class Container, class Row_container>
inline Heap_column<Master_matrix>::Heap_column(Index columnIndex,
                                               const Container& nonZeroRowIndices,
                                               Dimension dimension,
                                               Row_container* rowContainer,
                                               Column_settings* colSettings)
    : Dim_opt(dimension),
      Chain_opt([&] {
        if constexpr (Master_matrix::Option_list::is_z2) {
          return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : *std::prev(nonZeroRowIndices.end());
        } else {
          return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : std::prev(nonZeroRowIndices.end())->first;
        }
      }()),
      column_(nonZeroRowIndices.size(), nullptr),
      insertsSinceLastPrune_(0),
      operators_(nullptr),
      entryPool_(&(colSettings->entryConstructor))
{
  Index i = 0;
  if constexpr (Master_matrix::Option_list::is_z2) {
    for (ID_index id : nonZeroRowIndices) {
      column_[i++] = entryPool_->construct(id);
    }
  } else {
    operators_ = &(colSettings->operators);
    for (const auto& p : nonZeroRowIndices) {
      column_[i] = entryPool_->construct(p.first);
      column_[i++]->set_element(operators_->get_value(p.second));
    }
  }
  std::make_heap(column_.begin(), column_.end(), entryPointerComp_);
}

template <class Master_matrix>
template <class Row_container>
inline Heap_column<Master_matrix>::Heap_column(const Heap_column& column,
                                               Index columnIndex,
                                               Row_container* rowContainer,
                                               Column_settings* colSettings)
    : Dim_opt(static_cast<const Dim_opt&>(column)),
      Chain_opt(static_cast<const Chain_opt&>(column)),
      column_(column.column_.size(), nullptr),
      insertsSinceLastPrune_(0),
      operators_(colSettings == nullptr ? column.operators_ : nullptr),
      entryPool_(colSettings == nullptr ? column.entryPool_ : &(colSettings->entryConstructor))
{
  if constexpr (!Master_matrix::Option_list::is_z2) {
    if (colSettings != nullptr) operators_ = &(colSettings->operators);
  }

  Index i = 0;
  for (const Entry* entry : column.column_) {
    if constexpr (Master_matrix::Option_list::is_z2) {
      column_[i++] = entryPool_->construct(entry->get_row_index());
    } else {
      column_[i] = entryPool_->construct(entry->get_row_index());
      column_[i++]->set_element(entry->get_element());
    }
  }
  // column.column_ already ordered as a heap, so no need of make_heap.
}

template <class Master_matrix>
inline Heap_column<Master_matrix>::~Heap_column()
{
  for (auto* entry : column_) {
    entryPool_->destroy(entry);
  }
}

template <class Master_matrix>
inline std::vector<typename Heap_column<Master_matrix>::Field_element> Heap_column<Master_matrix>::get_content(
    int columnLength) const
{
  bool pivotLength = (columnLength < 0);
  if (columnLength < 0 && column_.size() > 0)
    columnLength = column_.front()->get_row_index() + 1;
  else if (columnLength < 0)
    return std::vector<Field_element>();

  std::vector<Field_element> container(columnLength, 0);
  for (auto it = column_.begin(); it != column_.end(); ++it) {
    if ((*it)->get_row_index() < static_cast<ID_index>(columnLength)) {
      if constexpr (Master_matrix::Option_list::is_z2) {
        container[(*it)->get_row_index()] = !container[(*it)->get_row_index()];
      } else {
        operators_->add_inplace(container[(*it)->get_row_index()], (*it)->get_element());
      }
    }
  }

  if (pivotLength) {
    while (!container.empty() && container.back() == 0u) container.pop_back();
  }

  return container;
}

template <class Master_matrix>
inline bool Heap_column<Master_matrix>::is_non_zero(ID_index rowIndex) const
{
  Field_element c(0);
  for (const Entry* entry : column_) {
    if (entry->get_row_index() == rowIndex) {
      if constexpr (Master_matrix::Option_list::is_z2)
        c = !c;
      else
        operators_->add_inplace(c, entry->get_element());
    }
  }
  return c != Field_operators::get_additive_identity();
}

template <class Master_matrix>
inline bool Heap_column<Master_matrix>::is_empty()
{
  Entry* pivot = _pop_pivot();
  if (pivot != nullptr) {
    column_.push_back(pivot);
    std::push_heap(column_.begin(), column_.end(), entryPointerComp_);
    return false;
  }
  return true;
}

template <class Master_matrix>
inline std::size_t Heap_column<Master_matrix>::size() const
{
  return column_.size();
}

template <class Master_matrix>
template <class Row_index_map>
inline void Heap_column<Master_matrix>::reorder(const Row_index_map& valueMap, [[maybe_unused]] Index columnIndex)
{
  static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type,
                "Method not available for chain columns.");

  Column_support tempCol;
  Entry* pivot = _pop_pivot();
  while (pivot != nullptr) {
    pivot->set_row_index(valueMap.at(pivot->get_row_index()));
    tempCol.push_back(pivot);
    pivot = _pop_pivot();
  }
  column_.swap(tempCol);
  std::make_heap(column_.begin(), column_.end(), entryPointerComp_);

  insertsSinceLastPrune_ = 0;
}

template <class Master_matrix>
inline void Heap_column<Master_matrix>::clear()
{
  static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type,
                "Method not available for chain columns as a base element should not be empty.");

  for (auto* entry : column_) {
    entryPool_->destroy(entry);
  }

  column_.clear();
  insertsSinceLastPrune_ = 0;
}

template <class Master_matrix>
inline void Heap_column<Master_matrix>::clear(ID_index rowIndex)
{
  static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type,
                "Method not available for chain columns.");

  Column_support tempCol;
  Entry* pivot = _pop_pivot();
  while (pivot != nullptr) {
    if (pivot->get_row_index() != rowIndex) {
      tempCol.push_back(pivot);
    } else {
      entryPool_->destroy(pivot);
    }
    pivot = _pop_pivot();
  }
  column_.swap(tempCol);

  insertsSinceLastPrune_ = 0;
}

template <class Master_matrix>
inline typename Heap_column<Master_matrix>::ID_index Heap_column<Master_matrix>::get_pivot()
{
  static_assert(Master_matrix::isNonBasic,
                "Method not available for base columns.");  // could technically be, but is the notion useful then?

  if constexpr (Master_matrix::Option_list::is_of_boundary_type) {
    Entry* pivot = _pop_pivot();
    if (pivot != nullptr) {
      column_.push_back(pivot);
      std::push_heap(column_.begin(), column_.end(), entryPointerComp_);
      return pivot->get_row_index();
    }
    return -1;
  } else {
    return Chain_opt::get_pivot();
  }
}

template <class Master_matrix>
inline typename Heap_column<Master_matrix>::Field_element Heap_column<Master_matrix>::get_pivot_value()
{
  static_assert(Master_matrix::isNonBasic,
                "Method not available for base columns.");  // could technically be, but is the notion useful then?

  if constexpr (Master_matrix::Option_list::is_z2) {
    return 1;
  } else {
    if constexpr (Master_matrix::Option_list::is_of_boundary_type) {
      Entry* pivot = _pop_pivot();
      if (pivot != nullptr) {
        column_.push_back(pivot);
        std::push_heap(column_.begin(), column_.end(), entryPointerComp_);
        return pivot->get_element();
      }
      return 0;
    } else {
      Field_element sum(0);
      if (Chain_opt::get_pivot() == static_cast<ID_index>(-1)) return sum;
      for (const Entry* entry : column_) {
        if (entry->get_row_index() == Chain_opt::get_pivot()) operators_->add_inplace(sum, entry->get_element());
      }
      return sum;  // should not be 0 if properly used.
    }
  }
}

template <class Master_matrix>
inline typename Heap_column<Master_matrix>::iterator Heap_column<Master_matrix>::begin() noexcept
{
  return column_.begin();
}

template <class Master_matrix>
inline typename Heap_column<Master_matrix>::const_iterator Heap_column<Master_matrix>::begin() const noexcept
{
  return column_.begin();
}

template <class Master_matrix>
inline typename Heap_column<Master_matrix>::iterator Heap_column<Master_matrix>::end() noexcept
{
  return column_.end();
}

template <class Master_matrix>
inline typename Heap_column<Master_matrix>::const_iterator Heap_column<Master_matrix>::end() const noexcept
{
  return column_.end();
}

template <class Master_matrix>
inline typename Heap_column<Master_matrix>::reverse_iterator Heap_column<Master_matrix>::rbegin() noexcept
{
  return column_.rbegin();
}

template <class Master_matrix>
inline typename Heap_column<Master_matrix>::const_reverse_iterator Heap_column<Master_matrix>::rbegin() const noexcept
{
  return column_.rbegin();
}

template <class Master_matrix>
inline typename Heap_column<Master_matrix>::reverse_iterator Heap_column<Master_matrix>::rend() noexcept
{
  return column_.rend();
}

template <class Master_matrix>
inline typename Heap_column<Master_matrix>::const_reverse_iterator Heap_column<Master_matrix>::rend() const noexcept
{
  return column_.rend();
}

template <class Master_matrix>
template <class Entry_range>
inline Heap_column<Master_matrix>& Heap_column<Master_matrix>::operator+=(const Entry_range& column)
{
  static_assert((!Master_matrix::isNonBasic || std::is_same_v<Entry_range, Heap_column>),
                "For boundary columns, the range has to be a column of same type to help ensure the validity of the "
                "base element.");  // could be removed, if we give the responsibility to the user.
  static_assert((!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type),
                "For chain columns, the given column cannot be constant.");

  _add(column);

  return *this;
}

template <class Master_matrix>
inline Heap_column<Master_matrix>& Heap_column<Master_matrix>::operator+=(Heap_column& column)
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
inline Heap_column<Master_matrix>& Heap_column<Master_matrix>::operator*=(unsigned int v)
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
    }
  }

  return *this;
}

template <class Master_matrix>
template <class Entry_range>
inline Heap_column<Master_matrix>& Heap_column<Master_matrix>::multiply_target_and_add(const Field_element& val,
                                                                                       const Entry_range& column)
{
  static_assert((!Master_matrix::isNonBasic || std::is_same_v<Entry_range, Heap_column>),
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
inline Heap_column<Master_matrix>& Heap_column<Master_matrix>::multiply_target_and_add(const Field_element& val,
                                                                                       Heap_column& column)
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
inline Heap_column<Master_matrix>& Heap_column<Master_matrix>::multiply_source_and_add(const Entry_range& column,
                                                                                       const Field_element& val)
{
  static_assert((!Master_matrix::isNonBasic || std::is_same_v<Entry_range, Heap_column>),
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
inline Heap_column<Master_matrix>& Heap_column<Master_matrix>::multiply_source_and_add(Heap_column& column,
                                                                                       const Field_element& val)
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
inline void Heap_column<Master_matrix>::push_back(const Entry& entry)
{
  static_assert(Master_matrix::Option_list::is_of_boundary_type, "`push_back` is not available for Chain matrices.");

  Entry* newEntry = entryPool_->construct(entry.get_row_index());
  if constexpr (!Master_matrix::Option_list::is_z2) {
    newEntry->set_element(operators_->get_value(entry.get_element()));
  }
  column_.push_back(newEntry);
  std::push_heap(column_.begin(), column_.end(), entryPointerComp_);
}

template <class Master_matrix>
inline Heap_column<Master_matrix>& Heap_column<Master_matrix>::operator=(const Heap_column& other)
{
  static_assert(!Master_matrix::Option_list::has_row_access, "= assignment not enabled with row access option.");

  Dim_opt::operator=(other);
  Chain_opt::operator=(other);

  while (column_.size() > other.column_.size()) {
    if (column_.back() != nullptr) entryPool_->destroy(column_.back());
    column_.pop_back();
  }

  column_.resize(other.column_.size(), nullptr);
  Index i = 0;
  for (const Entry* entry : other.column_) {
    if (column_[i] != nullptr) {
      entryPool_->destroy(column_[i]);
    }
    column_[i] = other.entryPool_->construct(entry->get_row_index());
    if constexpr (!Master_matrix::Option_list::is_z2) {
      column_[i]->set_element(entry->get_element());
    }
    ++i;
  }
  insertsSinceLastPrune_ = other.insertsSinceLastPrune_;
  operators_ = other.operators_;
  entryPool_ = other.entryPool_;

  return *this;
}

template <class Master_matrix>
inline std::size_t Heap_column<Master_matrix>::compute_hash_value()
{
  _prune();
  return hash_column(*this);
}

template <class Master_matrix>
inline void Heap_column<Master_matrix>::_prune()
{
  if (insertsSinceLastPrune_ == 0) return;

  Column_support tempCol;
  Entry* pivot = _pop_pivot();
  while (pivot != nullptr) {
    tempCol.push_back(pivot);
    pivot = _pop_pivot();
  }
  column_.swap(tempCol);

  insertsSinceLastPrune_ = 0;
}

template <class Master_matrix>
inline typename Heap_column<Master_matrix>::Entry* Heap_column<Master_matrix>::_pop_pivot()
{
  if (column_.empty()) {
    return nullptr;
  }

  Entry* pivot = column_.front();
  std::pop_heap(column_.begin(), column_.end(), entryPointerComp_);
  column_.pop_back();
  if constexpr (Master_matrix::Option_list::is_z2) {
    while (!column_.empty() && column_.front()->get_row_index() == pivot->get_row_index()) {
      std::pop_heap(column_.begin(), column_.end(), entryPointerComp_);
      entryPool_->destroy(column_.back());
      column_.pop_back();

      entryPool_->destroy(pivot);
      if (column_.empty()) {
        return nullptr;
      }
      pivot = column_.front();
      std::pop_heap(column_.begin(), column_.end(), entryPointerComp_);
      column_.pop_back();
    }
  } else {
    while (!column_.empty() && column_.front()->get_row_index() == pivot->get_row_index()) {
      operators_->add_inplace(pivot->get_element(), column_.front()->get_element());
      std::pop_heap(column_.begin(), column_.end(), entryPointerComp_);
      entryPool_->destroy(column_.back());
      column_.pop_back();
    }

    if (pivot->get_element() == Field_operators::get_additive_identity()) {
      entryPool_->destroy(pivot);
      return _pop_pivot();
    }
  }

  return pivot;
}

template <class Master_matrix>
template <class Entry_range>
inline bool Heap_column<Master_matrix>::_add(const Entry_range& column)
{
  if (column.begin() == column.end()) return false;
  if (column_.empty()) {  // chain should never enter here.
    column_.resize(column.size());
    Index i = 0;
    for (const Entry& entry : column) {
      column_[i] = entryPool_->construct(entry.get_row_index());
      if constexpr (!Master_matrix::Option_list::is_z2) {
        column_[i]->set_element(entry.get_element());
      }
      ++i;
    }
    insertsSinceLastPrune_ = column_.size();
    return true;
  }

  Field_element pivotVal(1);

  if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type)
    pivotVal = get_pivot_value();

  for (const Entry& entry : column) {
    ++insertsSinceLastPrune_;
    if constexpr (Master_matrix::Option_list::is_z2) {
      if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type) {
        if (entry.get_row_index() == Chain_opt::get_pivot()) pivotVal = !pivotVal;
      }
      column_.push_back(entryPool_->construct(entry.get_row_index()));
    } else {
      if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type) {
        if (entry.get_row_index() == Chain_opt::get_pivot()) operators_->add_inplace(pivotVal, entry.get_element());
      }
      column_.push_back(entryPool_->construct(entry.get_row_index()));
      column_.back()->set_element(entry.get_element());
    }
    std::push_heap(column_.begin(), column_.end(), entryPointerComp_);
  }

  if (2 * insertsSinceLastPrune_ > column_.size()) _prune();

  if constexpr (Master_matrix::Option_list::is_z2)
    return !pivotVal;
  else
    return pivotVal == Field_operators::get_additive_identity();
}

template <class Master_matrix>
template <class Entry_range>
inline bool Heap_column<Master_matrix>::_multiply_target_and_add(const Field_element& val, const Entry_range& column)
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
      column_[i] = entryPool_->construct(entry.get_row_index());
      column_[i++]->set_element(entry.get_element());
    }
    insertsSinceLastPrune_ = column_.size();
    return true;
  }

  Field_element pivotVal(0);

  for (Entry* entry : column_) {
    operators_->multiply_inplace(entry->get_element(), val);
    if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type) {
      if (entry->get_row_index() == Chain_opt::get_pivot()) operators_->add_inplace(pivotVal, entry->get_element());
    }
  }

  for (const Entry& entry : column) {
    ++insertsSinceLastPrune_;
    if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type) {
      if (entry.get_row_index() == Chain_opt::get_pivot()) operators_->add_inplace(pivotVal, entry.get_element());
    }
    column_.push_back(entryPool_->construct(entry.get_row_index()));
    column_.back()->set_element(entry.get_element());
    std::push_heap(column_.begin(), column_.end(), entryPointerComp_);
  }

  if (2 * insertsSinceLastPrune_ > column_.size()) _prune();

  if constexpr (Master_matrix::Option_list::is_z2)
    return !pivotVal;
  else
    return pivotVal == Field_operators::get_additive_identity();
}

template <class Master_matrix>
template <class Entry_range>
inline bool Heap_column<Master_matrix>::_multiply_source_and_add(const Entry_range& column, const Field_element& val)
{
  if (val == 0u || column.begin() == column.end()) {
    return false;
  }
  if (column_.empty()) {  // chain should never enter here.
    column_.resize(column.size());
    Index i = 0;
    for (const Entry& entry : column) {
      column_[i] = entryPool_->construct(entry.get_row_index());
      column_[i++]->set_element(entry.get_element());
    }
    insertsSinceLastPrune_ = column_.size();
    return true;
  }

  Field_element pivotVal(1);

  if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type)
    pivotVal = get_pivot_value();

  for (const Entry& entry : column) {
    ++insertsSinceLastPrune_;
    column_.push_back(entryPool_->construct(entry.get_row_index()));
    column_.back()->set_element(entry.get_element());
    operators_->multiply_inplace(column_.back()->get_element(), val);
    if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type) {
      if (entry.get_row_index() == Chain_opt::get_pivot()) {
        operators_->add_inplace(pivotVal, column_.back()->get_element());
      }
    }
    std::push_heap(column_.begin(), column_.end(), entryPointerComp_);
  }

  if (2 * insertsSinceLastPrune_ > column_.size()) _prune();

  return pivotVal == 0u;
}

}  // namespace persistence_matrix
}  // namespace Gudhi

/**
 * @ingroup persistence_matrix
 *
 * @brief Hash method for @ref Gudhi::persistence_matrix::Heap_column.
 *
 * @tparam Master_matrix Template parameter of @ref Gudhi::persistence_matrix::Heap_column.
 * @tparam Entry_constructor Template parameter of @ref Gudhi::persistence_matrix::Heap_column.
 */
template <class Master_matrix>
struct std::hash<Gudhi::persistence_matrix::Heap_column<Master_matrix> > {
  size_t operator()(Gudhi::persistence_matrix::Heap_column<Master_matrix>& column) const {
    return column.compute_hash_value();
  }
};

#endif  // PM_HEAP_COLUMN_H

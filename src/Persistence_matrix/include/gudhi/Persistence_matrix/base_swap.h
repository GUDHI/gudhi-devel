/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022-23 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file base_swap.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_matrix::Base_swap class and
 * @ref Gudhi::persistence_matrix::Dummy_base_swap structure.
 */

#ifndef PM_BASE_SWAP_H
#define PM_BASE_SWAP_H

#include <utility>    //std::swap, std::move & std::exchange
#include <algorithm>  //std::max

namespace Gudhi {
namespace persistence_matrix {

/**
 * @ingroup persistence_matrix
 *
 * @brief Empty structure.
 * Inherited instead of @ref Base_swap, when the column and row swaps are not enabled.
 */
struct Dummy_base_swap {
  friend void swap([[maybe_unused]] Dummy_base_swap& d1, [[maybe_unused]] Dummy_base_swap& d2) {}

  Dummy_base_swap([[maybe_unused]] unsigned int numberOfColumns = 0) {}
};

/**
 * @class Base_swap base_swap.h gudhi/Persistence_matrix/base_swap.h
 * @ingroup persistence_matrix
 *
 * @brief Class managing the column and row swaps in @ref Base_matrix and @ref Boundary_matrix.
 * 
 * @tparam Master_matrix An instantiation of @ref Matrix from which all types and options are deduced.
 * @tparam Base_matrix Either @ref Base_matrix or @ref Boundary_matrix.
 */
template <class Master_matrix, class Base_matrix>
class Base_swap {
 public:
  using Column_container = typename Master_matrix::Column_container;  /**< Column container type. */
  using Index = typename Master_matrix::Index;                        /**< Container index type. */
  using ID_index = typename Master_matrix::ID_index;                  /**< @ref IDIdx index type. */

  /**
   * @brief Default constructor.
   */
  Base_swap();
  /**
   * @brief As default constructor, but reserves spaces for @p numberOfColumns columns.
   * 
   * @param numberOfColumns Number of columns to reserve space for.
   */
  Base_swap(unsigned int numberOfColumns);
  /**
   * @brief Copy constructor.
   * 
   * @param matrixToCopy Matrix to copy.
   */
  Base_swap(const Base_swap& matrixToCopy) = default;
  /**
   * @brief Move constructor.
   * 
   * @param other Matrix to move.
   */
  Base_swap(Base_swap&& other) noexcept;

  /**
   * @brief Swaps the two columns at given indices in the column container. Does not updates the column index value,
   * potentially stored in the entries. This will be done when calling `_orderRows()`.
   * 
   * @param columnIndex1 First @ref MatIdx column index.
   * @param columnIndex2 Second @ref MatIdx column index.
   */
  void swap_columns(Index columnIndex1, Index columnIndex2);
  /**
   * @brief Swaps the two rows at the given indices, but in a lazy manner. That is, the swap is registered but not
   * executed. The reordering will be done when calling `_orderRows()`.
   * 
   * @param rowIndex1 First row index.
   * @param rowIndex2 Second row index.
   */
  void swap_rows(ID_index rowIndex1, ID_index rowIndex2);

  /**
   * @brief Assign operator.
   */
  Base_swap& operator=(Base_swap other);
  /**
   * @brief Swap operator.
   */
  friend void swap(Base_swap& base1, Base_swap& base2) {
    base1.indexToRow_.swap(base2.indexToRow_);
    base1.rowToIndex_.swap(base2.rowToIndex_);
    std::swap(base1.rowSwapped_, base2.rowSwapped_);
  }

 protected:
  using Index_dictionary = typename Master_matrix::template Dictionary<Index>;
  using Row_dictionary = typename Master_matrix::template Dictionary<ID_index>;

  Index_dictionary indexToRow_; /**< Map from row index to actual index in row container. */
  Row_dictionary rowToIndex_;   /**< Map from index in row container to "public" row index. */
  bool rowSwapped_;             /**< True if any rows were swapped since last call to `_orderRows()`. */

  void _orderRows();

  //access to inheriting matrix class
  constexpr Base_matrix* _matrix() { return static_cast<Base_matrix*>(this); }
  constexpr const Base_matrix* _matrix() const { return static_cast<const Base_matrix*>(this); }
};

template <class Master_matrix, class Base_matrix>
inline Base_swap<Master_matrix, Base_matrix>::Base_swap() : rowSwapped_(false)
{}

template <class Master_matrix, class Base_matrix>
inline Base_swap<Master_matrix, Base_matrix>::Base_swap(unsigned int numberOfColumns)
    : indexToRow_(numberOfColumns), rowToIndex_(numberOfColumns), rowSwapped_(false)
{
  for (Index i = 0; i < numberOfColumns; i++) {
    indexToRow_[i] = i;
    rowToIndex_[i] = i;
  }
}

template <class Master_matrix, class Base_matrix>
inline Base_swap<Master_matrix, Base_matrix>::Base_swap(Base_swap&& other) noexcept
    : indexToRow_(std::move(other.indexToRow_)),
      rowToIndex_(std::move(other.rowToIndex_)),
      rowSwapped_(std::exchange(other.rowSwapped_, 0))
{}

template <class Master_matrix, class Base_matrix>
inline void Base_swap<Master_matrix, Base_matrix>::swap_columns(Index columnIndex1, Index columnIndex2)
{
  swap(_matrix()->matrix_.at(columnIndex1), _matrix()->matrix_.at(columnIndex2));
  if constexpr (Master_matrix::Option_list::has_row_access) rowSwapped_ = true; //to update column index in entries.
}

template <class Master_matrix, class Base_matrix>
inline void Base_swap<Master_matrix, Base_matrix>::swap_rows(ID_index rowIndex1, ID_index rowIndex2)
{
  rowSwapped_ = true;

  if constexpr (Master_matrix::Option_list::has_map_column_container) {
    auto it1 = indexToRow_.find(rowIndex1);
    auto it2 = indexToRow_.find(rowIndex2);

    if (it1 == indexToRow_.end() && it2 == indexToRow_.end()) return;

    if (it1 == indexToRow_.end()) {
      indexToRow_.emplace(rowIndex1, it2->second);
      rowToIndex_.at(it2->second) = rowIndex1;
      indexToRow_.erase(it2->second);
      return;
    }

    if (it2 == indexToRow_.end()) {
      indexToRow_.emplace(rowIndex2, it1->second);
      rowToIndex_.at(it1->second) = rowIndex2;
      indexToRow_.erase(it1);
      return;
    }

    std::swap(rowToIndex_.at(it1->second), rowToIndex_.at(it2->second));
    std::swap(it1->second, it2->second);
  } else {
    for (auto i = indexToRow_.size(); i <= std::max(rowIndex1, rowIndex2); ++i) indexToRow_.push_back(i);

    std::swap(rowToIndex_[indexToRow_[rowIndex1]], rowToIndex_[indexToRow_[rowIndex2]]);
    std::swap(indexToRow_[rowIndex1], indexToRow_[rowIndex2]);
  }
}

template <class Master_matrix, class Base_matrix>
inline Base_swap<Master_matrix, Base_matrix>& Base_swap<Master_matrix, Base_matrix>::operator=(Base_swap other)
{
  indexToRow_.swap(other.indexToRow_);
  rowToIndex_.swap(other.rowToIndex_);
  std::swap(rowSwapped_, other.rowSwapped_);
  return *this;
}

template <class Master_matrix, class Base_matrix>
inline void Base_swap<Master_matrix, Base_matrix>::_orderRows()
{
  for (unsigned int i = 0; i < _matrix()->get_number_of_columns(); i++) {
    _matrix()->matrix_.at(i).reorder(rowToIndex_, i);
  }
  for (Index i = 0; i < _matrix()->get_number_of_columns(); i++) {
    indexToRow_[i] = i;
    rowToIndex_[i] = i;
  }
  rowSwapped_ = false;
}

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_BASE_SWAP_H

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
 * @file cell_types.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_matrix::Cell, @ref Gudhi::persistence_matrix::Cell_column_index and
 * @ref Gudhi::persistence_matrix::Cell_field_element classes, as well as the
 * @ref Gudhi::persistence_matrix::Dummy_cell_column_index_mixin and
 * @ref Gudhi::persistence_matrix::Dummy_cell_field_element_mixin structures.
 * Also defines the std::hash method for @ref Gudhi::persistence_matrix::Cell.
 */

#ifndef PM_MATRIX_CELL_H
#define PM_MATRIX_CELL_H

#include <utility>     //std::swap, std::exchange & std::move
#include <functional>  //std::hash

namespace Gudhi {
namespace persistence_matrix {

/**
 * @ingroup persistence_matrix
 *
 * @brief Empty structure.
 * Inherited instead of @ref Cell_column_index, when the row access is disabled.
 */
struct Dummy_cell_column_index_mixin 
{
  Dummy_cell_column_index_mixin() {}
  template <typename Index>
  Dummy_cell_column_index_mixin([[maybe_unused]] Index columnIndex) {}
};

/**
 * @ingroup persistence_matrix
 *
 * @brief Empty structure.
 * Inherited instead of @ref Cell_field_element, when @ref PersistenceMatrixOptions::is_z2 is true.
 */
struct Dummy_cell_field_element_mixin 
{
  Dummy_cell_field_element_mixin() {}
  template <class Field_element>
  Dummy_cell_field_element_mixin([[maybe_unused]] Field_element t) {}
};

/**
 * @ingroup persistence_matrix
 *
 * @brief Class managing the column index access of a cell.
 * 
 * @tparam Index @ref MatIdx index type.
 */
template <typename Index>
class Cell_column_index 
{
 public:
  /**
   * @brief Default constructor. Sets to the column index to -1.
   */
  Cell_column_index() : columnIndex_(-1){};
  /**
   * @brief Stores the given column index.
   * 
   * @param columnIndex Column index of the cell.
   */
  Cell_column_index(Index columnIndex) : columnIndex_(columnIndex){};
  /**
   * @brief Copy constructor.
   * 
   * @param cell Cell to copy.
   */
  Cell_column_index(const Cell_column_index& cell) : columnIndex_(cell.columnIndex_){};
  /**
   * @brief Move constructor.
   * 
   * @param cell Cell to move.
   */
  Cell_column_index(Cell_column_index&& cell) noexcept : columnIndex_(std::exchange(cell.columnIndex_, 0)){};

  /**
   * @brief Returns the @ref MatIdx column index stored in the cell.
   * 
   * @return Column index of the cell.
   */
  Index get_column_index() const { return columnIndex_; };
  /**
   * @brief Sets the column index to the given value.
   * 
   * @param columnIndex Column index of the cell.
   */
  void set_column_index(Index columnIndex) { columnIndex_ = columnIndex; }

  /**
   * @brief Assign operator.
   */
  Cell_column_index& operator=(Cell_column_index other) {
    std::swap(columnIndex_, other.columnIndex_);
    return *this;
  };

 private:
  Index columnIndex_;   /**< Column index. */
};

/**
 * @ingroup persistence_matrix
 *
 * @brief Class managing the value access of a cell.
 * 
 * @tparam Field_element Type of a cell value.
 */
template <class Field_element>
class Cell_field_element 
{
 public:
  /**
   * @brief Default constructor. Sets to the element to 0.
   */
  Cell_field_element() : element_(0){};
  /**
   * @brief Stores the given element.
   * 
   * @param element Value to store.
   */
  Cell_field_element(Field_element element) : element_(element){};
  /**
   * @brief Copy constructor.
   * 
   * @param cell Cell to copy.
   */
  Cell_field_element(const Cell_field_element& cell) : element_(cell.element_){};
  /**
   * @brief Move constructor.
   * 
   * @param cell Cell to move.
   */
  Cell_field_element(Cell_field_element&& cell) noexcept : element_(std::move(cell.element_)){};

  /**
   * @brief Returns the value stored in the cell.
   * 
   * @return Reference to the value of the cell.
   */
  Field_element& get_element() { return element_; };
  /**
   * @brief Returns the value stored in the cell.
   * 
   * @return Const reference to the value of the cell.
   */
  const Field_element& get_element() const { return element_; };
  /**
   * @brief Sets the value.
   * 
   * @param element Value to store.
   */
  void set_element(const Field_element& element) { element_ = element; }

  /**
   * @brief Assign operator.
   */
  Cell_field_element& operator=(Cell_field_element other) {
    std::swap(element_, other.element_);
    return *this;
  };

 private:
  Field_element element_;  /**< Value of the cell. */
};

/**
 * @class Cell cell_types.h gudhi/Persistence_matrix/columns/cell_types.h
 * @ingroup persistence_matrix
 *
 * @brief %Matrix cell class. Stores by default only the row index it belongs to, but can also store its
 * column index when the row access is enabled, as well as its value when they are different from only 0 and 1.
 * Zero-valued cells are never made explicit in the matrix.
 * 
 * @tparam Master_matrix An instantiation of @ref Matrix from which all types and options are deduced.
 */
template <class Master_matrix>
class Cell : public Master_matrix::Cell_column_index_option,
             public Master_matrix::Cell_field_element_option,
             public Master_matrix::Row_hook,
             public Master_matrix::Column_hook 
{
 private:
  using col_opt = typename Master_matrix::Cell_column_index_option;
  using field_opt = typename Master_matrix::Cell_field_element_option;

 public:
  using Master = Master_matrix;                           /**< Access to options from outside. */
  using Index = typename Master_matrix::Index;            /**< Column index type. */
  using ID_index = typename Master_matrix::ID_index;      /**< Row index type. */
  using Field_element = typename Master_matrix::Element;  /**< Value type. */

  /**
   * @brief Constructs an cell with all attributes at default values.
   */
  Cell(){};
  /**
   * @brief Constructs a cell with given row index. Other possible attributes are set at default values.
   * 
   * @param rowIndex @ref rowindex "Row index" of the cell.
   */
  Cell(ID_index rowIndex) : col_opt(), field_opt(), rowIndex_(rowIndex){};
  /**
   * @brief Constructs a cell with given row and column index. Other possible attributes are set at default values.
   * 
   * @param columnIndex Column index of the cell.
   * @param rowIndex @ref rowindex "Row index" of the cell.
   */
  Cell(Index columnIndex, ID_index rowIndex) : col_opt(columnIndex), field_opt(), rowIndex_(rowIndex){};
  /**
   * @brief Copy constructor.
   * 
   * @param cell Cell to copy.
   */
  Cell(const Cell& cell)
      : col_opt(static_cast<const col_opt&>(cell)),
        field_opt(static_cast<const field_opt&>(cell)),
        rowIndex_(cell.rowIndex_){};
  /**
   * @brief Move constructor.
   * 
   * @param cell Cell to move.
   */
  Cell(Cell&& cell) noexcept
      : col_opt(std::move(static_cast<col_opt&>(cell))),
        field_opt(std::move(static_cast<field_opt&>(cell))),
        rowIndex_(std::exchange(cell.rowIndex_, 0)){};

  /**
   * @brief Returns the row index stored in the cell.
   * 
   * @return @ref rowindex "Row index" of the cell.
   */
  ID_index get_row_index() const { return rowIndex_; };
  /**
   * @brief Sets the row index stored in the cell.
   * 
   * @param rowIndex @ref rowindex "Row index" of the cell.
   */
  void set_row_index(ID_index rowIndex) { rowIndex_ = rowIndex; };

  /**
   * @brief Assign operator.
   */
  Cell& operator=(Cell other) {
    col_opt::operator=(other);
    field_opt::operator=(other);
    std::swap(rowIndex_, other.rowIndex_);
    return *this;
  };

  /**
   * @brief Strictly smaller than comparator.
   * 
   * @param c1 First cell to compare.
   * @param c2 Second cell to compare.
   * @return true If the row index of the first cell is strictly smaller than the row index of the second cell.
   * @return false Otherwise.
   */
  friend bool operator<(const Cell& c1, const Cell& c2) { return c1.get_row_index() < c2.get_row_index(); }
  /**
   * @brief Equality comparator.
   * 
   * @param c1 First cell to compare.
   * @param c2 Second cell to compare.
   * @return true If the row index of the first cell is equal to the row index of the second cell.
   * @return false Otherwise.
   */
  friend bool operator==(const Cell& c1, const Cell& c2) { return c1.get_row_index() == c2.get_row_index(); }

  /**
   * @brief Converts the cell into a row index.
   * 
   * @return The row index of the cell.
   */
  operator ID_index() const { return rowIndex_; }
  /**
   * @brief Converts the cell into a pair of row index and cell value.
   * 
   * @return A std::pair with first element the row index and second element the value.
   */
  operator std::pair<ID_index, Field_element>() const {
    if constexpr (Master_matrix::Option_list::is_z2) {
      return {rowIndex_, 1};
    } else {
      return {rowIndex_, field_opt::element_};
    }
  }

 private:
  ID_index rowIndex_;   /**< Row index of the cell. */
};

}  // namespace persistence_matrix
}  // namespace Gudhi

/**
 * @ingroup persistence_matrix
 *
 * @brief Hash method for @ref Gudhi::persistence_matrix::Cell.
 *
 * The cells are differentiated by their row indices only. For example, two cells with the same row index
 * but different column indices have the same hash value.
 * 
 * @tparam Master_matrix Template parameter of @ref Gudhi::persistence_matrix::Cell.
 */
template <class Master_matrix>
struct std::hash<Gudhi::persistence_matrix::Cell<Master_matrix> > {
  std::size_t operator()(const Gudhi::persistence_matrix::Cell<Master_matrix>& cell) const {
    return std::hash<unsigned int>()(cell.get_row_index());
  }
};

#endif  // PM_MATRIX_CELL_H

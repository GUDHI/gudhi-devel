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
 * @file entry_types.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_matrix::Entry, @ref Gudhi::persistence_matrix::Entry_column_index and
 * @ref Gudhi::persistence_matrix::Entry_field_element classes, as well as the
 * @ref Gudhi::persistence_matrix::Dummy_entry_column_index_mixin and
 * @ref Gudhi::persistence_matrix::Dummy_entry_field_element_mixin structures.
 * Also defines the std::hash method for @ref Gudhi::persistence_matrix::Entry.
 */

#ifndef PM_MATRIX_ENTRY_H
#define PM_MATRIX_ENTRY_H

#include <utility>     //std::swap, std::exchange & std::move
#include <functional>  //std::hash

namespace Gudhi {
namespace persistence_matrix {

/**
 * @ingroup persistence_matrix
 *
 * @brief Empty structure.
 * Inherited instead of @ref Entry_column_index, when the row access is disabled.
 */
struct Dummy_entry_column_index_mixin 
{
  Dummy_entry_column_index_mixin() {}
  template <typename Index>
  Dummy_entry_column_index_mixin([[maybe_unused]] Index columnIndex) {}
};

/**
 * @ingroup persistence_matrix
 *
 * @brief Empty structure.
 * Inherited instead of @ref Entry_field_element, when @ref PersistenceMatrixOptions::is_z2 is true.
 */
struct Dummy_entry_field_element_mixin 
{
  Dummy_entry_field_element_mixin() {}
  template <class Field_element>
  Dummy_entry_field_element_mixin([[maybe_unused]] Field_element t) {}
};

/**
 * @ingroup persistence_matrix
 *
 * @brief Class managing the column index access of an entry.
 * 
 * @tparam Index @ref MatIdx index type.
 */
template <typename Index>
class Entry_column_index 
{
 public:
  /**
   * @brief Default constructor. Sets to the column index to -1.
   */
  Entry_column_index() : columnIndex_(-1){};
  /**
   * @brief Stores the given column index.
   * 
   * @param columnIndex Column index of the entry.
   */
  Entry_column_index(Index columnIndex) : columnIndex_(columnIndex){};
  /**
   * @brief Copy constructor.
   * 
   * @param entry Entry to copy.
   */
  Entry_column_index(const Entry_column_index& entry) : columnIndex_(entry.columnIndex_){};
  /**
   * @brief Move constructor.
   * 
   * @param entry Entry to move.
   */
  Entry_column_index(Entry_column_index&& entry) noexcept : columnIndex_(std::exchange(entry.columnIndex_, 0)){};

  /**
   * @brief Returns the @ref MatIdx column index stored in the entry.
   * 
   * @return Column index of the entry.
   */
  Index get_column_index() const { return columnIndex_; };
  /**
   * @brief Sets the column index to the given value.
   * 
   * @param columnIndex Column index of the entry.
   */
  void set_column_index(Index columnIndex) { columnIndex_ = columnIndex; }

  /**
   * @brief Assign operator.
   */
  Entry_column_index& operator=(Entry_column_index other) {
    std::swap(columnIndex_, other.columnIndex_);
    return *this;
  };

 private:
  Index columnIndex_;   /**< Column index. */
};

/**
 * @ingroup persistence_matrix
 *
 * @brief Class managing the value access of an entry.
 * 
 * @tparam Field_element Type of an entry value.
 */
template <class Field_element>
class Entry_field_element 
{
 public:
  /**
   * @brief Default constructor. Sets to the element to 0.
   */
  Entry_field_element() : element_(0){};
  /**
   * @brief Stores the given element.
   * 
   * @param element Value to store.
   */
  Entry_field_element(Field_element element) : element_(element){};
  /**
   * @brief Copy constructor.
   * 
   * @param entry Entry to copy.
   */
  Entry_field_element(const Entry_field_element& entry) : element_(entry.element_){};
  /**
   * @brief Move constructor.
   * 
   * @param entry Entry to move.
   */
  Entry_field_element(Entry_field_element&& entry) noexcept : element_(std::move(entry.element_)){};

  /**
   * @brief Returns the value stored in the entry.
   * 
   * @return Reference to the value of the entry.
   */
  Field_element& get_element() { return element_; };
  /**
   * @brief Returns the value stored in the entry.
   * 
   * @return Const reference to the value of the entry.
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
  Entry_field_element& operator=(Entry_field_element other) {
    std::swap(element_, other.element_);
    return *this;
  };

 private:
  Field_element element_;  /**< Value of the entry. */
};

/**
 * @class Entry entry_types.h gudhi/Persistence_matrix/columns/entry_types.h
 * @ingroup persistence_matrix
 *
 * @brief %Matrix entry class. Stores by default only the row index it belongs to, but can also store its
 * column index when the row access is enabled, as well as its value when they are different from only 0 and 1.
 * Zero-valued entries are never made explicit in the matrix.
 * 
 * @tparam Master_matrix An instantiation of @ref Matrix from which all types and options are deduced.
 */
template <class Master_matrix>
class Entry : public Master_matrix::Entry_column_index_option,
             public Master_matrix::Entry_field_element_option,
             public Master_matrix::Row_hook,
             public Master_matrix::Column_hook 
{
 private:
  using col_opt = typename Master_matrix::Entry_column_index_option;
  using field_opt = typename Master_matrix::Entry_field_element_option;

 public:
  using Master = Master_matrix;                           /**< Access to options from outside. */
  using Index = typename Master_matrix::Index;            /**< Column index type. */
  using ID_index = typename Master_matrix::ID_index;      /**< Row index type. */
  using Field_element = typename Master_matrix::Element;  /**< Value type. */

  /**
   * @brief Constructs an entry with all attributes at default values.
   */
  Entry(){};
  /**
   * @brief Constructs an entry with given row index. Other possible attributes are set at default values.
   * 
   * @param rowIndex @ref rowindex "Row index" of the entry.
   */
  Entry(ID_index rowIndex) : col_opt(), field_opt(), rowIndex_(rowIndex){};
  /**
   * @brief Constructs an entry with given row and column index. Other possible attributes are set at default values.
   * 
   * @param columnIndex Column index of the entry.
   * @param rowIndex @ref rowindex "Row index" of the entry.
   */
  Entry(Index columnIndex, ID_index rowIndex) : col_opt(columnIndex), field_opt(), rowIndex_(rowIndex){};
  /**
   * @brief Copy constructor.
   * 
   * @param entry Entry to copy.
   */
  Entry(const Entry& entry)
      : col_opt(static_cast<const col_opt&>(entry)),
        field_opt(static_cast<const field_opt&>(entry)),
        rowIndex_(entry.rowIndex_){};
  /**
   * @brief Move constructor.
   * 
   * @param entry Entry to move.
   */
  Entry(Entry&& entry) noexcept
      : col_opt(std::move(static_cast<col_opt&>(entry))),
        field_opt(std::move(static_cast<field_opt&>(entry))),
        rowIndex_(std::exchange(entry.rowIndex_, 0)){};

  /**
   * @brief Returns the row index stored in the entry.
   * 
   * @return @ref rowindex "Row index" of the entry.
   */
  ID_index get_row_index() const { return rowIndex_; };
  /**
   * @brief Sets the row index stored in the entry.
   * 
   * @param rowIndex @ref rowindex "Row index" of the entry.
   */
  void set_row_index(ID_index rowIndex) { rowIndex_ = rowIndex; };

  /**
   * @brief Assign operator.
   */
  Entry& operator=(Entry other) {
    col_opt::operator=(other);
    field_opt::operator=(other);
    std::swap(rowIndex_, other.rowIndex_);
    return *this;
  };

  /**
   * @brief Strictly smaller than comparator.
   * 
   * @param c1 First entry to compare.
   * @param c2 Second entry to compare.
   * @return true If the row index of the first entry is strictly smaller than the row index of the second entry.
   * @return false Otherwise.
   */
  friend bool operator<(const Entry& c1, const Entry& c2) { return c1.get_row_index() < c2.get_row_index(); }
  /**
   * @brief Equality comparator.
   * 
   * @param c1 First entry to compare.
   * @param c2 Second entry to compare.
   * @return true If the row index of the first entry is equal to the row index of the second entry.
   * @return false Otherwise.
   */
  friend bool operator==(const Entry& c1, const Entry& c2) { return c1.get_row_index() == c2.get_row_index(); }

  /**
   * @brief Converts the entry into a row index.
   * 
   * @return The row index of the entry.
   */
  operator ID_index() const { return rowIndex_; }
  /**
   * @brief Converts the entry into a pair of row index and entry value.
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
  ID_index rowIndex_;   /**< Row index of the entry. */
};

}  // namespace persistence_matrix
}  // namespace Gudhi

/**
 * @ingroup persistence_matrix
 *
 * @brief Hash method for @ref Gudhi::persistence_matrix::Entry.
 *
 * The entries are differentiated by their row indices only. For example, two entries with the same row index
 * but different column indices have the same hash value.
 * 
 * @tparam Master_matrix Template parameter of @ref Gudhi::persistence_matrix::Entry.
 */
template <class Master_matrix>
struct std::hash<Gudhi::persistence_matrix::Entry<Master_matrix> > {
  std::size_t operator()(const Gudhi::persistence_matrix::Entry<Master_matrix>& entry) const {
    return std::hash<unsigned int>()(entry.get_row_index());
  }
};

#endif  // PM_MATRIX_ENTRY_H

/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file PersistenceMatrixColumn.h
 * @author Hannah Schreiber
 * @brief Contains the concept for the matrix columns.
 */

namespace Gudhi {
namespace persistence_matrix {

/** @brief Concept of the column classes used by the @ref Matrix class.
 *
 * Implementations of this concept are @ref Heap_column, @ref List_column, @ref Vector_column, @ref Naive_vector_column
 * @ref Set_column, @ref Unordered_set_column, @ref PersistenceMatrixColumn and @ref Intrusive_set_column.
 */
class PersistenceMatrixColumn :
    /**
     * @brief If PersistenceMatrixOptions::has_row_access is true, then @ref Row_access.
     * Otherwise @ref Dummy_row_access. Can eventually be removed if the structure of the column does not allow
     * row access (as for @ref Heap_column), but then it needs to be notified in the documentation of @ref Column_types
     * and as static_assert in @ref Matrix::_assert_options. 
     */
    public Row_access_option,
    /**
     * @brief If PersistenceMatrixOptions::has_column_pairings ||
                 PersistenceMatrixOptions::has_vine_update ||
                 PersistenceMatrixOptions::can_retrieve_representative_cycles is true, then @ref
     Column_dimension_holder. Otherwise @ref Dummy_dimension_holder.
     *
     */
    public Column_dimension_option,
    /**
     * @brief If (PersistenceMatrixOptions::has_column_pairings ||
                  PersistenceMatrixOptions::has_vine_update ||
                  PersistenceMatrixOptions::can_retrieve_representative_cycles) &&
                  !PersistenceMatrixOptions::is_of_boundary_type is true, then @ref Chain_column_extra_properties.
       Otherwise @ref Dummy_chain_properties.
     */
    public Chain_column_option 
{
 public:
  using Master = unspecified;                 /**< Master matrix. */
  using Cell_constructor = unspecified;       /**< @ref Cell factory. */
  using Field_operators = unspecified;        /**< Follows the @ref FieldOperators concept. */
  using Field_element_type = unspecified;     /**< Type of a field element. */
  using index = unspecified;                  /**< Type of MatIdx index. */
  using id_index = unspecified;               /**< Type of IDIdx index. */
  using dimension_type = unspecified;         /**< Type for dimension value. */
  using Cell = unspecified;                   /**< @ref Cell. */
  using Column_type = unspecified;            /**< Type of cell container. */
  using iterator = unspecified;               /**< Column iterator type. */
  using const_iterator = unspecified;         /**< Column const_iterator type. */
  using reverse_iterator = unspecified;       /**< Column reverse_iterator type. */
  using const_reverse_iterator = unspecified; /**< Column const_reverse_iterator type. */

  /**
   * @brief Constructs an empty column. If @p cellConstructor is not specified or is set to nullptr, the column
   * can only be used as a dummy, i.e., no modifying method should be used or there will be a segmentation fault.
   * Same goes for @p operators if @ref PersistenceMatrixOptions::is_z2 is false.
   * 
   * @param operators Pointer to the field operators.
   * @param cellConstructor Pointer to the cell factory.
   */
  PersistenceMatrixColumn(Field_operators* operators = nullptr, Cell_constructor* cellConstructor = nullptr);
  /**
   * @brief Constructs a column from the given range of @ref Matrix::cell_rep_type. If the dimension is stored,
   * the face is assumed to be simplicial and its dimension to be @ref nonZeroRowIndices length - 1 or 0.
   * Otherwise, the dimension should be specified with another constructor.
   * 
   * @tparam Container_type Range of @ref Matrix::cell_rep_type. Assumed to have a begin(), end() and size() method.
   * @param nonZeroRowIndices Range of @ref Matrix::cell_rep_type representing all rows with non zero values.
   * @param operators Pointer to the field operators.
   * @param cellConstructor Pointer to the cell factory.
   */
  template <class Container_type = typename Master_matrix::boundary_type>
  PersistenceMatrixColumn(const Container_type& nonZeroRowIndices, 
                          Field_operators* operators,
                          Cell_constructor* cellConstructor);
  /**
   * @brief Constructs a column from the given range of @ref Matrix::cell_rep_type such that the rows can be accessed.
   * Each new cell in the column is also inserted in a row using @ref Row_access::insert_cell.
   * If the dimension is stored, the face is assumed to be simplicial and its dimension to be
   * @ref nonZeroRowIndices length - 1 or 0. Otherwise, the dimension should be specified with another constructor.
   * 
   * @tparam Container_type Range of @ref Matrix::cell_rep_type. Assumed to have a begin(), end() and size() method.
   * @tparam Row_container_type Either std::map if @ref PersistenceMatrixOptions::has_removable_rows is true or
   * std::vector<Row_type>.
   * @param columnIndex MatIdx column index that should be specified to the cells.
   * @param nonZeroRowIndices Range of @ref Matrix::cell_rep_type representing all rows with non zero values.
   * @param rowContainer Pointer to the row container that will be forwarded to @ref Row_access at construction.
   * @param operators Pointer to the field operators.
   * @param cellConstructor Pointer to the cell factory.
   */
  template <class Container_type = typename Master_matrix::boundary_type, class Row_container_type>
  PersistenceMatrixColumn(index columnIndex, 
                          const Container_type& nonZeroRowIndices, 
                          Row_container_type* rowContainer,
                          Field_operators* operators,
                          Cell_constructor* cellConstructor);
  /**
   * @brief Constructs a column from the given range of @ref Matrix::cell_rep_type and stores the given dimension
   * if @ref Column_dimension_option is not a dummy.
   * 
   * @tparam Container_type Range of @ref Matrix::cell_rep_type. Assumed to have a begin(), end() and size() method.
   * @param nonZeroChainRowIndices Range of @ref Matrix::cell_rep_type representing all rows with non zero values.
   * @param dimension Dimension of the column. Is ignored if the dimension is not stored.
   * @param operators Pointer to the field operators.
   * @param cellConstructor Pointer to the cell factory.
   */
  template <class Container_type = typename Master_matrix::boundary_type>
  PersistenceMatrixColumn(const Container_type& nonZeroChainRowIndices, 
                          dimension_type dimension,
                          Field_operators* operators,
                          Cell_constructor* cellConstructor);
  /**
   * @brief Constructs a column from the given range of @ref Matrix::cell_rep_type such that the rows can be accessed.
   * Each new cell in the column is also inserted in a row using @ref Row_access::insert_cell.
   * Stores the given dimension if @ref Column_dimension_option is not a dummy.
   * 
   * @tparam Container_type Range of @ref Matrix::cell_rep_type. Assumed to have a begin(), end() and size() method.
   * @tparam Row_container_type Either std::map if @ref PersistenceMatrixOptions::has_removable_rows is true or
   * std::vector<Row_type>.
   * @param columnIndex MatIdx column index that should be specified to the cells.
   * @param nonZeroRowIndices Range of @ref Matrix::cell_rep_type representing all rows with non zero values.
   * @param dimension Dimension of the column. Is ignored if the dimension is not stored.
   * @param rowContainer Pointer to the row container that will be forwarded to @ref Row_access at construction.
   * @param operators Pointer to the field operators.
   * @param cellConstructor Pointer to the cell factory.
   */
  template <class Container_type = typename Master_matrix::boundary_type, class Row_container_type>
  PersistenceMatrixColumn(index columnIndex, 
                          const Container_type& nonZeroChainRowIndices, 
                          dimension_type dimension,
                          Row_container_type* rowContainer, 
                          Field_operators* operators,
                          Cell_constructor* cellConstructor);
  /**
   * @brief Copy constructor. If @p operators or @p cellConstructor is not a null pointer, its value is kept
   * instead of the one in the copied column.
   * 
   * @param column Column to copy.
   * @param operators Pointer to the field operators.
   * If null pointer, the pointer in @p column is choosen instead.
   * @param cellConstructor Pointer to the cell factory.
   * If null pointer, the pointer in @p column is choosen instead.
   */
  PersistenceMatrixColumn(const PersistenceMatrixColumn& column, 
                          Field_operators* operators = nullptr,
                          Cell_constructor* cellConstructor = nullptr);
  /**
   * @brief Copy constructor with row access.
   * If @p operators or @p cellConstructor is not a null pointer, its value is kept
   * instead of the one in the copied column.
   * 
   * @tparam Row_container_type  Either std::map if @ref PersistenceMatrixOptions::has_removable_rows is true or
   * std::vector<Row_type>.
   * @param column Column to copy.
   * @param columnIndex MatIdx column index of the new column once copied.
   * @param rowContainer Pointer to the row container that will be forwarded to @ref Row_access.
   * @param operators  Pointer to the field operators.
   * If null pointer, the pointer in @p column is choosen instead.
   * @param cellConstructor Pointer to the cell factory.
   * If null pointer, the pointer in @p column is choosen instead.
   */
  template <class Row_container_type>
  PersistenceMatrixColumn(const PersistenceMatrixColumn& column, 
                          index columnIndex, 
                          Row_container_type* rowContainer,
                          Field_operators* operators = nullptr, 
                          Cell_constructor* cellConstructor = nullptr);
  /**
   * @brief Move constructor.
   * 
   * @param column Column to move.
   */
  PersistenceMatrixColumn(PersistenceMatrixColumn&& column) noexcept;
  /**
   * @brief Destructor.
   */
  ~PersistenceMatrixColumn();

  /**
   * @brief Returns the values of the column, zero values included.
   * 
   * @param columnLength Number of rows to be returned. If -1, the number of rows is fixed at the biggest 
   * row index with non zero value. Default value: -1.
   * @return Vector of @ref Field_element_type. At element \f$ i \f$ of the vector will be stored the value
   * at row \f$ i \f$ of the column.
   */
  std::vector<Field_element_type> get_content(int columnLength = -1) const;
  /**
   * @brief Indicates if the cell at given row index has value zero.
   * 
   * @param rowIndex Row index to look at.
   * @return true If the cell has value zero.
   * @return false Otherwise.
   */
  bool is_non_zero(id_index rowIndex) const;
  /**
   * @brief Indicates if the column is empty or has only zero values.
   * 
   * @return true If the column is empty or only has zero values.
   * @return false Otherwise.
   */
  bool is_empty();
  /**
   * @brief Returns the size of the underlying container.
   *
   * @warning Depending of the column type, the container does not have to contain only the non-zero cells.
   * Even if for most of the types, the size of the container will correspond to the number of non-zero values
   * in the column, it is not always the case. See description of the actual Column class for more details.
   * 
   * @return Size of the underlying container.
   */
  std::size_t size() const;

  /**
   * @brief Reorders the column with the given map of row indices. Also changes the column index stored in the
   * cells if row access is enabled and @p columnIndex is not -1.
   *
   * Only useful for base and boundary matrices using lazy swaps.
   * 
   * @tparam Map_type Map with an `at` method.
   * @param valueMap Map such that `valueMap.at(i)` indicates the new row index of the cell
   * at current row index `i`.
   * @param columnIndex New MatIdx column index of the column. If -1, the index does not change. Ignored if
   * the row access is not enabled. Default value: -1.
   */
  template <class Map_type>
  void reorder(const Map_type& valueMap, [[maybe_unused]] index columnIndex = -1);
  /**
   * @brief Zeros/empties the column.
   * 
   * Only useful for base and boundary matrices. Used in `zero_column` and in the reduction algorithm
   * for the persistence barcode.
   */
  void clear();
  /**
   * @brief Zeros the cell at given row index.
   * 
   * Only useful for base and boundary matrices. Used in `zero_cell` and during vine swaps.
   * 
   * @param rowIndex Row index of the cell to zero.
   */
  void clear(id_index rowIndex);

  /**
   * @brief Returns the row index of the pivot. If the column does not have a pivot, returns -1.
   *
   * Only useful for boundary and chain matrices.
   * 
   * @return Row index of the pivot or -1.
   */
  id_index get_pivot();
  /**
   * @brief Returns the value of the pivot. If the column does not have a pivot, returns 0.
   *
   * Has to have value 1 if \f$ Z_2 \f$ coefficients are used.
   *
   * Only useful for boundary and chain matrices.
   * 
   * @return The value of the pivot or 0.
   */
  Field_element_type get_pivot_value();

  /**
   * @brief Returns a begin @ref Cell iterator to iterate over all cells contained in the underlying container.
   *
   * @warning The iterators really just iterate over the underlying container. Depending of the column type,
   * neither the content nor the order is garanteed. See description of the actual Column class for more details.
   * 
   * @return @ref Cell iterator.
   */
  iterator begin() noexcept;
  /**
   * @brief Returns a begin @ref Cell const iterator to iterate over all cells contained in the underlying container.
   *
   * @warning The iterators really just iterate over the underlying container. Depending of the column type,
   * neither the content nor the order is garanteed. See description of the actual Column class for more details.
   * 
   * @return @ref Cell const iterator.
   */
  const_iterator begin() const noexcept;
  /**
   * @brief Returns a end @ref Cell iterator, iterating over all cells contained in the underlying container.
   *
   * @warning The iterators really just iterate over the underlying container. Depending of the column type,
   * neither the content nor the order is garanteed. See description of the actual Column class for more details.
   * 
   * @return @ref Cell iterator.
   */
  iterator end() noexcept;
  /**
   * @brief Returns a end @ref Cell const iterator, iterating over all cells contained in the underlying container.
   *
   * @warning The iterators really just iterate over the underlying container. Depending of the column type,
   * neither the content nor the order is garanteed. See description of the actual Column class for more details.
   * 
   * @return @ref Cell const iterator.
   */
  const_iterator end() const noexcept;
  /**
   * @brief Returns a begin @ref Cell reverse iterator to iterate over all cells contained in the underlying container.
   *
   * @warning The iterators really just iterate over the underlying container. Depending of the column type,
   * neither the content nor the order is garanteed. See description of the actual Column class for more details.
   * 
   * @return @ref Cell reverse iterator.
   */
  reverse_iterator rbegin() noexcept;
  /**
   * @brief Returns a begin @ref Cell const reverse iterator to iterate over all cells contained in the underlying
   * container.
   *
   * @warning The iterators really just iterate over the underlying container. Depending of the column type,
   * neither the content nor the order is garanteed. See description of the actual Column class for more details.
   * 
   * @return @ref Cell const reverse iterator.
   */
  const_reverse_iterator rbegin() const noexcept;
  /**
   * @brief Returns a end @ref Cell reverse iterator, iterating over all cells contained in the underlying container.
   *
   * @warning The iterators really just iterate over the underlying container. Depending of the column type,
   * neither the content nor the order is garanteed. See description of the actual Column class for more details.
   * 
   * @return @ref Cell reverse iterator.
   */
  reverse_iterator rend() noexcept;
  /**
   * @brief Returns a end @ref Cell const reverse iterator, iterating over all cells contained in the underlying
   * container.
   *
   * @warning The iterators really just iterate over the underlying container. Depending of the column type,
   * neither the content nor the order is garanteed. See description of the actual Column class for more details.
   * 
   * @return @ref Cell const reverse iterator.
   */
  const_reverse_iterator rend() const noexcept;

  /**
   * @brief Adds the given @ref Cell range onto the column.
   * 
   * @tparam Cell_range @ref Cell range with `begin` and `end` method.
   * @param column @ref Cell range. Every cell has to return the right value when using `get_row_index` and, 
   * if @ref PersistenceMatrixOptions::is_z2 is false, also when using `get_element`.
   * @return Reference to this column.
   */
  template <class Cell_range>
  PersistenceMatrixColumn& operator+=(const Cell_range& column);
  /**
   * @brief Adds the given column onto this column.
   * 
   * @param column Column to add.
   * @return Reference to this column.
   */
  PersistenceMatrixColumn& operator+=(PersistenceMatrixColumn& column);

  /**
   * @brief Multiplies all values in the column with the given value.
   * 
   * @param val Value to multiply.
   * @return Reference to this column.
   */
  PersistenceMatrixColumn& operator*=(const Field_element_type& val);

  /**
   * @brief @p this = @p val * @p this + @p column
   * 
   * @tparam Cell_range @ref Cell range with `begin` and `end` method.
   * @param val Value to multiply.
   * @param column @ref Cell range. Every cell has to return the right value when using `get_row_index` and, 
   * if @ref PersistenceMatrixOptions::is_z2 is false, also when using `get_element`.
   * @return Reference to this column.
   */
  template <class Cell_range>
  PersistenceMatrixColumn& multiply_and_add(const Field_element_type& val, const Cell_range& column);
  /**
   * @brief @p this = @p val * @p this + @p column
   * 
   * @param val Value to multiply.
   * @param column Column to add.
   * @return Reference to this column.
   */
  PersistenceMatrixColumn& multiply_and_add(const Field_element_type& val, PersistenceMatrixColumn& column);
  /**
   * @brief @p this = @p this + @p column * @p val
   * 
   * @tparam Cell_range @ref Cell range with `begin` and `end` method.
   * @param column @ref Cell range. Every cell has to return the right value when using `get_row_index` and, 
   * if @ref PersistenceMatrixOptions::is_z2 is false, also when using `get_element`.
   * @param val Value to multiply.
   * @return Reference to this column.
   */
  template <class Cell_range>
  PersistenceMatrixColumn& multiply_and_add(const Cell_range& column, const Field_element_type& val);
  /**
   * @brief @p this = @p this + @p column * @p val
   * 
   * @param column Column to add.
   * @param val Value to multiply.
   * @return Reference to this column.
   */
  PersistenceMatrixColumn& multiply_and_add(PersistenceMatrixColumn& column, const Field_element_type& val);

  /**
   * @brief Equality comparator. Equal in the sense that what is "supposed" to be contained in the columns is equal,
   * not what is actually stored in the underlying container. For exemple, the underlying container of 
   * @ref Vector_column can contain cells which were erased explicitely by @ref clear(index). Those cells should not
   * be taken into account while comparing.
   * 
   * @param c1 First column to compare.
   * @param c2 Second column to compare.
   * @return true If both column are equal.
   * @return false Otherwise.
   */
  friend bool operator==(const PersistenceMatrixColumn& c1, const PersistenceMatrixColumn& c2);
  /**
   * @brief "Strictly smaller than" comparator. Usually a lexicographical order, but what matters is that the
   * order is total. The order should apply on what is "supposed" to be contained in the columns,
   * not what is actually stored in the underlying container. For exemple, the underlying container of 
   * @ref Vector_column can contain cells which were erased explicitely by @ref clear(index). Those cells should not
   * be taken into account while comparing.
   * 
   * @param c1 First column to compare.
   * @param c2 Second column to compare.
   * @return true If the first column is strictly smaller than the second one.
   * @return false Otherwise.
   */
  friend bool operator<(const PersistenceMatrixColumn& c1, const PersistenceMatrixColumn& c2);

  /**
   * @brief Assign operator. Should be disabled when row access is enabled.
   */
  PersistenceMatrixColumn& operator=(const PersistenceMatrixColumn& other);
  /**
   * @brief Swap operator.
   */
  friend void swap(PersistenceMatrixColumn& col1, PersistenceMatrixColumn& col2);
};

}  // namespace persistence_matrix
}  // namespace Gudhi

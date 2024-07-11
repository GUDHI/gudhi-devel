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

/**
 * @ingroup persistence_matrix
 *
 * @brief If PersistenceMatrixOptions::has_row_access is true, then @ref Row_access. Otherwise @ref Dummy_row_access.
 * Can eventually be removed if the structure of the column does not allow row access (as for @ref Heap_column), but
 * then it needs to be notified in the documentation of @ref Column_types and as static_assert in
 * Matrix::_assert_options.
 */
using Row_access_option = Row_access;
/**
 * @ingroup persistence_matrix
 *
 * @brief If @ref PersistenceMatrixOptions::has_column_pairings or @ref PersistenceMatrixOptions::has_vine_update or
 * @ref PersistenceMatrixOptions::can_retrieve_representative_cycles is true, then @ref Column_dimension_holder.
 * Otherwise @ref Dummy_dimension_holder.
 */
using Column_dimension_option = Column_dimension_holder;
/**
 * @ingroup persistence_matrix
 *
 * @brief If @ref PersistenceMatrixOptions::is_of_boundary_type is false, and,
 * @ref PersistenceMatrixOptions::has_column_pairings or @ref PersistenceMatrixOptions::has_vine_update or
 * @ref PersistenceMatrixOptions::can_retrieve_representative_cycles is true, then @ref Chain_column_extra_properties.
 * Otherwise @ref Dummy_chain_properties.
 */
using Chain_column_option = Chain_column_extra_properties;

/** 
 * @ingroup persistence_matrix
 *
 * @brief Concept of the column classes used by the @ref Matrix class. The classes the columns inheritates from
 * are either real or dummy classes, see @ref Row_access_option, @ref Column_dimension_option, @ref Chain_column_option.
 * If used with column compression, the column type has to have its `std::hash` method.
 *
 * Implementations of this concept are @ref Heap_column, @ref List_column, @ref Vector_column, @ref Naive_vector_column
 * @ref Set_column, @ref Unordered_set_column, @ref Intrusive_list_column and @ref Intrusive_set_column.
 */
class PersistenceMatrixColumn :
    public Row_access_option,
    public Column_dimension_option,
    public Chain_column_option
{
 public:
  using Master = unspecified;                 /**< Master matrix, that is a templated @ref Matrix. */
  using index = unspecified;                  /**< Type of @ref MatIdx index. */
  using id_index = unspecified;               /**< Type of @ref IDIdx index. */
  using dimension_type = unspecified;         /**< Type for dimension value. */
  using Field_element_type = unspecified;     /**< Type of a field element. */
  using Cell = unspecified;                   /**< @ref Cell. */
  using Column_settings = unspecified;        /**< Structure giving access to external classes eventually necessary,
                                                   like a cell pool for example. */
  using iterator = unspecified;               /**< Column iterator type. */
  using const_iterator = unspecified;         /**< Column const_iterator type. */
  using reverse_iterator = unspecified;       /**< Column reverse_iterator type. */
  using const_reverse_iterator = unspecified; /**< Column const_reverse_iterator type. */

  /**
   * @brief Constructs an empty column. If @p cellConstructor is not specified or is set to `nullptr`, the column
   * can only be used as a dummy, i.e., no modifying method should be used or there will be a segmentation fault.
   * Same goes for @p operators if @ref PersistenceMatrixOptions::is_z2 is false.
   * 
   * @param colSettings Pointer to a setting structure or `nullptr`. The structure should contain all the necessary
   * classes specific to the column type, such as custom allocators. The specificities are this way hidden behind
   * a commun interface for all column types. If @p colSettings is not specified or is equal to `nullptr`, the column
   * should still be constructable eventhough not necessarily "usable".
   */
  PersistenceMatrixColumn(Column_settings* colSettings = nullptr);
  /**
   * @brief Constructs a column from the given range of @ref Matrix::cell_rep_type. If the dimension is stored,
   * the face is assumed to be simplicial and its dimension to be `nonZeroRowIndices length - 1` or `0`.
   * Otherwise, the dimension should be specified with another constructor.
   * 
   * @tparam Container_type Range of @ref Matrix::cell_rep_type. Assumed to have a %begin(), %end() and %size() method.
   * @param nonZeroRowIndices Range of @ref Matrix::cell_rep_type representing all rows with non zero values.
   * @param colSettings Pointer to an existing setting structure. The structure should contain all the necessary
   * classes specific to the column type, such as custom allocators. The specificities are this way hidden behind
   * a commun interface for all column types.
   */
  template <class Container_type = typename Master_matrix::boundary_type>
  PersistenceMatrixColumn(const Container_type& nonZeroRowIndices, 
                          Column_settings* colSettings);
  /**
   * @brief Constructs a column from the given range of @ref Matrix::cell_rep_type such that the rows can be accessed.
   * Each new cell in the column is also inserted in a row using @ref Row_access::insert_cell.
   * If the dimension is stored, the face is assumed to be simplicial and its dimension to be
   * `nonZeroRowIndices length - 1` or `0`. Otherwise, the dimension should be specified with another constructor.
   * 
   * @tparam Container_type Range of @ref Matrix::cell_rep_type. Assumed to have a %begin(), %end() and %size() method.
   * @tparam Row_container_type Either std::map if @ref PersistenceMatrixOptions::has_removable_rows is true or
   * std::vector<Row_type>.
   * @param columnIndex @ref MatIdx column index that should be specified to the cells.
   * @param nonZeroRowIndices Range of @ref Matrix::cell_rep_type representing all rows with non zero values.
   * @param rowContainer Pointer to the row container that will be forwarded to @ref Row_access at construction.
   * @param colSettings Pointer to an existing setting structure. The structure should contain all the necessary
   * classes specific to the column type, such as custom allocators. The specificities are this way hidden behind
   * a commun interface for all column types.
   */
  template <class Container_type = typename Master_matrix::boundary_type, class Row_container_type>
  PersistenceMatrixColumn(index columnIndex, 
                          const Container_type& nonZeroRowIndices, 
                          Row_container_type* rowContainer,
                          Column_settings* colSettings);
  /**
   * @brief Constructs a column from the given range of @ref Matrix::cell_rep_type and stores the given dimension
   * if @ref Column_dimension_option is not a dummy.
   * 
   * @tparam Container_type Range of @ref Matrix::cell_rep_type. Assumed to have a %begin(), %end() and %size() method.
   * @param nonZeroChainRowIndices Range of @ref Matrix::cell_rep_type representing all rows with non zero values.
   * @param dimension Dimension of the column. Is ignored if the dimension is not stored.
   * @param colSettings Pointer to an existing setting structure. The structure should contain all the necessary
   * classes specific to the column type, such as custom allocators. The specificities are this way hidden behind
   * a commun interface for all column types.
   */
  template <class Container_type = typename Master_matrix::boundary_type>
  PersistenceMatrixColumn(const Container_type& nonZeroChainRowIndices, 
                          dimension_type dimension,
                          Column_settings* colSettings);
  /**
   * @brief Constructs a column from the given range of @ref Matrix::cell_rep_type such that the rows can be accessed.
   * Each new cell in the column is also inserted in a row using @ref Row_access::insert_cell.
   * Stores the given dimension if @ref Column_dimension_option is not a dummy.
   * 
   * @tparam Container_type Range of @ref Matrix::cell_rep_type. Assumed to have a %begin(), %end() and %size() method.
   * @tparam Row_container_type Either std::map if @ref PersistenceMatrixOptions::has_removable_rows is true or
   * std::vector<Row_type>.
   * @param columnIndex @ref MatIdx column index that should be specified to the cells.
   * @param nonZeroChainRowIndices Range of @ref Matrix::cell_rep_type representing all rows with non zero values.
   * @param dimension Dimension of the column. Is ignored if the dimension is not stored.
   * @param rowContainer Pointer to the row container that will be forwarded to @ref Row_access at construction.
   * @param colSettings Pointer to an existing setting structure. The structure should contain all the necessary
   * classes specific to the column type, such as custom allocators. The specificities are this way hidden behind
   * a commun interface for all column types.
   */
  template <class Container_type = typename Master_matrix::boundary_type, class Row_container_type>
  PersistenceMatrixColumn(index columnIndex, 
                          const Container_type& nonZeroChainRowIndices, 
                          dimension_type dimension,
                          Row_container_type* rowContainer, 
                          Column_settings* colSettings);
  /**
   * @brief Copy constructor. If @p operators or @p cellConstructor is not a null pointer, its value is kept
   * instead of the one in the copied column.
   * 
   * @param column Column to copy.
   * @param colSettings Pointer to a setting structure or `nullptr`. The structure should contain all the necessary
   * classes specific to the column type, such as custom allocators. The specificities are this way hidden behind
   * a commun interface for all column types. If @p colSettings is not specified or is equal to `nullptr`, the structure
   * stored in @p column is used instead.
   */
  PersistenceMatrixColumn(const PersistenceMatrixColumn& column, 
                          Column_settings* colSettings = nullptr);
  /**
   * @brief Copy constructor with row access.
   * If @p operators or @p cellConstructor is not a null pointer, its value is kept
   * instead of the one in the copied column.
   * 
   * @tparam Row_container_type  Either std::map if @ref PersistenceMatrixOptions::has_removable_rows is true or
   * std::vector<Row_type>.
   * @param column Column to copy.
   * @param columnIndex @ref MatIdx column index of the new column once copied.
   * @param rowContainer Pointer to the row container that will be forwarded to @ref Row_access.
   * @param colSettings Pointer to a setting structure or `nullptr`. The structure should contain all the necessary
   * classes specific to the column type, such as custom allocators. The specificities are this way hidden behind
   * a commun interface for all column types. If @p colSettings is not specified or is equal to `nullptr`, the structure
   * stored in @p column is used instead.
   */
  template <class Row_container_type>
  PersistenceMatrixColumn(const PersistenceMatrixColumn& column, 
                          index columnIndex, 
                          Row_container_type* rowContainer,
                          Column_settings* colSettings = nullptr);
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
   * Only useful for @ref basematrix "base" and @ref boundarymatrix "boundary matrices" using lazy swaps.
   * 
   * @tparam Map_type Map with an %at() method.
   * @param valueMap Map such that `valueMap.at(i)` indicates the new row index of the cell
   * at current row index `i`.
   * @param columnIndex New @ref MatIdx column index of the column. If -1, the index does not change. Ignored if
   * the row access is not enabled. Default value: -1.
   */
  template <class Map_type>
  void reorder(const Map_type& valueMap, [[maybe_unused]] index columnIndex = -1);
  /**
   * @brief Zeros/empties the column.
   * 
   * Only useful for @ref basematrix "base" and @ref boundarymatrix "boundary matrices".
   * Used in @ref Matrix::zero_column and in the reduction algorithm for the persistence barcode.
   */
  void clear();
  /**
   * @brief Zeros the cell at given row index.
   * 
   * Only useful for @ref basematrix "base" and @ref boundarymatrix "boundary matrices".
   * Used in @ref Matrix::zero_cell and during vine swaps.
   *
   * @warning For @ref Vector_column, do not clear a cell that was already at zero or the results of @ref size and
   * @ref is_empty will be wrong.
   * 
   * @param rowIndex Row index of the cell to zero.
   */
  void clear(id_index rowIndex);

  /**
   * @brief Returns the row index of the pivot. If the column does not have a pivot, returns -1.
   *
   * Only useful for @ref boundarymatrix "boundary" and @ref chainmatrix "chain matrices".
   * 
   * @return Row index of the pivot or -1.
   */
  id_index get_pivot();
  /**
   * @brief Returns the value of the pivot. If the column does not have a pivot, returns 0.
   *
   * Has to have value 1 if \f$ Z_2 \f$ coefficients are used.
   *
   * Only useful for @ref boundarymatrix "boundary" and @ref chainmatrix "chain matrices".
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
   * @tparam Cell_range @ref Cell range with %begin() and %end() method.
   * Has to be ordered by row index if not specified otherwise.
   * @param column @ref Cell range. Only the stored row index and the stored element value
   * (if @ref PersistenceMatrixOptions::is_z2 is false) are token into account for this method.
   * Even if @ref PersistenceMatrixOptions::has_row_access is true, the column index does not need to be correct.
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
   * @brief `this = val * this + column`
   * 
   * @tparam Cell_range @ref Cell range with %begin() and %end() method.
   * Has to be ordered by row index if not specified otherwise.
   * @param val Value to multiply.
   * @param column @ref Cell range. Only the stored row index and the stored element value
   * (if @ref PersistenceMatrixOptions::is_z2 is false) are token into account for this method.
   * Even if @ref PersistenceMatrixOptions::has_row_access is true, the column index does not need to be correct.
   * @return Reference to this column.
   */
  template <class Cell_range>
  PersistenceMatrixColumn& multiply_target_and_add(const Field_element_type& val, const Cell_range& column);

  /**
   * @brief `this = this + column * val`
   * 
   * @tparam Cell_range @ref Cell range with %begin() and %end() method.
   * Has to be ordered by row index if not specified otherwise.
   * @param column @ref Cell range. Only the stored row index and the stored element value
   * (if @ref PersistenceMatrixOptions::is_z2 is false) are token into account for this method.
   * Even if @ref PersistenceMatrixOptions::has_row_access is true, the column index does not need to be correct.
   * @param val Value to multiply.
   * @return Reference to this column.
   */
  template <class Cell_range>
  PersistenceMatrixColumn& multiply_source_and_add(const Cell_range& column, const Field_element_type& val);

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

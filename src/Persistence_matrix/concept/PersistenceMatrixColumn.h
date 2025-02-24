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
 * @brief Concept of the column classes used by the @ref Matrix class. The classes the columns inherit from
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
  using Index = unspecified;                  /**< Type of @ref MatIdx index. */
  using ID_index = unspecified;               /**< Type of @ref IDIdx index. */
  using Dimension = unspecified;              /**< Type for dimension value. */
  using Field_element = unspecified;          /**< Type of a field element. */
  using Entry = unspecified;                  /**< @ref Entry. */
  using Column_settings = unspecified;        /**< Structure giving access to external classes eventually necessary,
                                                   like an entry pool for example. */
  using iterator = unspecified;               /**< Column iterator type. */
  using const_iterator = unspecified;         /**< Column const_iterator type. */
  using reverse_iterator = unspecified;       /**< Column reverse_iterator type. */
  using const_reverse_iterator = unspecified; /**< Column const_reverse_iterator type. */

  /**
   * @brief Constructs an empty column. If @p entryConstructor is not specified or is set to `nullptr`, the column
   * can only be used as a dummy, i.e., no modifying method should be used or there will be a segmentation fault.
   * Same goes for @p operators if @ref PersistenceMatrixOptions::is_z2 is false.
   *
   * @param colSettings Pointer to a setting structure or `nullptr`. The structure should contain all the necessary
   * classes specific to the column type, such as custom allocators. The specificities are this way hidden behind
   * a common interface for all column types. If @p colSettings is not specified or is equal to `nullptr`, the column
   * should still be constructable even though not necessarily "usable".
   */
  PersistenceMatrixColumn(Column_settings* colSettings = nullptr);
  /**
   * @brief Constructs a column from the given range of @ref Matrix::Entry_representative. If the dimension is stored,
   * the cell is assumed to be simplicial and its dimension to be `nonZeroRowIndices length - 1` or `0`.
   * Otherwise, the dimension should be specified with another constructor.
   *
   * @tparam Container Range of @ref Matrix::Entry_representative. Assumed to have a %begin(), %end() and %size()
   * method.
   * @param nonZeroRowIndices Range of @ref Matrix::Entry_representative representing all rows with non zero values.
   * @param colSettings Pointer to an existing setting structure. The structure should contain all the necessary
   * classes specific to the column type, such as custom allocators. The specificities are this way hidden behind
   * a common interface for all column types.
   */
  template <class Container = typename Master_matrix::Boundary>
  PersistenceMatrixColumn(const Container& nonZeroRowIndices,
                          Column_settings* colSettings);
  /**
   * @brief Constructs a column from the given range of @ref Matrix::Entry_representative such that the rows can be
   * accessed. Each new entry in the column is also inserted in a row using @ref Row_access::insert_entry.
   * If the dimension is stored, the cell is assumed to be simplicial and its dimension to be
   * `nonZeroRowIndices length - 1` or `0`. Otherwise, the dimension should be specified with another constructor.
   *
   * @tparam Container Range of @ref Matrix::Entry_representative. Assumed to have a %begin(), %end() and %size()
   * method.
   * @tparam Row_container Either std::map if @ref PersistenceMatrixOptions::has_removable_rows is true or
   * std::vector<Row>.
   * @param columnIndex @ref MatIdx column index that should be specified to the entries.
   * @param nonZeroRowIndices Range of @ref Matrix::Entry_representative representing all rows with non zero values.
   * @param rowContainer Pointer to the row container that will be forwarded to @ref Row_access at construction.
   * @param colSettings Pointer to an existing setting structure. The structure should contain all the necessary
   * classes specific to the column type, such as custom allocators. The specificities are this way hidden behind
   * a common interface for all column types.
   */
  template <class Container = typename Master_matrix::Boundary, class Row_container>
  PersistenceMatrixColumn(Index columnIndex,
                          const Container& nonZeroRowIndices,
                          Row_container* rowContainer,
                          Column_settings* colSettings);
  /**
   * @brief Constructs a column from the given range of @ref Matrix::Entry_representative and stores the given dimension
   * if @ref Column_dimension_option is not a dummy.
   *
   * @tparam Container Range of @ref Matrix::Entry_representative. Assumed to have a %begin(), %end() and %size()
   * method.
   * @param nonZeroChainRowIndices Range of @ref Matrix::Entry_representative representing all rows with non zero
   * values.
   * @param dimension Dimension of the column. Is ignored if the dimension is not stored.
   * @param colSettings Pointer to an existing setting structure. The structure should contain all the necessary
   * classes specific to the column type, such as custom allocators. The specificities are this way hidden behind
   * a common interface for all column types.
   */
  template <class Container = typename Master_matrix::Boundary>
  PersistenceMatrixColumn(const Container& nonZeroChainRowIndices,
                          Dimension dimension,
                          Column_settings* colSettings);
  /**
   * @brief Constructs a column from the given range of @ref Matrix::Entry_representative such that the rows can be
   * accessed. Each new entry in the column is also inserted in a row using @ref Row_access::insert_entry.
   * Stores the given dimension if @ref Column_dimension_option is not a dummy.
   *
   * @tparam Container Range of @ref Matrix::Entry_representative. Assumed to have a %begin(), %end() and %size()
   * method.
   * @tparam Row_container Either std::map if @ref PersistenceMatrixOptions::has_removable_rows is true or
   * std::vector<Row>.
   * @param columnIndex @ref MatIdx column index that should be specified to the entries.
   * @param nonZeroChainRowIndices Range of @ref Matrix::Entry_representative representing all rows with non zero
   * values.
   * @param dimension Dimension of the column. Is ignored if the dimension is not stored.
   * @param rowContainer Pointer to the row container that will be forwarded to @ref Row_access at construction.
   * @param colSettings Pointer to an existing setting structure. The structure should contain all the necessary
   * classes specific to the column type, such as custom allocators. The specificities are this way hidden behind
   * a common interface for all column types.
   */
  template <class Container = typename Master_matrix::Boundary, class Row_container>
  PersistenceMatrixColumn(Index columnIndex,
                          const Container& nonZeroChainRowIndices,
                          Dimension dimension,
                          Row_container* rowContainer,
                          Column_settings* colSettings);
  /**
   * @brief Copy constructor. If @p operators or @p entryConstructor is not a null pointer, its value is kept
   * instead of the one in the copied column.
   *
   * @param column Column to copy.
   * @param colSettings Pointer to a setting structure or `nullptr`. The structure should contain all the necessary
   * classes specific to the column type, such as custom allocators. The specificities are this way hidden behind
   * a common interface for all column types. If @p colSettings is not specified or is equal to `nullptr`, the structure
   * stored in @p column is used instead.
   */
  PersistenceMatrixColumn(const PersistenceMatrixColumn& column,
                          Column_settings* colSettings = nullptr);
  /**
   * @brief Copy constructor with row access.
   * If @p operators or @p entryConstructor stored in @p colSettings is not a null pointer, its value is kept
   * instead of the one in the copied column.
   *
   * @tparam Row_container  Either std::map if @ref PersistenceMatrixOptions::has_removable_rows is true or
   * std::vector<Row>.
   * @param column Column to copy.
   * @param columnIndex @ref MatIdx column index of the new column once copied.
   * @param rowContainer Pointer to the row container that will be forwarded to @ref Row_access.
   * @param colSettings Pointer to a setting structure or `nullptr`. The structure should contain all the necessary
   * classes specific to the column type, such as custom allocators. The specificities are this way hidden behind
   * a common interface for all column types. If @p colSettings is not specified or is equal to `nullptr`, the structure
   * stored in @p column is used instead.
   */
  template <class Row_container>
  PersistenceMatrixColumn(const PersistenceMatrixColumn& column,
                          Index columnIndex,
                          Row_container* rowContainer,
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
   * @return Vector of @ref Field_element. At element \f$ i \f$ of the vector will be stored the value
   * at row \f$ i \f$ of the column.
   */
  std::vector<Field_element> get_content(int columnLength = -1) const;
  /**
   * @brief Indicates if the entry at given row index has value zero.
   *
   * @param rowIndex Row index to look at.
   * @return true If the entry has value zero.
   * @return false Otherwise.
   */
  bool is_non_zero(ID_index rowIndex) const;
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
   * @warning Depending of the column type, the container does not have to contain only the non-zero entries.
   * Even if for most of the types, the size of the container will correspond to the number of non-zero values
   * in the column, it is not always the case. See description of the actual Column class for more details.
   *
   * @return Size of the underlying container.
   */
  std::size_t size() const;

  /**
   * @brief Reorders the column with the given map of row indices. Also changes the column index stored in the
   * entries if row access is enabled and @p columnIndex is not the @ref Matrix::get_null_value "null index".
   *
   * Only useful for @ref basematrix "base" and @ref boundarymatrix "boundary matrices" using lazy swaps.
   *
   * @tparam Row_index_map Map with an %at() method.
   * @param valueMap Map such that `valueMap.at(i)` indicates the new row index of the entry
   * at current row index `i`.
   * @param columnIndex New @ref MatIdx column index of the column. If @ref Matrix::get_null_value "null index",
   * the index does not change. Ignored if the row access is not enabled.
   * Default value: @ref Matrix::get_null_value "null index".
   */
  template <class Row_index_map>
  void reorder(const Row_index_map& valueMap, [[maybe_unused]] Index columnIndex = Matrix::get_null_value<Index>());
  /**
   * @brief Zeros/empties the column.
   *
   * Only useful for @ref basematrix "base" and @ref boundarymatrix "boundary matrices".
   * Used in @ref Matrix::zero_column and in the reduction algorithm for the persistence barcode.
   */
  void clear();
  /**
   * @brief Zeros the entry at given row index.
   *
   * Only useful for @ref basematrix "base" and @ref boundarymatrix "boundary matrices".
   * Used in @ref Matrix::zero_entry and during vine swaps.
   *
   * @warning For @ref Vector_column, do not clear an entry that was already at zero or the results of @ref size and
   * @ref is_empty will be wrong.
   *
   * @param rowIndex Row index of the entry to zero.
   */
  void clear(ID_index rowIndex);

  /**
   * @brief Returns the row index of the pivot. If the column does not have a pivot,
   * returns @ref Matrix::get_null_value "null index".
   *
   * Only useful for @ref boundarymatrix "boundary" and @ref chainmatrix "chain matrices".
   *
   * @return Row index of the pivot or @ref Matrix::get_null_value "null index".
   */
  ID_index get_pivot();
  /**
   * @brief Returns the value of the pivot. If the column does not have a pivot, returns 0.
   *
   * Has to have value 1 if \f$ Z_2 \f$ coefficients are used.
   *
   * Only useful for @ref boundarymatrix "boundary" and @ref chainmatrix "chain matrices".
   *
   * @return The value of the pivot or 0.
   */
  Field_element get_pivot_value();

  /**
   * @brief Returns a begin @ref Entry iterator to iterate over all entries contained in the underlying container.
   *
   * @warning The iterators really just iterate over the underlying container. Depending of the column type,
   * neither the content nor the order is guaranteed. See description of the actual Column class for more details.
   *
   * @return @ref Entry iterator.
   */
  iterator begin() noexcept;
  /**
   * @brief Returns a begin @ref Entry const iterator to iterate over all entries contained in the underlying container.
   *
   * @warning The iterators really just iterate over the underlying container. Depending of the column type,
   * neither the content nor the order is guaranteed. See description of the actual Column class for more details.
   *
   * @return @ref Entry const iterator.
   */
  const_iterator begin() const noexcept;
  /**
   * @brief Returns an end @ref Entry iterator, iterating over all entries contained in the underlying container.
   *
   * @warning The iterators really just iterate over the underlying container. Depending of the column type,
   * neither the content nor the order is guaranteed. See description of the actual Column class for more details.
   *
   * @return @ref Entry iterator.
   */
  iterator end() noexcept;
  /**
   * @brief Returns an end @ref Entry const iterator, iterating over all entries contained in the underlying container.
   *
   * @warning The iterators really just iterate over the underlying container. Depending of the column type,
   * neither the content nor the order is guaranteed. See description of the actual Column class for more details.
   *
   * @return @ref Entry const iterator.
   */
  const_iterator end() const noexcept;
  /**
   * @brief Returns a begin @ref Entry reverse iterator to iterate over all entries contained in the underlying
   * container.
   *
   * @warning The iterators really just iterate over the underlying container. Depending of the column type,
   * neither the content nor the order is guaranteed. See description of the actual Column class for more details.
   *
   * @return @ref Entry reverse iterator.
   */
  reverse_iterator rbegin() noexcept;
  /**
   * @brief Returns a begin @ref Entry const reverse iterator to iterate over all entries contained in the underlying
   * container.
   *
   * @warning The iterators really just iterate over the underlying container. Depending of the column type,
   * neither the content nor the order is guaranteed. See description of the actual Column class for more details.
   *
   * @return @ref Entry const reverse iterator.
   */
  const_reverse_iterator rbegin() const noexcept;
  /**
   * @brief Returns an end @ref Entry reverse iterator, iterating over all entries contained in the underlying
   * container.
   *
   * @warning The iterators really just iterate over the underlying container. Depending of the column type,
   * neither the content nor the order is guaranteed. See description of the actual Column class for more details.
   *
   * @return @ref Entry reverse iterator.
   */
  reverse_iterator rend() noexcept;
  /**
   * @brief Returns an end @ref Entry const reverse iterator, iterating over all entries contained in the underlying
   * container.
   *
   * @warning The iterators really just iterate over the underlying container. Depending of the column type,
   * neither the content nor the order is guaranteed. See description of the actual Column class for more details.
   *
   * @return @ref Entry const reverse iterator.
   */
  const_reverse_iterator rend() const noexcept;

  /**
   * @brief Adds the given @ref Entry range onto the column.
   *
   * @tparam Entry_range @ref Entry range with %begin() and %end() method.
   * Has to be ordered by row index if not specified otherwise.
   * @param column @ref Entry range. Only the stored row index and the stored element value
   * (if @ref PersistenceMatrixOptions::is_z2 is false) are token into account for this method.
   * Even if @ref PersistenceMatrixOptions::has_row_access is true, the column index does not need to be correct.
   * @return Reference to this column.
   */
  template <class Entry_range>
  PersistenceMatrixColumn& operator+=(const Entry_range& column);
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
  PersistenceMatrixColumn& operator*=(const Field_element& val);

  /**
   * @brief `this = val * this + column`
   *
   * @tparam Entry_range @ref Entry range with %begin() and %end() method.
   * Has to be ordered by row index if not specified otherwise.
   * @param val Value to multiply.
   * @param column @ref Entry range. Only the stored row index and the stored element value
   * (if @ref PersistenceMatrixOptions::is_z2 is false) are token into account for this method.
   * Even if @ref PersistenceMatrixOptions::has_row_access is true, the column index does not need to be correct.
   * @return Reference to this column.
   */
  template <class Entry_range>
  PersistenceMatrixColumn& multiply_target_and_add(const Field_element& val, const Entry_range& column);

  /**
   * @brief `this = this + column * val`
   *
   * @tparam Entry_range @ref Entry range with %begin() and %end() method.
   * Has to be ordered by row index if not specified otherwise.
   * @param column @ref Entry range. Only the stored row index and the stored element value
   * (if @ref PersistenceMatrixOptions::is_z2 is false) are token into account for this method.
   * Even if @ref PersistenceMatrixOptions::has_row_access is true, the column index does not need to be correct.
   * @param val Value to multiply.
   * @return Reference to this column.
   */
  template <class Entry_range>
  PersistenceMatrixColumn& multiply_source_and_add(const Entry_range& column, const Field_element& val);

  /**
   * @brief Adds a copy of the given entry at the end of the column. It is therefore assumed that the row index
   * of the entry is higher than the current pivot of the column. Not available for @ref chainmatrix "chain matrices"
   * and is only needed for @ref boundarymatrix "RU matrices".
   *
   * @param entry Entry to push back.
   */
  void push_back(const Entry& entry);

  /**
   * @brief Equality comparator. Equal in the sense that what is "supposed" to be contained in the columns is equal,
   * not what is actually stored in the underlying container. For example, the underlying container of
   * @ref Vector_column can contain entries which were erased explicitly by @ref clear(Index). Those entries should not
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
   * not what is actually stored in the underlying container. For example, the underlying container of
   * @ref Vector_column can contain entries which were erased explicitly by @ref clear(Index). Those entries should not
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

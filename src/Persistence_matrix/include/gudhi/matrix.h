/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef MASTER_MATRIX_H
#define MASTER_MATRIX_H

/** @file matrix.h
 * @brief Contains @ref Gudhi::persistence_matrix::Matrix class.
 */

#include <type_traits>
#include <vector>
#include <unordered_map>
#include <map>
#include <assert.h>
#include <initializer_list>

#include <boost/intrusive/list.hpp>

#include "persistence_matrix_options.h"

#include "Persistence_matrix/overlay_id_to_position_index.h"
#include "Persistence_matrix/overlay_position_to_id_index.h"

#include "Persistence_matrix/matrix_dimension_holders.h"
#include "Persistence_matrix/matrix_row_access.h"
#include "Persistence_matrix/base_swap.h"
#include "Persistence_matrix/base_pairing.h"
#include "Persistence_matrix/ru_pairing.h"
#include "Persistence_matrix/ru_vine_swap.h"
#include "Persistence_matrix/ru_rep_cycles.h"
#include "Persistence_matrix/chain_pairing.h"
#include "Persistence_matrix/chain_vine_swap.h"
#include "Persistence_matrix/chain_rep_cycles.h"

#include "Persistence_matrix/base_matrix.h"
#include "Persistence_matrix/base_matrix_with_column_compression.h"
#include "Persistence_matrix/boundary_matrix.h"
#include "Persistence_matrix/ru_matrix.h"
#include "Persistence_matrix/chain_matrix.h"

#include "Persistence_matrix/columns/cell_types.h"
#include "Persistence_matrix/columns/row_access.h"

#include "Persistence_matrix/columns/column_dimension_holder.h"
#include "Persistence_matrix/columns/chain_column_extra_properties.h"
#include "Persistence_matrix/columns/intrusive_list_column.h"
#include "Persistence_matrix/columns/intrusive_set_column.h"
#include "Persistence_matrix/columns/list_column.h"
#include "Persistence_matrix/columns/set_column.h"
#include "Persistence_matrix/columns/unordered_set_column.h"
#include "Persistence_matrix/columns/vector_column.h"
#include "Persistence_matrix/columns/naive_vector_column.h"
#include "Persistence_matrix/columns/heap_column.h"

namespace Gudhi {
namespace persistence_matrix {

/**
 * @brief Data structure for matrices, and in particular thought for matrices representing filtered complexes
 * in order to compute persistence and/or representative cycles.
 * 
 * The are roughly three types of matrices available and one is selected automatically depending on the options chosen:
 * - a basic matrix which can represent any matrix and therefore will not make any assumption on its content. 
 *   It is the only matrix type with the option of column compression (as it is the only one where it makes sense).
 *   This type is choosen by default when none of the homology related options are set to true: 
 *   @ref has_column_pairings, @ref has_vine_update and @ref can_retrieve_representative_cycles.
 * - a boundary matrix @$f B = R \cdot U @$f which either stores only @$f R @$f or the whole decomposition @$f R @$f 
 *   and @$f U @$f depending on the options. This type is selected if @ref is_of_boundary_type is set to true and 
 *   at least one of the following options: @ref has_column_pairings, @ref has_vine_update and 
 *   @ref can_retrieve_representative_cycles. If only @ref has_column_pairings is true, then only @$f R @$fn is stored,
 *   but if either @ref has_vine_update or @ref can_retrieve_representative_cycles is true, then @$f U @$f also needs 
 *   to be stored. Note that the option @ref is_indexed_by_position will produce a small overhead when set to **false**.
 * - a chain complex matrix representing a `compatible base` of a filtered chain complex (see TODO: cite Clément's zigzag paper here).
 *   This matrix is deduced from the boundary matrix and therefore encodes more or less the same information 
 *   but differently and can therefore be better suited for certain applications. This type can be used the same way 
 *   than the precedent type, only the option @ref is_of_boundary_type has to be set to false. So it is easy to switch
 *   from one representation to the other if one wants to test both. Just note that the option 
 *   @ref is_indexed_by_position will produce a small overhead when set to **true**.
 *
 * Indexation scheme:
 *
 * The indexation system for the different matrix types can be a bit tricky and different methods will not always take
 * the same type of index as input (for optimization purposes). So, to avoid confusion, we will name and define here the 
 * different possibilities, such that we can directly refer to it in the descriptions of the methods. 
 * Note that every column and row in a boundary or chain matrix is always associated to a single simplex/face, 
 * so in order to avoid repeating formulations like "of the simplex associated to the column" all the time, 
 * we will amalgamate both notions together.
 *
 * Let c be a column or a row.
 * - MatIdx: This will correspond to the position of c in the matrix, i.e., underlying_container[MatIdx] = c.
 *   This will be the only public indexing scheme for basic matrices (first of the list above).
 * - PosIdx: This will correspond to the position of c in the current filtration. For boundary matrices, PosIdx
 *   will always be equal to MatIdx, but this is not true for chain matrices, where PosIdx will correspond to
 *   the position of c in the original filtration (so PosIdx == MatIdx as long as no vine swaps are done).
 * - IDIdx: This will correspond to the ID of c in the complex. If an ID is not specified at the insertion of 
 *   a boundary, the default value of IDIdx will be the value of PosIdx at the insertion (so be careful, 
 *   as soon as vine swaps are performed, IDIdx and PosIdx will start to differ). 
 *   Otherwise, IDIdx will take the given value.
 * 
 * @tparam Options Structure encoding all the options of the matrix. 
 * See description of @ref Default_options for more details.
 */
template <class Options = Default_options<> >
class Matrix {
 public:
  using Option_list = Options;	//to make it accessible from the other classes
  using Field_type = typename Options::field_coeff_type;
  using index = typename Options::index_type;
  using dimension_type = typename Options::dimension_type;

  // TODO: move outside?
  /**
   * @brief Type for a bar in the computed barcode. Stores the birth, death and dimension of the bar.
   */
  struct Bar {
    Bar() : dim(-1), birth(-1), death(-1) {}

    Bar(dimension_type dim, int birth, int death) : dim(dim), birth(birth), death(death) {}

    dimension_type dim;
    int birth;
    int death;
  };

  //tags for boost to associate a row and a column to a same cell
  struct matrix_row_tag;
  struct matrix_column_tag;

  using base_hook_matrix_row =
      boost::intrusive::list_base_hook<boost::intrusive::tag<matrix_row_tag>,
                                       boost::intrusive::link_mode<boost::intrusive::auto_unlink> >;
  using base_hook_matrix_list_column =
      boost::intrusive::list_base_hook<boost::intrusive::tag<matrix_column_tag>,
                                       boost::intrusive::link_mode<boost::intrusive::safe_link> >;
  using base_hook_matrix_set_column =
      boost::intrusive::set_base_hook<boost::intrusive::tag<matrix_column_tag>,
                                      boost::intrusive::link_mode<boost::intrusive::safe_link> >;

  //Two dummies are necessary to avoid double inheritance as a cell can inherit both a row and a column hook.
  struct Dummy_row_hook {};
  struct Dummy_column_hook {};

  using row_hook_type = typename std::conditional<Options::has_row_access && Options::has_intrusive_rows,
                                                  base_hook_matrix_row, Dummy_row_hook
                                                 >::type;
  using column_hook_type =
      typename std::conditional<Options::column_type == Column_types::INTRUSIVE_LIST, 
                                base_hook_matrix_list_column,
                                typename std::conditional<Options::column_type == Column_types::INTRUSIVE_SET,
                                                          base_hook_matrix_set_column, 
                                                          Dummy_column_hook
                                                         >::type
                               >::type;

  //Option to store the column index within the cell (additionnaly to the row index). Necessary only with row access.
  using Cell_column_index_option =
      typename std::conditional<Options::has_row_access,
                                Cell_column_index<index>,
                                Dummy_cell_column_index_mixin
                               >::type;
  //Option to store the value of the cell. 
  //Unnecessary for values in Z_2 as there are always 1 (0-valued cells are never stored).
  using Cell_field_element_option =
      typename std::conditional<Options::is_z2,
                                Dummy_cell_field_element_mixin,
                                Cell_field_element<Field_type>
                               >::type;
  using Cell_type = Cell<Matrix<Options> >;

  using cell_rep_type = typename std::conditional<Options::is_z2,
                                                  index,
                                                  std::pair<index, Field_type>
                                                 >::type;

  /**
   * @brief Compaires two pairs, representing a cell (first = row index, second = value), 
   * by their position in the column and not their values. 
   * The two represented cells are therefore assumed to be in the same column.
   */
  struct CellPairComparator {
    bool operator()(const std::pair<index, Field_type>& p1, const std::pair<index, Field_type>& p2) const {
      return p1.first < p2.first;
    };
  };

  /**
   * @brief Compaires two cells by their position in the row. They are assume to be in the same row.
   */
  struct RowCellComp {
    bool operator()(const Cell_type& c1, const Cell_type& c2) const {
      return c1.get_column_index() < c2.get_column_index();
    }
  };

  using Row_type =
      typename std::conditional<Options::has_intrusive_rows,
                                boost::intrusive::list<Cell_type, 
                                                       boost::intrusive::constant_time_size<false>,
                                                       boost::intrusive::base_hook<base_hook_matrix_row>
                                                      >,
                                std::set<Cell_type, RowCellComp>
                               >::type;

  using row_container_type =
      typename std::conditional<Options::has_removable_rows,
                                std::map<index, Row_type>,
                                std::vector<Row_type>
                               >::type;

  //Row access at column level
  using Row_access_option =
      typename std::conditional<Options::has_row_access,
                                Row_access<Matrix<Options> >,
                                Dummy_row_access
                               >::type;
  //Row access at matrix level
  using Matrix_row_access_option =
      typename std::conditional<Options::has_row_access,
                                Matrix_row_access<Row_type, row_container_type, Options::has_removable_rows, index>,
                                Dummy_matrix_row_access
                               >::type;

  template <typename value_type>
  using dictionnary_type =
      typename std::conditional<Options::has_removable_columns,
                                std::unordered_map<unsigned int, value_type>,
                                std::vector<value_type>
                               >::type;

  static const bool isNonBasic =
      Options::has_column_pairings || Options::has_vine_update || Options::can_retrieve_representative_cycles;

  using Column_dimension_option =
      typename std::conditional<isNonBasic,
                                Column_dimension_holder<Matrix<Options> >,
                                Dummy_dimension_holder
                               >::type;
  //Extra information needed for a column when the matrix is a chain matrix. 
  using Chain_column_option =
      typename std::conditional<isNonBasic && !Options::is_of_boundary_type,
                                Chain_column_extra_properties<Matrix<Options> >,
                                Dummy_chain_properties
                               >::type;

  using Heap_column_type = Heap_column<Matrix<Options> >;
  using List_column_type = List_column<Matrix<Options> >;
  using Vector_column_type = Vector_column<Matrix<Options> >;
  using Naive_vector_column_type = Naive_vector_column<Matrix<Options> >;
  using Set_column_type = Set_column<Matrix<Options> >;
  using Unordered_set_column_type = Unordered_set_column<Matrix<Options> >;
  using Intrusive_list_column_type = Intrusive_list_column<Matrix<Options> >;
  using Intrusive_set_column_type = Intrusive_set_column<Matrix<Options> >;

  using Column_type = typename std::conditional<
        Options::column_type == Column_types::HEAP, 
        Heap_column_type,
        typename std::conditional<
             Options::column_type == Column_types::LIST, 
            List_column_type,
            typename std::conditional<
                Options::column_type == Column_types::SET, 
                Set_column_type,
                typename std::conditional<
                    Options::column_type == Column_types::UNORDERED_SET, 
                    Unordered_set_column_type,
                    typename std::conditional<
                        Options::column_type == Column_types::VECTOR, 
                        Vector_column_type,
                        typename std::conditional<
                            Options::column_type == Column_types::INTRUSIVE_LIST, 
                            Intrusive_list_column_type,
                            typename std::conditional<
                                Options::column_type == Column_types::NAIVE_VECTOR,
                                Naive_vector_column_type, 
                                Intrusive_set_column_type
                                >::type
                            >::type
                        >::type
                    >::type
                >::type
            >::type
        >::type;

  using column_container_type =
      typename std::conditional<Options::has_removable_columns, 
                                std::unordered_map<index, Column_type>,
                                std::vector<Column_type>
                               >::type;

  using barcode_type =
      typename std::conditional<Options::has_removable_columns, 
                                std::list<Bar>, 
                                std::vector<Bar>
                               >::type;
  using bar_dictionnary_type = typename std::conditional<Options::has_removable_columns,
                                                         std::unordered_map<index, typename barcode_type::iterator>,
                                                         std::vector<unsigned int>
                                                        >::type;

  using boundary_type = typename std::conditional<Options::is_z2, 
                                                  std::initializer_list<index>,
                                                  std::initializer_list<std::pair<index, Field_type> >
                                                 >::type;

  static const bool dimensionIsNeeded = Options::has_column_pairings && Options::is_of_boundary_type &&
                                        !Options::has_vine_update && !Options::can_retrieve_representative_cycles;

  using Matrix_dimension_option = typename std::conditional<
      Options::has_matrix_maximal_dimension_access || dimensionIsNeeded,
      typename std::conditional<Options::has_removable_columns, Matrix_all_dimension_holder<dimension_type>,
                                Matrix_max_dimension_holder<dimension_type>
                               >::type,
      Dummy_matrix_dimension_holder
    >::type;

  using Base_matrix_type =
      typename std::conditional<Options::has_column_compression, 
                                Base_matrix_with_column_compression<Matrix<Options> >,
                                Base_matrix<Matrix<Options> >
                               >::type;
  using Boundary_matrix_type = Boundary_matrix<Matrix<Options> >;
  using RU_matrix_type = RU_matrix<Matrix<Options> >;
  using Chain_matrix_type = Chain_matrix<Matrix<Options> >;

  template <class Base>
  using Base_swap_option = typename std::conditional<Options::has_vine_update || Options::has_column_and_row_swaps,
                                                     Base_swap<Matrix<Options>, Base>, 
                                                     Dummy_base_swap
                                                    >::type;
  using Base_pairing_option = typename std::conditional<Options::has_column_pairings && !Options::has_vine_update &&
                                                            !Options::can_retrieve_representative_cycles,
                                                        Base_pairing<Matrix<Options> >, 
                                                        Dummy_base_pairing
                                                       >::type;

  using RU_pairing_option = typename std::conditional<Options::has_column_pairings && !Options::has_vine_update,
                                                      RU_pairing<Matrix<Options> >, 
                                                      Dummy_ru_pairing
                                                     >::type;
  using RU_vine_swap_option =
      typename std::conditional<Options::has_vine_update,
                                RU_vine_swap<Matrix<Options> >,
                                Dummy_ru_vine_swap
                               >::type;
  using RU_representative_cycles_option =
      typename std::conditional<Options::can_retrieve_representative_cycles, 
                                RU_representative_cycles<Matrix<Options> >,
                                Dummy_ru_representative_cycles
                               >::type;

  using Chain_pairing_option = typename std::conditional<Options::has_column_pairings && !Options::has_vine_update,
                                                         Chain_pairing<Matrix<Options> >,
                                                         Dummy_chain_pairing
                                                        >::type;
  using Chain_vine_swap_option = typename std::conditional<Options::has_vine_update, 
                                                           Chain_vine_swap<Matrix<Options> >,
                                                           Dummy_chain_vine_swap
                                                          >::type;
  using Chain_representative_cycles_option =
      typename std::conditional<Options::can_retrieve_representative_cycles,
                                Chain_representative_cycles<Matrix<Options> >,
                                Dummy_chain_representative_cycles
                               >::type;

  using cycle_type = std::vector<index>;

  //Return types to factorize the corresponding methods
  using returned_column_type =
      typename std::conditional<!isNonBasic && Options::has_column_compression,
                                const Column_type,
                                Column_type
                               >::type;
  using returned_row_type =
      typename std::conditional<!isNonBasic && Options::has_column_compression,
                                const Row_type,
                                Row_type
                               >::type;
  using insertion_return_type =
      typename std::conditional<Options::is_of_boundary_type || !isNonBasic,
                                void,
                                std::vector<cell_rep_type>
                               >::type;

  /**
   * @brief Default constructor.
   */
  Matrix();
  /**
   * @brief Constructs a new matrix from the given matrix. 
   * If the columns are representing a boundary matrix, the indices of the simplices are assumed to be 
   * consecutifs and starting with 0.
   *
   * See options descriptions for futher details on how the given matrix is handled.
   * 
   * @tparam Container_type Range type for a column. Assumed to have a begin(), end() and size() method.
   * @param columns For a general/base matrix, the columns are copied as is. 
   * If options related to homology are activated, @p columns is interpreted as a boundary matrix of a 
   * **simplicial** complex. 
   * In this case, `columns[i]` should store the boundary of simplex `i` as an ordered list of indices 
   * of its facets (again those indices correspond to their respective position in the matrix). 
   * Therefore the indices of the simplices are assumed to be consecutifs and starting with 0 
   * (an empty boundary is interpreted as a vertex boundary and not as a non existing simplex). 
   * All dimensions up to the maximal dimension of interest have to be present. If only a higher dimension is of 
   * interest and not everything should be stored, then use the `insert_boundary` method instead (after creating 
   * the matrix with the `Matrix(int numberOfColumns)` constructor preferably).
   * If the persistence barcode has to be computed from this matrix, the simplices are also assumed to be ordered by 
   * appearance order in the filtration. Also, depending of the options, the matrix is eventually reduced on the fly 
   * or converted into a chain complex base, so the new matrix is not always identical to the old one.
   */
  template <class Container_type = boundary_type>
  Matrix(const std::vector<Container_type>& columns);
  /**
   * @brief Constructs a new empty matrix and reserves space for the given number of columns.
   * 
   * @param numberOfColumns Number of columns to reserve space for.
   */
  Matrix(int numberOfColumns);
  /**
   * @brief Constructs a new empty matrix with the given comparator functions. Only available when those comparators 
   * are necessary, i.e., when **all** following options have following values:
   *   - @ref is_of_boundary_type = false
   *   - @ref has_vine_update = true
   *   - @ref has_column_pairings = false
   * 
   * Those comparators are necesseray to distinguish cases in a vine update. When the matrix is of boundary type 
   * (i.e., a RU decomposition) or if the column pairing is activated (i.e., the barcode is stored), the comparators
   * can be easily deduced without overhead. If neither are true, we assume that one has additional information
   * outside of the matrix about the barcode to provide a better suited comparator adapted to the situation 
   * (as in the implementation of the Zigzag algorithm for example TODO: ref to zigzag module.)
   * 
   * @tparam BirthComparatorFunction Type of the birth comparator: (unsigned int, unsigned int) -> bool
   * @tparam DeathComparatorFunction Type of the death comparator: (unsigned int, unsigned int) -> bool
   * @param birthComparator Method taking two simplex indices as parameter and returns true if and only if the first 
   * simplex is associated to a bar with strictly smaller birth than the bar associated to the second one.
   * @param deathComparator Method taking two simplex indices as parameter and returns true if and only if the first 
   * simplex is associated to a bar with strictly smaller death than the bar associated to the second one.
   */
  template <typename BirthComparatorFunction, typename DeathComparatorFunction>
  Matrix(BirthComparatorFunction&& birthComparator, DeathComparatorFunction&& deathComparator);
  /**
   * @brief Constructs a new matrix from the given matrix with the given comparator functions. 
   * Only available when those comparators are necessary, i.e., when **all** following options have following values:
   *   - @ref is_of_boundary_type = false
   *   - @ref has_vine_update = true
   *   - @ref has_column_pairings = false
   *
   * See description of @ref Matrix(const std::vector<Container_type>& columns) for more information about 
   * @p orderedBoundaries and 
   * @ref Matrix(BirthComparatorFunction&& birthComparator, DeathComparatorFunction&& deathComparator) 
   * for more information about the comparators.
   * 
   * @tparam BirthComparatorFunction Type of the birth comparator: (unsigned int, unsigned int) -> bool
   * @tparam DeathComparatorFunction Type of the death comparator: (unsigned int, unsigned int) -> bool
   * @tparam Boundary_type Range type for a column. Assumed to have a begin(), end() and size() method.
   * @param orderedBoundaries Vector of boundaries in filtration order. Indexed continously starting at 0.
   * @param birthComparator Method taking two simplex indices as parameter and returns true if and only if the first 
   * simplex is associated to a bar with strictly smaller birth than the bar associated to the second one.
   * @param deathComparator Method taking two simplex indices as parameter and returns true if and only if the first 
   * simplex is associated to a bar with strictly smaller death than the bar associated to the second one.
   */
  template <typename BirthComparatorFunction, typename DeathComparatorFunction, class Boundary_type = boundary_type>
  Matrix(const std::vector<Boundary_type>& orderedBoundaries, 
         BirthComparatorFunction&& birthComparator,
         DeathComparatorFunction&& deathComparator);
  /**
   * @brief Constructs a new empty matrix and reserves space for the given number of columns.
   * Only available when those comparators are necessary, i.e., when **all** following options have following values:
   *   - @ref is_of_boundary_type = false
   *   - @ref has_vine_update = true
   *   - @ref has_column_pairings = false
   *
   * See description of 
   * @ref Matrix(BirthComparatorFunction&& birthComparator, DeathComparatorFunction&& deathComparator) 
   * for more information about the comparators.
   * 
   * @tparam BirthComparatorFunction Type of the birth comparator: (unsigned int, unsigned int) -> bool
   * @tparam DeathComparatorFunction Type of the death comparator: (unsigned int, unsigned int) -> bool
   * @param numberOfColumns Number of columns to reserve space for.
   * @param birthComparator Method taking two simplex indices as parameter and returns true if and only if the first 
   * simplex is associated to a bar with strictly smaller birth than the bar associated to the second one.
   * @param deathComparator Method taking two simplex indices as parameter and returns true if and only if the first 
   * simplex is associated to a bar with strictly smaller death than the bar associated to the second one.
   */
  template <typename BirthComparatorFunction, typename DeathComparatorFunction>
  Matrix(unsigned int numberOfColumns, 
         BirthComparatorFunction&& birthComparator,
         DeathComparatorFunction&& deathComparator);
  /**
   * @brief Copy constructor.
   * 
   * @param matrixToCopy Matrix to copy.
   */
  Matrix(const Matrix& matrixToCopy);
  /**
   * @brief Move constructor.
   * After the move, the given matrix will be empty.
   * 
   * @param other Matrix to move.
   */
  Matrix(Matrix&& other) noexcept;

  // (TODO: if there is no row access and the column type corresponds to the internal column type of the matrix, 
  // moving the column instead of copying it should be possible. Is it worth implementing it?)
  /**
   * @brief Inserts a new column at the end of the matrix by copying the given column. 
   * Only available when **all** of the following options are false:
   *   - @ref has_column_pairings
   *   - @ref has_vine_update
   *   - @ref can_retrieve_representative_cycles
   *
   * Otherwise use @ref insert_boundary which will deduce a new column from the boundary given.
   * 
   * @tparam Container_type Range type for a column. Assumed to have a begin(), end() and size() method.
   * @param column Column to be inserted.
   */
  template <class Container_type>
  void insert_column(const Container_type& column);
  /**
   * @brief Inserts a new column at the given index by copying the given column. There should not be any other column 
   * inserted at that index which was not explicitely removed before.
   * Only available if @ref has_removable_columns is true and **all** of the following options are false:
   *   - @ref has_column_pairings
   *   - @ref has_vine_update
   *   - @ref can_retrieve_representative_cycles
   *   - @ref has_column_compression
   * 
   * @tparam Container_type Range type for a column. Assumed to have a begin(), end() and size() method.
   * @param column Column to be inserted.
   * @param columnIndex Index to which the column has to be inserted.
   */
  template <class Container_type>
  void insert_column(const Container_type& column, int columnIndex);
  /**
   * @brief Inserts at the end of the matrix a new column corresponding to the given boundary. 
   * This means that we assume that the boundaries are inserted in the order of the filtration. 
   * The content of the new column will vary depending on the underlying type of the matrix.
   *
   * - If it is a basic matrix type, the boundary is copied as it is, i.e., the method is equivalent to 
   *   @ref insert_column.
   * - If it is a boundary type matrix and only R is stored, the boundary is also just copied. The column will only be 
   *   reduced later when the barcode is requested in order to apply some optimisations with the additional knowledge.
   *   Hence, the barcode will also not be updated, so call @ref get_current_barcode only when the matrix is complete.
   * - If it is a boundary type matrix and both R and U are stored, the new boundary is stored in its reduced form and
   *   the barcode, if active, is also updated.
   * - If it is a chain type matrix, the new column is of the form 
   *   `simplex index + linear combination of older columns`, where the combination is deduced while reducing the 
   *   given boundary. The simplex index will be new even if the same simplex was already inserted and then removed 
   *   once before. If the barcode is stored, it will also be updated.
   * 
   * @tparam Boundary_type Range type for a column. Assumed to have a begin(), end() and size() method.
   * @param boundary Boundary generating the new column.
   * @param dim Dimension of the face whose boundary is given. If the complex is simplicial, 
   * this parameter can be omitted as it can be deduced from the size of the boundary.
   * @return If it is a chain matrix, the method returns the indices of the unpaired chains used to reduce the boundary.
   * Otherwise, nothing.
   */
  template <class Boundary_type = boundary_type>
  insertion_return_type insert_boundary(const Boundary_type& boundary, dimension_type dim = -1);
  /**
   * @brief Only avalaible for matrices interfaced by IDIdx (in opposition to PosIdx --- see [TODO: ref to class description] 
   * for meaning and option @ref is_indexed_by_position). 
   * It does the same as the other version, but one chooses the used ID. Note that the given index needs to be new, 
   * the IDs can not "come back".
   * 
   * @tparam Boundary_type Range type for a column. Assumed to have a begin(), end() and size() method.
   * @param simplexIndex Index to be used to indentify the new face.
   * @param boundary Boundary generating the new column.
   * @param dim Dimension of the face whose boundary is given. If the complex is simplicial, 
   * this parameter can be omitted as it can be deduced from the size of the boundary.
   * @return If it is a chain matrix, the method returns the indices of the unpaired chains used to reduce the boundary.
   * Otherwise, nothing.
   */
  template <class Boundary_type>
  insertion_return_type insert_boundary(index simplexIndex, const Boundary_type& boundary, dimension_type dim = -1);

  // base
  // base comp: non const because of path compression in union-find
  // boundary
  // ru: inR = true forced
  // chain : id but not really id
  // id to pos
  // pos to id
  /**
   * @brief Returns the column at the given index for base matrices. For boundary
   * 
   * @param columnIndex 
   * @return returned_column_type& 
   */
  returned_column_type& get_column(index columnIndex);
  // chain
  // id to pos
  // pos to id
  const Column_type& get_column(index columnIndex) const;
  // ru
  const Column_type& get_column(index columnIndex, bool inR);
  // Warning: the get_column_index() function of the row cells returns not
  // the expected type of indices: for boundary matrices, it will return
  // the simplex number and for chain matrices, it will return the effectiv
  // column index, independently of the indexing chosen in the options.
  // get_row(rowIndex) --> simplex ID (=/= columnIndex)
  // base
  // boundary
  // ru: inR = true forced
  // chain
  // id to pos
  // pos to id
  returned_row_type& get_row(index rowIndex);
  // base comp
  // chain
  // id to pos
  // pos to id
  const Row_type& get_row(index rowIndex) const;
  // ru
  const Row_type& get_row(index rowIndex, bool inR);

  // base: (as direct link between column and row is not guaranteed)
  void erase_column(index columnIndex);
  // base: assumes the row is empty, just thought as an index cleanup
  // base comp: assumes the row is empty, just thought as an index cleanup
  // boundary: indirect, assumes the row is empty
  // ru: indirect, assumes the row is empty, only erase row in R, as U will never have an empty row
  // chain: indirect, assumes the row is empty
  // id to pos
  // pos to id
  void erase_row(index rowIndex);
  // boundary: update barcode if already computed, does not verify if it really was maximal
  // ru
  // chain
  // id to pos
  // pos to id
  void remove_maximal_simplex(index columnIndex);

  // boundary: indirect
  // ru
  // chain: indirect
  // id to pos
  // pos to id
  dimension_type get_max_dimension() const;
  // base
  // base comp
  // boundary
  // ru
  // chain
  // id to pos
  // pos to id
  unsigned int get_number_of_columns() const;
  // boundary
  // ru
  // chain
  // id to pos
  // pos to id
  dimension_type get_column_dimension(index columnIndex) const;

  //*** targetColumn += sourceColumn
  // base
  // base comp: modifies all similar columns to target together, not only target
  // boundary: avoid calling with pairing option or make it such that it makes sense for persistence
  // ru: avoid calling with specialized options or make it such that it makes sense for persistence
  // chain: avoid calling with specialized options or make it such that it makes sense for persistence
  // id to pos
  // pos to id
  template <typename Index_type>
  std::enable_if_t<std::is_integral_v<Index_type>> add_to(Index_type sourceColumnIndex, Index_type targetColumnIndex);
  // base
  // base comp
  template <class Cell_range>
  std::enable_if_t<!std::is_integral_v<Cell_range>> add_to(const Cell_range& sourceColumn, index targetColumnIndex);

  //*** targetColumn = (targetColumn * coefficient) + sourceColumn
  // base
  // base comp: modifies all similar columns to target together, not only target
  // boundary: avoid calling with pairing option or make it such that it makes sense for persistence
  // ru: avoid calling with specialized options or make it such that it makes sense for persistence
  // chain: avoid calling with specialized options or make it such that it makes sense for persistence
  // id to pos
  // pos to id
  template <typename Index_type>
  std::enable_if_t<std::is_integral_v<Index_type>> add_to(Index_type sourceColumnIndex,
                                                          const Field_type& coefficient,
                                                          Index_type targetColumnIndex);
  // base
  // base comp
  template <class Cell_range>
  std::enable_if_t<!std::is_integral_v<Cell_range>> add_to(const Cell_range& sourceColumn,
                                                           const Field_type& coefficient,
                                                           index targetColumnIndex);

  //*** targetColumn += (coefficient * sourceColumn)
  // base
  // base comp: modifies all similar columns to target together, not only target
  // boundary: avoid calling with pairing option or make it such that it makes sense for persistence
  // ru: avoid calling with specialized options or make it such that it makes sense for persistence
  // chain: avoid calling with specialized options or make it such that it makes sense for persistence
  // id to pos
  // pos to id
  template <typename Index_type>
  std::enable_if_t<std::is_integral_v<Index_type>> add_to(const Field_type& coefficient,
                                                          Index_type sourceColumnIndex,
                                                          Index_type targetColumnIndex);
  // base
  // base comp
  template <class Cell_range>
  std::enable_if_t<!std::is_integral_v<Cell_range>> add_to(const Field_type& coefficient,
                                                           const Cell_range& sourceColumn,
                                                           index targetColumnIndex);

  // base
  // boundary: avoid calling with pairing option or make it such that it makes sense for persistence
  // ru: inR = true forced, avoid calling with specialized options or make it such that it makes sense for persistence
  // id to pos
  void zero_cell(index columnIndex, index rowIndex);
  // ru
  void zero_cell(index columnIndex, index rowIndex, bool inR);
  // base
  // boundary: avoid calling with pairing option or make it such that it makes sense for persistence
  // ru: inR = true forced, avoid calling with specialized options or make it such that it makes sense for persistence
  // id to pos
  void zero_column(index columnIndex);
  // ru
  void zero_column(index columnIndex, bool inR);
  // base
  // base comp
  // boundary
  // ru: inR = true forced
  // chain
  // id to pos
  // pos to id
  bool is_zero_cell(index columnIndex, index rowIndex);
  // ru
  bool is_zero_cell(index columnIndex, index rowIndex, bool inR) const;
  // base
  // base comp
  // boundary
  // ru: inR = true forced
  // chain: just for sanity checks as a valid chain matrix never has an empty column.
  // id to pos
  // pos to id
  bool is_zero_column(index columnIndex);
  // ru
  bool is_zero_column(index columnIndex, bool inR);

  // ru: assumes that pivot exists
  // chain: assumes that pivot exists
  // id to pos
  // pos to id
  index get_column_with_pivot(index simplexIndex) const;
  // boundary
  // ru: only in R
  // chain
  // id to pos
  // pos to id
  int get_pivot(index columnIndex);

  Matrix& operator=(Matrix other);
  friend void swap(Matrix& matrix1, Matrix& matrix2) { swap(matrix1.matrix_, matrix2.matrix_); }

  void print();  // for debug

  //***** access to optionnal methods

  //*** Persistence diagram
  // boundary
  // ru: indirect through const version
  // chain: indirect through const version
  // id to pos
  // pos to id: indirect
  const barcode_type& get_current_barcode();
  // chain
  // ru
  // pos to id
  const barcode_type& get_current_barcode() const;

  //*** swap/vine
  // base
  // boundary: does not update barcode
  void swap_columns(index columnIndex1, index columnIndex2);
  // base
  // boundary: does not update barcode
  void swap_rows(index rowIndex1, index rowIndex2);
  // base
  // boundary: does not update barcode
  // id to pos: does not update barcode
  void swap_at_indices(index index1, index index2);
  // ru: returns true if barcode was changed
  // pos to id
  bool vine_swap_with_z_eq_1_case(index index);  // by column position with ordered column container
  // chain: returns index which was not modified, ie new i+1
  // id to pos
  index vine_swap_with_z_eq_1_case(index columnIndex1,
                                   index columnIndex2);  // by column id with potentielly unordered column container
  // ru: returns true if barcode was changed
  // pos to id
  bool vine_swap(index index);  // by column position with ordered column container
  // chain: returns index which was not modified, ie new i+1
  // id to pos
  index vine_swap(index columnIndex1, index columnIndex2);  // by column id with potentielly unordered column container

  //*** Representative cycles
  // chain
  // ru
  // id to pos
  // pos to id
  void update_representative_cycles();
  // chain
  // ru
  // id to pos
  // pos to id
  const std::vector<cycle_type>& get_representative_cycles();
  // chain
  // ru
  // id to pos
  // pos to id
  const cycle_type& get_representative_cycle(const Bar& bar);

 private:
  using matrix_type = typename std::conditional<
      isNonBasic,
      typename std::conditional<
          Options::is_of_boundary_type,
          typename std::conditional<
              Options::has_vine_update || Options::can_retrieve_representative_cycles,
              typename std::conditional<Options::is_indexed_by_position, 
                                        RU_matrix_type,
                                        Id_to_position_indexation_overlay<RU_matrix_type, Matrix<Options> >
                                       >::type,
              typename std::conditional<Options::is_indexed_by_position, 
                                        Boundary_matrix_type,
                                        Id_to_position_indexation_overlay<Boundary_matrix_type, Matrix<Options>>
                                       >::type
          >::type,
          typename std::conditional<Options::is_indexed_by_position,
                                    Position_to_id_indexation_overlay<Chain_matrix_type, Matrix<Options>>,
                                    Chain_matrix_type
                                   >::type
      >::type,
      Base_matrix_type
    >::type;

  matrix_type matrix_;

  static constexpr void _assert_options();
};

template <class Options>
inline Matrix<Options>::Matrix() 
{
  static_assert(
      Options::is_of_boundary_type || !Options::has_vine_update || Options::has_column_pairings,
      "When no barcode is recorded with vine swaps, comparaison functions for the columns have to be provided.");
  _assert_options();
}

template <class Options>
template <class Container_type>
inline Matrix<Options>::Matrix(const std::vector<Container_type>& columns) : matrix_(columns) 
{
  static_assert(Options::is_of_boundary_type || !Options::has_vine_update || Options::has_column_pairings,
                "When no barcode is recorded with vine swaps for chain matrices, comparaison functions for the columns "
                "have to be provided.");
  assert(Field_type::get_characteristic() != 0 &&
         "Columns cannot be initialized if the coefficient field characteristic is not specified. "
         "Use a compile-time characteristic initialized field type or use another constructor and call coefficient "
         "initializer of the chosen field class.");
  _assert_options();
}

template <class Options>
inline Matrix<Options>::Matrix(int numberOfColumns) : matrix_(numberOfColumns) 
{
  static_assert(Options::is_of_boundary_type || !Options::has_vine_update || Options::has_column_pairings,
                "When no barcode is recorded with vine swaps for chain matrices, comparaison functions for the columns "
                "have to be provided.");
  _assert_options();
}

template <class Options>
template <typename BirthComparatorFunction, typename DeathComparatorFunction>
inline Matrix<Options>::Matrix(BirthComparatorFunction&& birthComparator, DeathComparatorFunction&& deathComparator)
    : matrix_(birthComparator, deathComparator) 
{
  static_assert(
      !Options::is_of_boundary_type && Options::has_vine_update && !Options::has_column_pairings,
      "Constructor only available for chain matrices when vine swaps are enabled, but barcodes are not recorded.");
  _assert_options();
}

template <class Options>
template <typename BirthComparatorFunction, typename DeathComparatorFunction, class Boundary_type>
inline Matrix<Options>::Matrix(const std::vector<Boundary_type>& orderedBoundaries,
                               BirthComparatorFunction&& birthComparator, 
                               DeathComparatorFunction&& deathComparator)
    : matrix_(orderedBoundaries, birthComparator, deathComparator) 
{
  static_assert(
      !Options::is_of_boundary_type && Options::has_vine_update && !Options::has_column_pairings,
      "Constructor only available for chain matrices when vine swaps are enabled, but barcodes are not recorded.");
  assert(Field_type::get_characteristic() != 0 &&
         "Columns cannot be initialized if the coefficient field characteristic is not specified. "
         "Use a compile-time characteristic initialized field type or use another constructor and call coefficient "
         "initializer of the chosen field class.");
  _assert_options();
}

template <class Options>
template <typename BirthComparatorFunction, typename DeathComparatorFunction>
inline Matrix<Options>::Matrix(unsigned int numberOfColumns, 
                               BirthComparatorFunction&& birthComparator,
                               DeathComparatorFunction&& deathComparator)
    : matrix_(numberOfColumns, birthComparator, deathComparator) 
{
  static_assert(
      !Options::is_of_boundary_type && Options::has_vine_update && !Options::has_column_pairings,
      "Constructor only available for chain matrices when vine swaps are enabled, but barcodes are not recorded.");
  _assert_options();
}

template <class Options>
inline Matrix<Options>::Matrix(const Matrix& matrixToCopy) : matrix_(matrixToCopy.matrix_) {
  _assert_options();
}

template <class Options>
inline Matrix<Options>::Matrix(Matrix&& other) noexcept : matrix_(std::move(other.matrix_)) {
  _assert_options();
}

template <class Options>
template <class Container_type>
inline void Matrix<Options>::insert_column(const Container_type& column) {
  assert(Field_type::get_characteristic() != 0 &&
         "Columns cannot be initialized if the coefficient field characteristic is not specified. "
         "Use a compile-time characteristic initialized field type or call coefficient initializer of the chosen field "
         "class.");
  static_assert(
      !isNonBasic,
      "'insert_column' not available for the chosen options. The input has to be in the form of a simplex boundary.");
  matrix_.insert_column(column);
}

template <class Options>
template <class Container_type>
inline void Matrix<Options>::insert_column(const Container_type& column, int columnIndex) {
  assert(Field_type::get_characteristic() != 0 &&
         "Columns cannot be initialized if the coefficient field characteristic is not specified. "
         "Use a compile-time characteristic initialized field type or call coefficient initializer of the chosen field "
         "class.");
  static_assert(!isNonBasic && !Options::has_column_compression,
                "'insert_column' with those parameters is not available for the chosen options.");
  matrix_.insert_column(column, columnIndex);
}

template <class Options>
template <class Boundary_type>
inline typename Matrix<Options>::insertion_return_type Matrix<Options>::insert_boundary(const Boundary_type& boundary,
                                                                                        dimension_type dim) {
  assert(Field_type::get_characteristic() != 0 &&
         "Columns cannot be initialized if the coefficient field characteristic is not specified. "
         "Use a compile-time characteristic initialized field type or call coefficient initializer of the chosen field "
         "class.");
  return matrix_.insert_boundary(boundary, dim);
}

template <class Options>
template <class Boundary_type>
inline typename Matrix<Options>::insertion_return_type Matrix<Options>::insert_boundary(index simplexIndex,
                                                                                        const Boundary_type& boundary,
                                                                                        dimension_type dim) {
  assert(Field_type::get_characteristic() != 0 &&
         "Columns cannot be initialized if the coefficient field characteristic is not specified. "
         "Use a compile-time characteristic initialized field type or call coefficient initializer of the chosen field "
         "class.");
  static_assert(isNonBasic && !Options::is_indexed_by_position, "Only enabled for matrices index by simplex ID.");
  return matrix_.insert_boundary(simplexIndex, boundary, dim);
}

template <class Options>
inline typename Matrix<Options>::returned_column_type& Matrix<Options>::get_column(index columnIndex) {
  return matrix_.get_column(columnIndex);
}

template <class Options>
inline const typename Matrix<Options>::Column_type& Matrix<Options>::get_column(index columnIndex) const {
  return matrix_.get_column(columnIndex);
}

template <class Options>
inline const typename Matrix<Options>::Column_type& Matrix<Options>::get_column(index columnIndex, bool inR) {
  static_assert(isNonBasic && 
                    Options::is_of_boundary_type &&
                    (Options::has_vine_update || Options::can_retrieve_representative_cycles) &&
                    Options::is_indexed_by_position,
                "Only enabled for position indexed RU matrices.");

  return matrix_.get_column(columnIndex, inR);
}

template <class Options>
inline typename Matrix<Options>::returned_row_type& Matrix<Options>::get_row(index rowIndex) {
  static_assert(Options::has_row_access, "'get_row' is not available for the chosen options.");

  return matrix_.get_row(rowIndex);
}

template <class Options>
inline const typename Matrix<Options>::Row_type& Matrix<Options>::get_row(index rowIndex) const {
  static_assert(Options::has_row_access, "'get_row' is not available for the chosen options.");

  return matrix_.get_row(rowIndex);
}

template <class Options>
inline const typename Matrix<Options>::Row_type& Matrix<Options>::get_row(index rowIndex, bool inR) {
  static_assert(Options::has_row_access, "'get_row' is not available for the chosen options.");
  static_assert(isNonBasic && Options::is_of_boundary_type &&
                    (Options::has_vine_update || Options::can_retrieve_representative_cycles) &&
                    Options::is_indexed_by_position,
                "Only enabled for position indexed RU matrices.");

  return matrix_.get_row(rowIndex, inR);
}

template <class Options>
inline void Matrix<Options>::erase_column(index columnIndex) {
  static_assert(Options::has_removable_columns && !isNonBasic && !Options::has_column_compression,
                "'erase_column' is not available for the chosen options.");

  matrix_.erase_column(columnIndex);
}

template <class Options>
inline void Matrix<Options>::erase_row(index rowIndex) {
  static_assert(!isNonBasic || Options::has_removable_rows, "'erase_row' is not available for the chosen options.");

  matrix_.erase_row(rowIndex);
}

template <class Options>
inline void Matrix<Options>::remove_maximal_simplex(index columnIndex) {
  static_assert(Options::has_removable_columns && isNonBasic,
                "'remove_maximal_simplex' is not available for the chosen options.");

  matrix_.remove_maximal_simplex(columnIndex);
}

template <class Options>
inline typename Matrix<Options>::dimension_type Matrix<Options>::get_max_dimension() const {
  static_assert(isNonBasic, "'get_max_dimension' is not available for the chosen options.");

  return matrix_.get_max_dimension();
}

template <class Options>
inline unsigned int Matrix<Options>::get_number_of_columns() const {
  return matrix_.get_number_of_columns();
}

template <class Options>
inline typename Matrix<Options>::dimension_type Matrix<Options>::get_column_dimension(index columnIndex) const {
  static_assert(isNonBasic, "'get_column_dimension' is not available for the chosen options.");

  return matrix_.get_column_dimension(columnIndex);
}

template <class Options>
template <typename Index_type>
inline std::enable_if_t<std::is_integral_v<Index_type>> Matrix<Options>::add_to(Index_type sourceColumnIndex,
                                                                                Index_type targetColumnIndex) {
  return matrix_.add_to(sourceColumnIndex, targetColumnIndex);
}

template <class Options>
template <class Cell_range>
inline std::enable_if_t<!std::is_integral_v<Cell_range>> Matrix<Options>::add_to(const Cell_range& sourceColumn,
                                                                                 index targetColumnIndex) {
  static_assert(!isNonBasic,
                "For boundary or chain matrices, only additions with columns inside the matrix is allowed to maintain "
                "algebraic consistency.");

  return matrix_.add_to(sourceColumn, targetColumnIndex);
}

template <class Options>
template <typename Index_type>
inline std::enable_if_t<std::is_integral_v<Index_type>> Matrix<Options>::add_to(Index_type sourceColumnIndex,
                                                                                const Field_type& coefficient,
                                                                                Index_type targetColumnIndex) {
  return matrix_.add_to(sourceColumnIndex, coefficient, targetColumnIndex);
}

template <class Options>
template <class Cell_range>
inline std::enable_if_t<!std::is_integral_v<Cell_range>> Matrix<Options>::add_to(const Cell_range& sourceColumn,
                                                                                 const Field_type& coefficient,
                                                                                 index targetColumnIndex) {
  static_assert(!isNonBasic,
                "For boundary or chain matrices, only additions with columns inside the matrix is allowed to maintain "
                "algebraic consistency.");

  return matrix_.add_to(sourceColumn, coefficient, targetColumnIndex);
}

template <class Options>
template <typename Index_type>
inline std::enable_if_t<std::is_integral_v<Index_type>> Matrix<Options>::add_to(const Field_type& coefficient,
                                                                                Index_type sourceColumnIndex,
                                                                                Index_type targetColumnIndex) {
  return matrix_.add_to(coefficient, sourceColumnIndex, targetColumnIndex);
}

template <class Options>
template <class Cell_range>
inline std::enable_if_t<!std::is_integral_v<Cell_range>> Matrix<Options>::add_to(const Field_type& coefficient,
                                                                                 const Cell_range& sourceColumn,
                                                                                 index targetColumnIndex) {
  static_assert(!isNonBasic,
                "For boundary or chain matrices, only additions with columns inside the matrix is allowed to maintain "
                "algebraic consistency.");

  return matrix_.add_to(coefficient, sourceColumn, targetColumnIndex);
}

template <class Options>
inline void Matrix<Options>::zero_cell(index columnIndex, index rowIndex) {
  static_assert(Options::is_of_boundary_type && !Options::has_column_compression,
                "'zero_cell' is not available for the chosen options.");

  return matrix_.zero_cell(columnIndex, rowIndex);
}

template <class Options>
inline void Matrix<Options>::zero_cell(index columnIndex, index rowIndex, bool inR) {
  static_assert(isNonBasic && 
                    Options::is_of_boundary_type &&
                    (Options::has_vine_update || Options::can_retrieve_representative_cycles) &&
                    Options::is_indexed_by_position,
                "Only enabled for RU matrices.");

  return matrix_.zero_cell(columnIndex, rowIndex, inR);
}

template <class Options>
inline void Matrix<Options>::zero_column(index columnIndex) {
  static_assert(Options::is_of_boundary_type && !Options::has_column_compression,
                "'zero_column' is not available for the chosen options.");

  return matrix_.zero_column(columnIndex);
}

template <class Options>
inline void Matrix<Options>::zero_column(index columnIndex, bool inR) {
  static_assert(isNonBasic && Options::is_of_boundary_type &&
                    (Options::has_vine_update || Options::can_retrieve_representative_cycles) &&
                    Options::is_indexed_by_position,
                "Only enabled for RU matrices.");

  return matrix_.zero_column(columnIndex, inR);
}

template <class Options>
inline bool Matrix<Options>::is_zero_cell(index columnIndex, index rowIndex) {
  return matrix_.is_zero_cell(columnIndex, rowIndex);
}

template <class Options>
inline bool Matrix<Options>::is_zero_cell(index columnIndex, index rowIndex, bool inR) const {
  static_assert(isNonBasic && Options::is_of_boundary_type &&
                    (Options::has_vine_update || Options::can_retrieve_representative_cycles) &&
                    Options::is_indexed_by_position,
                "Only enabled for RU matrices.");

  return matrix_.is_zero_cell(columnIndex, rowIndex, inR);
}

template <class Options>
inline bool Matrix<Options>::is_zero_column(index columnIndex) {
  return matrix_.is_zero_column(columnIndex);
}

template <class Options>
inline bool Matrix<Options>::is_zero_column(index columnIndex, bool inR) {
  static_assert(isNonBasic && Options::is_of_boundary_type &&
                    (Options::has_vine_update || Options::can_retrieve_representative_cycles) &&
                    Options::is_indexed_by_position,
                "Only enabled for RU matrices.");

  return matrix_.is_zero_column(columnIndex, inR);
}

template <class Options>
inline typename Matrix<Options>::index Matrix<Options>::get_column_with_pivot(index simplexIndex) const {
  static_assert(isNonBasic && (!Options::is_of_boundary_type ||
                               (Options::has_vine_update || Options::can_retrieve_representative_cycles)),
                "'get_column_with_pivot' is not available for the chosen options.");

  return matrix_.get_column_with_pivot(simplexIndex);
}

template <class Options>
inline int Matrix<Options>::get_pivot(index columnIndex) {
  static_assert(isNonBasic, "'get_pivot' is not available for the chosen options.");

  return matrix_.get_pivot(columnIndex);
}

template <class Options>
inline Matrix<Options>& Matrix<Options>::operator=(Matrix other) {
  swap(matrix_, other.matrix_);

  return *this;
}

template <class Options>
inline void Matrix<Options>::print() {
  return matrix_.print();
}

template <class Options>
inline const typename Matrix<Options>::barcode_type& Matrix<Options>::get_current_barcode() {
  static_assert(Options::has_column_pairings, "This method was not enabled.");

  return matrix_.get_current_barcode();
}

template <class Options>
inline const typename Matrix<Options>::barcode_type& Matrix<Options>::get_current_barcode() const {
  static_assert(Options::has_column_pairings, "This method was not enabled.");
  static_assert(
      !Options::is_of_boundary_type || Options::has_vine_update || Options::can_retrieve_representative_cycles,
      "'get_current_barcode' is not const for boundary matrices as the barcode is only computed when explicitely "
      "asked.");

  return matrix_.get_current_barcode();
}

template <class Options>
inline void Matrix<Options>::swap_columns(index columnIndex1, index columnIndex2) {
  static_assert((!isNonBasic && !Options::has_column_compression) ||
                    (isNonBasic && Options::is_of_boundary_type && Options::is_indexed_by_position &&
                     !Options::has_vine_update && !Options::can_retrieve_representative_cycles),
                "This method was not enabled.");
  return matrix_.swap_columns(columnIndex1, columnIndex2);
}

template <class Options>
inline void Matrix<Options>::swap_rows(index rowIndex1, index rowIndex2) {
  static_assert((!isNonBasic && !Options::has_column_compression) ||
                    (isNonBasic && Options::is_of_boundary_type && Options::is_indexed_by_position &&
                     !Options::has_vine_update && !Options::can_retrieve_representative_cycles),
                "This method was not enabled.");
  return matrix_.swap_rows(rowIndex1, rowIndex2);
}

template <class Options>
inline void Matrix<Options>::swap_at_indices(index index1, index index2) {
  static_assert((!isNonBasic && !Options::has_column_compression) ||
                    (isNonBasic && Options::is_of_boundary_type && !Options::has_vine_update &&
                     !Options::can_retrieve_representative_cycles),
                "This method was not enabled.");
  return matrix_.swap_at_indices(index1, index2);
}

template <class Options>
inline bool Matrix<Options>::vine_swap_with_z_eq_1_case(index index) {
  static_assert(Options::has_vine_update && Options::is_indexed_by_position, "This method was not enabled.");
  return matrix_.vine_swap_with_z_eq_1_case(index);
}

template <class Options>
inline typename Matrix<Options>::index Matrix<Options>::vine_swap_with_z_eq_1_case(index columnIndex1,
                                                                                   index columnIndex2) {
  static_assert(Options::has_vine_update && !Options::is_indexed_by_position, "This method was not enabled.");
  return matrix_.vine_swap_with_z_eq_1_case(columnIndex1, columnIndex2);
}

template <class Options>
inline bool Matrix<Options>::vine_swap(index index) {
  static_assert(Options::has_vine_update && Options::is_indexed_by_position, "This method was not enabled.");
  return matrix_.vine_swap(index);
}

template <class Options>
inline typename Matrix<Options>::index Matrix<Options>::vine_swap(index columnIndex1, index columnIndex2) {
  static_assert(Options::has_vine_update && !Options::is_indexed_by_position, "This method was not enabled.");
  return matrix_.vine_swap(columnIndex1, columnIndex2);
}

template <class Options>
inline void Matrix<Options>::update_representative_cycles() {
  static_assert(Options::can_retrieve_representative_cycles, "This method was not enabled.");
  matrix_.update_representative_cycles();
}

template <class Options>
inline const std::vector<typename Matrix<Options>::cycle_type>& Matrix<Options>::get_representative_cycles() {
  static_assert(Options::can_retrieve_representative_cycles, "This method was not enabled.");
  return matrix_.get_representative_cycles();
}

template <class Options>
inline const typename Matrix<Options>::cycle_type& Matrix<Options>::get_representative_cycle(const Bar& bar) {
  static_assert(Options::can_retrieve_representative_cycles, "This method was not enabled.");
  return matrix_.get_representative_cycle(bar);
}

template <class Options>
inline constexpr void Matrix<Options>::_assert_options() {
  static_assert(Options::column_type != Column_types::HEAP || !Options::has_row_access,
                "Row access is not possible for heap columns.");
  static_assert(!Options::has_vine_update || Options::is_z2, "Vine update currently works only for Z_2 coefficients.");
  static_assert(Options::column_type != Column_types::HEAP || !Options::has_column_compression,
                "Column compression not compatible with heap columns.");

//   // This should be warnings instead, as Options::has_column_compression is just ignored in those cases and don't
//   // produces errors.
//   static_assert(!Options::has_column_compression || !Options::has_column_pairings,
//                 "Column compression not available to compute persistence homology (it would bring no advantages; "
//                 "use it for co-homology instead).");
//   static_assert(!Options::has_column_compression || !Options::has_vine_update,
//                 "Column compression not available for vineyards.");
//   static_assert(!Options::has_column_compression || !Options::can_retrieve_representative_cycles,
//                 "Column compression not available to retrieve representative cycles.");
//   // Would column removal while column compression be useful? If yes, should erase() remove a single column or the
//   // class of columns identical to the input?
//   // For a single column, I have an implementation for union-find (not the current one) which allows deleting a
//   // single element in constant time, but find becomes log n in worst case.
//   // For a column class, we either just ignore the removed class (constant time), but this results in memory
//   // residues, or, we have an erase method which is at least linear in the size of the class.
//   static_assert(!Options::has_column_compression || !Options::has_removable_columns,
//                 "When column compression is used, the removal of columns is not implemented yet.");
}

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // MASTER_MATRIX_H

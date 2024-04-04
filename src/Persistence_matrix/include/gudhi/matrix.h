/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022-24 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/** @file matrix.h
 * @author Hannah Schreiber
 * @brief Contains @ref Gudhi::persistence_matrix::Matrix class.
 */

#ifndef MASTER_MATRIX_H
#define MASTER_MATRIX_H

#include <type_traits>
#include <vector>
#include <unordered_map>
#include <map>
#include <assert.h>
#include <initializer_list>

#include <boost/intrusive/list.hpp>

#include <gudhi/persistence_matrix_options.h>

#include <gudhi/Persistence_matrix/overlay_ididx_to_matidx.h>
#include <gudhi/Persistence_matrix/overlay_posidx_to_matidx.h>

#include <gudhi/Persistence_matrix/matrix_dimension_holders.h>
#include <gudhi/Persistence_matrix/matrix_row_access.h>
#include <gudhi/Persistence_matrix/base_swap.h>
#include <gudhi/Persistence_matrix/base_pairing.h>
#include <gudhi/Persistence_matrix/ru_pairing.h>
#include <gudhi/Persistence_matrix/ru_vine_swap.h>
#include <gudhi/Persistence_matrix/ru_rep_cycles.h>
#include <gudhi/Persistence_matrix/chain_pairing.h>
#include <gudhi/Persistence_matrix/chain_vine_swap.h>
#include <gudhi/Persistence_matrix/chain_rep_cycles.h>

#include <gudhi/Persistence_matrix/base_matrix.h>
#include <gudhi/Persistence_matrix/base_matrix_with_column_compression.h>
#include <gudhi/Persistence_matrix/boundary_matrix.h>
#include <gudhi/Persistence_matrix/ru_matrix.h>
#include <gudhi/Persistence_matrix/chain_matrix.h>

#include <gudhi/Persistence_matrix/columns/cell_constructors.h>
#include <gudhi/Persistence_matrix/columns/cell_types.h>
#include <gudhi/Persistence_matrix/columns/row_access.h>

#include <gudhi/Persistence_matrix/columns/column_dimension_holder.h>
#include <gudhi/Persistence_matrix/columns/chain_column_extra_properties.h>
#include <gudhi/Persistence_matrix/columns/intrusive_list_column.h>
#include <gudhi/Persistence_matrix/columns/intrusive_set_column.h>
#include <gudhi/Persistence_matrix/columns/list_column.h>
#include <gudhi/Persistence_matrix/columns/set_column.h>
#include <gudhi/Persistence_matrix/columns/unordered_set_column.h>
#include <gudhi/Persistence_matrix/columns/vector_column.h>
#include <gudhi/Persistence_matrix/columns/naive_vector_column.h>
#include <gudhi/Persistence_matrix/columns/heap_column.h>

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
 * - a boundary matrix @f$ B = R \cdot U @f$ which either stores only @f$ R @f$ or the whole decomposition @f$ R @f$ 
 *   and @f$ U @f$ depending on the options. This type is selected if @ref is_of_boundary_type is set to true and 
 *   at least one of the following options: @ref has_column_pairings, @ref has_vine_update and 
 *   @ref can_retrieve_representative_cycles. If only @ref has_column_pairings is true, then only @f$ R @f$ is stored,
 *   but if either @ref has_vine_update or @ref can_retrieve_representative_cycles is true, then @f$ U @f$ also needs 
 *   to be stored. Note that the option @ref column_indexation_type will produce a small overhead when set to @ref Column_indexation_types::IDENTIFIER.
 * - a chain complex matrix representing a `compatible base` of a filtered chain complex (see TODO: cite Clément's zigzag paper here).
 *   This matrix is deduced from the boundary matrix and therefore encodes more or less the same information 
 *   but differently and can therefore be better suited for certain applications. This type can be used the same way 
 *   than the precedent type, only the option @ref is_of_boundary_type has to be set to false. Note that the option 
 *   @ref column_indexation_type will produce a small overhead when set to  @ref Column_indexation_types::POSITION or  @ref Column_indexation_types::IDENTIFIER.
 *
 * __Indexation scheme:__
 *
 * The indexation system for columns of the different matrix types can be a bit tricky and different methods will not always take
 * the same type of index as input (for optimization purposes). So, to avoid confusion, we will name and define here the 
 * different possibilities, such that we can directly refer to it in the descriptions of the methods. 
 * Note that every column and row in a boundary or chain matrix is always associated to a single simplex/face, 
 * so in order to avoid repeating formulations like "of the simplex associated to the column" all the time, 
 * we will amalgamate both notions together.
 *
 * Let c be a column.
 * - MatIdx: This will correspond to the position of c in the matrix, i.e., underlying_container[MatIdx] = c.
 *   This will be the only public indexing scheme for basic matrices (first of the list above).
 * - PosIdx: This will correspond to the relative position of c in the current filtration compared to the other columns,
 *   starting the count at 0. For boundary matrices, PosIdx will always be equal to MatIdx, but this is not true for
 *   chain matrices, where MatIdx will correspond to the relative position of c in the original filtration
 *   (so PosIdx == MatIdx as long as no vine swaps are done).
 * - IDIdx: This will correspond to the ID of c in the complex used to identify it in the boundaries.
 *   If at the insertion of c, its ID was not specified and it was the nth insertion, it is assumed that the ID is n
 *   (which means that IDIdx and PosIdx will only start to differ when swaps or removals are performed).
 *   If an ID is specified at the insertion of c, the ID is stored as the IDIdx of c. IDs can be freely choosed with 
 *   the only restriction that they have to be strictly increasing in the order of the filtration.
 *
 * In conclusion, with default values, if no vine swaps or removals occurs, all three indexing schemes are the same.
 *
 * Let r be a row. Rows are indexed in two ways depending only if the matrix is a chain matrix or not.
 * If the matrix is a chain matrix, r is always indexed by its ID, so it correspond to the IDIdx indexing scheme. 
 * If the matrix is not a chain matrix, r will originaly also be indexed by the ID, but when a swap occurs,
 * the rows also swap IDs and the new ID has to be used to access r. This means that when the default IDIdx scheme 
 * is used, the indexation of the rows correspond to PosIdx.
 * 
 * @tparam Options Structure encoding all the options of the matrix. 
 * See description of @ref Default_options for more details.
 */
template <class PersistenceMatrixOptions = Default_options<> >
class Matrix {
 public:
  using Option_list = PersistenceMatrixOptions; //to make it accessible from the other classes
  using index = typename PersistenceMatrixOptions::index_type;                 /**< Type of MatIdx index. */
  using id_index = typename PersistenceMatrixOptions::index_type;              /**< Type of IDIdx index. */
  using pos_index = typename PersistenceMatrixOptions::index_type;             /**< Type of PosIdx index. */
  using dimension_type = typename PersistenceMatrixOptions::dimension_type;    /**< Type for dimension value. */

  struct Dummy_field_operators{
    using element_type = unsigned int;
    using characteristic_type = element_type;

    Dummy_field_operators([[maybe_unused]] characteristic_type characteristic = 0){}

    friend void swap([[maybe_unused]] Dummy_field_operators& d1, [[maybe_unused]] Dummy_field_operators& d2){}

    static constexpr characteristic_type get_characteristic() { return 2; }
  };

  using Field_operators =
      typename std::conditional<PersistenceMatrixOptions::is_z2, 
                                Dummy_field_operators, 
                                typename PersistenceMatrixOptions::Field_coeff_operators
                               >::type;                       /**< Coefficiants field type. */
  using element_type = typename std::conditional<PersistenceMatrixOptions::is_z2, 
                                                 bool, 
                                                 typename Field_operators::element_type
                                                >::type;
  using characteristic_type =
      typename std::conditional<PersistenceMatrixOptions::is_z2, 
                                unsigned int, 
                                typename Field_operators::characteristic_type
                               >::type;

  // TODO: move outside? unify with other bar types in Gudhi?
  /**
   * @brief Type for a bar in the computed barcode. Stores the birth, death and dimension of the bar.
   */
  struct Bar {
    Bar() : dim(-1), birth(-1), death(-1) {}

    Bar(dimension_type dim, pos_index birth, pos_index death) : dim(dim), birth(birth), death(death) {}

    dimension_type dim; /**< Dimension of the bar.*/
    pos_index birth;    /**< Birth index in the current filtration. */
    pos_index death;    /**< Death index in the current filtration. */

    inline friend std::ostream &operator<<(std::ostream &stream, const Bar& bar) {
      stream << "[" << bar.dim << "] ";
      stream << bar.birth << ", " << bar.death;
      stream << "\n";
      return stream;
    }
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

  using row_hook_type = typename std::conditional<PersistenceMatrixOptions::has_row_access &&
                                                      PersistenceMatrixOptions::has_intrusive_rows,
                                                  base_hook_matrix_row, 
                                                  Dummy_row_hook
                                                 >::type;
  using column_hook_type = typename std::conditional<
      PersistenceMatrixOptions::column_type == Column_types::INTRUSIVE_LIST, 
      base_hook_matrix_list_column,
      typename std::conditional<PersistenceMatrixOptions::column_type == Column_types::INTRUSIVE_SET,
                                base_hook_matrix_set_column, 
                                Dummy_column_hook
                               >::type
    >::type;

  //Option to store the column index within the cell (additionnaly to the row index). Necessary only with row access.
  using Cell_column_index_option =
      typename std::conditional<PersistenceMatrixOptions::has_row_access,
                                Cell_column_index<index>,
                                Dummy_cell_column_index_mixin
                               >::type;
  //Option to store the value of the cell. 
  //Unnecessary for values in Z_2 as there are always 1 (0-valued cells are never stored).
  using Cell_field_element_option =
      typename std::conditional<PersistenceMatrixOptions::is_z2,
                                Dummy_cell_field_element_mixin,
                                Cell_field_element<element_type>
                               >::type;
  /**
   * @brief Type of a matrix cell. See @ref Cell for a more detailed description.
   */
  using Cell_type = Cell<Matrix<PersistenceMatrixOptions> >;

  /**
   * @brief Default cell constructor/destructor, using classic new and delete.
   * For now, only used as default value for columns constructed independently outside of the matrix by the user.
   * Could be used in the futur when parallel options are implemented, as usual pools are not thread safe.
   */
  inline static New_cell_constructor<Cell_type> defaultCellConstructor;
  /**
   * @brief Cell constructor/destructor used by the matrix. Uses a pool of cells to accelerate memory management,
   * as cells are constructed and destroyed a lot during reduction, swaps or additions.
   */
  using Cell_constructor = Pool_cell_constructor<Cell_type>;

  /**
   * @brief Type used to identify a cell, for exemple when inserting a boundary. If @ref is_z2 is true, the type is
   * an IDIdx and corresponds to the row index of the cell (the cell value is assumed to be 1). 
   * If @ref is_z2 is false, the type is a pair whose first element is the row index of the cell and 
   * the second element is the value of the cell (which again is assumed to be non-zero). The column index of the row
   * is always deduced from the context in which the type is used.
   */
  using cell_rep_type = typename std::conditional<PersistenceMatrixOptions::is_z2,
                                                  id_index,
                                                  std::pair<id_index, element_type>
                                                 >::type;

  /**
   * @brief Compaires two pairs, representing a cell (first = row index, second = value), 
   * by their position in the column and not their values. 
   * The two represented cells are therefore assumed to be in the same column.
   */
  struct CellPairComparator {
    bool operator()(const std::pair<id_index, element_type>& p1, const std::pair<id_index, element_type>& p2) const {
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

  /**
   * @brief Type of the rows storted in the matrix. Is either an intrusive list of @ref Cell_type (not ordered) if 
   * @ref has_intrusive_rows is true, or a set of @ref Cell_type (ordered by @ref get_column_index) otherwise.
   */
  using Row_type =
      typename std::conditional<PersistenceMatrixOptions::has_intrusive_rows,
                                boost::intrusive::list<Cell_type, 
                                                       boost::intrusive::constant_time_size<false>,
                                                       boost::intrusive::base_hook<base_hook_matrix_row>
                                                      >,
                                std::set<Cell_type, RowCellComp>
                               >::type;

  using row_container_type =
      typename std::conditional<PersistenceMatrixOptions::has_removable_rows,
                                std::map<id_index, Row_type>,
                                std::vector<Row_type>
                               >::type;

  //Row access at column level
  using Row_access_option =
      typename std::conditional<PersistenceMatrixOptions::has_row_access,
                                Row_access<Matrix<PersistenceMatrixOptions> >,
                                Dummy_row_access
                               >::type;
  //Row access at matrix level
  using Matrix_row_access_option =
      typename std::conditional<PersistenceMatrixOptions::has_row_access,
                                Matrix_row_access<Row_type, 
                                                  row_container_type, 
                                                  PersistenceMatrixOptions::has_removable_rows, 
                                                  id_index>,
                                Dummy_matrix_row_access
                               >::type;

  template <typename value_type>
  using dictionnary_type =
      typename std::conditional<PersistenceMatrixOptions::has_map_column_container,
                                std::unordered_map<unsigned int, value_type>,
                                std::vector<value_type>
                               >::type;

  static const bool isNonBasic = PersistenceMatrixOptions::has_column_pairings ||
                                 PersistenceMatrixOptions::has_vine_update ||
                                 PersistenceMatrixOptions::can_retrieve_representative_cycles;

  using Column_dimension_option =
      typename std::conditional<isNonBasic,
                                Column_dimension_holder<Matrix<PersistenceMatrixOptions> >,
                                Dummy_dimension_holder
                               >::type;
  //Extra information needed for a column when the matrix is a chain matrix. 
  using Chain_column_option =
      typename std::conditional<isNonBasic && !PersistenceMatrixOptions::is_of_boundary_type,
                                Chain_column_extra_properties<Matrix<PersistenceMatrixOptions> >,
                                Dummy_chain_properties
                               >::type;

  using Heap_column_type = Heap_column<Matrix<PersistenceMatrixOptions>, Cell_constructor>;
  using List_column_type = List_column<Matrix<PersistenceMatrixOptions>, Cell_constructor>;
  using Vector_column_type = Vector_column<Matrix<PersistenceMatrixOptions>, Cell_constructor>;
  using Naive_vector_column_type = Naive_vector_column<Matrix<PersistenceMatrixOptions>, Cell_constructor>;
  using Set_column_type = Set_column<Matrix<PersistenceMatrixOptions>, Cell_constructor>;
  using Unordered_set_column_type = Unordered_set_column<Matrix<PersistenceMatrixOptions>, Cell_constructor>;
  using Intrusive_list_column_type = Intrusive_list_column<Matrix<PersistenceMatrixOptions>, Cell_constructor>;
  using Intrusive_set_column_type = Intrusive_set_column<Matrix<PersistenceMatrixOptions>, Cell_constructor>;

  /**
   * @brief Type of the columns stored in the matrix. The type depends on the value of @ref column_type defined
   * in the given options. See @ref Column_types for a more detailed description. All columns follow the
   * @ref PersistenceMatrixColumn concept.
   */
  using Column_type = typename std::conditional<
        PersistenceMatrixOptions::column_type == Column_types::HEAP, 
        Heap_column_type,
        typename std::conditional<
            PersistenceMatrixOptions::column_type == Column_types::LIST, 
            List_column_type,
            typename std::conditional<
                PersistenceMatrixOptions::column_type == Column_types::SET, 
                Set_column_type,
                typename std::conditional<
                    PersistenceMatrixOptions::column_type == Column_types::UNORDERED_SET, 
                    Unordered_set_column_type,
                    typename std::conditional<
                        PersistenceMatrixOptions::column_type == Column_types::VECTOR, 
                        Vector_column_type,
                        typename std::conditional<
                            PersistenceMatrixOptions::column_type == Column_types::INTRUSIVE_LIST, 
                            Intrusive_list_column_type,
                            typename std::conditional<
                                PersistenceMatrixOptions::column_type == Column_types::NAIVE_VECTOR,
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
      typename std::conditional<PersistenceMatrixOptions::has_map_column_container, 
                                std::unordered_map<index, Column_type>,
                                std::vector<Column_type>
                               >::type;

  static const bool hasFixedBarcode = Option_list::is_of_boundary_type && !PersistenceMatrixOptions::has_vine_update;
  /**
   * @brief Type of the computed barcode. If @ref has_map_column_container is true, it is a list of @ref Bar,
   * otherwise it is a vector of @ref Bar.
   */
  using barcode_type =
      typename std::conditional<hasFixedBarcode, 
                                std::vector<Bar>, 
                                typename std::conditional<PersistenceMatrixOptions::has_removable_columns, 
                                  std::list<Bar>, 
                                  std::vector<Bar>
                                >::type
                               >::type;
  using bar_dictionnary_type = 
      typename std::conditional<hasFixedBarcode,
                                typename std::conditional<PersistenceMatrixOptions::can_retrieve_representative_cycles,
                                  std::vector<index>,                   //RU
                                  std::unordered_map<pos_index, index>  //boundary
                                >::type,
                                typename std::conditional<PersistenceMatrixOptions::has_removable_columns,
                                  std::unordered_map<pos_index, typename barcode_type::iterator>,
                                  std::vector<index>
                                >::type
                               >::type;

  //default type for boundaries to permit list initialization directly in function parameters
  using boundary_type = typename std::conditional<PersistenceMatrixOptions::is_z2, 
                                                  std::initializer_list<id_index>,
                                                  std::initializer_list<std::pair<id_index, element_type> >
                                                 >::type;

  //i.e. is simple boundary matrix. Also, only needed because of the reduction algorithm. 
  //TODO: remove the necessity and recalculate when needed or keep it like that?
  static const bool maxDimensionIsNeeded =
      PersistenceMatrixOptions::has_column_pairings && PersistenceMatrixOptions::is_of_boundary_type &&
      !PersistenceMatrixOptions::has_vine_update && !PersistenceMatrixOptions::can_retrieve_representative_cycles;

  using Matrix_dimension_option = typename std::conditional<
      PersistenceMatrixOptions::has_matrix_maximal_dimension_access || maxDimensionIsNeeded,
      typename std::conditional<PersistenceMatrixOptions::has_removable_columns, 
                                Matrix_all_dimension_holder<dimension_type>,
                                Matrix_max_dimension_holder<dimension_type>
                               >::type,
      Dummy_matrix_dimension_holder
    >::type;

  using Base_matrix_type =
      typename std::conditional<PersistenceMatrixOptions::has_column_compression, 
                                Base_matrix_with_column_compression<Matrix<PersistenceMatrixOptions> >,
                                Base_matrix<Matrix<PersistenceMatrixOptions> >
                               >::type;
  using Boundary_matrix_type = Boundary_matrix<Matrix<PersistenceMatrixOptions> >;
  using RU_matrix_type = RU_matrix<Matrix<PersistenceMatrixOptions> >;
  using Chain_matrix_type = Chain_matrix<Matrix<PersistenceMatrixOptions> >;

  template <class Base>
  using Base_swap_option =
      typename std::conditional<PersistenceMatrixOptions::has_vine_update ||
                                    PersistenceMatrixOptions::has_column_and_row_swaps,
                                Base_swap<Matrix<PersistenceMatrixOptions>, Base>, 
                                Dummy_base_swap
                               >::type;
  using Base_pairing_option =
      typename std::conditional<PersistenceMatrixOptions::has_column_pairings &&
                                    !PersistenceMatrixOptions::has_vine_update &&
                                    !PersistenceMatrixOptions::can_retrieve_representative_cycles,
                                Base_pairing<Matrix<PersistenceMatrixOptions> >, 
                                Dummy_base_pairing
                               >::type;

  using RU_pairing_option =
      typename std::conditional<PersistenceMatrixOptions::has_column_pairings &&
                                    !PersistenceMatrixOptions::has_vine_update,
                                RU_pairing<Matrix<PersistenceMatrixOptions> >, 
                                Dummy_ru_pairing
                               >::type;
  using RU_vine_swap_option =
      typename std::conditional<PersistenceMatrixOptions::has_vine_update,
                                RU_vine_swap<Matrix<PersistenceMatrixOptions> >,
                                Dummy_ru_vine_swap
                               >::type;
  using RU_representative_cycles_option =
      typename std::conditional<PersistenceMatrixOptions::can_retrieve_representative_cycles, 
                                RU_representative_cycles<Matrix<PersistenceMatrixOptions> >,
                                Dummy_ru_representative_cycles
                               >::type;

  using Chain_pairing_option =
      typename std::conditional<PersistenceMatrixOptions::has_column_pairings &&
                                    !PersistenceMatrixOptions::has_vine_update,
                                Chain_pairing<Matrix<PersistenceMatrixOptions> >, 
                                Dummy_chain_pairing
                               >::type;
  using Chain_vine_swap_option = typename std::conditional<PersistenceMatrixOptions::has_vine_update, 
                                                           Chain_vine_swap<Matrix<PersistenceMatrixOptions> >,
                                                           Dummy_chain_vine_swap
                                                          >::type;
  using Chain_representative_cycles_option =
      typename std::conditional<PersistenceMatrixOptions::can_retrieve_representative_cycles,
                                Chain_representative_cycles<Matrix<PersistenceMatrixOptions> >,
                                Dummy_chain_representative_cycles
                               >::type;

  /**
   * @brief Type of a representative cycle. Vector of row indices, see [TODO: index paragraph].
   */
  using cycle_type = std::vector<id_index>; //TODO: add coefficients

  //Return types to factorize the corresponding methods

  //The returned column is `const` if the matrix uses column compression
  using returned_column_type =
      typename std::conditional<!isNonBasic && PersistenceMatrixOptions::has_column_compression,
                                const Column_type,
                                Column_type
                               >::type;
  //The returned row is `const` if the matrix uses column compression
  using returned_row_type =
      typename std::conditional<!isNonBasic && PersistenceMatrixOptions::has_column_compression,
                                const Row_type,
                                Row_type
                               >::type;
  //If the matrix is a chain matrix, the insertion method returns the pivots of its unpaired columns used to reduce
  //the inserted boundary. Otherwise, void.
  using insertion_return_type =
      typename std::conditional<PersistenceMatrixOptions::is_of_boundary_type || !isNonBasic ||
                                    PersistenceMatrixOptions::column_indexation_type ==
                                        Column_indexation_types::POSITION,
                                void, 
                                std::vector<cell_rep_type>
                               >::type;

  /**
   * @brief Default constructor. Initializes an empty matrix.
   */
  Matrix();
  /**
   * @brief Constructs a new matrix from the given ranges of @ref cell_rep_type. Each range corresponds to a column 
   * (the order of the ranges are preserved). The content of the ranges is assumed to be sorted by increasing IDs.
   * If the columns are representing a boundary matrix, the IDs of the simplices are also assumed to be 
   * consecutifs, ordered by filtration value, starting with 0. 
   *
   * See options descriptions for futher details on how the given matrix is handled.
   * 
   * @tparam Container_type Range type for @ref cell_rep_type ranges. Assumed to have a begin(), end() and size() method.
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
   * @param characteristic Characteristic of the coefficient field. Has to be specified if @ref is_z2 is false.
   * Default value is 11. Ignored if @ref is_z2 is true.
   */
  template <class Container_type = boundary_type>
  Matrix(const std::vector<Container_type>& columns, characteristic_type characteristic = 11);
  /**
   * @brief Constructs a new empty matrix and reserves space for the given number of columns.
   * 
   * @param numberOfColumns Number of columns to reserve space for.
   * @param characteristic Characteristic of the coefficient field. If not specified and @ref is_z2 is false,
   * the characteristic has to be set later with the use of `set_characteristic()` before calling for the first time
   * a method needing it. Ignored if @ref is_z2 is true.
   */
  Matrix(int numberOfColumns, characteristic_type characteristic = 0);
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
   * @tparam EventComparatorFunction Type of the birth or death comparator: (@ref pos_index, @ref pos_index) -> bool
   * @param birthComparator Method taking two PosIdx indices as parameter and returns true if and only if the first 
   * face is associated to a bar with strictly smaller birth than the bar associated to the second one.
   * @param deathComparator Method taking two PosIdx indices as parameter and returns true if and only if the first 
   * face is associated to a bar with strictly smaller death than the bar associated to the second one.
   */
  template <typename EventComparatorFunction>
  Matrix(EventComparatorFunction&& birthComparator, 
         EventComparatorFunction&& deathComparator);
  /**
   * @brief Constructs a new matrix from the given ranges with the given comparator functions. 
   * Only available when those comparators are necessary, i.e., when **all** following options have following values:
   *   - @ref is_of_boundary_type = false
   *   - @ref has_vine_update = true
   *   - @ref has_column_pairings = false
   *
   * See description of @ref Matrix(const std::vector<Container_type>& columns) for more information about 
   * @p orderedBoundaries and 
   * @ref Matrix(EventComparatorFunction&& birthComparator, EventComparatorFunction&& deathComparator) 
   * for more information about the comparators.
   * 
   * @tparam EventComparatorFunction Type of the birth or death comparator: (@ref pos_index, @ref pos_index) -> bool
   * @tparam Boundary_type Range type for a column. Assumed to have a begin(), end() and size() method.
   * @param orderedBoundaries Vector of ordered boundaries in filtration order. Indexed continously starting at 0.
   * @param birthComparator Method taking two PosIdx indices as parameter and returns true if and only if the first 
   * face is associated to a bar with strictly smaller birth than the bar associated to the second one.
   * @param deathComparator Method taking two PosIdx indices as parameter and returns true if and only if the first 
   * face is associated to a bar with strictly smaller death than the bar associated to the second one.
   * @param characteristic Characteristic of the coefficient field. Has to be specified if @ref is_z2 is false.
   * Default value is 11. Ignored if @ref is_z2 is true.
   */
  template <typename EventComparatorFunction, class Boundary_type = boundary_type>
  Matrix(const std::vector<Boundary_type>& orderedBoundaries, 
         EventComparatorFunction&& birthComparator,
         EventComparatorFunction&& deathComparator, 
         characteristic_type characteristic = 11);
  /**
   * @brief Constructs a new empty matrix and reserves space for the given number of columns.
   * Only available when those comparators are necessary, i.e., when **all** following options have following values:
   *   - @ref is_of_boundary_type = false
   *   - @ref has_vine_update = true
   *   - @ref has_column_pairings = false
   *
   * See description of 
   * @ref Matrix(EventComparatorFunction&& birthComparator, EventComparatorFunction&& deathComparator) 
   * for more information about the comparators.
   * 
   * @tparam EventComparatorFunction Type of the birth or death comparator: (@ref pos_index, @ref pos_index) -> bool
   * @param numberOfColumns Number of columns to reserve space for.
   * @param birthComparator Method taking two PosIdx indices as parameter and returns true if and only if the first 
   * face is associated to a bar with strictly smaller birth than the bar associated to the second one.
   * @param deathComparator Method taking two PosIdx indices as parameter and returns true if and only if the first 
   * face is associated to a bar with strictly smaller death than the bar associated to the second one.
   * @param characteristic Characteristic of the coefficient field. If not specified and @ref is_z2 is false,
   * the characteristic has to be set later with the use of `set_characteristic()` before calling for the first time
   * a method needing it. Ignored if @ref is_z2 is true.
   */
  template <typename EventComparatorFunction>
  Matrix(unsigned int numberOfColumns, 
         EventComparatorFunction&& birthComparator,
         EventComparatorFunction&& deathComparator, 
         characteristic_type characteristic = 0);
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

  ~Matrix();

  //TODO: compatibily with multi fields:
  //  - set_characteristic(characteristic_type min, characteristic_type max)
  //  - readapt reduction?
  /**
   * @brief Sets the characteristic of the coefficient field if @ref is_z2 is false, does nothing otherwise.
   * Should be used if no characteristic could be specified at the creation of the empty matrix.
   * Do not change the value of the characteristic once used.
   *
   * @warning The coefficient values stored in the matrix are stored after computing the corresponding modulo.
   * Therefore, changing the characteristic after is very likely to invalidate all cell values.
   * 
   * @param characteristic The characteristic to set.
   */
  void set_characteristic(characteristic_type characteristic);

  // (TODO: if there is no row access and the column type corresponds to the internal column type of the matrix, 
  // moving the column instead of copying it should be possible. Is it worth implementing it?)
  /**
   * @brief Inserts a new ordered column at the end of the matrix by copying the given range of @ref cell_rep_type.
   * The content of the range is assumed to be sorted by increasing ID value. 
   *
   * Only available for base matrices, see [TODO: description].
   * Otherwise use @ref insert_boundary which will deduce a new column from the boundary given.
   * 
   * @tparam Container_type Range of @ref cell_rep_type. Assumed to have a begin(), end() and size() method.
   * @param column Column to be inserted.
   */
  template <class Container_type>
  void insert_column(const Container_type& column);
  /**
   * @brief Inserts a new ordered column at the given index by copying the given range of @ref cell_rep_type.
   * There should not be any other column inserted at that index which was not explicitely removed before.
   * The content of the range is assumed to be sorted by increasing ID value. 
   *
   * Only available for base matrices without column compression and without row access, see [TODO: description and Options].
   * 
   * @tparam Container_type Range of @ref cell_rep_type. Assumed to have a begin(), end() and size() method.
   * @param column Column to be inserted.
   * @param columnIndex MatIdx index to which the column has to be inserted.
   */
  template <class Container_type>
  void insert_column(const Container_type& column, index columnIndex);
  //TODO: for simple boundary matrices, add an index pointing to the first column inserted after the last call of 
  //get_current_barcode to enable several calls to get_current_barcode
  /**
   * @brief Inserts at the end of the matrix a new ordered column corresponding to the given boundary. 
   * This means that it is assumed that this method is called on boundaries in the order of the filtration. 
   * It also assumes that the faces in the given boundary are identified by their relative position in the filtration, 
   * starting at 0. If it is not the case, use the other `insert_boundary` instead by indicating the face ID
   * used in the boundaries when the face is inserted.
   *
   * Different to the constructor, the boundaries do not have to come from a simplicial complex, but also from
   * a more general cell complex. This includes cubical complexes or Morse complexes for example.
   *
   * The content of the new column will vary depending on the underlying type of the matrix (see [TODO: description]):
   * - If it is a basic matrix type, the boundary is copied as it is, i.e., the method is equivalent to 
   *   @ref insert_column.
   * - If it is a boundary type matrix and only R is stored, the boundary is also just copied. The column will only be 
   *   reduced later when the barcode is requested in order to apply some optimisations with the additional knowledge.
   *   Hence, the barcode will also not be updated, so call @ref get_current_barcode only when the matrix is complete.
   * - If it is a boundary type matrix and both R and U are stored, the new boundary is stored in its reduced form and
   *   the barcode, if active, is also updated.
   * - If it is a chain type matrix, the new column is of the form 
   *   `IDIdx + linear combination of older column IDIdxs`, where the combination is deduced while reducing the 
   *   given boundary. If the barcode is stored, it will also be updated.
   * 
   * @tparam Boundary_type Range of @ref cell_rep_type. Assumed to have a begin(), end() and size() method.
   * @param boundary Boundary generating the new column. The content should be ordered by ID.
   * @param dim Dimension of the face whose boundary is given. If the complex is simplicial, 
   * this parameter can be omitted as it can be deduced from the size of the boundary.
   * @return If it is a chain matrix, the method returns the MatIdx indices of the unpaired chains used to reduce
   * the boundary. Otherwise, nothing.
   */
  template <class Boundary_type = boundary_type>
  insertion_return_type insert_boundary(const Boundary_type& boundary, dimension_type dim = -1);
  /**
   * @brief Only avalaible for non basic matrices.
   * It does the same as the other version, but permits the boundary faces to be identified without restrictions
   * except that all IDs have to be strictly increasing in the order of filtration. Note that you should avoid then
   * to use the other insertion method to avoid overwriting IDs.
   *
   * As a face has to be inserted before one of its cofaces in a valid filtration (recall that it is assumed that
   * for non basic matrices, the faces are inserted by order of filtration), it is sufficient to indicate the ID
   * of the face being inserted.
   * 
   * @tparam Boundary_type Range of @ref cell_rep_type. Assumed to have a begin(), end() and size() method.
   * @param faceIndex IDIdx index to use to indentify the new face.
   * @param boundary Boundary generating the new column. The indices of the boundary have to correspond to the 
   * @p faceIndex values of precedent calls of the method for the corresponding faces and should be ordered in 
   * increasing order.
   * @param dim Dimension of the face whose boundary is given. If the complex is simplicial, 
   * this parameter can be omitted as it can be deduced from the size of the boundary.
   * @return If it is a chain matrix, the method returns the MatIdx indices of the unpaired chains used to reduce the boundary.
   * Otherwise, nothing.
   */
  template <class Boundary_type>
  insertion_return_type insert_boundary(id_index faceIndex, const Boundary_type& boundary, dimension_type dim = -1);

  /**
   * @brief Returns the column at the given MatIdx index. 
   * For RU matrices, is equivalent to `get_column(columnIndex, true)`.
   * The type of the column depends on the choosen options, see @ref PersistenceMatrixOptions::column_type.
   * 
   * @param columnIndex MatIdx index of the column to return.
   * @return Reference to the column. Is `const` if the matrix has column compression.
   */
  returned_column_type& get_column(index columnIndex);
  /**
   * @brief Only available for chain matrices. Returns the column at the given MatIdx index.
   * The type of the column depends on the choosen options, see @ref PersistenceMatrixOptions::column_type.
   * 
   * @param columnIndex MatIdx index of the column to return.
   * @return Const reference to the column.
   */
  const Column_type& get_column(index columnIndex) const;
  //TODO: there is no particular reason that this method is not available for identifier indexing,
  // just has to be added to the interface...
  /**
   * @brief Only available for RU matrices without @ref Column_indexation_types::IDENTIFIER indexing. 
   * Returns the column at the given MatIdx index in R if @p inR is true and in U if @p inR is false.
   * The type of the column depends on the choosen options, see @ref PersistenceMatrixOptions::column_type.
   * 
   * @param columnIndex MatIdx index of the column to return.
   * @param inR If true, returns the column in R, if false, returns the column in U.
   * @return Const reference to the column.
   */
  const Column_type& get_column(index columnIndex, bool inR);

  //TODO: update column indices when reordering rows (after lazy swap) such that always MatIdx are returned.
  /**
   * @brief Only available if @ref has_row_access is true. Returns the row at the given row index, see [TODO: description].
   * For RU matrices, is equivalent to `get_row(columnIndex, true)`.
   * The type of the row depends on the choosen options, see @ref PersistenceMatrixOptions::has_intrusive_rows.
   * 
   * @param rowIndex Row index of the row to return: IDIdx for chain matrices or updated IDIdx for boundary matrices
   * if swaps occured, see [TODO: description].
   * @return Reference to the row. Is `const` if the matrix has column compression.
   */
  returned_row_type& get_row(id_index rowIndex);
  /**
   * @brief Only available for chain matrices and matrices with column compression.
   * Returns the row at the given row index, see [TODO: description].
   * The type of the row depends on the choosen options, see @ref PersistenceMatrixOptions::has_intrusive_rows.
   * 
   * @param rowIndex Row index of the row to return: IDIdx for chain matrices or updated IDIdx for boundary matrices
   * if swaps occured, see [TODO: description].
   * @return Const reference to the row.
   */
  const Row_type& get_row(id_index rowIndex) const;
  //TODO: there is no particular reason that this method is not available for identifier indexing,
  // just has to be added to the interface...
  /**
   * @brief Only available for RU matrices without @ref Column_indexation_types::IDENTIFIER indexing. 
   * Returns the row at the given row index (see [TODO: description]) in R if @p inR is true and in U if @p inR is false.
   * The type of the row depends on the choosen options, see @ref PersistenceMatrixOptions::has_intrusive_rows.
   * 
   * @param rowIndex Row index of the row to return: updated IDIdx if swaps occured, see [TODO: description].
   * @param inR If true, returns the row in R, if false, returns the row in U.
   * @return Const reference to the row.
   */
  const Row_type& get_row(id_index rowIndex, bool inR);

  /**
   * @brief Only available for base matrices without column compression and if @ref has_map_column_container is true.
   * Otherwise, see @ref remove_last.
   * Erases the given column from the matrix.
   * If @ref has_row_access is also true, the deleted column cells are also automatically removed from their 
   * respective rows.
   * 
   * @param columnIndex MatIdx index of the column to remove.
   */
  void remove_column(index columnIndex);
  //TODO: rename method to be less confusing.
  /**
   * @brief The effect varies depending on the matrices and the options:
   * - Base matrix and boundary matrix:
   *    - @ref has_map_column_container and @ref has_column_and_row_swaps are true:
   *      cleans up maps used for the lazy row swaps.
   *    - @ref has_row_access and @ref has_removable_rows are true: assumes that the row is empty and removes it. 
   *    - Otherwise, does nothing.
   * - Boundary matrix with U stored: only R is affected by the above. If properly used, U will never have empty rows.
   * - Chain matrix: only available if @ref has_row_access and @ref has_removable_rows are true.
   *   Assumes that the row is empty and removes it. 
   *
   * @warning The removed rows are always assumed to be empty. If it is not the case, the deleted row cells are not
   * removed from their columns. And in the case of intrusive rows, this will generate a segmentation fault when 
   * the column cells are destroyed later. The row access is just meant as a "read only" access to the rows and the
   * `erase_row` method just as a way to specify that a row is empty and can therefore be removed from dictionnaries.
   * This allows to avoid testing the emptiness of a row at each column cell removal, what can be quite frequent. 
   * 
   * @param rowIndex Row index of the empty row to remove, see [TODO: description].
   */
  void erase_row(id_index rowIndex);
  //TODO: for chain matrices, replace IDIdx input with MatIdx input to homogenise.
  /**
   * @brief Only available for RU and chain matrices and if @ref has_removable_columns and @ref has_vine_update are true.
   * For chain matrices, @ref has_map_column_container and @ref has_column_pairings also need to be true.
   * Assumes that the face is maximal in the current complex and removes it such that the matrix remains consistent
   * (i.e., RU is still an upper triangular decomposition of the boundary matrix and chain is still a compatible
   * bases of the chain complex in the sense of @cite [TODO: zigzag paper]).
   * The maximality of the face is not verified.
   * Also updates the barcode if it was computed.
   *
   * For chain matrices, using the other version of the method could perform better depending on how the data is 
   * maintained on the side of the user. Then, @ref has_column_pairings also do not need to be true.
   *
   * See also @ref remove_last and @ref remove_column.
   * 
   * @param columnIndex If boundary matrix, MatIdx index of the face to remove, otherwise the IDIdx index.
   */
  void remove_maximal_face(index columnIndex);
  //TODO: See if it would be better to use something more general than a vector for columnsToSwap, such that
  // the user do not have to construct the vector from scratch. Like passing iterators instead. But it would be nice,
  // to still be able to do (face, {})...
  /**
   * @brief Only available for chain matrices and if @ref has_removable_columns, @ref has_vine_update 
   * and @ref has_map_column_container are true.
   * Assumes that the face is maximal in the current complex and removes it such that the matrix remains consistent
   * (i.e., it is still a compatible bases of the chain complex in the sense of @cite [TODO: zigzag paper]).
   * The maximality of the face is not verified.
   * Also updates the barcode if it was computed.
   *
   * To maintain the compatibility, vine swaps are done to move the face up to the end of the filtration. Once at 
   * the end, the removal is trivial. But for chain matrices, swaps do not actually swap the position of the column
   * every time, so the faces appearing after @p faceIndex in the filtration have to be searched first within the matrix.
   * If the user has an easy access to the IDIdx of the faces in the order of filtration, passing them by argument with
   * @p columnsToSwap allows to skip a linear search process. Typically, if the user knows that the face he wants to
   * remove is already the last face of the filtration, calling `remove_maximal_face(faceIndex, {})` will be faster
   * than `remove_last()`.
   *
   * See also @ref remove_last.
   * 
   * @param faceIndex IDIdx index of the face to remove
   * @param columnsToSwap Vector of IDIdx indices of the faces coming after @p faceIndex in the filtration.
   */
  void remove_maximal_face(id_index faceIndex, const std::vector<id_index>& columnsToSwap);
  /**
   * @brief Removes the last inserted column/face from the matrix.
   * If the matrix is non basic, @ref has_removable_columns has to be true for the method to be available.
   * Additionnaly, if the matrix is a chain matrix, either @ref has_map_column_container has to be true or
   * @ref has_vine_update has to be false. And if the matrix is a base matrix it should be without column compression.
   *
   * See also @ref remove_maximal_face and @ref remove_column.
   *
   * For chain matrices, if @ref has_vine_update is true, the last face does not have to be at the end of the matrix
   * and therefore has to be searched first. In this case, if the user already knows the IDIdx of the last face,
   * calling `remove_maximal_face(faceID, {})` instead allows to skip the search.
   */
  void remove_last();

  /**
   * @brief Returns the maximal dimension of a face stored in the matrix. Only available for non basic matrices and
   * if @ref has_matrix_maximal_dimension_access is true.
   * 
   * @return The maximal dimension.
   */
  dimension_type get_max_dimension() const;
  /**
   * @brief Returns the current number of columns in the matrix.
   * 
   * @return The number of columns.
   */
  index get_number_of_columns() const;
  /**
   * @brief Returns the dimension of the given face. Only available for non basic matrices.
   * 
   * @param columnIndex MatIdx index of the column representing the face.
   * @return Dimension of the face.
   */
  dimension_type get_column_dimension(index columnIndex) const;

  /**
   * @brief Adds column at @p sourceColumnIndex onto the column at @p targetColumnIndex in the matrix. Is available
   * for every matrix type, but should be used with care with non basic matrices, as they will be no verification
   * to ensure that the addition makes sense for the meaning of the underlying object. For example, a right-to-left 
   * addition could corrupt the computation of the barcode or the representative cycles if done blindly.
   *
   * For basic matrices with column compression, the representatives are summed together, which means that
   * all column compressed together with the target column are affected by the change, not only the target.
   * 
   * @tparam Index_type Any signed or unsigned integer type.
   * @param sourceColumnIndex MatIdx index of the column to add.
   * @param targetColumnIndex MatIdx index of the target column.
   */
  template <typename Index_type>
  std::enable_if_t<std::is_integral_v<Index_type> > add_to(Index_type sourceColumnIndex, Index_type targetColumnIndex);
  /**
   * @brief Adds the given range of @ref Cell onto the column at @p targetColumnIndex in the matrix. Only available 
   * for basic matrices.
   *
   * For basic matrices with column compression, the range is summed onto the representative, which means that
   * all column compressed together with the target column are affected by the change, not only the target.
   * 
   * @tparam Cell_range Range of @ref Cell. Needs a begin() and end() method. A column index does not need to be
   * stored in the cells, even if @ref has_row_access is true.
   * @param sourceColumn Source cell range.
   * @param targetColumnIndex MatIdx index of the target column.
   */
  template <class Cell_range>
  std::enable_if_t<!std::is_integral_v<Cell_range> > add_to(const Cell_range& sourceColumn, index targetColumnIndex);

  /**
   * @brief Multiplies the target column with the coefficiant and then adds the source column to it.
   * That is: targetColumn = (targetColumn * coefficient) + sourceColumn.
   * Is available for every matrix type, but should be used with care with non basic matrices, as they will be no
   * verification to ensure that the addition makes sense for the meaning of the underlying object.
   * For example, a right-to-left addition could corrupt the computation of the barcode or the representative cycles
   * if done blindly.
   *
   * For basic matrices with column compression, the representatives are summed together, which means that
   * all column compressed together with the target column are affected by the change, not only the target.
   * 
   * @tparam Index_type Any signed or unsigned integer type.
   * @param sourceColumnIndex MatIdx index of the column to add.
   * @param coefficient Value to multiply.
   * @param targetColumnIndex MatIdx index of the target column.
   */
  template <typename Index_type>
  std::enable_if_t<std::is_integral_v<Index_type> > multiply_target_and_add_to(Index_type sourceColumnIndex,
                                                                               int coefficient,
                                                                               Index_type targetColumnIndex);
  /**
   * @brief Multiplies the target column with the coefficiant and then adds the given range of @ref Cell to it.
   * That is: targetColumn = (targetColumn * coefficient) + sourceColumn. Only available for basic matrices.
   *
   * For basic matrices with column compression, the range is summed onto the representative, which means that
   * all column compressed together with the target column are affected by the change, not only the target.
   * 
   * @tparam Cell_range Range of @ref Cell. Needs a begin() and end() method. A column index does not need to be
   * stored in the cells, even if @ref has_row_access is true.
   * @param sourceColumn Source cell range.
   * @param coefficient Value to multiply.
   * @param targetColumnIndex MatIdx index of the target column.
   */
  template <class Cell_range>
  std::enable_if_t<!std::is_integral_v<Cell_range> > multiply_target_and_add_to(const Cell_range& sourceColumn,
                                                                                int coefficient,
                                                                                index targetColumnIndex);

  /**
   * @brief Multiplies the source column with the coefficiant before adding it to the target column.
   * That is: targetColumn += (coefficient * sourceColumn). The source column will **not** be modified.
   * Is available for every matrix type, but should be used with care with non basic matrices, as they will be no
   * verification to ensure that the addition makes sense for the meaning of the underlying object.
   * For example, a right-to-left addition could corrupt the computation of the barcode or the representative cycles
   * if done blindly.
   *
   * For basic matrices with column compression, the representatives are summed together, which means that
   * all column compressed together with the target column are affected by the change, not only the target.
   * 
   * @tparam Index_type Any signed or unsigned integer type.
   * @param coefficient Value to multiply.
   * @param sourceColumnIndex MatIdx index of the column to add.
   * @param targetColumnIndex MatIdx index of the target column.
   */
  template <typename Index_type>
  std::enable_if_t<std::is_integral_v<Index_type> > multiply_source_and_add_to(int coefficient,
                                                                               Index_type sourceColumnIndex,
                                                                               Index_type targetColumnIndex);
  /**
   * @brief Multiplies the source column with the coefficiant before adding it to the target column.
   * That is: targetColumn += (coefficient * sourceColumn). The source column will **not** be modified.
   * Only available for basic matrices.
   *
   * For basic matrices with column compression, the range is summed onto the representative, which means that
   * all column compressed together with the target column are affected by the change, not only the target.
   * 
   * @tparam Cell_range Range of @ref Cell. Needs a begin() and end() method. A column index does not need to be
   * stored in the cells, even if @ref has_row_access is true.
   * @param coefficient Value to multiply.
   * @param sourceColumn Source cell range.
   * @param targetColumnIndex MatIdx index of the target column.
   */
  template <class Cell_range>
  std::enable_if_t<!std::is_integral_v<Cell_range> > multiply_source_and_add_to(int coefficient,
                                                                                const Cell_range& sourceColumn,
                                                                                index targetColumnIndex);

  /**
   * @brief Zeroes the cell at the given coordinates. Not available for chain matrices and for base matrices with 
   * column compression. In general, should be used with care with non basic matrices to not destroy the validity 
   * of the persistence related properties of the matrix.
   *
   * For RU matrices, equivalent to `zero_cell(columnIndex, rowIndex, true)`.
   * 
   * @param columnIndex MatIdx index of the column of the cell.
   * @param rowIndex Row index of the row of the cell.
   */
  void zero_cell(index columnIndex, id_index rowIndex);
  /**
   * @brief Only available for RU matrices. Zeroes the cell at the given coordinates in R if @p inR is true or in
   * U if @p inR is false. Should be used with care to not destroy the validity of the persistence related properties
   * of the matrix.
   * 
   * @param columnIndex MatIdx index of the column of the cell.
   * @param rowIndex Row index of the row of the cell.
   * @param inR Boolean indicating in which matrix to zero: if true in R and if false in U.
   */
  void zero_cell(index columnIndex, id_index rowIndex, bool inR);
  /**
   * @brief Zeroes the column at the given index. Not available for chain matrices and for base matrices with 
   * column compression. In general, should be used with care with non basic matrices to not destroy the validity 
   * of the persistence related properties of the matrix.
   *
   * For RU matrices, equivalent to `zero_column(columnIndex, true)`.
   * 
   * @param columnIndex MatIdx index of the column to zero.
   */
  void zero_column(index columnIndex);
  /**
   * @brief Only available for RU matrices. Zeroes the column at the given index in R if @p inR is true or in
   * U if @p inR is false. Should be used with care to not destroy the validity of the persistence related properties
   * of the matrix.
   * 
   * @param columnIndex MatIdx index of the column to zero.
   * @param inR Boolean indicating in which matrix to zero: if true in R and if false in U.
   */
  void zero_column(index columnIndex, bool inR);
  /**
   * @brief Indicates if the cell at given coordinates has value zero.
   *
   * For RU matrices, equivalent to `is_zero_cell(columnIndex, rowIndex, true)`.
   * 
   * @param columnIndex MatIdx index of the column of the cell.
   * @param rowIndex Row index of the row of the cell.
   * @return true If the cell has value zero.
   * @return false Otherwise.
   */
  bool is_zero_cell(index columnIndex, id_index rowIndex);
  /**
   * @brief Only available for RU matrices. Indicates if the cell at given coordinates has value zero in R if
   * @p inR is true or in U if @p inR is false.
   * 
   * @param columnIndex MatIdx index of the column of the cell.
   * @param rowIndex Row index of the row of the cell.
   * @param inR Boolean indicating in which matrix to look: if true in R and if false in U.
   * @return true If the cell has value zero.
   * @return false Otherwise.
   */
  bool is_zero_cell(index columnIndex, id_index rowIndex, bool inR) const;
  /**
   * @brief Indicates if the column at given index has value zero.
   *
   * For RU matrices, equivalent to `is_zero_column(columnIndex, true)`.
   *
   * Note that for chain matrices, this method should always return false, as a valid chain matrix never has
   * empty columns.
   * 
   * @param columnIndex MatIdx index of the column.
   * @return true If the column has value zero.
   * @return false Otherwise.
   */
  bool is_zero_column(index columnIndex);
  /**
   * @brief Only available for RU matrices. Indicates if the column at given index has value zero in R if
   * @p inR is true or in U if @p inR is false.
   *
   * Note that if @p inR is false, this method should usually return false.
   * 
   * @param columnIndex MatIdx index of the column.
   * @param inR Boolean indicating in which matrix to look: if true in R and if false in U.
   * @return true If the column has value zero.
   * @return false Otherwise.
   */
  bool is_zero_column(index columnIndex, bool inR);

  /**
   * @brief Returns the MatIdx index of the column which has the given row index as pivot. Only available for 
   * RU and chain matrices. Assumes that the pivot exists. For RU matrices, the column is returned from R.
   *
   * Recall that the row indices for chain matrices correspond to the IDIdx indices and that the row indices
   * for a RU matrix correspond to the updated IDIdx indices which got potentially swapped by a vine swap.
   * 
   * @param faceIndex Row index of the pivot.
   * @return MatIdx index of the column with the given pivot.
   */
  index get_column_with_pivot(id_index faceIndex) const;
  /**
   * @brief Returns the row index of the pivot of the given column. Only available for non basic matrices.
   * 
   * @param columnIndex MatIdx index of the column
   * @return The row index of the pivot.
   */
  id_index get_pivot(index columnIndex);

  /**
   * @brief Assign operator.
   * 
   * @param other Matrix to copy
   * @return Reference to this object.
   */
  Matrix& operator=(Matrix other);
  /**
   * @brief Swap operator for two matrices.
   * 
   * @param matrix1 First matrix to swap.
   * @param matrix2 Second matrix to swap.
   */
  friend void swap(Matrix& matrix1, Matrix& matrix2) { 
    swap(matrix1.matrix_, matrix2.matrix_);
    std::swap(matrix1.operators_, matrix2.operators_);
    std::swap(matrix1.cellPool_, matrix2.cellPool_);
  }

  void print();  // for debug

  //TODO: change the behaviour for boundary matrices.
  /**
   * @brief Returns the current barcode of the matrix. Available only if @ref has_column_pairings is true.
   *
   * Recall that we assume that the boundaries were inserted in the order of filtration for the barcode to be valid.
   *
   * @warning For simple boundary matrices (only storing R), we assume that `get_current_barcode` is only called 
   * once, when the matrix is completed.
   * 
   * @return A reference to the barcode. The barcode is a vector of @ref Bar. A bar stores three informations:
   * the PosIdx birth index, the PosIdx death index and the dimension of the bar.
   */
  const barcode_type& get_current_barcode();
  /**
   * @brief Returns the current barcode of the matrix. Available only if @ref has_column_pairings is true.
   *
   * Recall that we assume that the boundaries were inserted in the order of filtration for the barcode to be valid.
   *
   * @warning For simple boundary matrices (only storing R), we assume that `get_current_barcode` is only called 
   * once, when the matrix is completed.
   * 
   * @return A reference to the barcode. The barcode is a vector of @ref Bar. A bar stores three informations:
   * the PosIdx birth index, the PosIdx death index and the dimension of the bar.
   */
  const barcode_type& get_current_barcode() const;

  /**
   * @brief Only available for base matrices without column compression and simple boundary matrices (only storing R)
   * and if @ref has_column_and_row_swaps is true.
   * Swaps the two given columns. Note for boundary matrices, that it really just swaps two columns and do not updates
   * anything else, nor performs additions to maintain some properties on the matrix.
   * 
   * @param columnIndex1 First column MatIdx index to swap.
   * @param columnIndex2 Second column MatIdx index to swap.
   */
  void swap_columns(index columnIndex1, index columnIndex2);
  /**
   * @brief Only available for base matrices without column compression and simple boundary matrices (only storing R)
   * and if @ref has_column_and_row_swaps is true.
   * Swaps the two given rows. Note for boundary matrices, that it really just swaps two rows and do not updates
   * anything else, nor performs additions to maintain some properties on the matrix.
   * 
   * @param rowIndex1 First row index to swap.
   * @param rowIndex2 Second row index to swap.
   */
  void swap_rows(index rowIndex1, index rowIndex2);
  //TODO: find better name. And benchmark also to verify if it is really worth it to have this extra version in addition
  //to vine_swap.
  /**
   * @brief Only available if @ref has_vine_update is true and if it is either a bounary matrix or
   * @ref column_indexation_type is set to @ref Column_indexation_types::POSITION.
   * Does the same than @ref vine_swap, but assumes that the swap is non trivial and
   * therefore skips a part of the case study.
   * 
   * @param index PosIdx index of the first face to swap. The second one has to be at (@p index + 1). Recall that
   * for boundary matrices, PosIdx == MatIdx.
   * @return true If the barcode changed from the swap.
   * @return false Otherwise.
   */
  bool vine_swap_with_z_eq_1_case(pos_index index);
  /**
   * @brief Only available if @ref has_vine_update is true and if it is either a chain matrix or
   * @ref column_indexation_type is set to @ref Column_indexation_types::IDENTIFIER.
   * Does the same than @ref vine_swap, but assumes that the swap is non trivial and
   * therefore skips a part of the case study.
   * 
   * @param columnIndex1 MatIdx index of the first face.
   * @param columnIndex2 MatIdx index of the second face. It is assumed that the PosIdx of both only differs by one.
   * @return Let pos1 be the PosIdx index of @p columnIndex1 and pos2 be the PosIdx index of @p columnIndex2.
   * The method returns the MatIdx of the column which has now, after the swap, the PosIdx max(pos1, pos2).
   */
  index vine_swap_with_z_eq_1_case(index columnIndex1, index columnIndex2);
  /**
   * @brief Only available if @ref has_vine_update is true and if it is either a bounary matrix or
   * @ref column_indexation_type is set to @ref Column_indexation_types::POSITION.
   * Does a vine swap between two faces which are consecutives in the filtration. Roughly, if \f$ F \f$ is the current
   * filtration represented by the matrix, the method modifies the matrix such that the new state corresponds to 
   * a valid state for the filtration \f$ F' \f$ equal to \f$ F \f$ but with the two faces at position @p index
   * and @p index + 1 swapped. Of course, the two faces should not have a face/coface relation which each other ;
   * \f$ F' \f$ has to be a valid filtration.
   * See @cite [TODO: vineyard paper] for more information about vine and vineyards.
   * 
   * @param index PosIdx index of the first face to swap. The second one has to be at (@p index + 1). Recall that
   * for boundary matrices, PosIdx == MatIdx.
   * @return true If the barcode changed from the swap.
   * @return false Otherwise.
   */
  bool vine_swap(pos_index index);
  /**
   * @brief Only available if @ref has_vine_update is true and if it is either a chain matrix or
   * @ref column_indexation_type is set to @ref Column_indexation_types::IDENTIFIER.
   * Does a vine swap between two faces which are consecutives in the filtration. Roughly, if \f$ F \f$ is
   * the current filtration represented by the matrix, the method modifies the matrix such that the new state
   * corresponds to a valid state for the filtration \f$ F' \f$ equal to \f$ F \f$ but with the two given faces
   * at swapped positions. Of course, the two faces should not have a face/coface relation which each other ;
   * \f$ F' \f$ has to be a valid filtration.
   * See @cite [TODO: vineyard paper] for more information about vine and vineyards.
   * 
   * @param columnIndex1 MatIdx index of the first face.
   * @param columnIndex2 MatIdx index of the second face. It is assumed that the PosIdx of both only differs by one.
   * @return Let pos1 be the PosIdx index of @p columnIndex1 and pos2 be the PosIdx index of @p columnIndex2.
   * The method returns the MatIdx of the column which has now, after the swap, the PosIdx max(pos1, pos2).
   */
  index vine_swap(index columnIndex1, index columnIndex2);

  //TODO: Rethink the interface for representative cycles
  /**
   * @brief Only available if @ref can_retrieve_representative_cycles is true. Precomputes the representative cycles
   * of the current state of the filtration represented by the matrix.
   * It does not need to be called before `get_representative_cycles` is called for the first time, but needs to be
   * called before calling `get_representative_cycles` again if the matrix was modified in between. Otherwise the
   * old cycles will be returned.
   */
  void update_representative_cycles();
  /**
   * @brief Only available if @ref can_retrieve_representative_cycles is true.
   * Returns all representative cycles of the current filtration.
   * 
   * @return A const reference to the vector of representative cycles.
   */
  const std::vector<cycle_type>& get_representative_cycles();
  /**
   * @brief Only available if @ref can_retrieve_representative_cycles is true.
   * Returns the cycle representing the given bar.
   * 
   * @param bar A bar from the current barcode.
   * @return A const reference to the cycle representing @p bar.
   */
  const cycle_type& get_representative_cycle(const Bar& bar);

 private:
  using matrix_type = 
    typename std::conditional<
        isNonBasic,
        typename std::conditional<
            PersistenceMatrixOptions::is_of_boundary_type,
            typename std::conditional<
                PersistenceMatrixOptions::has_vine_update || 
                  PersistenceMatrixOptions::can_retrieve_representative_cycles,
                typename std::conditional<
                    PersistenceMatrixOptions::column_indexation_type == Column_indexation_types::CONTAINER || 
                      PersistenceMatrixOptions::column_indexation_type == Column_indexation_types::POSITION,
                    RU_matrix_type, Id_to_index_overlay<RU_matrix_type, Matrix<PersistenceMatrixOptions> >
                >::type,
                typename std::conditional<
                    PersistenceMatrixOptions::column_indexation_type == Column_indexation_types::CONTAINER || 
                      PersistenceMatrixOptions::column_indexation_type == Column_indexation_types::POSITION,
                    Boundary_matrix_type,
                    Id_to_index_overlay<Boundary_matrix_type, Matrix<PersistenceMatrixOptions> >
                >::type
            >::type,
            typename std::conditional<
                PersistenceMatrixOptions::column_indexation_type == Column_indexation_types::CONTAINER, 
                Chain_matrix_type,
                typename std::conditional<
                    PersistenceMatrixOptions::column_indexation_type == Column_indexation_types::POSITION,
                    Position_to_index_overlay<Chain_matrix_type, Matrix<PersistenceMatrixOptions> >,
                    Id_to_index_overlay<Chain_matrix_type, Matrix<PersistenceMatrixOptions> >
                >::type
            >::type
        >::type,
        Base_matrix_type
    >::type;

  Field_operators* operators_;
  Cell_constructor* cellPool_;
  matrix_type matrix_;

  static constexpr void _assert_options();
};

template <class Options>
inline Matrix<Options>::Matrix()
    : operators_(new Field_operators()), cellPool_(new Cell_constructor()), matrix_(operators_, cellPool_) 
{
  static_assert(
      Options::is_of_boundary_type || !Options::has_vine_update || Options::has_column_pairings,
      "When no barcode is recorded with vine swaps, comparaison functions for the columns have to be provided.");
  _assert_options();
}

template <class Options>
template <class Container_type>
inline Matrix<Options>::Matrix(const std::vector<Container_type>& columns, characteristic_type characteristic)
    : operators_(new Field_operators(characteristic)),
      cellPool_(new Cell_constructor()),
      matrix_(columns, operators_, cellPool_) 
{
  static_assert(Options::is_of_boundary_type || !Options::has_vine_update || Options::has_column_pairings,
                "When no barcode is recorded with vine swaps for chain matrices, comparaison functions for the columns "
                "have to be provided.");
  _assert_options();
}

template <class Options>
inline Matrix<Options>::Matrix(int numberOfColumns, characteristic_type characteristic)
    : operators_(new Field_operators(characteristic)),
      cellPool_(new Cell_constructor()),
      matrix_(numberOfColumns, operators_, cellPool_) 
{
  static_assert(Options::is_of_boundary_type || !Options::has_vine_update || Options::has_column_pairings,
                "When no barcode is recorded with vine swaps for chain matrices, comparaison functions for the columns "
                "have to be provided.");
  _assert_options();
}

template <class Options>
template <typename EventComparatorFunction>
inline Matrix<Options>::Matrix(EventComparatorFunction&& birthComparator, EventComparatorFunction&& deathComparator)
    : operators_(new Field_operators()),
      cellPool_(new Cell_constructor()),
      matrix_(operators_, cellPool_, birthComparator, deathComparator) 
{
  static_assert(
      !Options::is_of_boundary_type && Options::has_vine_update && !Options::has_column_pairings,
      "Constructor only available for chain matrices when vine swaps are enabled, but barcodes are not recorded.");
  _assert_options();
}

template <class Options>
template <typename EventComparatorFunction, class Boundary_type>
inline Matrix<Options>::Matrix(const std::vector<Boundary_type>& orderedBoundaries,
                               EventComparatorFunction&& birthComparator, 
                               EventComparatorFunction&& deathComparator,
                               characteristic_type characteristic)
    : operators_(new Field_operators(characteristic)),
      cellPool_(new Cell_constructor()),
      matrix_(orderedBoundaries, operators_, cellPool_, birthComparator, deathComparator) 
{
  static_assert(
      !Options::is_of_boundary_type && Options::has_vine_update && !Options::has_column_pairings,
      "Constructor only available for chain matrices when vine swaps are enabled, but barcodes are not recorded.");
  _assert_options();
}

template <class Options>
template <typename EventComparatorFunction>
inline Matrix<Options>::Matrix(unsigned int numberOfColumns, 
                               EventComparatorFunction&& birthComparator,
                               EventComparatorFunction&& deathComparator, 
                               characteristic_type characteristic)
    : operators_(new Field_operators(characteristic)),
      cellPool_(new Cell_constructor()),
      matrix_(numberOfColumns, operators_, cellPool_, birthComparator, deathComparator) 
{
  static_assert(
      !Options::is_of_boundary_type && Options::has_vine_update && !Options::has_column_pairings,
      "Constructor only available for chain matrices when vine swaps are enabled, but barcodes are not recorded.");
  _assert_options();
}

template <class Options>
inline Matrix<Options>::Matrix(const Matrix& matrixToCopy)
    : operators_(new Field_operators(matrixToCopy.operators_->get_characteristic())),
      cellPool_(new Cell_constructor()),
      matrix_(matrixToCopy.matrix_, operators_, cellPool_) 
{
  _assert_options();
}

template <class Options>
inline Matrix<Options>::Matrix(Matrix&& other) noexcept
    : operators_(std::exchange(other.operators_, nullptr)),
      cellPool_(std::exchange(other.cellPool_, nullptr)),
      matrix_(std::move(other.matrix_)) 
{
  // TODO: verify that the address of operators_ == address of other.operators_ after move
  // and that therefore the addresses stored in matrix_ are correct.
  _assert_options();
}

template <class Options>
inline Matrix<Options>::~Matrix() 
{
  matrix_.reset(operators_, cellPool_);
  delete cellPool_;
  delete operators_;
}

template <class Options>
inline void Matrix<Options>::set_characteristic(characteristic_type characteristic) 
{
  if constexpr (!Options::is_z2) {
    if (operators_->get_characteristic() != 0) {
      std::cerr << "Warning: Characteristic already initialised. Changing it could lead to incoherences in the matrice "
                   "as the modulo was already applied to values in existing columns.";
    }

    operators_->set_characteristic(characteristic);
  }
}

template <class Options>
template <class Container_type>
inline void Matrix<Options>::insert_column(const Container_type& column) 
{
  assert(operators_->get_characteristic() != 0 &&
         "Columns cannot be initialized if the coefficient field characteristic is not specified.");
  static_assert(
      !isNonBasic,
      "'insert_column' not available for the chosen options. The input has to be in the form of a face boundary.");
  matrix_.insert_column(column);
}

template <class Options>
template <class Container_type>
inline void Matrix<Options>::insert_column(const Container_type& column, index columnIndex) 
{
  assert(operators_->get_characteristic() != 0 &&
         "Columns cannot be initialized if the coefficient field characteristic is not specified.");
  static_assert(!isNonBasic && !Options::has_column_compression,
                "'insert_column' with those parameters is not available for the chosen options.");
  static_assert(!Options::has_row_access,
                "Columns have to be inserted at the end of the matrix when row access is enabled.");
  matrix_.insert_column(column, columnIndex);
}

template <class Options>
template <class Boundary_type>
inline typename Matrix<Options>::insertion_return_type Matrix<Options>::insert_boundary(const Boundary_type& boundary,
                                                                                        dimension_type dim) 
{
  assert(operators_->get_characteristic() != 0 &&
         "Columns cannot be initialized if the coefficient field characteristic is not specified.");
  if constexpr (isNonBasic && !Options::is_of_boundary_type &&
                Options::column_indexation_type == Column_indexation_types::CONTAINER)
    return matrix_.insert_boundary(boundary, dim);
  else
    matrix_.insert_boundary(boundary, dim);
}

template <class Options>
template <class Boundary_type>
inline typename Matrix<Options>::insertion_return_type Matrix<Options>::insert_boundary(id_index faceIndex,
                                                                                        const Boundary_type& boundary,
                                                                                        dimension_type dim) 
{
  assert(operators_->get_characteristic() != 0 &&
         "Columns cannot be initialized if the coefficient field characteristic is not specified.");
  static_assert(isNonBasic, "Only enabled for non-basic matrices.");
  if constexpr (!Options::is_of_boundary_type &&
                Options::column_indexation_type == Column_indexation_types::CONTAINER)
    return matrix_.insert_boundary(faceIndex, boundary, dim);
  else
    matrix_.insert_boundary(faceIndex, boundary, dim);
}

template <class Options>
inline typename Matrix<Options>::returned_column_type& Matrix<Options>::get_column(index columnIndex) 
{
  return matrix_.get_column(columnIndex);
}

template <class Options>
inline const typename Matrix<Options>::Column_type& Matrix<Options>::get_column(index columnIndex) const 
{
  return matrix_.get_column(columnIndex);
}

template <class Options>
inline const typename Matrix<Options>::Column_type& Matrix<Options>::get_column(index columnIndex, bool inR) 
{
  // TODO: I don't think there is a particular reason why the indexation is forced, should be removed.
  static_assert(isNonBasic && Options::is_of_boundary_type &&
                    (Options::has_vine_update || Options::can_retrieve_representative_cycles) &&
                    Options::column_indexation_type != Column_indexation_types::IDENTIFIER,
                "Only enabled for position indexed RU matrices.");

  return matrix_.get_column(columnIndex, inR);
}

template <class Options>
inline typename Matrix<Options>::returned_row_type& Matrix<Options>::get_row(id_index rowIndex) 
{
  static_assert(Options::has_row_access, "'get_row' is not available for the chosen options.");

  return matrix_.get_row(rowIndex);
}

template <class Options>
inline const typename Matrix<Options>::Row_type& Matrix<Options>::get_row(id_index rowIndex) const 
{
  static_assert(Options::has_row_access, "'get_row' is not available for the chosen options.");

  return matrix_.get_row(rowIndex);
}

template <class Options>
inline const typename Matrix<Options>::Row_type& Matrix<Options>::get_row(id_index rowIndex, bool inR) 
{
  static_assert(Options::has_row_access, "'get_row' is not available for the chosen options.");
  // TODO: I don't think there is a particular reason why the indexation is forced, should be removed.
  static_assert(isNonBasic && Options::is_of_boundary_type &&
                    (Options::has_vine_update || Options::can_retrieve_representative_cycles) &&
                    Options::column_indexation_type != Column_indexation_types::IDENTIFIER,
                "Only enabled for position indexed RU matrices.");

  return matrix_.get_row(rowIndex, inR);
}

template <class Options>
inline void Matrix<Options>::remove_column(index columnIndex) 
{
  static_assert(Options::has_map_column_container && !isNonBasic && !Options::has_column_compression,
                "'remove_column' is not available for the chosen options.");

  matrix_.remove_column(columnIndex);
}

template <class Options>
inline void Matrix<Options>::erase_row(id_index rowIndex) 
{
  static_assert(!isNonBasic || Options::is_of_boundary_type || Options::has_removable_rows,
                "'erase_row' is not available for the chosen options.");

  matrix_.erase_row(rowIndex);
}

template <class Options>
inline void Matrix<Options>::remove_maximal_face(index columnIndex) 
{
  static_assert(Options::has_removable_columns,
                "'remove_maximal_face(id_index)' is not available for the chosen options.");
  static_assert(isNonBasic && Options::has_vine_update,
                "'remove_maximal_face(id_index)' is not available for the chosen options.");
  static_assert(Options::is_of_boundary_type || (Options::has_map_column_container && Options::has_column_pairings),
                "'remove_maximal_face(id_index)' is not available for the chosen options.");

  matrix_.remove_maximal_face(columnIndex);
}

template <class Options>
inline void Matrix<Options>::remove_maximal_face(id_index faceIndex, const std::vector<id_index>& columnsToSwap) 
{
  static_assert(Options::has_removable_columns,
                "'remove_maximal_face(id_index,const std::vector<index>&)' is not available for the chosen options.");
  static_assert(isNonBasic && !Options::is_of_boundary_type,
                "'remove_maximal_face(id_index,const std::vector<index>&)' is not available for the chosen options.");
  static_assert(Options::has_map_column_container && Options::has_vine_update,
                "'remove_maximal_face(id_index,const std::vector<index>&)' is not available for the chosen options.");

  matrix_.remove_maximal_face(faceIndex, columnsToSwap);
}

template <class Options>
inline void Matrix<Options>::remove_last() 
{
  static_assert(Options::has_removable_columns || !isNonBasic,
                "'remove_last' is not available for the chosen options.");
  static_assert(!Options::has_column_compression || isNonBasic,
                "'remove_last' is not available for the chosen options.");
  static_assert(
      !isNonBasic || Options::is_of_boundary_type || Options::has_map_column_container || !Options::has_vine_update,
      "'remove_last' is not available for the chosen options.");

  matrix_.remove_last();
}

template <class Options>
inline typename Matrix<Options>::dimension_type Matrix<Options>::get_max_dimension() const 
{
  static_assert(isNonBasic, "'get_max_dimension' is not available for the chosen options.");

  return matrix_.get_max_dimension();
}

template <class Options>
inline typename Matrix<Options>::index Matrix<Options>::get_number_of_columns() const 
{
  return matrix_.get_number_of_columns();
}

template <class Options>
inline typename Matrix<Options>::dimension_type Matrix<Options>::get_column_dimension(index columnIndex) const 
{
  static_assert(isNonBasic, "'get_column_dimension' is not available for the chosen options.");

  return matrix_.get_column_dimension(columnIndex);
}

template <class Options>
template <typename Index_type>
inline std::enable_if_t<std::is_integral_v<Index_type> > Matrix<Options>::add_to(Index_type sourceColumnIndex,
                                                                                Index_type targetColumnIndex) 
{
  matrix_.add_to(sourceColumnIndex, targetColumnIndex);
}

template <class Options>
template <class Cell_range>
inline std::enable_if_t<!std::is_integral_v<Cell_range> > Matrix<Options>::add_to(const Cell_range& sourceColumn,
                                                                                 index targetColumnIndex) 
{
  static_assert(!isNonBasic,
                "For boundary or chain matrices, only additions with columns inside the matrix is allowed to maintain "
                "algebraic consistency.");

  matrix_.add_to(sourceColumn, targetColumnIndex);
}

template <class Options>
template <typename Index_type>
inline std::enable_if_t<std::is_integral_v<Index_type> > Matrix<Options>::multiply_target_and_add_to(
    Index_type sourceColumnIndex, int coefficient, Index_type targetColumnIndex) 
{
  if constexpr (Options::is_z2) {
    // coef will be converted to bool, because of element_type
    matrix_.multiply_target_and_add_to(sourceColumnIndex, coefficient % 2, targetColumnIndex);
  } else {
    matrix_.multiply_target_and_add_to(sourceColumnIndex, operators_->get_value(coefficient), targetColumnIndex);
  }
}

template <class Options>
template <class Cell_range>
inline std::enable_if_t<!std::is_integral_v<Cell_range> > Matrix<Options>::multiply_target_and_add_to(
    const Cell_range& sourceColumn, int coefficient, index targetColumnIndex) 
{
  static_assert(!isNonBasic,
                "For boundary or chain matrices, only additions with columns inside the matrix is allowed to maintain "
                "algebraic consistency.");

  if constexpr (Options::is_z2) {
    // coef will be converted to bool, because of element_type
    matrix_.multiply_target_and_add_to(sourceColumn, coefficient % 2, targetColumnIndex);
  } else {
    matrix_.multiply_target_and_add_to(sourceColumn, operators_->get_value(coefficient), targetColumnIndex);
  }
}

template <class Options>
template <typename Index_type>
inline std::enable_if_t<std::is_integral_v<Index_type> > Matrix<Options>::multiply_source_and_add_to(
    int coefficient, Index_type sourceColumnIndex, Index_type targetColumnIndex) 
{
  if constexpr (Options::is_z2) {
    // coef will be converted to bool, because of element_type
    matrix_.multiply_source_and_add_to(coefficient % 2, sourceColumnIndex, targetColumnIndex);
  } else {
    matrix_.multiply_source_and_add_to(operators_->get_value(coefficient), sourceColumnIndex, targetColumnIndex);
  }
}

template <class Options>
template <class Cell_range>
inline std::enable_if_t<!std::is_integral_v<Cell_range> > Matrix<Options>::multiply_source_and_add_to(
    int coefficient, const Cell_range& sourceColumn, index targetColumnIndex) 
{
  static_assert(!isNonBasic,
                "For boundary or chain matrices, only additions with columns inside the matrix is allowed to maintain "
                "algebraic consistency.");

  if constexpr (Options::is_z2) {
    // coef will be converted to bool, because of element_type
    matrix_.multiply_source_and_add_to(coefficient % 2, sourceColumn, targetColumnIndex);
  } else {
    matrix_.multiply_source_and_add_to(operators_->get_value(coefficient), sourceColumn, targetColumnIndex);
  }
}

template <class Options>
inline void Matrix<Options>::zero_cell(index columnIndex, id_index rowIndex) 
{
  static_assert(Options::is_of_boundary_type && !Options::has_column_compression,
                "'zero_cell' is not available for the chosen options.");

  return matrix_.zero_cell(columnIndex, rowIndex);
}

template <class Options>
inline void Matrix<Options>::zero_cell(index columnIndex, id_index rowIndex, bool inR) 
{
  // TODO: I don't think there is a particular reason why the indexation is forced, should be removed.
  static_assert(isNonBasic && Options::is_of_boundary_type &&
                    (Options::has_vine_update || Options::can_retrieve_representative_cycles) &&
                    Options::column_indexation_type != Column_indexation_types::IDENTIFIER,
                "Only enabled for RU matrices.");

  return matrix_.zero_cell(columnIndex, rowIndex, inR);
}

template <class Options>
inline void Matrix<Options>::zero_column(index columnIndex) 
{
  static_assert(Options::is_of_boundary_type && !Options::has_column_compression,
                "'zero_column' is not available for the chosen options.");

  return matrix_.zero_column(columnIndex);
}

template <class Options>
inline void Matrix<Options>::zero_column(index columnIndex, bool inR) 
{
  // TODO: I don't think there is a particular reason why the indexation is forced, should be removed.
  static_assert(isNonBasic && Options::is_of_boundary_type &&
                    (Options::has_vine_update || Options::can_retrieve_representative_cycles) &&
                    Options::column_indexation_type != Column_indexation_types::IDENTIFIER,
                "Only enabled for RU matrices.");

  return matrix_.zero_column(columnIndex, inR);
}

template <class Options>
inline bool Matrix<Options>::is_zero_cell(index columnIndex, id_index rowIndex) 
{
  return matrix_.is_zero_cell(columnIndex, rowIndex);
}

template <class Options>
inline bool Matrix<Options>::is_zero_cell(index columnIndex, id_index rowIndex, bool inR) const 
{
  // TODO: I don't think there is a particular reason why the indexation is forced, should be removed.
  static_assert(isNonBasic && Options::is_of_boundary_type &&
                    (Options::has_vine_update || Options::can_retrieve_representative_cycles) &&
                    Options::column_indexation_type != Column_indexation_types::IDENTIFIER,
                "Only enabled for RU matrices.");

  return matrix_.is_zero_cell(columnIndex, rowIndex, inR);
}

template <class Options>
inline bool Matrix<Options>::is_zero_column(index columnIndex) 
{
  return matrix_.is_zero_column(columnIndex);
}

template <class Options>
inline bool Matrix<Options>::is_zero_column(index columnIndex, bool inR) 
{
  // TODO: I don't think there is a particular reason why the indexation is forced, should be removed.
  static_assert(isNonBasic && Options::is_of_boundary_type &&
                    (Options::has_vine_update || Options::can_retrieve_representative_cycles) &&
                    Options::column_indexation_type != Column_indexation_types::IDENTIFIER,
                "Only enabled for RU matrices.");

  return matrix_.is_zero_column(columnIndex, inR);
}

template <class Options>
inline typename Matrix<Options>::index Matrix<Options>::get_column_with_pivot(id_index faceIndex) const 
{
  static_assert(isNonBasic && (!Options::is_of_boundary_type ||
                               (Options::has_vine_update || Options::can_retrieve_representative_cycles)),
                "'get_column_with_pivot' is not available for the chosen options.");

  return matrix_.get_column_with_pivot(faceIndex);
}

template <class Options>
inline typename Matrix<Options>::id_index Matrix<Options>::get_pivot(index columnIndex) 
{
  static_assert(isNonBasic, "'get_pivot' is not available for the chosen options.");

  return matrix_.get_pivot(columnIndex);
}

template <class Options>
inline Matrix<Options>& Matrix<Options>::operator=(Matrix other) 
{
  swap(matrix_, other.matrix_);
  std::swap(operators_, other.operators_);
  std::swap(cellPool_, other.cellPool_);
  // if constexpr (!Options::is_z2) matrix_.set_operators(&operators_);

  return *this;
}

template <class Options>
inline void Matrix<Options>::print() 
{
  return matrix_.print();
}

template <class Options>
inline const typename Matrix<Options>::barcode_type& Matrix<Options>::get_current_barcode() 
{
  static_assert(Options::has_column_pairings, "This method was not enabled.");

  return matrix_.get_current_barcode();
}

template <class Options>
inline const typename Matrix<Options>::barcode_type& Matrix<Options>::get_current_barcode() const 
{
  static_assert(Options::has_column_pairings, "This method was not enabled.");
  static_assert(
      !Options::is_of_boundary_type || Options::has_vine_update || Options::can_retrieve_representative_cycles,
      "'get_current_barcode' is not const for boundary matrices as the barcode is only computed when explicitely "
      "asked.");

  return matrix_.get_current_barcode();
}

template <class Options>
inline void Matrix<Options>::swap_columns(index columnIndex1, index columnIndex2) 
{
  static_assert((!isNonBasic && !Options::has_column_compression) ||
                    (isNonBasic && Options::is_of_boundary_type && !Options::has_vine_update &&
                     !Options::can_retrieve_representative_cycles),
                "This method was not enabled.");
  return matrix_.swap_columns(columnIndex1, columnIndex2);
}

template <class Options>
inline void Matrix<Options>::swap_rows(index rowIndex1, index rowIndex2) 
{
  static_assert((!isNonBasic && !Options::has_column_compression) ||
                    (isNonBasic && Options::is_of_boundary_type && !Options::has_vine_update &&
                     !Options::can_retrieve_representative_cycles),
                "This method was not enabled.");
  return matrix_.swap_rows(rowIndex1, rowIndex2);
}

template <class Options>
inline bool Matrix<Options>::vine_swap_with_z_eq_1_case(pos_index index) 
{
  static_assert(Options::has_vine_update, "This method was not enabled.");
  static_assert(
      Options::column_indexation_type == Column_indexation_types::POSITION ||
          (Options::is_of_boundary_type && Options::column_indexation_type == Column_indexation_types::CONTAINER),
      "This method was not enabled.");
  return matrix_.vine_swap_with_z_eq_1_case(index);
}

template <class Options>
inline typename Matrix<Options>::index Matrix<Options>::vine_swap_with_z_eq_1_case(index columnIndex1,
                                                                                   index columnIndex2) 
{
  static_assert(Options::has_vine_update, "This method was not enabled.");
  static_assert(
      Options::column_indexation_type == Column_indexation_types::IDENTIFIER ||
          (!Options::is_of_boundary_type && Options::column_indexation_type == Column_indexation_types::CONTAINER),
      "This method was not enabled.");
  return matrix_.vine_swap_with_z_eq_1_case(columnIndex1, columnIndex2);
}

template <class Options>
inline bool Matrix<Options>::vine_swap(pos_index index) 
{
  static_assert(Options::has_vine_update, "This method was not enabled.");
  static_assert(
      Options::column_indexation_type == Column_indexation_types::POSITION ||
          (Options::is_of_boundary_type && Options::column_indexation_type == Column_indexation_types::CONTAINER),
      "This method was not enabled.");
  return matrix_.vine_swap(index);
}

template <class Options>
inline typename Matrix<Options>::index Matrix<Options>::vine_swap(index columnIndex1, index columnIndex2) 
{
  static_assert(Options::has_vine_update, "This method was not enabled.");
  static_assert(
      Options::column_indexation_type == Column_indexation_types::IDENTIFIER ||
          (!Options::is_of_boundary_type && Options::column_indexation_type == Column_indexation_types::CONTAINER),
      "This method was not enabled.");
  return matrix_.vine_swap(columnIndex1, columnIndex2);
}

template <class Options>
inline void Matrix<Options>::update_representative_cycles() 
{
  static_assert(Options::can_retrieve_representative_cycles, "This method was not enabled.");
  matrix_.update_representative_cycles();
}

template <class Options>
inline const std::vector<typename Matrix<Options>::cycle_type>& Matrix<Options>::get_representative_cycles() 
{
  static_assert(Options::can_retrieve_representative_cycles, "This method was not enabled.");
  return matrix_.get_representative_cycles();
}

template <class Options>
inline const typename Matrix<Options>::cycle_type& Matrix<Options>::get_representative_cycle(const Bar& bar) 
{
  static_assert(Options::can_retrieve_representative_cycles, "This method was not enabled.");
  return matrix_.get_representative_cycle(bar);
}

template <class Options>
inline constexpr void Matrix<Options>::_assert_options() 
{
  static_assert(Options::column_type != Column_types::HEAP || !Options::has_row_access,
                "Row access is not possible for heap columns.");
  static_assert(!Options::has_vine_update || Options::is_z2, "Vine update currently works only for Z_2 coefficients.");
  // static_assert(!Options::can_retrieve_representative_cycles || Options::is_z2,
  //               "Representaive cycles can currently only be computed with Z_2 coefficients.");
  static_assert(Options::column_type != Column_types::HEAP || !Options::has_column_compression,
                "Column compression not compatible with heap columns.");

//   // This should be warnings instead, as Options::has_column_compression is just ignored in those cases and don't
//   // produces errors as long as the corresponding methods are not called.
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
//   static_assert(!Options::has_column_compression || !Options::has_map_column_container,
//                 "When column compression is used, the removal of columns is not implemented yet.");
}

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // MASTER_MATRIX_H

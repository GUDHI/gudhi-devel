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
 *   to be stored. Note that the option @ref is_indexed_by_position will produce a small overhead when set to **false**.
 * - a chain complex matrix representing a `compatible base` of a filtered chain complex (see TODO: cite Clément's zigzag paper here).
 *   This matrix is deduced from the boundary matrix and therefore encodes more or less the same information 
 *   but differently and can therefore be better suited for certain applications. This type can be used the same way 
 *   than the precedent type, only the option @ref is_of_boundary_type has to be set to false. So it is easy to switch
 *   from one representation to the other if one wants to test both. Just note that the option 
 *   @ref is_indexed_by_position will produce a small overhead when set to **true**.
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
 *
 * In conclusion, with default values, if no vine swaps or removals occurs, all three indexing schemes are the same.
 * 
 * Different from columns, rows are always indexed by ID, i.e., the values used in column/boundary for insertions.
 * So, for row access, it is usualy necessary to remember the indexing scheme used.
 * 
 * @tparam Options Structure encoding all the options of the matrix. 
 * See description of @ref Default_options for more details.
 */
template <class Options = Default_options<> >
class Matrix {
 public:
  using Option_list = Options;	//to make it accessible from the other classes
  using index = typename Options::index_type;                 /**< Type of MatIdx index. */
  using id_index = typename Options::id_type;                 /**< Type of IDIdx index or row index. */
  using pos_index = typename Options::pos_type;               /**< Type of PosIdx index. */
  using dimension_type = typename Options::dimension_type;    /**< Type for dimension. */

  struct Dummy_field_operators{
    using element_type = unsigned int;
    using characteristic_type = element_type;

    Dummy_field_operators([[maybe_unused]] characteristic_type characteristic = 0){}

    friend void swap([[maybe_unused]] Dummy_field_operators& d1, [[maybe_unused]] Dummy_field_operators& d2){}

    static constexpr characteristic_type get_characteristic() { return 2; }
  };

  using Field_operators =
      typename std::conditional<Options::is_z2, 
                                Dummy_field_operators, 
                                typename Options::field_coeff_operators
                               >::type;                       /**< Coefficiants field type. */
  using element_type = typename std::conditional<Options::is_z2, 
                                                 bool, 
                                                 typename Field_operators::element_type
                                                >::type;
  using characteristic_type =
      typename std::conditional<Options::is_z2, 
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
                                Cell_field_element<element_type>
                               >::type;
  /**
   * @brief Type of a matrix cell. See @ref Cell for a more detailed description.
   */
  using Cell_type = Cell<Matrix<Options> >;
  inline static New_cell_constructor<Cell_type> defaultCellConstructor;
  using Cell_constructor = Pool_cell_constructor<Cell_type>;

  using cell_rep_type = typename std::conditional<Options::is_z2,
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
   * @brief Type of the rows storted in the matrix. Is either an intrsuive list of @ref Cell_type (not ordered) if 
   * @ref has_intrusive_rows is true, or a set of @ref Cell_type (ordered by @ref get_column_index) otherwise.
   */
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
                                std::map<id_index, Row_type>,
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
                                Matrix_row_access<Row_type, row_container_type, Options::has_removable_rows, id_index>,
                                Dummy_matrix_row_access
                               >::type;

  template <typename value_type>
  using dictionnary_type =
      typename std::conditional<Options::has_map_column_container,
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

  using Heap_column_type = Heap_column<Matrix<Options>, Cell_constructor>;
  using List_column_type = List_column<Matrix<Options>, Cell_constructor>;
  using Vector_column_type = Vector_column<Matrix<Options>, Cell_constructor>;
  using Naive_vector_column_type = Naive_vector_column<Matrix<Options>, Cell_constructor>;
  using Set_column_type = Set_column<Matrix<Options>, Cell_constructor>;
  using Unordered_set_column_type = Unordered_set_column<Matrix<Options>, Cell_constructor>;
  using Intrusive_list_column_type = Intrusive_list_column<Matrix<Options>, Cell_constructor>;
  using Intrusive_set_column_type = Intrusive_set_column<Matrix<Options>, Cell_constructor>;

  /**
   * @brief Type of the columns stored in the matrix. The type depends on the value of @ref column_type defined
   * in the given options. See @ref Column_types for a more detailed description.
   */
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
      typename std::conditional<Options::has_map_column_container, 
                                std::unordered_map<index, Column_type>,
                                std::vector<Column_type>
                               >::type;

  static const bool hasFixedBarcode = Option_list::is_of_boundary_type && !Options::has_vine_update;
  /**
   * @brief Type of the computed barcode. If @ref has_map_column_container is true, it is a list of @ref Bar,
   * otherwise it is a vector of @ref Bar.
   */
  using barcode_type =
      typename std::conditional<hasFixedBarcode, 
                                std::vector<Bar>, 
                                typename std::conditional<Options::has_removable_columns, 
                                  std::list<Bar>, 
                                  std::vector<Bar>
                                >::type
                               >::type;
  using bar_dictionnary_type = 
      typename std::conditional<hasFixedBarcode,
                                typename std::conditional<Options::can_retrieve_representative_cycles,
                                  std::vector<index>,                   //RU
                                  std::unordered_map<pos_index, index>  //boundary
                                >::type,
                                typename std::conditional<Options::has_removable_columns,
                                  std::unordered_map<pos_index, typename barcode_type::iterator>,
                                  std::vector<index>
                                >::type
                               >::type;

  //default type for boundaries to permit list initialization directly in function parameters
  using boundary_type = typename std::conditional<Options::is_z2, 
                                                  std::initializer_list<id_index>,
                                                  std::initializer_list<std::pair<id_index, element_type> >
                                                 >::type;

  static const bool dimensionIsNeeded = Options::has_column_pairings && Options::is_of_boundary_type &&
                                        !Options::has_vine_update && !Options::can_retrieve_representative_cycles;

  using Matrix_dimension_option = typename std::conditional<
      Options::has_matrix_maximal_dimension_access || dimensionIsNeeded,
      typename std::conditional<Options::has_removable_columns, 
                                Matrix_all_dimension_holder<dimension_type>,
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

  /**
   * @brief Type of a representative cycle. Vector of PosIdx indices for boundary matrices and vector of IDIdx
   * indices for chain matrices.
   */
  using cycle_type = std::vector<id_index>;	//TODO: add coefficients

  //Return types to factorize the corresponding methods

  //The returned column is `const` if the matrix uses column compression
  using returned_column_type =
      typename std::conditional<!isNonBasic && Options::has_column_compression,
                                const Column_type,
                                Column_type
                               >::type;
  //The returned row is `const` if the matrix uses column compression
  using returned_row_type =
      typename std::conditional<!isNonBasic && Options::has_column_compression,
                                const Row_type,
                                Row_type
                               >::type;
  //If the matrix is a chain matrix, the insertion method returns the pivots of its unpaired columns used to reduce
  //the inserted boundary. Otherwise, void.
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
  Matrix(const std::vector<Container_type>& columns, characteristic_type characteristic = 11);
  /**
   * @brief Constructs a new empty matrix and reserves space for the given number of columns.
   * 
   * @param numberOfColumns Number of columns to reserve space for.
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
   * @tparam EventComparatorFunction Type of the birth comparator: (unsigned int, unsigned int) -> bool
   * @tparam EventComparatorFunction Type of the death comparator: (unsigned int, unsigned int) -> bool
   * @param birthComparator Method taking two IDIdx indices as parameter and returns true if and only if the first 
   * face is associated to a bar with strictly smaller birth than the bar associated to the second one.
   * @param deathComparator Method taking two IDIdx indices as parameter and returns true if and only if the first 
   * face is associated to a bar with strictly smaller death than the bar associated to the second one.
   */
  template <typename EventComparatorFunction>
  Matrix(EventComparatorFunction&& birthComparator, 
         EventComparatorFunction&& deathComparator);
  /**
   * @brief Constructs a new matrix from the given matrix with the given comparator functions. 
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
   * @tparam EventComparatorFunction Type of the birth comparator: (unsigned int, unsigned int) -> bool
   * @tparam EventComparatorFunction Type of the death comparator: (unsigned int, unsigned int) -> bool
   * @tparam Boundary_type Range type for a column. Assumed to have a begin(), end() and size() method.
   * @param orderedBoundaries Vector of ordered boundaries in filtration order. Indexed continously starting at 0.
   * @param birthComparator Method taking two IDIdx indices as parameter and returns true if and only if the first 
   * face is associated to a bar with strictly smaller birth than the bar associated to the second one.
   * @param deathComparator Method taking two IDIdx indices as parameter and returns true if and only if the first 
   * face is associated to a bar with strictly smaller death than the bar associated to the second one.
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
   * @tparam EventComparatorFunction Type of the birth comparator: (unsigned int, unsigned int) -> bool
   * @tparam EventComparatorFunction Type of the death comparator: (unsigned int, unsigned int) -> bool
   * @param numberOfColumns Number of columns to reserve space for.
   * @param birthComparator Method taking two IDIdx indices as parameter and returns true if and only if the first 
   * face is associated to a bar with strictly smaller birth than the bar associated to the second one.
   * @param deathComparator Method taking two IDIdx indices as parameter and returns true if and only if the first 
   * face is associated to a bar with strictly smaller death than the bar associated to the second one.
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
  void set_characteristic(characteristic_type characteristic);

  // (TODO: if there is no row access and the column type corresponds to the internal column type of the matrix, 
  // moving the column instead of copying it should be possible. Is it worth implementing it?)
  /**
   * @brief Inserts a new ordered column at the end of the matrix by copying the given column. 
   * Only available when **all** of the following options are **false**:
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
   * @brief Inserts a new ordered column at the given index by copying the given column.
   * There should not be any other column inserted at that index which was not explicitely removed before.
   * Only available if @ref has_map_column_container is true and **all** of the following options are **false**:
   *   - @ref has_column_pairings
   *   - @ref has_vine_update
   *   - @ref can_retrieve_representative_cycles
   *   - @ref has_column_compression
   * 
   * @tparam Container_type Range type for a column. Assumed to have a begin(), end() and size() method.
   * @param column Column to be inserted.
   * @param columnIndex MatIdx index to which the column has to be inserted.
   */
  template <class Container_type>
  void insert_column(const Container_type& column, index columnIndex);
  /**
   * @brief Inserts at the end of the matrix a new ordered column corresponding to the given boundary. 
   * This means that we assume that the boundaries are inserted in the order of the filtration. 
   * For chain matrices or ID indexed matrices, we also assume that the faces in the given boundary are identified by
   * their position in the filtration (not counting removals in the case of zigzag), starting at 0. If it is not the
   * case, use the other `insert_boundary` instead by indicating the face ID used in the boundaries when the face is
   * inserted.
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
   *   `IDIdx + linear combination of older column IDIdxs`, where the combination is deduced while reducing the 
   *   given boundary. If the barcode is stored, it will also be updated.
   * 
   * @tparam Boundary_type Range type for a column. Assumed to have a begin(), end() and size() method.
   * @param boundary Boundary generating the new column.
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
   * @tparam Boundary_type Range type for a column. Assumed to have a begin(), end() and size() method.
   * @param faceIndex IDIdx index to be used to indentify the new face.
   * @param boundary Boundary generating the new column. The indices of the boundary have to correspond to the 
   * @p faceIndex values of precedent calls of the method for the corresponding faces.
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
   * 
   * @param columnIndex MatIdx index of the column to return.
   * @return Reference to the column. Is `const` if the matrix has column compression.
   */
  returned_column_type& get_column(index columnIndex);
  /**
   * @brief Only available for chain matrices. Returns the column at the given MatIdx index.
   * 
   * @param columnIndex MatIdx index of the column to return.
   * @return Const reference to the column.
   */
  const Column_type& get_column(index columnIndex) const;
  /**
   * @brief Only available for RU matrices with position indexing. 
   * Returns the column at the given MatIdx index in R if @p inR is true and in U if @p inR is false.
   * 
   * @param columnIndex MatIdx index of the column to return.
   * @param inR If true, returns the column in R, if false, returns the column in U.
   * @return Const reference to the column.
   */
  const Column_type& get_column(index columnIndex, bool inR);

  //TODO: update column indices when reordering rows (after lazy swap) such that always MatIdx are returned.
  /**
   * @brief Only available if @ref has_row_access is true. Returns the row at the given IDIdx index.
   * For RU matrices, is equivalent to `get_row(columnIndex, true)`.
   *
   * @warning The @ref get_column_index method of the row cells returns IDIdx indices for boundary matrices and
   * MatIdx indices for chain matrices.
   * 
   * @param rowIndex IDIdx index of the row to return.
   * @return Reference to the row. Is `const` if the matrix has column compression.
   */
  returned_row_type& get_row(id_index rowIndex);
  /**
   * @brief Only available for chain matrices and matrices with column compression.
   * Returns the row at the given IDIdx index.
   *
   * @warning The @ref get_column_index method of the row cells returns IDIdx indices for boundary matrices and
   * MatIdx indices for chain matrices.
   * 
   * @param rowIndex IDIdx index of the row to return.
   * @return Const reference to the row.
   */
  const Row_type& get_row(id_index rowIndex) const;
  /**
   * @brief Only available for RU matrices with position indexing. 
   * Returns the row at the given IDIdx index in R if @p inR is true and in U if @p inR is false.
   *
   * @warning The @ref get_column_index method of the row cells returns IDIdx indices for boundary matrices and
   * MatIdx indices for chain matrices.
   * 
   * @param rowIndex IDIdx index of the row to return.
   * @param inR If true, returns the row in R, if false, returns the row in U.
   * @return Const reference to the row.
   */
  const Row_type& get_row(id_index rowIndex, bool inR);

  /**
   * @brief Only available for base matrices and if @ref has_map_column_container is true.
   * For other matrices, see @ref remove_maximal_face.
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
   *    - @ref has_map_column_container and has_column_and_row_swaps are true: cleans up maps used for the lazy row swaps.
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
   * The emptiness of a row is therefore not tested at each column cell removal. 
   * 
   * @param rowIndex IDIdx index of the row to remove.
   */
  void erase_row(id_index rowIndex);
  // boundary: update barcode if already computed, does not verify if it really was maximal
  // ru
  // chain
  // id to pos
  // pos to id
  //TODO: for chain matrices, replace IDIdx input with MatIdx input to homogenise.
  /**
   * @brief Only available for boundary and chain matrices and if @ref has_map_column_container is true.
   * For base matrices, see @ref remove_column.
   * Assumes that the face is maximal in the current complex and removes it. The maximality of the face is not verified.
   * Also updates the barcode if it was computed.
   * 
   * @param columnIndex If boundary matrix, MatIdx index of the face to remove, otherwise the IDIdx index.
   */
  void remove_maximal_face(index columnIndex);
  void remove_maximal_face(id_index faceIndex, const std::vector<index>& columnsToSwap);
  void remove_last();

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
  index get_number_of_columns() const;
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
  std::enable_if_t<std::is_integral_v<Index_type>> multiply_target_and_add_to(Index_type sourceColumnIndex,
                                                                              int coefficient,
                                                                              Index_type targetColumnIndex);
  // base
  // base comp
  template <class Cell_range>
  std::enable_if_t<!std::is_integral_v<Cell_range>> multiply_target_and_add_to(const Cell_range& sourceColumn,
                                                                               int coefficient,
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
  std::enable_if_t<std::is_integral_v<Index_type>> multiply_source_and_add_to(int coefficient,
                                                                              Index_type sourceColumnIndex,
                                                                              Index_type targetColumnIndex);
  // base
  // base comp
  template <class Cell_range>
  std::enable_if_t<!std::is_integral_v<Cell_range>> multiply_source_and_add_to(int coefficient,
                                                                               const Cell_range& sourceColumn,
                                                                               index targetColumnIndex);

  // base
  // boundary: avoid calling with pairing option or make it such that it makes sense for persistence
  // ru: inR = true forced, avoid calling with specialized options or make it such that it makes sense for persistence
  // id to pos
  void zero_cell(index columnIndex, id_index rowIndex);
  // ru
  void zero_cell(index columnIndex, id_index rowIndex, bool inR);
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
  bool is_zero_cell(index columnIndex, id_index rowIndex);
  // ru
  bool is_zero_cell(index columnIndex, id_index rowIndex, bool inR) const;
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
  index get_column_with_pivot(id_index faceIndex) const;
  // boundary
  // ru: only in R
  // chain
  // id to pos
  // pos to id
  id_index get_pivot(index columnIndex);

  Matrix& operator=(Matrix other);
  friend void swap(Matrix& matrix1, Matrix& matrix2) { 
    swap(matrix1.matrix_, matrix2.matrix_);
    std::swap(matrix1.operators_, matrix2.operators_);
    std::swap(matrix1.cellPool_, matrix2.cellPool_);
    // matrix1.matrix_.set_operators(&matrix1.operators_);
    // matrix2.matrix_.set_operators(&matrix2.operators_);
  }

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
  // ru: returns true if barcode was changed
  // pos to id
  bool vine_swap_with_z_eq_1_case(pos_index index);  // by column position with ordered column container
  // chain: returns index which was not modified, ie new i+1
  // id to pos
  index vine_swap_with_z_eq_1_case(index columnIndex1,
                                   index columnIndex2);  // by column id with potentielly unordered column container
  // ru: returns true if barcode was changed
  // pos to id
  bool vine_swap(pos_index index);  // by column position with ordered column container
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
  using matrix_type = 
    typename std::conditional<
        isNonBasic,
        typename std::conditional<
            Options::is_of_boundary_type,
            typename std::conditional<
                Options::has_vine_update || Options::can_retrieve_representative_cycles,
                typename std::conditional<
                    Options::column_indexation_type == Column_indexation_types::CONTAINER || 
                      Options::column_indexation_type == Column_indexation_types::POSITION,
                    RU_matrix_type, Id_to_index_overlay<RU_matrix_type, Matrix<Options> >
                >::type,
                typename std::conditional<
                    Options::column_indexation_type == Column_indexation_types::CONTAINER || 
                      Options::column_indexation_type == Column_indexation_types::POSITION,
                    Boundary_matrix_type,
                    Id_to_index_overlay<Boundary_matrix_type, Matrix<Options> >
                >::type
            >::type,
            typename std::conditional<
                Options::column_indexation_type == Column_indexation_types::CONTAINER, 
                Chain_matrix_type,
                typename std::conditional<
                    Options::column_indexation_type == Column_indexation_types::POSITION,
                    Position_to_index_overlay<Chain_matrix_type, Matrix<Options> >,
                    Id_to_index_overlay<Chain_matrix_type, Matrix<Options> >
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
inline Matrix<Options>::Matrix() : operators_(new Field_operators()), cellPool_(new Cell_constructor()), matrix_(operators_, cellPool_)
{
  static_assert(
      Options::is_of_boundary_type || !Options::has_vine_update || Options::has_column_pairings,
      "When no barcode is recorded with vine swaps, comparaison functions for the columns have to be provided.");
  _assert_options();
}

template <class Options>
template <class Container_type>
inline Matrix<Options>::Matrix(const std::vector<Container_type>& columns, characteristic_type characteristic)
    : operators_(new Field_operators(characteristic)), cellPool_(new Cell_constructor()), matrix_(columns, operators_, cellPool_)
{
  static_assert(Options::is_of_boundary_type || !Options::has_vine_update || Options::has_column_pairings,
                "When no barcode is recorded with vine swaps for chain matrices, comparaison functions for the columns "
                "have to be provided.");
  _assert_options();
}

template <class Options>
inline Matrix<Options>::Matrix(int numberOfColumns, characteristic_type characteristic)
    : operators_(new Field_operators(characteristic)), cellPool_(new Cell_constructor()), matrix_(numberOfColumns, operators_, cellPool_)
{
  static_assert(Options::is_of_boundary_type || !Options::has_vine_update || Options::has_column_pairings,
                "When no barcode is recorded with vine swaps for chain matrices, comparaison functions for the columns "
                "have to be provided.");
  _assert_options();
}

template <class Options>
template <typename EventComparatorFunction>
inline Matrix<Options>::Matrix(EventComparatorFunction&& birthComparator, 
                               EventComparatorFunction&& deathComparator)
    : operators_(new Field_operators()), cellPool_(new Cell_constructor()), matrix_(operators_, cellPool_, birthComparator, deathComparator)
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
    : operators_(new Field_operators(characteristic)), cellPool_(new Cell_constructor()), matrix_(orderedBoundaries, operators_, cellPool_, birthComparator, deathComparator)
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
    : operators_(new Field_operators(characteristic)), cellPool_(new Cell_constructor()), matrix_(numberOfColumns, operators_, cellPool_, birthComparator, deathComparator)
{
  static_assert(
      !Options::is_of_boundary_type && Options::has_vine_update && !Options::has_column_pairings,
      "Constructor only available for chain matrices when vine swaps are enabled, but barcodes are not recorded.");
  _assert_options();
}

template <class Options>
inline Matrix<Options>::Matrix(const Matrix& matrixToCopy) 
    : operators_(new Field_operators(matrixToCopy.operators_->get_characteristic())), cellPool_(new Cell_constructor()), matrix_(matrixToCopy.matrix_, operators_, cellPool_)
{
  _assert_options();
}

template <class Options>
inline Matrix<Options>::Matrix(Matrix&& other) noexcept 
    : operators_(std::exchange(other.operators_, nullptr)), cellPool_(std::exchange(other.cellPool_, nullptr)), matrix_(std::move(other.matrix_))
{
  //TODO: verify that the address of operators_ == address of other.operators_ after move
  //and that therefore the addresses stored in matrix_ are correct.
  _assert_options();
}

template <class Options>
inline Matrix<Options>::~Matrix(){
  matrix_.reset(operators_, cellPool_);
  delete cellPool_;
  delete operators_;
}

template <class Options>
inline void Matrix<Options>::set_characteristic(characteristic_type characteristic){
  static_assert(!Options::is_z2, "The characteristic is definitely set to 2.");

  if (operators_->get_characteristic() != 0) {
    std::cerr << "Warning: Characteristic already initialised. Changing it could lead to incoherences in the matrice "
                 "as the modulo was already applied to values in existing columns.";
  }

  operators_->set_characteristic(characteristic);
}

template <class Options>
template <class Container_type>
inline void Matrix<Options>::insert_column(const Container_type& column) {
  assert(operators_->get_characteristic() != 0 &&
         "Columns cannot be initialized if the coefficient field characteristic is not specified.");
  static_assert(
      !isNonBasic,
      "'insert_column' not available for the chosen options. The input has to be in the form of a face boundary.");
  matrix_.insert_column(column);
}

template <class Options>
template <class Container_type>
inline void Matrix<Options>::insert_column(const Container_type& column, index columnIndex) {
 assert(operators_->get_characteristic() != 0 &&
         "Columns cannot be initialized if the coefficient field characteristic is not specified.");
  static_assert(!isNonBasic && !Options::has_column_compression,
                "'insert_column' with those parameters is not available for the chosen options.");
  matrix_.insert_column(column, columnIndex);
}

template <class Options>
template <class Boundary_type>
inline typename Matrix<Options>::insertion_return_type Matrix<Options>::insert_boundary(const Boundary_type& boundary,
                                                                                        dimension_type dim) {
  assert(operators_->get_characteristic() != 0 &&
         "Columns cannot be initialized if the coefficient field characteristic is not specified.");
  if constexpr (isNonBasic && !Options::is_of_boundary_type) return matrix_.insert_boundary(boundary, dim);
  else matrix_.insert_boundary(boundary, dim);
}

template <class Options>
template <class Boundary_type>
inline typename Matrix<Options>::insertion_return_type Matrix<Options>::insert_boundary(id_index faceIndex,
                                                                                        const Boundary_type& boundary,
                                                                                        dimension_type dim) {
  assert(operators_->get_characteristic() != 0 &&
         "Columns cannot be initialized if the coefficient field characteristic is not specified.");
  static_assert(isNonBasic, "Only enabled for non-basic matrices.");
  if constexpr (!Options::is_of_boundary_type) return matrix_.insert_boundary(faceIndex, boundary, dim);
  else matrix_.insert_boundary(faceIndex, boundary, dim);
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
  //TODO: I don't think there is a particular reason why the indexation is forced, should be removed.
  static_assert(isNonBasic && 
                    Options::is_of_boundary_type &&
                    (Options::has_vine_update || Options::can_retrieve_representative_cycles) &&
                    Options::column_indexation_type != Column_indexation_types::IDENTIFIER,
                "Only enabled for position indexed RU matrices.");

  return matrix_.get_column(columnIndex, inR);
}

template <class Options>
inline typename Matrix<Options>::returned_row_type& Matrix<Options>::get_row(id_index rowIndex) {
  static_assert(Options::has_row_access, "'get_row' is not available for the chosen options.");

  return matrix_.get_row(rowIndex);
}

template <class Options>
inline const typename Matrix<Options>::Row_type& Matrix<Options>::get_row(id_index rowIndex) const {
  static_assert(Options::has_row_access, "'get_row' is not available for the chosen options.");

  return matrix_.get_row(rowIndex);
}

template <class Options>
inline const typename Matrix<Options>::Row_type& Matrix<Options>::get_row(id_index rowIndex, bool inR) {
  static_assert(Options::has_row_access, "'get_row' is not available for the chosen options.");
  //TODO: I don't think there is a particular reason why the indexation is forced, should be removed.
  static_assert(isNonBasic && Options::is_of_boundary_type &&
                    (Options::has_vine_update || Options::can_retrieve_representative_cycles) &&
                    Options::column_indexation_type != Column_indexation_types::IDENTIFIER,
                "Only enabled for position indexed RU matrices.");

  return matrix_.get_row(rowIndex, inR);
}

template <class Options>
inline void Matrix<Options>::remove_column(index columnIndex) {
  static_assert(Options::has_map_column_container && !isNonBasic && !Options::has_column_compression,
                "'remove_column' is not available for the chosen options.");

  matrix_.remove_column(columnIndex);
}

template <class Options>
inline void Matrix<Options>::erase_row(id_index rowIndex) {
  static_assert(!isNonBasic || Options::has_removable_rows, "'erase_row' is not available for the chosen options.");

  matrix_.erase_row(rowIndex);
}

template <class Options>
inline void Matrix<Options>::remove_maximal_face(index columnIndex) {
  static_assert(Options::has_removable_columns, 
                "'remove_maximal_face(id_index)' is not available for the chosen options.");
  static_assert(isNonBasic && Options::has_vine_update, 
                "'remove_maximal_face(id_index)' is not available for the chosen options.");
  static_assert(Options::is_of_boundary_type || (Options::has_map_column_container && Options::has_column_pairings),
                "'remove_maximal_face(id_index)' is not available for the chosen options.");

  matrix_.remove_maximal_face(columnIndex);
}

template <class Options>
inline void Matrix<Options>::remove_maximal_face(id_index faceIndex, const std::vector<index>& columnsToSwap) {
  static_assert(Options::has_removable_columns, 
                "'remove_maximal_face(id_index,const std::vector<index>&)' is not available for the chosen options.");
  static_assert(isNonBasic && !Options::is_of_boundary_type, 
                "'remove_maximal_face(id_index,const std::vector<index>&)' is not available for the chosen options.");
  static_assert(Options::has_map_column_container && Options::has_vine_update, 
                "'remove_maximal_face(id_index,const std::vector<index>&)' is not available for the chosen options.");

  matrix_.remove_maximal_face(faceIndex, columnsToSwap);
}

template <class Options>
inline void Matrix<Options>::remove_last() {
  static_assert(Options::has_removable_columns || !isNonBasic, 
                "'remove_last' is not available for the chosen options.");
  static_assert(!Options::has_column_compression || isNonBasic,
                "'remove_last' is not available for the chosen options.");
  static_assert(Options::is_of_boundary_type || Options::has_map_column_container || !Options::has_vine_update,
                "'remove_last' is not available for the chosen options.");

  matrix_.remove_last();
}

template <class Options>
inline typename Matrix<Options>::dimension_type Matrix<Options>::get_max_dimension() const {
  static_assert(isNonBasic, "'get_max_dimension' is not available for the chosen options.");

  return matrix_.get_max_dimension();
}

template <class Options>
inline typename Matrix<Options>::index Matrix<Options>::get_number_of_columns() const {
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
inline std::enable_if_t<std::is_integral_v<Index_type>> Matrix<Options>::multiply_target_and_add_to(
    Index_type sourceColumnIndex, int coefficient, Index_type targetColumnIndex) 
{
  if constexpr (Options::is_z2){
    return matrix_.multiply_target_and_add_to(sourceColumnIndex, coefficient % 2, targetColumnIndex); //will be converted to bool
  } else {
    return matrix_.multiply_target_and_add_to(sourceColumnIndex, operators_->get_value(coefficient), targetColumnIndex);
  }
}

template <class Options>
template <class Cell_range>
inline std::enable_if_t<!std::is_integral_v<Cell_range>> Matrix<Options>::multiply_target_and_add_to(
    const Cell_range& sourceColumn, int coefficient, index targetColumnIndex) 
{
  static_assert(!isNonBasic,
                "For boundary or chain matrices, only additions with columns inside the matrix is allowed to maintain "
                "algebraic consistency.");

  if constexpr (Options::is_z2){
    //coef will be converted to bool, because of element_type
    return matrix_.multiply_target_and_add_to(sourceColumn, coefficient % 2, targetColumnIndex);
  } else {
    return matrix_.multiply_target_and_add_to(sourceColumn, operators_->get_value(coefficient), targetColumnIndex);
  }
}

template <class Options>
template <typename Index_type>
inline std::enable_if_t<std::is_integral_v<Index_type>> Matrix<Options>::multiply_source_and_add_to(
    int coefficient, Index_type sourceColumnIndex, Index_type targetColumnIndex) 
{
  if constexpr (Options::is_z2){
    //coef will be converted to bool, because of element_type
    return matrix_.multiply_source_and_add_to(coefficient % 2, sourceColumnIndex, targetColumnIndex);
  } else {
    return matrix_.multiply_source_and_add_to(operators_->get_value(coefficient), sourceColumnIndex, targetColumnIndex);
  }
}

template <class Options>
template <class Cell_range>
inline std::enable_if_t<!std::is_integral_v<Cell_range>> Matrix<Options>::multiply_source_and_add_to(
    int coefficient, const Cell_range& sourceColumn, index targetColumnIndex) 
{
  static_assert(!isNonBasic,
                "For boundary or chain matrices, only additions with columns inside the matrix is allowed to maintain "
                "algebraic consistency.");

  if constexpr (Options::is_z2){
    //coef will be converted to bool, because of element_type
    return matrix_.multiply_source_and_add_to(coefficient % 2, sourceColumn, targetColumnIndex);
  } else {
    return matrix_.multiply_source_and_add_to(operators_->get_value(coefficient), sourceColumn, targetColumnIndex);
  }
}

template <class Options>
inline void Matrix<Options>::zero_cell(index columnIndex, id_index rowIndex) {
  static_assert(Options::is_of_boundary_type && !Options::has_column_compression,
                "'zero_cell' is not available for the chosen options.");

  return matrix_.zero_cell(columnIndex, rowIndex);
}

template <class Options>
inline void Matrix<Options>::zero_cell(index columnIndex, id_index rowIndex, bool inR) {
  //TODO: I don't think there is a particular reason why the indexation is forced, should be removed.
  static_assert(isNonBasic && 
                    Options::is_of_boundary_type &&
                    (Options::has_vine_update || Options::can_retrieve_representative_cycles) &&
                    Options::column_indexation_type != Column_indexation_types::IDENTIFIER,
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
  //TODO: I don't think there is a particular reason why the indexation is forced, should be removed.
  static_assert(isNonBasic && Options::is_of_boundary_type &&
                    (Options::has_vine_update || Options::can_retrieve_representative_cycles) &&
                    Options::column_indexation_type != Column_indexation_types::IDENTIFIER,
                "Only enabled for RU matrices.");

  return matrix_.zero_column(columnIndex, inR);
}

template <class Options>
inline bool Matrix<Options>::is_zero_cell(index columnIndex, id_index rowIndex) {
  return matrix_.is_zero_cell(columnIndex, rowIndex);
}

template <class Options>
inline bool Matrix<Options>::is_zero_cell(index columnIndex, id_index rowIndex, bool inR) const {
  //TODO: I don't think there is a particular reason why the indexation is forced, should be removed.
  static_assert(isNonBasic && Options::is_of_boundary_type &&
                    (Options::has_vine_update || Options::can_retrieve_representative_cycles) &&
                    Options::column_indexation_type != Column_indexation_types::IDENTIFIER,
                "Only enabled for RU matrices.");

  return matrix_.is_zero_cell(columnIndex, rowIndex, inR);
}

template <class Options>
inline bool Matrix<Options>::is_zero_column(index columnIndex) {
  return matrix_.is_zero_column(columnIndex);
}

template <class Options>
inline bool Matrix<Options>::is_zero_column(index columnIndex, bool inR) {
  //TODO: I don't think there is a particular reason why the indexation is forced, should be removed.
  static_assert(isNonBasic && Options::is_of_boundary_type &&
                    (Options::has_vine_update || Options::can_retrieve_representative_cycles) &&
                    Options::column_indexation_type != Column_indexation_types::IDENTIFIER,
                "Only enabled for RU matrices.");

  return matrix_.is_zero_column(columnIndex, inR);
}

template <class Options>
inline typename Matrix<Options>::index Matrix<Options>::get_column_with_pivot(id_index faceIndex) const {
  static_assert(isNonBasic && (!Options::is_of_boundary_type ||
                               (Options::has_vine_update || Options::can_retrieve_representative_cycles)),
                "'get_column_with_pivot' is not available for the chosen options.");

  return matrix_.get_column_with_pivot(faceIndex);
}

template <class Options>
inline typename Matrix<Options>::id_index Matrix<Options>::get_pivot(index columnIndex) {
  static_assert(isNonBasic, "'get_pivot' is not available for the chosen options.");

  return matrix_.get_pivot(columnIndex);
}

template <class Options>
inline Matrix<Options>& Matrix<Options>::operator=(Matrix other) {
  swap(matrix_, other.matrix_);
  std::swap(operators_, other.operators_);
  std::swap(cellPool_, other.cellPool_);
//   if constexpr (!Options::is_z2) matrix_.set_operators(&operators_);

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
                    (isNonBasic && Options::is_of_boundary_type && 
                     !Options::has_vine_update && !Options::can_retrieve_representative_cycles),
                "This method was not enabled.");
  return matrix_.swap_columns(columnIndex1, columnIndex2);
}

template <class Options>
inline void Matrix<Options>::swap_rows(index rowIndex1, index rowIndex2) {
  static_assert((!isNonBasic && !Options::has_column_compression) ||
                    (isNonBasic && Options::is_of_boundary_type && 
                     !Options::has_vine_update && !Options::can_retrieve_representative_cycles),
                "This method was not enabled.");
  return matrix_.swap_rows(rowIndex1, rowIndex2);
}

template <class Options>
inline bool Matrix<Options>::vine_swap_with_z_eq_1_case(pos_index index) {
  static_assert(Options::has_vine_update, 
                "This method was not enabled.");
  static_assert(Options::column_indexation_type == Column_indexation_types::POSITION || 
                (Options::is_of_boundary_type && 
                Options::column_indexation_type == Column_indexation_types::CONTAINER), 
                "This method was not enabled.");
  return matrix_.vine_swap_with_z_eq_1_case(index);
}

template <class Options>
inline typename Matrix<Options>::index Matrix<Options>::vine_swap_with_z_eq_1_case(index columnIndex1,
                                                                                   index columnIndex2) {
  static_assert(Options::has_vine_update, 
                "This method was not enabled.");
  static_assert(Options::column_indexation_type == Column_indexation_types::IDENTIFIER || 
                (!Options::is_of_boundary_type && 
                Options::column_indexation_type == Column_indexation_types::CONTAINER), 
                "This method was not enabled.");
  return matrix_.vine_swap_with_z_eq_1_case(columnIndex1, columnIndex2);
}

template <class Options>
inline bool Matrix<Options>::vine_swap(pos_index index) {
  static_assert(Options::has_vine_update, 
                "This method was not enabled.");
  static_assert(Options::column_indexation_type == Column_indexation_types::POSITION || 
                (Options::is_of_boundary_type && 
                Options::column_indexation_type == Column_indexation_types::CONTAINER), 
                "This method was not enabled.");
  return matrix_.vine_swap(index);
}

template <class Options>
inline typename Matrix<Options>::index Matrix<Options>::vine_swap(index columnIndex1, index columnIndex2) {
  static_assert(Options::has_vine_update, 
                "This method was not enabled.");
  static_assert(Options::column_indexation_type == Column_indexation_types::IDENTIFIER || 
                (!Options::is_of_boundary_type && 
                Options::column_indexation_type == Column_indexation_types::CONTAINER), 
                "This method was not enabled.");
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
//   static_assert(!Options::has_column_compression || !Options::has_map_column_container,
//                 "When column compression is used, the removal of columns is not implemented yet.");
}

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // MASTER_MATRIX_H

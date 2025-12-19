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
 * @file Matrix.h
 * @author Hannah Schreiber
 * @brief Contains @ref Gudhi::persistence_matrix::Matrix class.
 */

#ifndef MASTER_MATRIX_H
#define MASTER_MATRIX_H

#include <type_traits>
#include <vector>
#include <unordered_map>
#include <map>
#include <initializer_list>

#include <boost/intrusive/list.hpp>

#include <gudhi/Debug_utils.h>

#include <gudhi/persistence_matrix_options.h>
#include <gudhi/persistence_interval.h>

#include <gudhi/Fields/Z2_field_operators.h>

#include <gudhi/Persistence_matrix/Id_to_index_overlay.h>
#include <gudhi/Persistence_matrix/Position_to_index_overlay.h>

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

#include <gudhi/Persistence_matrix/Base_matrix.h>
#include <gudhi/Persistence_matrix/Base_matrix_with_column_compression.h>
#include <gudhi/Persistence_matrix/Boundary_matrix.h>
#include <gudhi/Persistence_matrix/RU_matrix.h>
#include <gudhi/Persistence_matrix/Chain_matrix.h>

#include <gudhi/Persistence_matrix/allocators/entry_constructors.h>
#include <gudhi/Persistence_matrix/columns/entry_types.h>
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

/// Gudhi namespace.
namespace Gudhi {
/// Persistence matrix namespace.
namespace persistence_matrix {

/**
 * @class Matrix matrix.h gudhi/matrix.h
 * @ingroup persistence_matrix
 *
 * @brief Data structure for matrices, and in particular thought for matrices representing filtered complexes
 * in order to compute persistence and/or representative cycles.
 *
 * @anchor mp_matrices __%Matrix types:__
 *
 * The are roughly three types of matrices available and one is selected automatically depending on the options chosen:
 * - @anchor basematrix a @ref basematrix "basic matrix" which can represent any matrix and therefore will not make any
 *   assumption on its content. It is the only matrix type with the option of column compression (as it is the only one
 *   where it makes sense). This type is chosen by default when none of the homology related options are set to true:
 *   @ref PersistenceMatrixOptions::has_column_pairings, @ref PersistenceMatrixOptions::has_vine_update and
 *   @ref PersistenceMatrixOptions::can_retrieve_representative_cycles.
 * - @anchor boundarymatrix a @ref boundarymatrix "boundary matrix" @f$ B = R \cdot U @f$ which either stores only
 *   @f$ R @f$ or the whole decomposition @f$ R @f$ and @f$ U @f$ depending on the options. This type is selected if
 *   @ref PersistenceMatrixOptions::is_of_boundary_type is set to true and at least one of the following options:
 *   @ref PersistenceMatrixOptions::has_column_pairings, @ref PersistenceMatrixOptions::has_vine_update and
 *   @ref PersistenceMatrixOptions::can_retrieve_representative_cycles. If only
 *   @ref PersistenceMatrixOptions::has_column_pairings is true, then only @f$ R @f$ is stored, but if either
 *   @ref PersistenceMatrixOptions::has_vine_update or @ref PersistenceMatrixOptions::can_retrieve_representative_cycles
 *   is true, then @f$ U @f$ also needs to be stored. Note that the option
 *   @ref PersistenceMatrixOptions::column_indexation_type will produce a small overhead when set to
 *   @ref Column_indexation_types::IDENTIFIER.
 * - @anchor chainmatrix a @ref chainmatrix "chain complex matrix" representing a "compatible base" of a filtered chain
 *   complex (see @cite zigzag). This matrix is deduced from the boundary matrix and therefore encodes more or less the
 *   same information but differently and can therefore be better suited for certain applications. This type can be used
 *   the same way than the precedent type, only the option @ref PersistenceMatrixOptions::is_of_boundary_type has to be
 *   set to false. Note that the option @ref PersistenceMatrixOptions::column_indexation_type will produce a small
 *   overhead when set to  @ref Column_indexation_types::POSITION or  @ref Column_indexation_types::IDENTIFIER.
 *
 * @anchor mp_indexation __Indexation scheme:__
 *
 * The indexation system for columns of the different matrix types can be a bit tricky and different methods will not
 * always take the same type of index as input (for optimization purposes). So, to avoid confusion, we will name and
 * define here the different possibilities, such that we can directly refer to it in the descriptions of the methods.
 * Note that every column and row in a @ref boundarymatrix "boundary" or @ref chainmatrix "chain matrix" is always
 * associated to a single simplex/cell, so in order to avoid repeating formulations like "of the simplex associated to
 * the column" all the time, we will amalgamate both notions together.
 *
 * Let @f$ c @f$ be a column.
 * - @anchor MatIdx @ref MatIdx "MatIdx": This will correspond to the position of @f$ c @f$ in the matrix, i.e.,
 *   @f$ underlying\_container[MatIdx] = c @f$. This will be the only public indexing scheme for
 *   @ref basematrix "basic matrices".
 * - @anchor PosIdx @ref PosIdx "PosIdx": This will correspond to the relative position of @f$ c @f$ in the current
 *   filtration compared to the other columns, starting the count at 0. For @ref boundarymatrix "boundary matrices",
 *   @ref PosIdx will always be equal to @ref MatIdx, but this is not true for @ref chainmatrix "chain matrices" when
 *   swaps or removals were performed.
 * - @anchor IDIdx @ref IDIdx "IDIdx": This will correspond to the ID of @f$ c @f$ in the complex used to identify it in
 *   the boundaries. If at the insertion of @f$ c @f$, its ID was not specified and it was the @f$ n^{th} @f$ insertion,
 *   it is assumed that the ID is @f$ n @f$ (which means that @ref IDIdx and @ref PosIdx will only start to differ when
 *   swaps or removals are performed). If an ID is specified at the insertion of @f$ c @f$, the ID is stored as the
 *   @ref IDIdx of @f$ c @f$. IDs can be freely chosen with the only restriction that they have to be strictly
 *   increasing in the order of the filtration at initialisation.
 *
 * In conclusion, with default values, if no vine swaps or removals occurs, all three indexing schemes are the same.
 *
 * @anchor rowindex Let @f$ r @f$ be a row. Rows are indexed in two ways depending only if the matrix is a
 * @ref chainmatrix "chain matrix" or not. If the matrix is a @ref chainmatrix "chain matrix", @f$ r @f$ is always
 * indexed by its ID, so it correspond to the @ref IDIdx indexing scheme. If the matrix is not a
 * @ref chainmatrix "chain matrix", @f$ r @f$ will originally also be indexed by the ID, but when a swap occurs,
 * the rows also swap IDs and the new ID has to be used to access @f$ r @f$. This means that when the default
 * @ref IDIdx scheme is used (the cells are numerated in order of appearance in the filtration starting at 0),
 * the indexation of the rows correspond to @ref PosIdx.
 *
 * @tparam PersistenceMatrixOptions Structure encoding all the options of the matrix.
 * See description of @ref PersistenceMatrixOptions for more details.
 */
template <class PersistenceMatrixOptions = Default_options<> >
class Matrix
{
 public:
  using Option_list = PersistenceMatrixOptions;                   // to make it accessible from the other classes
  using Index = typename PersistenceMatrixOptions::Index;         /**< Type of @ref MatIdx index. */
  using ID_index = typename PersistenceMatrixOptions::Index;      /**< Type of @ref IDIdx index. */
  using Pos_index = typename PersistenceMatrixOptions::Index;     /**< Type of @ref PosIdx index. */
  using Dimension = typename PersistenceMatrixOptions::Dimension; /**< Type for dimension value. */

  /**
   * @brief coefficients field type.
   */
  using Field_operators = std::conditional_t<PersistenceMatrixOptions::is_z2,
                                             Gudhi::persistence_fields::Z2_field_operators,
                                             typename PersistenceMatrixOptions::Field_coeff_operators>;
  /**
   * @brief Type of a field element.
   */
  using Element = typename Field_operators::Element;
  using Characteristic = typename Field_operators::Characteristic;

  /**
   * @brief Returns value from a type when not set.
   */
  template <typename T>
  static constexpr T get_null_value()
  {
    return -1;
  }

  /**
   * @brief Type for a bar in the computed barcode. Stores the birth, death and dimension of the bar.
   */
  using Bar = Persistence_interval<Dimension, Pos_index>;

  // tags for boost to associate a row and a column to a same entry
  struct Matrix_row_tag;
  struct Matrix_column_tag;

  using Base_hook_matrix_row =
      boost::intrusive::list_base_hook<boost::intrusive::tag<Matrix_row_tag>,
                                       boost::intrusive::link_mode<boost::intrusive::auto_unlink> >;
  using Base_hook_matrix_list_column =
      boost::intrusive::list_base_hook<boost::intrusive::tag<Matrix_column_tag>,
                                       boost::intrusive::link_mode<boost::intrusive::safe_link> >;
  using Base_hook_matrix_set_column =
      boost::intrusive::set_base_hook<boost::intrusive::tag<Matrix_column_tag>,
                                      boost::intrusive::link_mode<boost::intrusive::safe_link> >;

  // Two dummies are necessary to avoid double inheritance as an entry can inherit both a row and a column hook.
  struct Dummy_row_hook {};

  struct Dummy_column_hook {};

  using Row_hook =
      std::conditional_t<PersistenceMatrixOptions::has_row_access && PersistenceMatrixOptions::has_intrusive_rows,
                         Base_hook_matrix_row,
                         Dummy_row_hook>;
  using Column_hook =
      std::conditional_t<PersistenceMatrixOptions::column_type == Column_types::INTRUSIVE_LIST,
                         Base_hook_matrix_list_column,
                         std::conditional_t<PersistenceMatrixOptions::column_type == Column_types::INTRUSIVE_SET,
                                            Base_hook_matrix_set_column,
                                            Dummy_column_hook> >;

  // Option to store the column index within the entry (additionally to the row index). Necessary only with row access.
  using Entry_column_index_option = std::conditional_t<PersistenceMatrixOptions::has_row_access,
                                                       Entry_column_index<Index>,
                                                       Dummy_entry_column_index_mixin>;
  // Option to store the value of the entry.
  // Unnecessary for values in Z_2 as there are always 1 (0-valued entries are never stored).
  using Entry_field_element_option = std::
      conditional_t<PersistenceMatrixOptions::is_z2, Dummy_entry_field_element_mixin, Entry_field_element<Element> >;
  /**
   * @brief Type of a matrix entry. See @ref Entry for a more detailed description.
   */
  using Matrix_entry = Entry<Matrix<PersistenceMatrixOptions> >;

  /**
   * @brief Default entry constructor/destructor, using classic new and delete.
   * For now, only used as default value for columns constructed independently outside of the matrix by the user.
   * Could be used in the future when parallel options are implemented, as usual pools are not thread safe.
   */
  inline static New_entry_constructor<Matrix_entry> defaultEntryConstructor;
  /**
   * @brief Entry constructor/destructor used by the matrix. Uses a pool of entries to accelerate memory management,
   * as entries are constructed and destroyed a lot during reduction, swaps or additions.
   */
  using Entry_constructor = Pool_entry_constructor<Matrix_entry>;

  /**
   * @brief Type used to identify an entry, for example when inserting a boundary.
   * If @ref PersistenceMatrixOptions::is_z2 is true, the type is an @ref IDIdx and corresponds to the row index of the
   * entry (the entry value is assumed to be 1). If @ref PersistenceMatrixOptions::is_z2 is false, the type is a pair
   * whose first element is the row index of the entry and the second element is the value of the entry (which again is
   * assumed to be non-zero). The column index of the row is always deduced from the context in which the type is used.
   */
  using Entry_representative =
      std::conditional_t<PersistenceMatrixOptions::is_z2, ID_index, std::pair<ID_index, Element> >;

  /**
   * @private
   */
  static ID_index get_row_index(const Matrix_entry& e)
  {
    return e.get_row_index();
  }
  /**
   * @private
   */
  static ID_index get_row_index(const Entry_representative& e)
  {
    if constexpr (PersistenceMatrixOptions::is_z2) {
      return e;
    } else {
      return e.first;
    }
  }
  /**
   * @private
   */
  static ID_index& get_row_index(Entry_representative& e)
  {
    if constexpr (PersistenceMatrixOptions::is_z2) {
      return e;
    } else {
      return e.first;
    }
  }
  /**
   * @private
   */
  static Element get_element(const Matrix_entry& e)
  {
    return e.get_element();
  }
  /**
   * @private
   */
  static Element get_element(const Entry_representative& e)
  {
    if constexpr (PersistenceMatrixOptions::is_z2) {
      return Field_operators::get_multiplicative_identity();
    } else {
      return e.second;
    }
  }

  /**
   * @brief Compares two entries by their position in the row. They are assume to be in the same row.
   */
  struct RowEntryComp {
    bool operator()(const Matrix_entry& c1, const Matrix_entry& c2) const
    {
      return c1.get_column_index() < c2.get_column_index();
    }
  };

  /**
   * @brief Type of the rows stored in the matrix. Is either an intrusive list of @ref Matrix_entry (not ordered) if
   * @ref PersistenceMatrixOptions::has_intrusive_rows is true, or a set of @ref Matrix_entry (ordered by
   * column index) otherwise.
   */
  using Row = std::conditional_t<PersistenceMatrixOptions::has_intrusive_rows,
                                 boost::intrusive::list<Matrix_entry,
                                                        boost::intrusive::constant_time_size<false>,
                                                        boost::intrusive::base_hook<Base_hook_matrix_row> >,
                                 std::set<Matrix_entry, RowEntryComp> >;

  using Row_container =
      std::conditional_t<PersistenceMatrixOptions::has_removable_rows, std::map<ID_index, Row>, std::vector<Row> >;

  // Row access at column level
  using Row_access_option = std::conditional_t<PersistenceMatrixOptions::has_row_access,
                                               Row_access<Matrix<PersistenceMatrixOptions> >,
                                               Dummy_row_access>;
  // Row access at matrix level
  using Matrix_row_access_option =
      std::conditional_t<PersistenceMatrixOptions::has_row_access,
                         Matrix_row_access<Row, Row_container, PersistenceMatrixOptions::has_removable_rows, ID_index>,
                         Dummy_matrix_row_access>;

  template <typename value_type>
  using Dictionary = std::conditional_t<PersistenceMatrixOptions::has_map_column_container,
                                        std::unordered_map<unsigned int, value_type>,
                                        std::vector<value_type> >;

  static const bool isNonBasic = PersistenceMatrixOptions::has_column_pairings ||
                                 PersistenceMatrixOptions::has_vine_update ||
                                 PersistenceMatrixOptions::can_retrieve_representative_cycles;

  using Column_dimension_option = std::
      conditional_t<isNonBasic, Column_dimension_holder<Matrix<PersistenceMatrixOptions> >, Dummy_dimension_holder>;
  // Extra information needed for a column when the matrix is a @ref chainmatrix "chain matrix".
  using Chain_column_option = std::conditional_t<isNonBasic && !PersistenceMatrixOptions::is_of_boundary_type,
                                                 Chain_column_extra_properties<Matrix<PersistenceMatrixOptions> >,
                                                 Dummy_chain_properties>;

  using Matrix_heap_column = Heap_column<Matrix<PersistenceMatrixOptions> >;
  using Matrix_list_column = List_column<Matrix<PersistenceMatrixOptions> >;
  using Matrix_vector_column = Vector_column<Matrix<PersistenceMatrixOptions> >;
  using Matrix_naive_vector_column = Naive_std_vector_column<Matrix<PersistenceMatrixOptions> >;
  using Matrix_small_vector_column = Naive_small_vector_column<Matrix<PersistenceMatrixOptions> >;
  using Matrix_set_column = Set_column<Matrix<PersistenceMatrixOptions> >;
  using Matrix_unordered_set_column = Unordered_set_column<Matrix<PersistenceMatrixOptions> >;
  using Matrix_intrusive_list_column = Intrusive_list_column<Matrix<PersistenceMatrixOptions> >;
  using Matrix_intrusive_set_column = Intrusive_set_column<Matrix<PersistenceMatrixOptions> >;

  /**
   * @brief Type of the columns stored in the matrix. The type depends on the value of
   * @ref PersistenceMatrixOptions::column_type defined in the given options. See @ref Column_types for a more detailed
   * description. All columns follow the @ref PersistenceMatrixColumn concept.
   */
  using Column = std::conditional_t<
      PersistenceMatrixOptions::column_type == Column_types::HEAP,
      Matrix_heap_column,
      std::conditional_t<
          PersistenceMatrixOptions::column_type == Column_types::LIST,
          Matrix_list_column,
          std::conditional_t<
              PersistenceMatrixOptions::column_type == Column_types::SET,
              Matrix_set_column,
              std::conditional_t<
                  PersistenceMatrixOptions::column_type == Column_types::UNORDERED_SET,
                  Matrix_unordered_set_column,
                  std::conditional_t<
                      PersistenceMatrixOptions::column_type == Column_types::VECTOR,
                      Matrix_vector_column,
                      std::conditional_t<
                          PersistenceMatrixOptions::column_type == Column_types::INTRUSIVE_LIST,
                          Matrix_intrusive_list_column,
                          std::conditional_t<
                              PersistenceMatrixOptions::column_type == Column_types::NAIVE_VECTOR,
                              Matrix_naive_vector_column,
                              std::conditional_t<PersistenceMatrixOptions::column_type == Column_types::SMALL_VECTOR,
                                                 Matrix_small_vector_column,
                                                 Matrix_intrusive_set_column>
                              > > > > > > >;

  struct Column_z2_settings {
    Column_z2_settings() : entryConstructor() {}
    Column_z2_settings([[maybe_unused]] Characteristic characteristic) : entryConstructor() {}
    Column_z2_settings([[maybe_unused]] const Column_z2_settings& toCopy) : entryConstructor() {}
    Column_z2_settings([[maybe_unused]] Column_z2_settings&& toMove) noexcept : entryConstructor() {}
    ~Column_z2_settings() = default;
    Column_z2_settings& operator=([[maybe_unused]] const Column_z2_settings& other) { return *this; }
    Column_z2_settings& operator=([[maybe_unused]] Column_z2_settings&& other) noexcept { return *this; }

    Entry_constructor entryConstructor;  // will be replaced by more specific allocators depending on the column type.
  };

  struct Column_zp_settings {
    Column_zp_settings() : operators(), entryConstructor() {}
    Column_zp_settings(Characteristic characteristic) : operators(characteristic), entryConstructor() {}
    Column_zp_settings(const Column_zp_settings& toCopy) : operators(toCopy.operators), entryConstructor() {}
    Column_zp_settings(Column_zp_settings&& toMove) noexcept
        : operators(std::move(toMove.operators)), entryConstructor() {}
    ~Column_zp_settings() = default;
    Column_zp_settings& operator=(const Column_zp_settings& other)
    {
      operators = other.operators;
      return *this;
    }
    Column_zp_settings& operator=(Column_zp_settings&& other) noexcept
    {
      if (this == &other) return *this;

      operators = std::move(other.operators);
      return *this;
    }

    Field_operators operators;
    Entry_constructor entryConstructor;  // will be replaced by more specific allocators depending on the column type.
  };

  // struct Column_z2_with_rows_settings {
  //   Column_z2_with_rows_settings() : entryConstructor(), rows(nullptr) {}
  //   Column_z2_with_rows_settings([[maybe_unused]] Characteristic characteristic)
  //       : entryConstructor(), rows(nullptr) {}

  //   Entry_constructor entryConstructor;
  //   Row_container* rows;
  // };

  // struct Column_zp_with_rows_settings {
  //   Column_zp_with_rows_settings() : operators(), entryConstructor(), rows(nullptr) {}
  //   Column_zp_with_rows_settings(Characteristic characteristic)
  //       : operators(characteristic), entryConstructor(), rows(nullptr) {}

  //   Field_operators operators;
  //   Entry_constructor entryConstructor;
  //   Row_container* rows;
  // };

  // To prepare a more flexible use of the column types later (custom allocators depending on the column type etc.)
  using Column_settings = std::conditional_t<PersistenceMatrixOptions::is_z2, Column_z2_settings, Column_zp_settings>;

  /**
   * @private
   */
  static Field_operators const* get_operator_ptr(Column_settings const* colSettings)
  {
    if constexpr (PersistenceMatrixOptions::is_z2) {
      return nullptr;
    } else {
      if (colSettings == nullptr) return nullptr; // used for dummy columns
      return &(colSettings->operators);
    }
  }
  /**
   * @private
   */
  template <typename T>
  static Element get_coefficient_value(T v, [[maybe_unused]] Field_operators const* operators)
  {
    if constexpr (PersistenceMatrixOptions::is_z2) {
      return Field_operators::get_value(v);
    } else {
      return operators->get_value(v);
    }
  }
  /**
   * @private
   */
  static Element get_coefficient_value(bool v, [[maybe_unused]] Field_operators const* operators) { return v; }

  // using Column_settings = typename std::conditional<
  //     PersistenceMatrixOptions::is_z2,
  //     typename std::conditional<PersistenceMatrixOptions::has_row_access,
  //                               Column_z2_with_rows_settings,
  //                               Column_z2_settings
  //                              >::type,
  //     typename std::conditional<PersistenceMatrixOptions::has_row_access,
  //                               Column_zp_with_rows_settings,
  //                               Column_zp_settings
  //                              >::type
  //   >::type;

  using Column_container = std::conditional_t<PersistenceMatrixOptions::has_map_column_container,
                                              std::unordered_map<Index, Column>,
                                              std::vector<Column> >;

  static const bool hasFixedBarcode = Option_list::is_of_boundary_type && !PersistenceMatrixOptions::has_vine_update;
  /**
   * @brief Type of the computed barcode. It is either a list of @ref Matrix::Bar or a vector of @ref Matrix::Bar,
   * depending if bars need to be removed from the container as some point or not.
   */
  using Barcode = std::conditional_t<
      hasFixedBarcode,
      std::vector<Bar>,
      std::conditional_t<PersistenceMatrixOptions::has_removable_columns, std::list<Bar>, std::vector<Bar> > >;
  using Bar_dictionary =
      std::conditional_t<hasFixedBarcode,
                         std::conditional_t<PersistenceMatrixOptions::can_retrieve_representative_cycles,
                                            std::vector<Index>,                   // RU
                                            std::unordered_map<Pos_index, Index>  // boundary
                                            >,
                         std::conditional_t<PersistenceMatrixOptions::has_removable_columns,
                                            std::unordered_map<Pos_index, typename Barcode::iterator>,
                                            std::vector<Index> > >;

  // default type for boundaries to permit list initialization directly in function parameters
  using Boundary = std::conditional_t<PersistenceMatrixOptions::is_z2,
                                      std::initializer_list<ID_index>,
                                      std::initializer_list<std::pair<ID_index, Element> > >;

  // i.e. is simple @ref boundarymatrix "boundary matrix". Also, only needed because of the reduction algorithm.
  // TODO: remove the necessity and recalculate when needed or keep it like that?
  static const bool maxDimensionIsNeeded =
      PersistenceMatrixOptions::has_column_pairings && PersistenceMatrixOptions::is_of_boundary_type &&
      !PersistenceMatrixOptions::has_vine_update && !PersistenceMatrixOptions::can_retrieve_representative_cycles;

  using Matrix_dimension_option =
      std::conditional_t<PersistenceMatrixOptions::has_matrix_maximal_dimension_access || maxDimensionIsNeeded,
                         std::conditional_t<PersistenceMatrixOptions::has_removable_columns,
                                            Matrix_all_dimension_holder<Dimension>,
                                            Matrix_max_dimension_holder<Dimension> >,
                         Dummy_matrix_dimension_holder>;

  using Master_base_matrix = std::conditional_t<PersistenceMatrixOptions::has_column_compression,
                                                Base_matrix_with_column_compression<Matrix<PersistenceMatrixOptions> >,
                                                Base_matrix<Matrix<PersistenceMatrixOptions> > >;
  using Master_boundary_matrix = Boundary_matrix<Matrix<PersistenceMatrixOptions> >;
  using Master_RU_matrix = RU_matrix<Matrix<PersistenceMatrixOptions> >;
  using Master_chain_matrix = Chain_matrix<Matrix<PersistenceMatrixOptions> >;

  template <class Base>
  using Base_swap_option = std::conditional_t<PersistenceMatrixOptions::has_vine_update ||
                                                  PersistenceMatrixOptions::has_column_and_row_swaps,
                                              Base_swap<Matrix<PersistenceMatrixOptions>, Base>,
                                              Dummy_base_swap>;
  using Base_pairing_option =
      std::conditional_t<PersistenceMatrixOptions::has_column_pairings && !PersistenceMatrixOptions::has_vine_update &&
                             !PersistenceMatrixOptions::can_retrieve_representative_cycles,
                         Base_pairing<Matrix<PersistenceMatrixOptions> >,
                         Dummy_base_pairing>;

  using RU_pairing_option = std::conditional_t<PersistenceMatrixOptions::has_column_pairings,
                                               std::conditional_t<PersistenceMatrixOptions::has_vine_update,
                                                                  RU_barcode_swap<Matrix<PersistenceMatrixOptions> >,
                                                                  RU_pairing<Matrix<PersistenceMatrixOptions> > >,
                                               Dummy_ru_pairing>;
  using RU_vine_swap_option = std::conditional_t<PersistenceMatrixOptions::has_vine_update,
                                                 RU_vine_swap<Matrix<PersistenceMatrixOptions> >,
                                                 Dummy_ru_vine_swap>;
  using RU_representative_cycles_option =
      std::conditional_t<PersistenceMatrixOptions::can_retrieve_representative_cycles,
                         RU_representative_cycles<Matrix<PersistenceMatrixOptions> >,
                         Dummy_ru_representative_cycles>;

  using Chain_pairing_option =
      std::conditional_t<PersistenceMatrixOptions::has_column_pairings,
                         std::conditional_t<PersistenceMatrixOptions::has_vine_update,
                                            Chain_barcode_swap<Matrix<PersistenceMatrixOptions> >,
                                            Chain_pairing<Matrix<PersistenceMatrixOptions> > >,
                         Dummy_chain_pairing>;
  using Chain_vine_swap_option = std::conditional_t<PersistenceMatrixOptions::has_vine_update,
                                                    Chain_vine_swap<Matrix<PersistenceMatrixOptions> >,
                                                    Dummy_chain_vine_swap>;
  using Chain_representative_cycles_option =
      std::conditional_t<PersistenceMatrixOptions::can_retrieve_representative_cycles,
                         Chain_representative_cycles<Matrix<PersistenceMatrixOptions> >,
                         Dummy_chain_representative_cycles>;

  /**
   * @brief Type of a representative cycle. Vector of @ref rowindex "row indices".
   */
  using Cycle = std::vector<Entry_representative>;

  /**
   * @private
   */
  template <class EntryRange>
  static Cycle build_cycle_from_range(const EntryRange& entries)
  {
    Cycle cycle;
    if constexpr (RangeTraits<EntryRange>::has_size) {
      cycle.reserve(entries.size());
    }
    for (const auto& c : entries) {
      if constexpr (PersistenceMatrixOptions::is_z2) {
        cycle.push_back(c.get_row_index());
      } else {
        cycle.push_back({c.get_row_index(), c.get_element()});
      }
    }
    return cycle;
  }

  // Return types to factorize the corresponding methods

  // The returned column is `const` if the matrix uses column compression
  using Returned_column =
      std::conditional_t<!isNonBasic && PersistenceMatrixOptions::has_column_compression, const Column, Column>;
  // The returned row is `const` if the matrix uses column compression
  using Returned_row =
      std::conditional_t<!isNonBasic && PersistenceMatrixOptions::has_column_compression, const Row, Row>;
  // If the matrix is a chain matrix, the insertion method returns the pivots of its unpaired columns used to reduce the
  // inserted boundary. Otherwise, void.
  using Insertion_return =
      std::conditional_t<PersistenceMatrixOptions::is_of_boundary_type || !isNonBasic ||
                             PersistenceMatrixOptions::column_indexation_type == Column_indexation_types::POSITION,
                         void,
                         std::vector<Entry_representative> >;

  /**
   * @brief Default constructor. Initializes an empty matrix.
   */
  Matrix();
  /**
   * @brief Constructs a new matrix from the given ranges of @ref Entry_representative. Each range corresponds to a
   * column (the order of the ranges are preserved). The content of the ranges is assumed to be sorted by increasing
   * IDs. If the columns are representing a boundary matrix, the IDs of the simplices are also assumed to be
   * consecutive, ordered by filtration value, starting with 0.
   *
   * See @ref mp_matrices "matrix descriptions" for further details on how the given matrix is handled.
   *
   * @tparam Container Range type for @ref Entry_representative ranges. Assumed to have a begin(), end() and size()
   * method.
   * @param columns For a @ref basematrix "base matrix", the columns are copied as is. If options related to homology
   * are activated, @p columns is interpreted as a boundary matrix of a **simplicial** complex. In this case,
   * `columns[i]` should store the boundary of simplex `i` as an ordered list of indices of its facets (again those
   * indices correspond to their respective position in the matrix). Therefore the indices of the simplices are assumed
   * to be consecutive and starting with 0 (an empty boundary is interpreted as a vertex boundary and not as a non
   * existing simplex). All dimensions up to the maximal dimension of interest have to be present. If only a higher
   * dimension is of interest and not everything should be stored, then use the @ref insert_boundary method instead
   * (after creating the matrix with the @ref Matrix(unsigned int numberOfColumns, Characteristic characteristic)
   * constructor preferably). If the persistence barcode has to be computed from this matrix, the simplices are also
   * assumed to be ordered by appearance order in the filtration. Also, depending of the options, the matrix is
   * eventually reduced on the fly or converted into a chain complex base, so the new matrix is not always identical to
   * the old one.
   * @param characteristic Characteristic of the coefficient field. Has to be specified if
   * @ref PersistenceMatrixOptions::is_z2 is false. Default value is 11. Ignored if
   * @ref PersistenceMatrixOptions::is_z2 is true.
   */
  template <class Container = Boundary>
  Matrix(const std::vector<Container>& columns, Characteristic characteristic = 11);
  /**
   * @brief Constructs a new empty matrix and reserves space for the given number of columns.
   *
   * @param numberOfColumns Number of columns to reserve space for.
   * @param characteristic Characteristic of the coefficient field. If not specified and
   * @ref PersistenceMatrixOptions::is_z2 is false, the characteristic has to be set later with the use of
   * @ref set_characteristic before calling for the first time a method needing it. Ignored if
   * @ref PersistenceMatrixOptions::is_z2 is true.
   */
  Matrix(unsigned int numberOfColumns, Characteristic characteristic = Field_operators::nullCharacteristic);
  /**
   * @brief Constructs a new empty matrix with the given comparator functions. Only available when those comparators
   * are necessary.
   *
   * That is, when **all** following options have following values:
   *   - @ref PersistenceMatrixOptions::is_of_boundary_type = false
   *   - @ref PersistenceMatrixOptions::has_vine_update = true
   *   - @ref PersistenceMatrixOptions::has_column_pairings = false
   *
   * Those comparators are necessary to distinguish cases in a vine update. When the matrix is of
   * @ref boundarymatrix "boundary type" or if the column pairing is activated (i.e., the barcode is stored),
   * the comparators can be easily deduced without overhead. If neither are true, we assume that one has additional
   * information outside of the matrix about the barcode to provide a better suited comparator adapted to the situation
   * (as in the implementation of the Zigzag algorithm @cite zigzag for example.)
   *
   * @param birthComparator Method taking two @ref PosIdx indices as parameter and returns true if and only if the first
   * cell is associated to a bar with strictly smaller birth than the bar associated to the second one.
   * @param deathComparator Method taking two @ref PosIdx indices as parameter and returns true if and only if the first
   * cell is associated to a bar with strictly smaller death than the bar associated to the second one.
   */
  Matrix(const std::function<bool(Pos_index, Pos_index)>& birthComparator,
         const std::function<bool(Pos_index, Pos_index)>& deathComparator);
  /**
   * @brief Constructs a new matrix from the given ranges with the given comparator functions.
   * Only available when those comparators are necessary.
   *
   * That is, when **all** following options have following values:
   *   - @ref PersistenceMatrixOptions::is_of_boundary_type = false
   *   - @ref PersistenceMatrixOptions::has_vine_update = true
   *   - @ref PersistenceMatrixOptions::has_column_pairings = false
   *
   * See description of @ref Matrix(const std::vector<Container>& columns, Characteristic characteristic)
   * for more information about  @p orderedBoundaries and
   * @ref Matrix(const std::function<bool(Pos_index,Pos_index)>&, const std::function<bool(Pos_index,Pos_index)>&)
   * for more information about the comparators.
   *
   * @tparam Boundary_range Range type for @ref Entry_representative ranges. Assumed to have a begin(), end() and size()
   * method.
   * @param orderedBoundaries Vector of ordered boundaries in filtration order. Indexed continuously starting at 0.
   * @param birthComparator Method taking two @ref PosIdx indices as parameter and returns true if and only if the first
   * cell is associated to a bar with strictly smaller birth than the bar associated to the second one.
   * @param deathComparator Method taking two @ref PosIdx indices as parameter and returns true if and only if the first
   * cell is associated to a bar with strictly smaller death than the bar associated to the second one.
   * @param characteristic Characteristic of the coefficient field. Has to be specified if
   * @ref PersistenceMatrixOptions::is_z2 is false. Default value is 11.
   * Ignored if @ref PersistenceMatrixOptions::is_z2 is true.
   */
  template <class Boundary_range = Boundary>
  Matrix(const std::vector<Boundary_range>& orderedBoundaries,
         const std::function<bool(Pos_index, Pos_index)>& birthComparator,
         const std::function<bool(Pos_index, Pos_index)>& deathComparator,
         Characteristic characteristic = 11);
  /**
   * @brief Constructs a new empty matrix and reserves space for the given number of columns.
   * Only available when those comparators are necessary.
   *
   * That is, when **all** following options have following values:
   *   - @ref PersistenceMatrixOptions::is_of_boundary_type = false
   *   - @ref PersistenceMatrixOptions::has_vine_update = true
   *   - @ref PersistenceMatrixOptions::has_column_pairings = false
   *
   * See description of
   * @ref Matrix(const std::function<bool(Pos_index,Pos_index)>&, const std::function<bool(Pos_index,Pos_index)>&)
   * for more information about the comparators.
   *
   * @param numberOfColumns Number of columns to reserve space for.
   * @param birthComparator Method taking two @ref PosIdx indices as parameter and returns true if and only if the first
   * cell is associated to a bar with strictly smaller birth than the bar associated to the second one.
   * @param deathComparator Method taking two @ref PosIdx indices as parameter and returns true if and only if the first
   * cell is associated to a bar with strictly smaller death than the bar associated to the second one.
   * @param characteristic Characteristic of the coefficient field. If not specified and
   * @ref PersistenceMatrixOptions::is_z2 is false, the characteristic has to be set later with the use of
   * @ref set_characteristic before calling for the first time a method needing it.
   * Ignored if @ref PersistenceMatrixOptions::is_z2 is true.
   */
  Matrix(unsigned int numberOfColumns,
         const std::function<bool(Pos_index, Pos_index)>& birthComparator,
         const std::function<bool(Pos_index, Pos_index)>& deathComparator,
         Characteristic characteristic = Field_operators::nullCharacteristic);
  /**
   * @brief Copy constructor.
   *
   * @param matrixToCopy %Matrix to copy.
   */
  Matrix(const Matrix& matrixToCopy);
  /**
   * @brief Move constructor.
   * After the move, the given matrix will be empty.
   *
   * @param other %Matrix to move.
   */
  Matrix(Matrix&& other) noexcept;

  ~Matrix();

  // TODO: compatibility with multi fields:
  //   - set_characteristic(Characteristic min, Characteristic max)
  //   - readapt reduction?
  /**
   * @brief Sets the characteristic of the coefficient field if @ref PersistenceMatrixOptions::is_z2 is false,
   * does nothing otherwise.
   * Should be used if no characteristic could be specified at the creation of the empty matrix.
   * Do not change the value of the characteristic once used.
   *
   * @warning The coefficient values stored in the matrix are stored after computing the corresponding modulo.
   * Therefore, changing the characteristic after is very likely to invalidate all entry values.
   *
   * @param characteristic The characteristic to set.
   */
  void set_characteristic(Characteristic characteristic);

  // (TODO: if there is no row access and the column type corresponds to the internal column type of the matrix,
  // moving the column instead of copying it should be possible. Is it worth implementing it?)
  /**
   * @brief Inserts a new ordered column at the end of the matrix by copying the given range of
   * @ref Entry_representative. The content of the range is assumed to be sorted by increasing ID value.
   *
   * Only available for @ref basematrix "base matrices".
   * Otherwise use @ref insert_boundary which will deduce a new column from the boundary given.
   *
   * @tparam Container Range of @ref Entry_representative. Assumed to have a begin(), end() and size() method.
   * @param column Column to be inserted.
   */
  template <class Container, class = std::enable_if_t<!std::is_arithmetic_v<Container> > >
  void insert_column(const Container& column);
  /**
   * @brief Inserts a new ordered column at the given index by copying the given range of @ref Entry_representative.
   * There should not be any other column inserted at that index which was not explicitly removed before.
   * The content of the range is assumed to be sorted by increasing ID value.
   *
   * Only available for @ref basematrix "base matrices" without column compression and without row access.
   *
   * @tparam Container Range of @ref Entry_representative. Assumed to have a begin(), end() and size() method.
   * @param column Column to be inserted.
   * @param columnIndex @ref MatIdx index to which the column has to be inserted.
   */
  template <class Container, class = std::enable_if_t<!std::is_arithmetic_v<Container> > >
  void insert_column(const Container& column, Index columnIndex);
  /**
   * @brief Inserts a new column at the end of the matrix. The column will consist of the given index only.
   * 
   * @param idx Entry ID.
   * @param e Entry coefficient. Ignored if the coefficient field is Z2. Default value: 1.
   */
  void insert_column(ID_index idx, Element e = 1U);
  // TODO: for simple boundary matrices, add an index pointing to the first column inserted after the last call of
  // get_current_barcode to enable several calls to get_current_barcode
  /**
   * @brief Inserts at the end of the matrix a new ordered column corresponding to the given boundary.
   * This means that it is assumed that this method is called on boundaries in the order of the filtration.
   * It also assumes that the cells in the given boundary are identified by their relative position in the filtration,
   * starting at 0. If it is not the case, use the other
   * @ref insert_boundary(ID_index cellIndex, const Boundary_range& boundary, Dimension dim) "insert_boundary"
   * instead by indicating the cell ID used in the boundaries when the cell is inserted.
   *
   * Different to the constructor, the boundaries do not have to come from a simplicial complex, but also from
   * a more general entry complex. This includes cubical complexes or Morse complexes for example.
   *
   * The content of the new column will vary depending on the underlying @ref mp_matrices "type of the matrix":
   * - If it is a @ref basematrix "basic matrix" type, the boundary is copied as it is, i.e., the method is equivalent
   * to @ref insert_column.
   * - If it is a @ref boundarymatrix "boundary type matrix" and only \f$ R \f$ is stored, the boundary is also just
   * copied. The column will only be reduced later when the barcode is requested in order to apply some optimizations
   * with the additional knowledge. Hence, the barcode will also not be updated, so call @ref get_current_barcode only
   * when the matrix is complete.
   * - If it is a @ref boundarymatrix "boundary type matrix" and both \f$ R \f$ and \f$ U \f$ are stored, the new
   * boundary is stored in its reduced form and the barcode, if active, is also updated.
   * - If it is a @ref chainmatrix "chain type matrix", the new column is of the form `IDIdx + linear combination of
   * older column IDIdxs`, where the combination is deduced while reducing the given boundary. If the barcode is stored,
   * it will also be updated.
   *
   * @tparam Boundary_range Range of @ref Entry_representative. Assumed to have a begin(), end() and size() method.
   * @param boundary Boundary generating the new column. The content should be ordered by ID.
   * @param dim Dimension of the cell whose boundary is given. If the complex is simplicial,
   * this parameter can be omitted as it can be deduced from the size of the boundary.
   * @return If it is a @ref chainmatrix "chain matrix", the method returns the @ref MatIdx indices of the unpaired
   * chains used to reduce the boundary. Otherwise, nothing.
   */
  template <class Boundary_range = Boundary>
  Insertion_return insert_boundary(const Boundary_range& boundary, Dimension dim = Matrix::get_null_value<Dimension>());
  /**
   * @brief Only available for @ref mp_matrices "non-basic matrices".
   * It does the same as the other version, but allows the boundary cells to be identified without restrictions
   * except that all IDs have to be strictly increasing in the order of filtration. Note that you should avoid then
   * to use the other insertion method to avoid overwriting IDs.
   *
   * As a cell has to be inserted before one of its cofaces in a valid filtration (recall that it is assumed that
   * for @ref mp_matrices "non-basic matrices", the cells are inserted by order of filtration), it is sufficient to
   * indicate the ID of the cell being inserted.
   *
   * @tparam Boundary_range Range of @ref Entry_representative. Assumed to have a begin(), end() and size() method.
   * @param cellIndex @ref IDIdx index to use to identify the new cell.
   * @param boundary Boundary generating the new column. The indices of the boundary have to correspond to the
   * @p cellIndex values of precedent calls of the method for the corresponding cells and should be ordered in
   * increasing order.
   * @param dim Dimension of the cell whose boundary is given. If the complex is simplicial,
   * this parameter can be omitted as it can be deduced from the size of the boundary.
   * @return If it is a @ref chainmatrix "chain matrix", the method returns the @ref MatIdx indices of the unpaired
   * chains used to reduce the boundary. Otherwise, nothing.
   */
  template <class Boundary_range = Boundary>
  Insertion_return insert_boundary(ID_index cellIndex,
                                   const Boundary_range& boundary,
                                   Dimension dim = Matrix::get_null_value<Dimension>());

  /**
   * @brief Returns the column at the given @ref MatIdx index.
   * For @ref boundarymatrix "RU matrices", is equivalent to
   * @ref get_column(Index columnIndex, bool inR) "get_column(columnIndex, true)".
   * The type of the column depends on the chosen options, see @ref PersistenceMatrixOptions::column_type.
   *
   * @param columnIndex @ref MatIdx index of the column to return.
   * @return Reference to the column. Is `const` if the matrix has column compression.
   */
  Returned_column& get_column(Index columnIndex);
  /**
   * @brief Only available for @ref chainmatrix "chain matrices". Returns the column at the given @ref MatIdx index.
   * The type of the column depends on the chosen options, see @ref PersistenceMatrixOptions::column_type.
   *
   * @param columnIndex @ref MatIdx index of the column to return.
   * @return Const reference to the column.
   */
  const Column& get_column(Index columnIndex) const;
  // TODO: there is no particular reason that this method is not available for identifier indexing,
  //  just has to be added to the interface...
  /**
   * @brief Only available for @ref boundarymatrix "RU matrices" without @ref Column_indexation_types::IDENTIFIER
   * indexing. Returns the column at the given @ref MatIdx index in \f$ R \f$ if @p inR is true and in \f$ U \f$ if
   * @p inR is false. The type of the column depends on the chosen options,
   * see @ref PersistenceMatrixOptions::column_type.
   *
   * @param columnIndex @ref MatIdx index of the column to return.
   * @param inR If true, returns the column in \f$ R \f$, if false, returns the column in \f$ U \f$.
   * @return Const reference to the column.
   */
  const Column& get_column(Index columnIndex, bool inR);

  // TODO: update column indices when reordering rows (after lazy swap) such that always MatIdx are returned.
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::has_row_access is true. Returns the row at the given
   * @ref rowindex "row index". For @ref boundarymatrix "RU matrices", is equivalent to
   * @ref get_row(ID_index rowIndex, bool inR) "get_row(columnIndex, true)". The type of the row depends on the
   * chosen options, see @ref PersistenceMatrixOptions::has_intrusive_rows.
   *
   * @param rowIndex @ref rowindex "Row index" of the row to return: @ref IDIdx for @ref chainmatrix "chain matrices" or
   * updated @ref IDIdx for @ref boundarymatrix "boundary matrices" if swaps occurred.
   * @return Reference to the row. Is `const` if the matrix has column compression.
   */
  Returned_row& get_row(ID_index rowIndex);
  /**
   * @brief Only available for @ref chainmatrix "chain matrices" and matrices with column compression.
   * Returns the row at the given @ref rowindex "row index".
   * The type of the row depends on the chosen options, see @ref PersistenceMatrixOptions::has_intrusive_rows.
   *
   * @param rowIndex @ref rowindex "Row index" of the row to return: @ref IDIdx for @ref chainmatrix "chain matrices"
   * or updated @ref IDIdx for @ref boundarymatrix "boundary matrices" if swaps occurred.
   * @return Const reference to the row.
   */
  const Row& get_row(ID_index rowIndex) const;
  // TODO: there is no particular reason that this method is not available for identifier indexing,
  //  just has to be added to the interface...
  /**
   * @brief Only available for @ref boundarymatrix "RU matrices" without @ref Column_indexation_types::IDENTIFIER
   * indexing. Returns the row at the given @ref rowindex "row index" in \f$ R \f$ if @p inR is true and in \f$ U \f$ if
   * @p inR is false. The type of the row depends on the chosen options, see
   * @ref PersistenceMatrixOptions::has_intrusive_rows.
   *
   * @param rowIndex @ref rowindex "Row index" of the row to return: updated @ref IDIdx if swaps occurred.
   * @param inR If true, returns the row in \f$ R \f$, if false, returns the row in \f$ U \f$.
   * @return Const reference to the row.
   */
  const Row& get_row(ID_index rowIndex, bool inR);

  /**
   * @brief Only available for @ref basematrix "base matrices" without column compression and if
   * @ref PersistenceMatrixOptions::has_map_column_container is true. Otherwise, see @ref remove_last. Erases the given
   * column from the matrix. If @ref PersistenceMatrixOptions::has_row_access is also true, the deleted column entries
   * are also automatically removed from their respective rows.
   *
   * @param columnIndex @ref MatIdx index of the column to remove.
   */
  void remove_column(Index columnIndex);
  // TODO: rename method to be less confusing.
  /**
   * @brief The effect varies depending on the matrices and the options:
   * - @ref basematrix "base matrix" and @ref boundarymatrix "boundary matrix":
   *    - @ref PersistenceMatrixOptions::has_map_column_container and @ref
   *      PersistenceMatrixOptions::has_column_and_row_swaps are true: cleans up maps used for the lazy row swaps.
   *    - @ref PersistenceMatrixOptions::has_row_access and @ref PersistenceMatrixOptions::has_removable_rows are true:
   *      assumes that the row is empty and removes it.
   *    - Otherwise, does nothing.
   * - @ref boundarymatrix "boundary matrix" with \f$ U \f$ stored: only \f$ R \f$ is affected by the above. If properly
   *   used, \f$ U \f$ will never have empty rows.
   * - @ref chainmatrix "chain matrix": only available if @ref PersistenceMatrixOptions::has_row_access and
   *   @ref PersistenceMatrixOptions::has_removable_rows are true. Assumes that the row is empty and removes it.
   *
   * @warning The removed rows are always assumed to be empty. If it is not the case, the deleted row entries are not
   * removed from their columns. And in the case of intrusive rows, this will generate a segmentation fault when
   * the column entries are destroyed later. The row access is just meant as a "read only" access to the rows and the
   * @ref erase_empty_row method just as a way to specify that a row is empty and can therefore be removed from
   * dictionaries. This allows to avoid testing the emptiness of a row at each column entry removal, what can be quite
   * frequent.
   *
   * @param rowIndex @ref rowindex "Row index" of the empty row to remove.
   */
  void erase_empty_row(ID_index rowIndex);
  // TODO: for chain matrices, replace IDIdx input with MatIdx input to homogenize.
  /**
   * @brief Only available for @ref boundarymatrix "RU" and @ref chainmatrix "chain matrices" and if
   * @ref PersistenceMatrixOptions::has_removable_columns and @ref PersistenceMatrixOptions::has_vine_update are true.
   * For @ref chainmatrix "chain matrices", @ref PersistenceMatrixOptions::has_map_column_container and
   * @ref PersistenceMatrixOptions::has_column_pairings also need to be true. Assumes that the cell is maximal in the
   * current complex and removes it such that the matrix remains consistent (i.e., RU is still an upper triangular
   * decomposition of the boundary matrix and chain is still a compatible bases of the chain complex in the sense
   * of @cite zigzag). The maximality of the cell is not verified. Also updates the barcode if it was computed.
   *
   * For @ref chainmatrix "chain matrices", using the other version of the method could perform better depending on how
   * the data is maintained on the side of the user. Then, @ref PersistenceMatrixOptions::has_column_pairings also do
   * not need to be true.
   *
   * See also @ref remove_last and @ref remove_column.
   *
   * @param columnIndex If @ref boundarymatrix "boundary matrix", @ref MatIdx index of the cell to remove, otherwise the
   * @ref IDIdx index.
   */
  void remove_maximal_cell(Index columnIndex);
  // TODO: See if it would be better to use something more general than a vector for columnsToSwap, such that
  //  the user do not have to construct the vector from scratch. Like passing iterators instead. But it would be nice,
  //  to still be able to do (cell, {})...
  /**
   * @brief Only available for @ref chainmatrix "chain matrices" and if
   * @ref PersistenceMatrixOptions::has_removable_columns, @ref PersistenceMatrixOptions::has_vine_update and
   * @ref PersistenceMatrixOptions::has_map_column_container are true. Assumes that the cell is maximal in the current
   * complex and removes it such that the matrix remains consistent (i.e., it is still a compatible bases of the chain
   * complex in the sense of @cite zigzag). The maximality of the cell is not verified. Also updates the barcode if it
   * was computed.
   *
   * To maintain the compatibility, vine swaps are done to move the cell up to the end of the filtration. Once at the
   * end, the removal is trivial. But for @ref chainmatrix "chain matrices", swaps do not actually swap the position of
   * the column every time, so the cells appearing after @p cellIndex in the filtration have to be searched first within
   * the matrix. If the user has an easy access to the @ref IDIdx of the cells in the order of filtration, passing them
   * by argument with @p columnsToSwap allows to skip a linear search process. Typically, if the user knows that the
   * cell he wants to remove is already the last cell of the filtration, calling
   * @ref remove_maximal_cell(ID_index cellIndex, const std::vector<ID_index>& columnsToSwap)
   * "remove_maximal_cell(cellID, {})" will be faster than @ref remove_last().
   *
   * See also @ref remove_last.
   *
   * @param cellIndex @ref IDIdx index of the cell to remove
   * @param columnsToSwap Vector of @ref IDIdx indices of the cells coming after @p cellIndex in the filtration.
   */
  void remove_maximal_cell(ID_index cellIndex, const std::vector<ID_index>& columnsToSwap);
  /**
   * @brief Removes the last inserted column/cell from the matrix.
   * If the matrix is @ref mp_matrices "non basic", @ref PersistenceMatrixOptions::has_removable_columns has to be true
   * for the method to be available. Additionally, if the matrix is a @ref chainmatrix "chain matrix", either
   * @ref PersistenceMatrixOptions::has_map_column_container has to be true or
   * @ref PersistenceMatrixOptions::has_vine_update has to be false. And if the matrix is a
   * @ref basematrix "base matrix" it should be without column compression.
   *
   * See also @ref remove_maximal_cell and @ref remove_column.
   *
   * For @ref chainmatrix "chain matrices", if @ref PersistenceMatrixOptions::has_vine_update is true, the last cell
   * does not have to be at the end of the matrix and therefore has to be searched first. In this case, if the user
   * already knows the @ref IDIdx of the last cell, calling
   * @ref remove_maximal_cell(ID_index cellIndex, const std::vector<ID_index>& columnsToSwap)
   * "remove_maximal_cell(cellID, {})" instead allows to skip the search.
   */
  void remove_last();

  /**
   * @brief Returns the maximal dimension of a cell stored in the matrix. Only available for
   * @ref mp_matrices "non-basic matrices" and if @ref PersistenceMatrixOptions::has_matrix_maximal_dimension_access
   * is true.
   *
   * @return The maximal dimension.
   */
  Dimension get_max_dimension() const;
  /**
   * @brief Returns the current number of columns in the matrix.
   *
   * @return The number of columns.
   */
  Index get_number_of_columns() const;
  /**
   * @brief Returns the dimension of the given cell. Only available for @ref mp_matrices "non-basic matrices".
   *
   * @param columnIndex @ref MatIdx index of the column representing the cell.
   * @return Dimension of the cell.
   */
  Dimension get_column_dimension(Index columnIndex) const;

  /**
   * @brief Adds column at @p sourceColumnIndex onto the column at @p targetColumnIndex in the matrix. Is available
   * for every matrix type, but should be used with care with @ref mp_matrices "non-basic matrices", as they will be
   * no verification to ensure that the addition makes sense for the meaning of the underlying object. For example,
   * a right-to-left addition could corrupt the computation of the barcode or the representative cycles if done blindly.
   *
   * For @ref basematrix "basic matrices" with column compression, the representatives are summed together, which means
   * that all column compressed together with the target column are affected by the change, not only the target.
   *
   * @tparam Integer_index Any signed or unsigned integer type.
   * @param sourceColumnIndex @ref MatIdx index of the column to add.
   * @param targetColumnIndex @ref MatIdx index of the target column.
   */
  template <typename Integer_index>
  std::enable_if_t<std::is_integral_v<Integer_index> > add_to(Integer_index sourceColumnIndex,
                                                              Integer_index targetColumnIndex);
  /**
   * @brief Adds the given range of @ref Entry onto the column at @p targetColumnIndex in the matrix. Only available
   * for @ref basematrix "basic matrices".
   *
   * For @ref basematrix "basic matrices" with column compression, the range is summed onto the representative, which
   * means that all column compressed together with the target column are affected by the change, not only the target.
   *
   * @tparam Entry_range Range of @ref Entry. Needs a begin() and end() method. A column index does not need to be
   * stored in the entries, even if @ref PersistenceMatrixOptions::has_row_access is true.
   * @param sourceColumn Source @ref Entry range.
   * @param targetColumnIndex @ref MatIdx index of the target column.
   */
  template <class Entry_range>
  std::enable_if_t<!std::is_integral_v<Entry_range> > add_to(const Entry_range& sourceColumn, Index targetColumnIndex);

  /**
   * @brief Multiplies the target column with the coefficient and then adds the source column to it.
   * That is: `targetColumn = (targetColumn * coefficient) + sourceColumn`.
   * Is available for every matrix type, but should be used with care with @ref mp_matrices "non-basic matrices",
   * as they will be no verification to ensure that the addition makes sense for the meaning of the underlying object.
   * For example, a right-to-left addition could corrupt the computation of the barcode or the representative cycles
   * if done blindly.
   *
   * For @ref basematrix "basic matrices" with column compression, the representatives are summed together, which means
   * that all column compressed together with the target column are affected by the change, not only the target.
   *
   * @tparam Integer_index Any signed or unsigned integer type.
   * @param sourceColumnIndex @ref MatIdx index of the column to add.
   * @param coefficient Value to multiply.
   * @param targetColumnIndex @ref MatIdx index of the target column.
   */
  template <typename Integer_index>
  std::enable_if_t<std::is_integral_v<Integer_index> > multiply_target_and_add_to(Integer_index sourceColumnIndex,
                                                                                  int coefficient,
                                                                                  Integer_index targetColumnIndex);
  /**
   * @brief Multiplies the target column with the coefficient and then adds the given range of @ref Entry to it.
   * That is: `targetColumn = (targetColumn * coefficient) + sourceColumn`. Only available for
   * @ref basematrix "basic matrices".
   *
   * For @ref basematrix "basic matrices" with column compression, the range is summed onto the representative, which
   * means that all column compressed together with the target column are affected by the change, not only the target.
   *
   * @tparam Entry_range Range of @ref Entry. Needs a begin() and end() method. A column index does not need to be
   * stored in the entries, even if @ref PersistenceMatrixOptions::has_row_access is true.
   * @param sourceColumn Source @ref Entry range.
   * @param coefficient Value to multiply.
   * @param targetColumnIndex @ref MatIdx index of the target column.
   */
  template <class Entry_range>
  std::enable_if_t<!std::is_integral_v<Entry_range> > multiply_target_and_add_to(const Entry_range& sourceColumn,
                                                                                 int coefficient,
                                                                                 Index targetColumnIndex);

  /**
   * @brief Multiplies the source column with the coefficient before adding it to the target column.
   * That is: `targetColumn += (coefficient * sourceColumn)`. The source column will **not** be modified.
   * Is available for every matrix type, but should be used with care with @ref mp_matrices "non-basic matrices", as
   * they will be no verification to ensure that the addition makes sense for the meaning of the underlying object.
   * For example, a right-to-left addition could corrupt the computation of the barcode or the representative cycles
   * if done blindly.
   *
   * For @ref basematrix "basic matrices" with column compression, the representatives are summed together, which means
   * that all column compressed together with the target column are affected by the change, not only the target.
   *
   * @tparam Integer_index Any signed or unsigned integer type.
   * @param coefficient Value to multiply.
   * @param sourceColumnIndex @ref MatIdx index of the column to add.
   * @param targetColumnIndex @ref MatIdx index of the target column.
   */
  template <typename Integer_index>
  std::enable_if_t<std::is_integral_v<Integer_index> > multiply_source_and_add_to(int coefficient,
                                                                                  Integer_index sourceColumnIndex,
                                                                                  Integer_index targetColumnIndex);
  /**
   * @brief Multiplies the source column with the coefficient before adding it to the target column.
   * That is: `targetColumn += (coefficient * sourceColumn)`. The source column will **not** be modified.
   * Only available for @ref basematrix "basic matrices".
   *
   * For @ref basematrix "basic matrices" with column compression, the range is summed onto the representative, which
   * means that all column compressed together with the target column are affected by the change, not only the target.
   *
   * @tparam Entry_range Range of @ref Entry. Needs a begin() and end() method. A column index does not need to be
   * stored in the entries, even if @ref PersistenceMatrixOptions::has_row_access is true.
   * @param coefficient Value to multiply.
   * @param sourceColumn Source @ref Entry range.
   * @param targetColumnIndex @ref MatIdx index of the target column.
   */
  template <class Entry_range>
  std::enable_if_t<!std::is_integral_v<Entry_range> > multiply_source_and_add_to(int coefficient,
                                                                                 const Entry_range& sourceColumn,
                                                                                 Index targetColumnIndex);

  /**
   * @brief Zeroes the entry at the given coordinates. Not available for @ref chainmatrix "chain matrices" and for
   * @ref basematrix "base matrices" with column compression. In general, should be used with care with
   * @ref mp_matrices "non-basic matrices" to not destroy the validity of the persistence related properties of the
   * matrix.
   *
   * For @ref boundarymatrix "RU matrices", equivalent to
   * @ref zero_entry(Index columnIndex, ID_index rowIndex, bool inR) "zero_entry(columnIndex, rowIndex, true)".
   *
   * @param columnIndex @ref MatIdx index of the column of the entry.
   * @param rowIndex @ref rowindex "Row index" of the row of the entry.
   */
  void zero_entry(Index columnIndex, ID_index rowIndex);
  /**
   * @brief Only available for @ref boundarymatrix "RU matrices". Zeroes the entry at the given coordinates in \f$ R \f$
   * if @p inR is true or in \f$ U \f$ if @p inR is false. Should be used with care to not destroy the validity of the
   * persistence related properties of the matrix.
   *
   * @param columnIndex @ref MatIdx index of the column of the entry.
   * @param rowIndex @ref rowindex "Row index" of the row of the entry.
   * @param inR Boolean indicating in which matrix to zero: if true in \f$ R \f$ and if false in \f$ U \f$.
   */
  void zero_entry(Index columnIndex, ID_index rowIndex, bool inR);
  /**
   * @brief Zeroes the column at the given index. Not available for @ref chainmatrix "chain matrices" and for
   * @ref basematrix "base matrices" with column compression. In general, should be used with care with
   * @ref mp_matrices "non-basic matrices" to not destroy the validity of the persistence related properties of the
   * matrix.
   *
   * For @ref boundarymatrix "RU matrices", equivalent to
   * @ref zero_column(Index columnIndex, bool inR) "zero_column(columnIndex, true)".
   *
   * @param columnIndex @ref MatIdx index of the column to zero.
   */
  void zero_column(Index columnIndex);
  /**
   * @brief Only available for @ref boundarymatrix "RU matrices". Zeroes the column at the given index in \f$ R \f$ if
   * @p inR is true or in \f$ U \f$ if @p inR is false. Should be used with care to not destroy the validity of the
   * persistence related properties of the matrix.
   *
   * @param columnIndex @ref MatIdx index of the column to zero.
   * @param inR Boolean indicating in which matrix to zero: if true in \f$ R \f$ and if false in \f$ U \f$.
   */
  void zero_column(Index columnIndex, bool inR);
  /**
   * @brief Indicates if the entry at given coordinates has value zero.
   *
   * For @ref boundarymatrix "RU matrices", equivalent to
   * @ref is_zero_entry(Index columnIndex, ID_index rowIndex, bool inR) const
   * "is_zero_entry(columnIndex, rowIndex, true)".
   *
   * @param columnIndex @ref MatIdx index of the column of the entry.
   * @param rowIndex @ref rowindex "Row index" of the row of the entry.
   * @return true If the entry has value zero.
   * @return false Otherwise.
   */
  bool is_zero_entry(Index columnIndex, ID_index rowIndex);
  /**
   * @brief Only available for @ref boundarymatrix "RU matrices". Indicates if the entry at given coordinates has value
   * zero in \f$ R \f$ if @p inR is true or in \f$ U \f$ if @p inR is false.
   *
   * @param columnIndex @ref MatIdx index of the column of the entry.
   * @param rowIndex @ref rowindex "Row index" of the row of the entry.
   * @param inR Boolean indicating in which matrix to look: if true in \f$ R \f$ and if false in \f$ U \f$.
   * @return true If the entry has value zero.
   * @return false Otherwise.
   */
  bool is_zero_entry(Index columnIndex, ID_index rowIndex, bool inR) const;
  /**
   * @brief Indicates if the column at given index has value zero.
   *
   * For @ref boundarymatrix "RU matrices", equivalent to
   * @ref is_zero_column(Index columnIndex, bool inR) "is_zero_column(columnIndex, true)".
   *
   * Note that for @ref chainmatrix "chain matrices", this method should always return false, as a valid
   * @ref chainmatrix "chain matrix" never has empty columns.
   *
   * @param columnIndex @ref MatIdx index of the column.
   * @return true If the column has value zero.
   * @return false Otherwise.
   */
  bool is_zero_column(Index columnIndex);
  /**
   * @brief Only available for @ref boundarymatrix "RU matrices". Indicates if the column at given index has value zero
   * in \f$ R \f$ if @p inR is true or in \f$ U \f$ if @p inR is false.
   *
   * Note that if @p inR is false, this method should usually return false.
   *
   * @param columnIndex @ref MatIdx index of the column.
   * @param inR Boolean indicating in which matrix to look: if true in \f$ R \f$ and if false in \f$ U \f$.
   * @return true If the column has value zero.
   * @return false Otherwise.
   */
  bool is_zero_column(Index columnIndex, bool inR);

  /**
   * @brief Returns the @ref MatIdx index of the column which has the given @ref rowindex "row index" as pivot. Only
   * available for @ref boundarymatrix "RU" and @ref chainmatrix "chain matrices". Assumes that the pivot exists. For
   * @ref boundarymatrix "RU matrices", the column is returned from \f$ R \f$.
   *
   * Recall that the @ref rowindex "row indices" for @ref chainmatrix "chain matrices" correspond to the @ref IDIdx
   * indices and that the @ref rowindex "row indices" for a @ref boundarymatrix "RU matrix" correspond to the updated
   * @ref IDIdx indices which got potentially swapped by a vine swap.
   *
   * @param cellIndex @ref rowindex "Row index" of the pivot.
   * @return @ref MatIdx index of the column with the given pivot.
   */
  Index get_column_with_pivot(ID_index cellIndex) const;
  /**
   * @brief Returns the @ref rowindex "row index" of the pivot of the given column. Only available for
   * @ref mp_matrices "non-basic matrices".
   *
   * @param columnIndex @ref MatIdx index of the column
   * @return The @ref rowindex "row index" of the pivot.
   */
  ID_index get_pivot(Index columnIndex);

  /**
   * @brief Assign operator.
   *
   * @param other %Matrix to copy
   * @return Reference to this object.
   */
  Matrix& operator=(Matrix other) &;
  /**
   * @brief Assign operator.
   *
   * @param other %Matrix to move
   * @return Reference to this object.
   */
  Matrix& operator=(Matrix&& other) && noexcept;

  /**
   * @brief Swap operator for two matrices.
   *
   * @param matrix1 First matrix to swap.
   * @param matrix2 Second matrix to swap.
   */
  friend void swap(Matrix& matrix1, Matrix& matrix2) noexcept
  {
    swap(matrix1.matrix_, matrix2.matrix_);
    std::swap(matrix1.colSettings_, matrix2.colSettings_);
  }

  void print();  // for debug

  // TODO: change the behaviour for boundary matrices.
  /**
   * @brief Returns the current barcode of the matrix. Available only if
   * @ref PersistenceMatrixOptions::has_column_pairings is true.
   *
   * Recall that we assume that the boundaries were inserted in the order of filtration for the barcode to be valid.
   *
   * @warning For simple @ref boundarymatrix "boundary matrices" (only storing \f$ R \f$), we assume that
   * @ref get_current_barcode is only called once the matrix is completed and won't be modified again.
   *
   * @return A reference to the barcode. The barcode is a vector of @ref Matrix::Bar. A bar stores three information:
   * the @ref PosIdx birth index, the @ref PosIdx death index and the dimension of the bar.
   */
  const Barcode& get_current_barcode();
  /**
   * @brief Returns the current barcode of the matrix. Available only if
   * @ref PersistenceMatrixOptions::has_column_pairings is true.
   *
   * Recall that we assume that the boundaries were inserted in the order of filtration for the barcode to be valid.
   *
   * @warning For simple @ref boundarymatrix "boundary matrices" (only storing \f$ R \f$), we assume that
   * @ref get_current_barcode is only called once the matrix is completed and won't be modified again.
   *
   * @return A const reference to the barcode. The barcode is a vector of @ref Matrix::Bar. A bar stores three
   * information: the @ref PosIdx birth index, the @ref PosIdx death index and the dimension of the bar.
   */
  const Barcode& get_current_barcode() const;

  /**
   * @brief Only available for @ref basematrix "base matrices" without column compression and simple
   * @ref boundarymatrix "boundary matrices" (only storing \f$ R \f$) and if
   * @ref PersistenceMatrixOptions::has_column_and_row_swaps is true. Swaps the two given columns. Note for
   * @ref boundarymatrix "boundary matrices", that it really just swaps two columns and does not update anything else,
   * nor performs additions to maintain some properties on the matrix.
   *
   * @param columnIndex1 First column @ref MatIdx index to swap.
   * @param columnIndex2 Second column @ref MatIdx index to swap.
   */
  void swap_columns(Index columnIndex1, Index columnIndex2);
  /**
   * @brief Only available for @ref basematrix "base matrices" without column compression and simple
   * @ref boundarymatrix "boundary matrices" (only storing \f$ R \f$) and if
   * @ref PersistenceMatrixOptions::has_column_and_row_swaps is true. Swaps the two given rows. Note for
   * @ref boundarymatrix "boundary matrices", that it really just swaps two rows and does not update anything else,
   * nor performs additions to maintain some properties on the matrix.
   *
   * @param rowIndex1 First @ref rowindex "row index" to swap.
   * @param rowIndex2 Second @ref rowindex "row index" to swap.
   */
  void swap_rows(Index rowIndex1, Index rowIndex2);
  // TODO: find better name. And benchmark also to verify if it is really worth it to have this extra version in
  // addition to vine_swap.
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::has_vine_update is true and if it is either a boundary
   * matrix or @ref PersistenceMatrixOptions::column_indexation_type is set to @ref Column_indexation_types::POSITION.
   * Does the same than @ref vine_swap, but assumes that the swap is non trivial and therefore skips a part of the case
   * study.
   *
   * @param index @ref PosIdx index of the first cell to swap. The second one has to be at `index + 1`. Recall that for
   * @ref boundarymatrix "boundary matrices", @ref PosIdx == @ref MatIdx.
   * @return true If the barcode changed from the swap.
   * @return false Otherwise.
   */
  bool vine_swap_with_z_eq_1_case(Pos_index index);
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::has_vine_update is true and if it is either a
   * @ref chainmatrix "chain matrix" or @ref PersistenceMatrixOptions::column_indexation_type is set to
   * @ref Column_indexation_types::IDENTIFIER. Does the same than @ref vine_swap, but assumes that the swap is
   * non-trivial and therefore skips a part of the case study.
   *
   * @param columnIndex1 @ref MatIdx index of the first cell.
   * @param columnIndex2 @ref MatIdx index of the second cell. It is assumed that the @ref PosIdx of both only differs
   * by one.
   * @return Let \f$ pos1 \f$ be the @ref PosIdx index of @p columnIndex1 and \f$ pos2 \f$ be the @ref PosIdx index of
   * @p columnIndex2. The method returns the @ref MatIdx of the column which has now, after the swap, the @ref PosIdx
   * \f$ max(pos1, pos2) \f$.
   */
  Index vine_swap_with_z_eq_1_case(Index columnIndex1, Index columnIndex2);
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::has_vine_update is true and if it is either a
   * @ref boundarymatrix "boundary matrix" or @ref PersistenceMatrixOptions::column_indexation_type is set to
   * @ref Column_indexation_types::POSITION. Does a vine swap between two cells which are consecutive in the
   * filtration. Roughly, if \f$ F \f$ is the current filtration represented by the matrix, the method modifies the
   * matrix such that the new state corresponds to a valid state for the filtration \f$ F' \f$ equal to \f$ F \f$ but
   * with the two cells at position `index` and `index + 1` swapped. Of course, the two cells should not have a
   * face/coface relation which each other ; \f$ F' \f$ has to be a valid filtration.
   * See @cite vineyards for more information about vine and vineyards.
   *
   * @param index @ref PosIdx index of the first cell to swap. The second one has to be at `index + 1`. Recall that for
   * @ref boundarymatrix "boundary matrices", @ref PosIdx == @ref MatIdx.
   * @return true If the barcode changed from the swap.
   * @return false Otherwise.
   */
  bool vine_swap(Pos_index index);
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::has_vine_update is true and if it is either a
   * @ref chainmatrix "chain matrix" or @ref PersistenceMatrixOptions::column_indexation_type is set to
   * @ref Column_indexation_types::IDENTIFIER. Does a vine swap between two cells which are consecutive in the
   * filtration. Roughly, if \f$ F \f$ is the current filtration represented by the matrix, the method modifies the
   * matrix such that the new state corresponds to a valid state for the filtration \f$ F' \f$ equal to \f$ F \f$ but
   * with the two given cells at swapped positions. Of course, the two cells should not have a face/coface relation
   * which each other ; \f$ F' \f$ has to be a valid filtration.
   * See @cite vineyards for more information about vine and vineyards.
   *
   * @param columnIndex1 @ref MatIdx index of the first cell.
   * @param columnIndex2 @ref MatIdx index of the second cell. It is assumed that the @ref PosIdx of both only differs
   * by one.
   * @return Let \f$ pos1 \f$ be the @ref PosIdx index of @p columnIndex1 and \f$ pos2 \f$ be the @ref PosIdx index of
   * @p columnIndex2. The method returns the @ref MatIdx of the column which has now, after the swap, the @ref PosIdx
   * \f$ max(pos1, pos2) \f$.
   */
  Index vine_swap(Index columnIndex1, Index columnIndex2);

  /**
   * @brief Only available if @ref PersistenceMatrixOptions::can_retrieve_representative_cycles is true. Pre-computes
   * the representative cycles of the current state of the filtration represented by the matrix. It needs to be called
   * before calling @ref get_representative_cycles if the matrix was modified since last call. Otherwise the old cycles
   * will be returned.
   *
   * @param dim If different from default value, only the cycles of the given dimension are updated.
   * All others are erased.
   */
  void update_representative_cycles(Dimension dim = get_null_value<Dimension>());
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::can_retrieve_representative_cycles is true. Pre-computes
   * the representative cycle in the current matrix state of the given bar. It needs to be called
   * before calling @ref get_representative_cycle if the matrix was modified since last call. Otherwise the old cycle
   * will be returned.
   *
   * @param bar Bar corresponding to the wanted representative cycle.
   */
  void update_representative_cycle(const Bar& bar);
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::can_retrieve_representative_cycles is true.
   * Returns all representative cycles of the current filtration. @ref update_representative_cycles has to be called
   * first if a modification to the matrix has to be token into account since last call.
   *
   * @return A const reference to the vector of representative cycles.
   */
  const std::vector<Cycle>& get_representative_cycles();
  /**
   * @brief Only available if @ref PersistenceMatrixOptions::can_retrieve_representative_cycles is true.
   * Returns the cycle representing the given bar. @ref update_representative_cycles or
   * @ref update_representative_cycle have to be called first if a modification to the matrix has to be token into
   * account since last call.
   *
   * @param bar A bar from the current barcode.
   * @return A const reference to the cycle representing @p bar.
   */
  const Cycle& get_representative_cycle(const Bar& bar);

 private:
  using Underlying_matrix = std::conditional_t<
      isNonBasic,
      std::conditional_t<
          PersistenceMatrixOptions::is_of_boundary_type,
          std::conditional_t<
              PersistenceMatrixOptions::has_vine_update || PersistenceMatrixOptions::can_retrieve_representative_cycles,
              std::conditional_t<
                  PersistenceMatrixOptions::column_indexation_type == Column_indexation_types::CONTAINER ||
                      PersistenceMatrixOptions::column_indexation_type == Column_indexation_types::POSITION,
                  Master_RU_matrix,
                  Id_to_index_overlay<Master_RU_matrix, Matrix<PersistenceMatrixOptions>
              > >,
              std::conditional_t<
                  PersistenceMatrixOptions::column_indexation_type == Column_indexation_types::CONTAINER ||
                      PersistenceMatrixOptions::column_indexation_type == Column_indexation_types::POSITION,
                  Master_boundary_matrix,
                  Id_to_index_overlay<Master_boundary_matrix, Matrix<PersistenceMatrixOptions>
              > > >,
          std::conditional_t<
              PersistenceMatrixOptions::column_indexation_type == Column_indexation_types::CONTAINER,
              Master_chain_matrix,
              std::conditional_t<PersistenceMatrixOptions::column_indexation_type == Column_indexation_types::POSITION,
                                 Position_to_index_overlay<Master_chain_matrix, Matrix<PersistenceMatrixOptions> >,
                                 Id_to_index_overlay<Master_chain_matrix, Matrix<PersistenceMatrixOptions>
                                > > > >,
      Master_base_matrix>;

  // Field_operators* operators_;
  // Entry_constructor* entryPool_;
  Column_settings* colSettings_;  // pointer because the of swap operator on matrix_ which also stores the pointer
  Underlying_matrix matrix_;

  static constexpr void _assert_options();

  Element _get_value(int coefficient) const;
};

template <class PersistenceMatrixOptions>
inline Matrix<PersistenceMatrixOptions>::Matrix() : colSettings_(new Column_settings()), matrix_(colSettings_)
{
  static_assert(
      PersistenceMatrixOptions::is_of_boundary_type || !PersistenceMatrixOptions::has_vine_update ||
          PersistenceMatrixOptions::has_column_pairings,
      "When no barcode is recorded with vine swaps, comparison functions for the columns have to be provided.");
  _assert_options();
}

template <class PersistenceMatrixOptions>
template <class Container>
inline Matrix<PersistenceMatrixOptions>::Matrix(const std::vector<Container>& columns, Characteristic characteristic)
    : colSettings_(new Column_settings(characteristic)), matrix_(columns, colSettings_)
{
  static_assert(PersistenceMatrixOptions::is_of_boundary_type || !PersistenceMatrixOptions::has_vine_update ||
                    PersistenceMatrixOptions::has_column_pairings,
                "When no barcode is recorded with vine swaps for chain matrices, comparison functions for the columns "
                "have to be provided.");
  _assert_options();
}

template <class PersistenceMatrixOptions>
inline Matrix<PersistenceMatrixOptions>::Matrix(unsigned int numberOfColumns, Characteristic characteristic)
    : colSettings_(new Column_settings(characteristic)), matrix_(numberOfColumns, colSettings_)
{
  static_assert(PersistenceMatrixOptions::is_of_boundary_type || !PersistenceMatrixOptions::has_vine_update ||
                    PersistenceMatrixOptions::has_column_pairings,
                "When no barcode is recorded with vine swaps for chain matrices, comparison functions for the columns "
                "have to be provided.");
  _assert_options();
}

template <class PersistenceMatrixOptions>
inline Matrix<PersistenceMatrixOptions>::Matrix(const std::function<bool(Pos_index, Pos_index)>& birthComparator,
                                                const std::function<bool(Pos_index, Pos_index)>& deathComparator)
    : colSettings_(new Column_settings()), matrix_(colSettings_, birthComparator, deathComparator)
{
  static_assert(
      !PersistenceMatrixOptions::is_of_boundary_type && PersistenceMatrixOptions::has_vine_update &&
          !PersistenceMatrixOptions::has_column_pairings,
      "Constructor only available for chain matrices when vine swaps are enabled, but the barcode is not recorded.");
  _assert_options();
}

template <class PersistenceMatrixOptions>
template <class Boundary_range>
inline Matrix<PersistenceMatrixOptions>::Matrix(const std::vector<Boundary_range>& orderedBoundaries,
                                                const std::function<bool(Pos_index, Pos_index)>& birthComparator,
                                                const std::function<bool(Pos_index, Pos_index)>& deathComparator,
                                                Characteristic characteristic)
    : colSettings_(new Column_settings(characteristic)),
      matrix_(orderedBoundaries, colSettings_, birthComparator, deathComparator)
{
  static_assert(
      !PersistenceMatrixOptions::is_of_boundary_type && PersistenceMatrixOptions::has_vine_update &&
          !PersistenceMatrixOptions::has_column_pairings,
      "Constructor only available for chain matrices when vine swaps are enabled, but the barcode is not recorded.");
  _assert_options();
}

template <class PersistenceMatrixOptions>
inline Matrix<PersistenceMatrixOptions>::Matrix(unsigned int numberOfColumns,
                                                const std::function<bool(Pos_index, Pos_index)>& birthComparator,
                                                const std::function<bool(Pos_index, Pos_index)>& deathComparator,
                                                Characteristic characteristic)
    : colSettings_(new Column_settings(characteristic)),
      matrix_(numberOfColumns, colSettings_, birthComparator, deathComparator)
{
  static_assert(
      !PersistenceMatrixOptions::is_of_boundary_type && PersistenceMatrixOptions::has_vine_update &&
          !PersistenceMatrixOptions::has_column_pairings,
      "Constructor only available for chain matrices when vine swaps are enabled, but the barcode is not recorded.");
  _assert_options();
}

template <class PersistenceMatrixOptions>
inline Matrix<PersistenceMatrixOptions>::Matrix(const Matrix& matrixToCopy)
    : colSettings_(new Column_settings(*matrixToCopy.colSettings_)), matrix_(matrixToCopy.matrix_, colSettings_)
{
  _assert_options();
}

template <class PersistenceMatrixOptions>
inline Matrix<PersistenceMatrixOptions>::Matrix(Matrix&& other) noexcept
    : colSettings_(std::exchange(other.colSettings_, nullptr)), matrix_(std::move(other.matrix_))
{
  _assert_options();
}

template <class PersistenceMatrixOptions>
inline Matrix<PersistenceMatrixOptions>::~Matrix()
{
  matrix_.reset(colSettings_);  // to avoid crashes at destruction, all columns have to be destroyed first
  delete colSettings_;
}

template <class PersistenceMatrixOptions>
inline void Matrix<PersistenceMatrixOptions>::set_characteristic(Characteristic characteristic)
{
  if constexpr (!PersistenceMatrixOptions::is_z2) {
    if (colSettings_->operators.get_characteristic() != Field_operators::nullCharacteristic) {
      std::cerr << "Warning: Characteristic already initialised. Changing it could lead to incoherences in the matrix "
                   "as the modulo was already applied to values in existing columns.";
    }

    colSettings_->operators.set_characteristic(characteristic);
  }
}

template <class PersistenceMatrixOptions>
template <class Container, class>
inline void Matrix<PersistenceMatrixOptions>::insert_column(const Container& column)
{
  if constexpr (!PersistenceMatrixOptions::is_z2) {
    GUDHI_CHECK(colSettings_->operators.get_characteristic() != Field_operators::nullCharacteristic,
                std::logic_error("Matrix::insert_column - Columns cannot be initialized if the coefficient field "
                                 "characteristic is not specified."));
  }

  static_assert(
      !isNonBasic,
      "'insert_column' not available for the chosen options. The input has to be in the form of a cell boundary.");
  matrix_.insert_column(column);
}

template <class PersistenceMatrixOptions>
template <class Container, class>
inline void Matrix<PersistenceMatrixOptions>::insert_column(const Container& column, Index columnIndex)
{
  if constexpr (!PersistenceMatrixOptions::is_z2) {
    GUDHI_CHECK(colSettings_->operators.get_characteristic() != Field_operators::nullCharacteristic,
                std::logic_error("Matrix::insert_column - Columns cannot be initialized if the coefficient field "
                                 "characteristic is not specified."));
  }

  static_assert(!isNonBasic && !PersistenceMatrixOptions::has_column_compression,
                "'insert_column' with those parameters is not available for the chosen options.");
  static_assert(!PersistenceMatrixOptions::has_row_access,
                "Columns have to be inserted at the end of the matrix when row access is enabled.");
  matrix_.insert_column(column, columnIndex);
}

template <class PersistenceMatrixOptions>
inline void Matrix<PersistenceMatrixOptions>::insert_column(ID_index idx, [[maybe_unused]] Element e)
{
  static_assert(
      !isNonBasic,
      "'insert_column' not available for the chosen options. The input has to be in the form of a cell boundary.");

  if constexpr (PersistenceMatrixOptions::is_z2) {
    matrix_.insert_column(idx);
  } else {
    GUDHI_CHECK(colSettings_->operators.get_characteristic() != Field_operators::nullCharacteristic,
                std::logic_error("Matrix::insert_column - Columns cannot be initialized if the coefficient field "
                                 "characteristic is not specified."));
    matrix_.insert_column(idx, e);
  }
}

template <class PersistenceMatrixOptions>
template <class Boundary_range>
inline typename Matrix<PersistenceMatrixOptions>::Insertion_return Matrix<PersistenceMatrixOptions>::insert_boundary(
    const Boundary_range& boundary,
    Dimension dim)
{
  if constexpr (!PersistenceMatrixOptions::is_z2) {
    GUDHI_CHECK(colSettings_->operators.get_characteristic() != Field_operators::nullCharacteristic,
                std::logic_error("Matrix::insert_boundary - Columns cannot be initialized if the coefficient field "
                                 "characteristic is not specified."));
  }

  if constexpr (isNonBasic && !PersistenceMatrixOptions::is_of_boundary_type &&
                PersistenceMatrixOptions::column_indexation_type == Column_indexation_types::CONTAINER)
    return matrix_.insert_boundary(boundary, dim);
  else
    matrix_.insert_boundary(boundary, dim);
}

template <class PersistenceMatrixOptions>
template <class Boundary_range>
inline typename Matrix<PersistenceMatrixOptions>::Insertion_return
Matrix<PersistenceMatrixOptions>::insert_boundary(ID_index cellIndex, const Boundary_range& boundary, Dimension dim)
{
  if constexpr (!PersistenceMatrixOptions::is_z2) {
    GUDHI_CHECK(colSettings_->operators.get_characteristic() != Field_operators::nullCharacteristic,
                std::logic_error("Matrix::insert_boundary - Columns cannot be initialized if the coefficient field "
                                 "characteristic is not specified."));
  }

  static_assert(isNonBasic, "Only enabled for non-basic matrices.");
  if constexpr (!PersistenceMatrixOptions::is_of_boundary_type &&
                PersistenceMatrixOptions::column_indexation_type == Column_indexation_types::CONTAINER)
    return matrix_.insert_boundary(cellIndex, boundary, dim);
  else
    matrix_.insert_boundary(cellIndex, boundary, dim);
}

template <class PersistenceMatrixOptions>
inline typename Matrix<PersistenceMatrixOptions>::Returned_column& Matrix<PersistenceMatrixOptions>::get_column(
    Index columnIndex)
{
  return matrix_.get_column(columnIndex);
}

template <class PersistenceMatrixOptions>
inline const typename Matrix<PersistenceMatrixOptions>::Column& Matrix<PersistenceMatrixOptions>::get_column(
    Index columnIndex) const
{
  return matrix_.get_column(columnIndex);
}

template <class PersistenceMatrixOptions>
inline const typename Matrix<PersistenceMatrixOptions>::Column& Matrix<PersistenceMatrixOptions>::get_column(
    Index columnIndex,
    bool inR)
{
  // TODO: I don't think there is a particular reason why the indexation is forced, should be removed.
  static_assert(
      isNonBasic && PersistenceMatrixOptions::is_of_boundary_type &&
          (PersistenceMatrixOptions::has_vine_update || PersistenceMatrixOptions::can_retrieve_representative_cycles) &&
          PersistenceMatrixOptions::column_indexation_type != Column_indexation_types::IDENTIFIER,
      "Only enabled for position indexed RU matrices.");

  return matrix_.get_column(columnIndex, inR);
}

template <class PersistenceMatrixOptions>
inline typename Matrix<PersistenceMatrixOptions>::Returned_row& Matrix<PersistenceMatrixOptions>::get_row(
    ID_index rowIndex)
{
  static_assert(PersistenceMatrixOptions::has_row_access, "'get_row' is not available for the chosen options.");

  return matrix_.get_row(rowIndex);
}

template <class PersistenceMatrixOptions>
inline const typename Matrix<PersistenceMatrixOptions>::Row& Matrix<PersistenceMatrixOptions>::get_row(
    ID_index rowIndex) const
{
  static_assert(PersistenceMatrixOptions::has_row_access, "'get_row' is not available for the chosen options.");

  return matrix_.get_row(rowIndex);
}

template <class PersistenceMatrixOptions>
inline const typename Matrix<PersistenceMatrixOptions>::Row& Matrix<PersistenceMatrixOptions>::get_row(
    ID_index rowIndex,
    bool inR)
{
  static_assert(PersistenceMatrixOptions::has_row_access, "'get_row' is not available for the chosen options.");
  // TODO: I don't think there is a particular reason why the indexation is forced, should be removed.
  static_assert(
      isNonBasic && PersistenceMatrixOptions::is_of_boundary_type &&
          (PersistenceMatrixOptions::has_vine_update || PersistenceMatrixOptions::can_retrieve_representative_cycles) &&
          PersistenceMatrixOptions::column_indexation_type != Column_indexation_types::IDENTIFIER,
      "Only enabled for position indexed RU matrices.");

  return matrix_.get_row(rowIndex, inR);
}

template <class PersistenceMatrixOptions>
inline void Matrix<PersistenceMatrixOptions>::remove_column(Index columnIndex)
{
  static_assert(PersistenceMatrixOptions::has_map_column_container && !isNonBasic &&
                    !PersistenceMatrixOptions::has_column_compression,
                "'remove_column' is not available for the chosen options.");

  matrix_.remove_column(columnIndex);
}

template <class PersistenceMatrixOptions>
inline void Matrix<PersistenceMatrixOptions>::erase_empty_row(ID_index rowIndex)
{
  static_assert(
      !isNonBasic || PersistenceMatrixOptions::is_of_boundary_type || PersistenceMatrixOptions::has_removable_rows,
      "'erase_empty_row' is not available for the chosen options.");

  matrix_.erase_empty_row(rowIndex);
}

template <class PersistenceMatrixOptions>
inline void Matrix<PersistenceMatrixOptions>::remove_maximal_cell(Index columnIndex)
{
  static_assert(PersistenceMatrixOptions::has_removable_columns,
                "'remove_maximal_cell(ID_index)' is not available for the chosen options.");
  static_assert(isNonBasic && PersistenceMatrixOptions::has_vine_update,
                "'remove_maximal_cell(ID_index)' is not available for the chosen options.");
  static_assert(PersistenceMatrixOptions::is_of_boundary_type || (PersistenceMatrixOptions::has_map_column_container &&
                                                                  PersistenceMatrixOptions::has_column_pairings),
                "'remove_maximal_cell(ID_index)' is not available for the chosen options.");

  matrix_.remove_maximal_cell(columnIndex);
}

template <class PersistenceMatrixOptions>
inline void Matrix<PersistenceMatrixOptions>::remove_maximal_cell(ID_index cellIndex,
                                                                  const std::vector<ID_index>& columnsToSwap)
{
  static_assert(PersistenceMatrixOptions::has_removable_columns,
                "'remove_maximal_cell(ID_index,const std::vector<Index>&)' is not available for the chosen options.");
  static_assert(isNonBasic && !PersistenceMatrixOptions::is_of_boundary_type,
                "'remove_maximal_cell(ID_index,const std::vector<Index>&)' is not available for the chosen options.");
  static_assert(PersistenceMatrixOptions::has_map_column_container && PersistenceMatrixOptions::has_vine_update,
                "'remove_maximal_cell(ID_index,const std::vector<Index>&)' is not available for the chosen options.");

  matrix_.remove_maximal_cell(cellIndex, columnsToSwap);
}

template <class PersistenceMatrixOptions>
inline void Matrix<PersistenceMatrixOptions>::remove_last()
{
  static_assert(PersistenceMatrixOptions::has_removable_columns || !isNonBasic,
                "'remove_last' is not available for the chosen options.");
  static_assert(!PersistenceMatrixOptions::has_column_compression || isNonBasic,
                "'remove_last' is not available for the chosen options.");
  static_assert(!isNonBasic || PersistenceMatrixOptions::is_of_boundary_type ||
                    PersistenceMatrixOptions::has_map_column_container || !PersistenceMatrixOptions::has_vine_update,
                "'remove_last' is not available for the chosen options.");

  matrix_.remove_last();
}

template <class PersistenceMatrixOptions>
inline typename Matrix<PersistenceMatrixOptions>::Dimension Matrix<PersistenceMatrixOptions>::get_max_dimension() const
{
  static_assert(isNonBasic, "'get_max_dimension' is not available for the chosen options.");

  return matrix_.get_max_dimension();
}

template <class PersistenceMatrixOptions>
inline typename Matrix<PersistenceMatrixOptions>::Index Matrix<PersistenceMatrixOptions>::get_number_of_columns() const
{
  return matrix_.get_number_of_columns();
}

template <class PersistenceMatrixOptions>
inline typename Matrix<PersistenceMatrixOptions>::Dimension Matrix<PersistenceMatrixOptions>::get_column_dimension(
    Index columnIndex) const
{
  static_assert(isNonBasic, "'get_column_dimension' is not available for the chosen options.");

  return matrix_.get_column_dimension(columnIndex);
}

template <class PersistenceMatrixOptions>
template <typename Integer_index>
inline std::enable_if_t<std::is_integral_v<Integer_index> > Matrix<PersistenceMatrixOptions>::add_to(
    Integer_index sourceColumnIndex,
    Integer_index targetColumnIndex)
{
  matrix_.add_to(sourceColumnIndex, targetColumnIndex);
}

template <class PersistenceMatrixOptions>
template <class Entry_range>
inline std::enable_if_t<!std::is_integral_v<Entry_range> > Matrix<PersistenceMatrixOptions>::add_to(
    const Entry_range& sourceColumn,
    Index targetColumnIndex)
{
  static_assert(!isNonBasic,
                "For boundary or chain matrices, only additions with columns inside the matrix is allowed to maintain "
                "algebraic consistency.");

  matrix_.add_to(sourceColumn, targetColumnIndex);
}

template <class PersistenceMatrixOptions>
template <typename Integer_index>
inline std::enable_if_t<std::is_integral_v<Integer_index> >
Matrix<PersistenceMatrixOptions>::multiply_target_and_add_to(Integer_index sourceColumnIndex,
                                                             int coefficient,
                                                             Integer_index targetColumnIndex)
{
  matrix_.multiply_target_and_add_to(sourceColumnIndex, _get_value(coefficient), targetColumnIndex);
}

template <class PersistenceMatrixOptions>
template <class Entry_range>
inline std::enable_if_t<!std::is_integral_v<Entry_range> > Matrix<PersistenceMatrixOptions>::multiply_target_and_add_to(
    const Entry_range& sourceColumn,
    int coefficient,
    Index targetColumnIndex)
{
  static_assert(!isNonBasic,
                "For boundary or chain matrices, only additions with columns inside the matrix is allowed to maintain "
                "algebraic consistency.");

  matrix_.multiply_target_and_add_to(sourceColumn, _get_value(coefficient), targetColumnIndex);
}

template <class PersistenceMatrixOptions>
template <typename Integer_index>
inline std::enable_if_t<std::is_integral_v<Integer_index> >
Matrix<PersistenceMatrixOptions>::multiply_source_and_add_to(int coefficient,
                                                             Integer_index sourceColumnIndex,
                                                             Integer_index targetColumnIndex)
{
  matrix_.multiply_source_and_add_to(_get_value(coefficient), sourceColumnIndex, targetColumnIndex);
}

template <class PersistenceMatrixOptions>
template <class Entry_range>
inline std::enable_if_t<!std::is_integral_v<Entry_range> > Matrix<PersistenceMatrixOptions>::multiply_source_and_add_to(
    int coefficient,
    const Entry_range& sourceColumn,
    Index targetColumnIndex)
{
  static_assert(!isNonBasic,
                "For boundary or chain matrices, only additions with columns inside the matrix is allowed to maintain "
                "algebraic consistency.");

  matrix_.multiply_source_and_add_to(_get_value(coefficient), sourceColumn, targetColumnIndex);
}

template <class PersistenceMatrixOptions>
inline void Matrix<PersistenceMatrixOptions>::zero_entry(Index columnIndex, ID_index rowIndex)
{
  static_assert(PersistenceMatrixOptions::is_of_boundary_type && !PersistenceMatrixOptions::has_column_compression,
                "'zero_entry' is not available for the chosen options.");

  return matrix_.zero_entry(columnIndex, rowIndex);
}

template <class PersistenceMatrixOptions>
inline void Matrix<PersistenceMatrixOptions>::zero_entry(Index columnIndex, ID_index rowIndex, bool inR)
{
  // TODO: I don't think there is a particular reason why the indexation is forced, should be removed.
  static_assert(
      isNonBasic && PersistenceMatrixOptions::is_of_boundary_type &&
          (PersistenceMatrixOptions::has_vine_update || PersistenceMatrixOptions::can_retrieve_representative_cycles) &&
          PersistenceMatrixOptions::column_indexation_type != Column_indexation_types::IDENTIFIER,
      "Only enabled for RU matrices.");

  return matrix_.zero_entry(columnIndex, rowIndex, inR);
}

template <class PersistenceMatrixOptions>
inline void Matrix<PersistenceMatrixOptions>::zero_column(Index columnIndex)
{
  static_assert(PersistenceMatrixOptions::is_of_boundary_type && !PersistenceMatrixOptions::has_column_compression,
                "'zero_column' is not available for the chosen options.");

  return matrix_.zero_column(columnIndex);
}

template <class PersistenceMatrixOptions>
inline void Matrix<PersistenceMatrixOptions>::zero_column(Index columnIndex, bool inR)
{
  // TODO: I don't think there is a particular reason why the indexation is forced, should be removed.
  static_assert(
      isNonBasic && PersistenceMatrixOptions::is_of_boundary_type &&
          (PersistenceMatrixOptions::has_vine_update || PersistenceMatrixOptions::can_retrieve_representative_cycles) &&
          PersistenceMatrixOptions::column_indexation_type != Column_indexation_types::IDENTIFIER,
      "Only enabled for RU matrices.");

  return matrix_.zero_column(columnIndex, inR);
}

template <class PersistenceMatrixOptions>
inline bool Matrix<PersistenceMatrixOptions>::is_zero_entry(Index columnIndex, ID_index rowIndex)
{
  return matrix_.is_zero_entry(columnIndex, rowIndex);
}

template <class PersistenceMatrixOptions>
inline bool Matrix<PersistenceMatrixOptions>::is_zero_entry(Index columnIndex, ID_index rowIndex, bool inR) const
{
  // TODO: I don't think there is a particular reason why the indexation is forced, should be removed.
  static_assert(
      isNonBasic && PersistenceMatrixOptions::is_of_boundary_type &&
          (PersistenceMatrixOptions::has_vine_update || PersistenceMatrixOptions::can_retrieve_representative_cycles) &&
          PersistenceMatrixOptions::column_indexation_type != Column_indexation_types::IDENTIFIER,
      "Only enabled for RU matrices.");

  return matrix_.is_zero_entry(columnIndex, rowIndex, inR);
}

template <class PersistenceMatrixOptions>
inline bool Matrix<PersistenceMatrixOptions>::is_zero_column(Index columnIndex)
{
  return matrix_.is_zero_column(columnIndex);
}

template <class PersistenceMatrixOptions>
inline bool Matrix<PersistenceMatrixOptions>::is_zero_column(Index columnIndex, bool inR)
{
  // TODO: I don't think there is a particular reason why the indexation is forced, should be removed.
  static_assert(
      isNonBasic && PersistenceMatrixOptions::is_of_boundary_type &&
          (PersistenceMatrixOptions::has_vine_update || PersistenceMatrixOptions::can_retrieve_representative_cycles) &&
          PersistenceMatrixOptions::column_indexation_type != Column_indexation_types::IDENTIFIER,
      "Only enabled for RU matrices.");

  return matrix_.is_zero_column(columnIndex, inR);
}

template <class PersistenceMatrixOptions>
inline typename Matrix<PersistenceMatrixOptions>::Index Matrix<PersistenceMatrixOptions>::get_column_with_pivot(
    ID_index cellIndex) const
{
  static_assert(isNonBasic && (!PersistenceMatrixOptions::is_of_boundary_type ||
                               (PersistenceMatrixOptions::has_vine_update ||
                                PersistenceMatrixOptions::can_retrieve_representative_cycles)),
                "'get_column_with_pivot' is not available for the chosen options.");

  return matrix_.get_column_with_pivot(cellIndex);
}

template <class PersistenceMatrixOptions>
inline typename Matrix<PersistenceMatrixOptions>::ID_index Matrix<PersistenceMatrixOptions>::get_pivot(
    Index columnIndex)
{
  static_assert(isNonBasic, "'get_pivot' is not available for the chosen options.");

  return matrix_.get_pivot(columnIndex);
}

template <class PersistenceMatrixOptions>
inline Matrix<PersistenceMatrixOptions>& Matrix<PersistenceMatrixOptions>::operator=(Matrix other) &
{
  swap(matrix_, other.matrix_);
  std::swap(colSettings_, other.colSettings_);

  return *this;
}

template <class PersistenceMatrixOptions>
inline Matrix<PersistenceMatrixOptions>& Matrix<PersistenceMatrixOptions>::operator=(Matrix&& other) && noexcept
{
  if (this == &other) return *this;

  matrix_ = std::move(other.matrix_);
  colSettings_ = std::exchange(other.colSettings_, nullptr);

  return *this;
}

template <class PersistenceMatrixOptions>
inline void Matrix<PersistenceMatrixOptions>::print()
{
  return matrix_.print();
}

template <class PersistenceMatrixOptions>
inline const typename Matrix<PersistenceMatrixOptions>::Barcode& Matrix<PersistenceMatrixOptions>::get_current_barcode()
{
  static_assert(PersistenceMatrixOptions::has_column_pairings, "This method was not enabled.");

  return matrix_.get_current_barcode();
}

template <class PersistenceMatrixOptions>
inline const typename Matrix<PersistenceMatrixOptions>::Barcode& Matrix<PersistenceMatrixOptions>::get_current_barcode()
    const
{
  static_assert(PersistenceMatrixOptions::has_column_pairings, "This method was not enabled.");
  static_assert(
      !PersistenceMatrixOptions::is_of_boundary_type || PersistenceMatrixOptions::has_vine_update ||
          PersistenceMatrixOptions::can_retrieve_representative_cycles,
      "'get_current_barcode' is not const for boundary matrices as the barcode is only computed when explicitly "
      "asked.");

  return matrix_.get_current_barcode();
}

template <class PersistenceMatrixOptions>
inline void Matrix<PersistenceMatrixOptions>::swap_columns(Index columnIndex1, Index columnIndex2)
{
  static_assert(
      (!isNonBasic && !PersistenceMatrixOptions::has_column_compression) ||
          (isNonBasic && PersistenceMatrixOptions::is_of_boundary_type && !PersistenceMatrixOptions::has_vine_update &&
           !PersistenceMatrixOptions::can_retrieve_representative_cycles),
      "This method was not enabled.");
  return matrix_.swap_columns(columnIndex1, columnIndex2);
}

template <class PersistenceMatrixOptions>
inline void Matrix<PersistenceMatrixOptions>::swap_rows(Index rowIndex1, Index rowIndex2)
{
  static_assert(
      (!isNonBasic && !PersistenceMatrixOptions::has_column_compression) ||
          (isNonBasic && PersistenceMatrixOptions::is_of_boundary_type && !PersistenceMatrixOptions::has_vine_update &&
           !PersistenceMatrixOptions::can_retrieve_representative_cycles),
      "This method was not enabled.");
  return matrix_.swap_rows(rowIndex1, rowIndex2);
}

template <class PersistenceMatrixOptions>
inline bool Matrix<PersistenceMatrixOptions>::vine_swap_with_z_eq_1_case(Pos_index index)
{
  static_assert(PersistenceMatrixOptions::has_vine_update, "This method was not enabled.");
  static_assert(PersistenceMatrixOptions::column_indexation_type == Column_indexation_types::POSITION ||
                    (PersistenceMatrixOptions::is_of_boundary_type &&
                     PersistenceMatrixOptions::column_indexation_type == Column_indexation_types::CONTAINER),
                "This method was not enabled.");
  return matrix_.vine_swap_with_z_eq_1_case(index);
}

template <class PersistenceMatrixOptions>
inline typename Matrix<PersistenceMatrixOptions>::Index Matrix<PersistenceMatrixOptions>::vine_swap_with_z_eq_1_case(
    Index columnIndex1,
    Index columnIndex2)
{
  static_assert(PersistenceMatrixOptions::has_vine_update, "This method was not enabled.");
  static_assert(PersistenceMatrixOptions::column_indexation_type == Column_indexation_types::IDENTIFIER ||
                    (!PersistenceMatrixOptions::is_of_boundary_type &&
                     PersistenceMatrixOptions::column_indexation_type == Column_indexation_types::CONTAINER),
                "This method was not enabled.");

  return matrix_.vine_swap_with_z_eq_1_case(columnIndex1, columnIndex2);
}

template <class PersistenceMatrixOptions>
inline bool Matrix<PersistenceMatrixOptions>::vine_swap(Pos_index index)
{
  static_assert(PersistenceMatrixOptions::has_vine_update, "This method was not enabled.");
  static_assert(PersistenceMatrixOptions::column_indexation_type == Column_indexation_types::POSITION ||
                    (PersistenceMatrixOptions::is_of_boundary_type &&
                     PersistenceMatrixOptions::column_indexation_type == Column_indexation_types::CONTAINER),
                "This method was not enabled.");
  return matrix_.vine_swap(index);
}

template <class PersistenceMatrixOptions>
inline typename Matrix<PersistenceMatrixOptions>::Index Matrix<PersistenceMatrixOptions>::vine_swap(Index columnIndex1,
                                                                                                    Index columnIndex2)
{
  static_assert(PersistenceMatrixOptions::has_vine_update, "This method was not enabled.");
  static_assert(PersistenceMatrixOptions::column_indexation_type == Column_indexation_types::IDENTIFIER ||
                    (!PersistenceMatrixOptions::is_of_boundary_type &&
                     PersistenceMatrixOptions::column_indexation_type == Column_indexation_types::CONTAINER),
                "This method was not enabled.");
  return matrix_.vine_swap(columnIndex1, columnIndex2);
}

template <class PersistenceMatrixOptions>
inline void Matrix<PersistenceMatrixOptions>::update_representative_cycles(Dimension dim)
{
  static_assert(PersistenceMatrixOptions::can_retrieve_representative_cycles, "This method was not enabled.");
  matrix_.update_representative_cycles(dim);
}

template <class PersistenceMatrixOptions>
inline void Matrix<PersistenceMatrixOptions>::update_representative_cycle(const Bar& bar)
{
  static_assert(PersistenceMatrixOptions::can_retrieve_representative_cycles, "This method was not enabled.");
  matrix_.update_representative_cycle(bar);
}

template <class PersistenceMatrixOptions>
inline const std::vector<typename Matrix<PersistenceMatrixOptions>::Cycle>&
Matrix<PersistenceMatrixOptions>::get_representative_cycles()
{
  static_assert(PersistenceMatrixOptions::can_retrieve_representative_cycles, "This method was not enabled.");
  return matrix_.get_representative_cycles();
}

template <class PersistenceMatrixOptions>
inline const typename Matrix<PersistenceMatrixOptions>::Cycle&
Matrix<PersistenceMatrixOptions>::get_representative_cycle(const Bar& bar)
{
  static_assert(PersistenceMatrixOptions::can_retrieve_representative_cycles, "This method was not enabled.");
  return matrix_.get_representative_cycle(bar);
}

template <class PersistenceMatrixOptions>
constexpr void Matrix<PersistenceMatrixOptions>::_assert_options()
{
  static_assert(
      PersistenceMatrixOptions::column_type != Column_types::HEAP || !PersistenceMatrixOptions::has_row_access,
      "Row access is not possible for heap columns.");
  static_assert(!PersistenceMatrixOptions::has_vine_update || PersistenceMatrixOptions::is_z2,
                "Vine update currently works only for Z_2 coefficients.");
  // static_assert(!PersistenceMatrixOptions::can_retrieve_representative_cycles || PersistenceMatrixOptions::is_z2,
  //               "Representative cycles can currently only be computed with Z_2 coefficients.");
  static_assert(
      PersistenceMatrixOptions::column_type != Column_types::HEAP || !PersistenceMatrixOptions::has_column_compression,
      "Column compression not compatible with heap columns.");

  // // This should be warnings instead, as PersistenceMatrixOptions::has_column_compression is just ignored in those
  // // cases and don't produces errors as long as the corresponding methods are not called.
  // static_assert(!PersistenceMatrixOptions::has_column_compression || !PersistenceMatrixOptions::has_column_pairings,
  //               "Column compression not available to compute persistence homology (it would bring no advantages; "
  //               "use it for co-homology instead).");
  // static_assert(!PersistenceMatrixOptions::has_column_compression || !PersistenceMatrixOptions::has_vine_update,
  //               "Column compression not available for vineyards.");
  // static_assert(!PersistenceMatrixOptions::has_column_compression ||
  //                   !PersistenceMatrixOptions::can_retrieve_representative_cycles,
  //               "Column compression not available to retrieve representative cycles.");
  // // Would column removal while column compression be useful? If yes, should erase() remove a single column or the
  // // class of columns identical to the input?
  // // For a single column, I have an implementation for union-find (not the current one) which allows deleting a
  // // single element in constant time, but find becomes log n in worst case.
  // // For a column class, we either just ignore the removed class (constant time), but this results in memory
  // // residues, or, we have an erase method which is at least linear in the size of the class.
  // static_assert(
  //     !PersistenceMatrixOptions::has_column_compression || !PersistenceMatrixOptions::has_map_column_container,
  //     "When column compression is used, the removal of columns is not implemented yet.");
}

template <class PersistenceMatrixOptions>
inline typename Matrix<PersistenceMatrixOptions>::Element Matrix<PersistenceMatrixOptions>::_get_value(
    int coefficient) const
{
  return get_coefficient_value(coefficient, get_operator_ptr(colSettings_));
}

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // MASTER_MATRIX_H

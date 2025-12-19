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
 * @file persistence_matrix_options.h
 * @author Hannah Schreiber
 * @brief Contains the options for the matrix template.
 */

#ifndef PM_OPTIONS_INCLUDED
#define PM_OPTIONS_INCLUDED

#include <cstdint>  // std::uint8_t

#include <gudhi/Fields/Zp_field_operators.h>

namespace Gudhi {
namespace persistence_matrix {

/**
 * @ingroup persistence_matrix
 *
 * @brief List of column types.
 */
enum class Column_types : std::uint8_t {
  LIST,           /**< @ref List_column "": Underlying container is a std::list<@ref Entry*>. */
  SET,            /**< @ref Set_column "": Underlying container is a std::set<@ref Entry*>. */
  HEAP,           /**< @ref Heap_column "": Underlying container is a std::vector<@ref Entry*> ordered as a heap.
                       Is not compatible with row access and column compression. */
  VECTOR,         /**< @ref Vector_column "": Underlying container is a std::vector<@ref Entry*>
                       with a lazy removal method. */
  NAIVE_VECTOR,   /**< @ref Naive_vector_column "": Underlying container is a std::vector<@ref Entry*>. */
  SMALL_VECTOR,   /**< @ref Naive_vector_column "": Underlying container is a
                       boost::container::small_vector<@ref Entry*, 8>. */
  UNORDERED_SET,  /**< @ref Unordered_set_column "": Underlying container is a std::unordered_set<@ref Entry*>. */
  INTRUSIVE_LIST, /**< @ref Intrusive_list_column "": Underlying container is a boost::intrusive::list<@ref Entry>. */
  INTRUSIVE_SET   /**< @ref Intrusive_set_column "": Underlying container is a boost::intrusive::set<@ref Entry>. */
};

/**
 * @ingroup persistence_matrix
 *
 * @brief List if indexation schemes. See @ref mp_indexation "description of indexation schemes" for more details
 * about the meaning of the indexation types.
 */
enum class Column_indexation_types : std::uint8_t {
  CONTAINER, /**< Default use of @ref MatIdx indices. */
  POSITION,  /**< All input and output @ref MatIdx indices are replaced with @ref PosIdx indices. */
  IDENTIFIER /**< All input and output @ref MatIdx indices are replaced with @ref IDIdx indices. */
};

/**
 * @struct Default_options persistence_matrix_options.h gudhi/persistence_matrix_options.h
 * @ingroup persistence_matrix
 *
 * @brief Default option structure for @ref Matrix class.
 * See the @ref PersistenceMatrixOptions concept for a more detailed description of the fields.
 * Produces a @ref basematrix "base matrix" with no enabled option.
 *
 * To create other matrix types, the easiest is to simply inherit from this structure and overwrite only the options
 * one is interested in.
 *
 * @tparam col_type Column type for the matrix. Default value: @ref Column_types::INTRUSIVE_SET
 * @tparam is_z2_only Flag indicating if only \f$Z_2\f$ coefficient will be used with the matrix. Set to true if it
 * is the case, false otherwise. Default value: true.
 * @tparam FieldOperators Field operators used by the matrix, see FieldOperators concept.
 * Only necessary if @p is_z2_only is false.
 * Default value: @ref Gudhi::persistence_fields::Zp_field_operators<>.
 */
template <Column_types col_type = Column_types::INTRUSIVE_SET,
          bool is_z2_only = true,
          class FieldOperators = persistence_fields::Zp_field_operators<> >
struct Default_options {
  using Field_coeff_operators = FieldOperators;
  using Dimension = int;
  using Index = unsigned int;

  static const bool is_z2 = is_z2_only;
  static const Column_types column_type = col_type;

  static const Column_indexation_types column_indexation_type = Column_indexation_types::CONTAINER;

  static const bool has_column_compression = false;
  static const bool has_column_and_row_swaps = false;

  static const bool has_map_column_container = false;
  static const bool has_removable_columns = false;

  static const bool has_row_access = false;
  static const bool has_intrusive_rows = true;
  static const bool has_removable_rows = false;

  static const bool is_of_boundary_type = true;

  static const bool has_matrix_maximal_dimension_access = false;
  static const bool has_column_pairings = false;
  static const bool has_vine_update = false;
  static const bool can_retrieve_representative_cycles = false;
};

// TODO: The following structures are the one used by the other modules or debug tests.
//  They will probably be removed once the module was properly integrated.

/**
 * @brief Options used for the Zigzag persistence module.
 *
 * @tparam column_type Column type for the matrix.
 */
template <Column_types column_type = Column_types::INTRUSIVE_LIST>
struct Zigzag_options : Default_options<column_type, true> {
  static const bool has_row_access = true;
  static const bool has_column_pairings = false;
  static const bool has_vine_update = true;
  static const bool is_of_boundary_type = false;
  static const bool has_map_column_container = true;
  static const bool has_removable_columns = true;
  static const bool has_removable_rows = true;
};

/**
 * @brief Options needed to use the representative cycles.
 *
 * @tparam col_type Column type for the matrix.
 */
template <Column_types col_type = Column_types::INTRUSIVE_SET>
struct Representative_cycles_options : Default_options<col_type, true> {
  static const bool has_column_pairings = true;
  static const bool can_retrieve_representative_cycles = true;
};

/**
 * @brief Options used by the Multipersistence module.
 *
 * @tparam column_type Column type for the matrix.
 */
template <Column_types column_type = Column_types::INTRUSIVE_SET>
struct Multi_persistence_options : Default_options<column_type, true> {
  static const bool has_column_pairings = true;
  static const bool has_vine_update = true;
};

/**
 * @brief Options used by the cohomology module.
 *
 * @tparam column_type Column type for the matrix.
 * @tparam is_z2_only True if Z2.
 * @tparam FieldOperators Field operator.
 */
template <Column_types column_type = Column_types::INTRUSIVE_LIST,
          bool is_z2_only = true,
          class FieldOperators = persistence_fields::Zp_field_operators<> >
struct Cohomology_persistence_options : Default_options<column_type, is_z2_only, FieldOperators> {
  static const bool has_row_access = true;
  static const bool has_column_compression = true;
  static const bool has_removable_rows = true;
};

/**
 * @private
 */
template <typename T>
class RangeTraits
{
 private:
  static auto check_begin(...) -> std::false_type;
  template <typename U>
  static auto check_begin(const U& x) -> decltype(x.begin(), std::true_type{});

  static auto check_size(...) -> std::false_type;
  template <typename U>
  static auto check_size(const U& x) -> decltype(x.size(), std::true_type{});

 public:
  static constexpr bool has_begin = decltype(check_begin(std::declval<T>()))::value;
  static constexpr bool has_size = decltype(check_size(std::declval<T>()))::value;
};

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_OPTIONS_INCLUDED

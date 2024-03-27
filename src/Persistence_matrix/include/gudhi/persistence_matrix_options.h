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

#include "Fields/Zp_field_operators.h"

namespace Gudhi {
namespace persistence_matrix {

/**
 * @brief List of column types.
 */
enum Column_types { 
  LIST,           /**< Underlying container is a std::list. */
  SET,            /**< Underlying container is a std::set. */
  HEAP,           /**< Underlying container is a std::vector ordered as a heap.
                       Is not compatible with row access and column compression */
  VECTOR,         /**< Underlying container is a std::vector with a lazy removal method. */
  NAIVE_VECTOR,   /**< Underlying container is a std::vector. */
  UNORDERED_SET,  /**< Underlying container is a std::unordered_set. */
  INTRUSIVE_LIST, /**< Underlying container is a boost::intrusive::list. */
  INTRUSIVE_SET   /**< Underlying container is a boost::intrusive::set. */
};

/**
 * @brief List if indexation schemes. See [TODO: ref to introduction] for more details about the meaning
 * of the indexation types.
 */
enum Column_indexation_types { 
  CONTAINER,  /**< Default use of MatIdx indices. */
  POSITION,   /**< All input and output MatIdx indices are replaced with PosIdx indices. */
  IDENTIFIER  /**< All input and output MatIdx indices are replaced with IDIdx indices. */
};

/**
 * @brief Default option structure for @ref Matrix class.
 * See [TODO: concept] for a more detailed description of the fields.
 * Produces a base matrix with no enabled option.
 *
 * To create other matrix types, the easiest is to simply inherit from this structure and overwrite only the options
 * one is interested in.
 * 
 * @tparam col_type Column type for the matrix. Default value: @ref Column_types::INTRUSIVE_SET
 * @tparam is_z2_only Flag indicating if only \f$Z_2\f$ coefficient will be used with the matrix. Set to true if it
 * is the case, false otherwise. Default value: true.
 * @tparam Field_type Field operators used by the matrix, see [TODO: concept]. Only necessary if @p is_z2_only is false. 
 * Default value: @ref Zp_field_operators<>.
 */
template <Column_types col_type = Column_types::INTRUSIVE_SET, 
          bool is_z2_only = true,
          class Field_type = Zp_field_operators<> >
struct Default_options 
{
  using field_coeff_operators = Field_type;
  // if not signed, max value should not be used.
  using dimension_type = int;
  // type used for the different indexations, if unsigned, the max value is reserved.
  using index_type = unsigned int;

  static const bool is_z2 = is_z2_only;
  static const Column_types column_type = col_type;

  static const Column_indexation_types column_indexation_type = Column_indexation_types::CONTAINER;

  // can be enabled only if no specialised method is enabled: has_column_pairings, has_vine_update,
  // can_retrieve_representative_cycles, has_map_column_container
  static const bool has_column_compression = false;
  // ignored if has_vine_update or can_retrieve_representative_cycles is true.
  static const bool has_column_and_row_swaps = false;

  static const bool has_map_column_container = false;
  static const bool has_removable_columns = false;

  static const bool has_row_access = false;
  // ignored if has_row_access = false
  static const bool has_intrusive_rows = true;
  // ignored if has_row_access = false
  static const bool has_removable_rows = false;

  // ignored if not at least one specialised method is enabled: has_column_pairings, has_vine_update,
  // can_retrieve_representative_cycles
  static const bool is_of_boundary_type = true;

  // ignored and put to false if matrix is not specialised, as the notion of dimension makes no sense for free columns.
  // Also ignored but set to true for base `has_column_pairings`.
  static const bool has_matrix_maximal_dimension_access = false;
  static const bool has_column_pairings = false;
  static const bool has_vine_update = false;
  static const bool can_retrieve_representative_cycles = false;

  // not implemented yet
  //  static const bool is_separated_by_dimension = false;
  //  static const bool is_parallelizable = false;
};

//TODO: The following structures are the one used by the other modules or debug tests.
// They will probably be removed once the module was properly integrated.

/**
 * @brief Options used for the Zigzag persistence module.
 * 
 * @tparam column_type Column type for the matrix.
 */
template <Column_types column_type = Column_types::INTRUSIVE_LIST>
struct Zigzag_options : Default_options<column_type, true> 
{
  static const bool has_row_access = true;
  static const bool has_column_pairings = false;
  static const bool has_vine_update = true;
  static const bool is_of_boundary_type = false;
  static const bool has_map_column_container = true;
  static const bool has_removable_rows = true;
};

/**
 * @brief Options needed to use the representative cycles.
 * 
 * @tparam col_type Column type for the matrix.
 */
template <Column_types col_type = Column_types::INTRUSIVE_SET>
struct Representative_cycles_options : Default_options<col_type, true> 
{
  static const bool has_column_pairings = true;
  static const bool can_retrieve_representative_cycles = true;
};

/**
 * @brief Options used by the Multipersistence module.
 * 
 * @tparam column_type Column type for the matrix.
 */
template <Column_types column_type = Column_types::INTRUSIVE_SET>
struct Multi_persistence_options : Default_options<column_type, true> 
{
  static const bool has_column_pairings = true;
  static const bool has_vine_update = true;
};

/**
 * @brief Options used by the cohomology module.
 * 
 * @tparam column_type Column type for the matrix.
 * @tparam is_z2_only True if Z2.
 * @tparam Field_type Field operator.
 */
template <Column_types column_type = Column_types::INTRUSIVE_LIST, 
          bool is_z2_only = true,
          class Field_type = Zp_field_operators<> >
struct Cohomology_persistence_options : Default_options<column_type, is_z2_only, Field_type> 
{
  static const bool has_row_access = true;
  static const bool has_column_compression = true;
  static const bool has_removable_rows = true;
};

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_OPTIONS_INCLUDED

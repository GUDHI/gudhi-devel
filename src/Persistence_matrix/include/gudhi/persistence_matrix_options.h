/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */
/**
 * @file options.h
 * @author Hannah Schreiber
 * @brief Contains the options for the matrix template.
 */

#ifndef OPTIONS_INCLUDED
#define OPTIONS_INCLUDED

#include "Fields/Zp_field_operators.h"

namespace Gudhi {
namespace persistence_matrix {

enum Column_types {
	LIST,
	SET,
	HEAP,
	VECTOR,
	NAIVE_VECTOR,
	UNORDERED_SET,
	INTRUSIVE_LIST,
	INTRUSIVE_SET
};

enum Column_indexation_types {
	CONTAINER,
	POSITION,
	IDENTIFIER
};

template<Column_types col_type = Column_types::INTRUSIVE_SET, bool is_z2_only = true, class Field_type = Zp_field_operators<> >
struct Default_options{
	using field_coeff_operators = Field_type;
	using dimension_type = int;			//if not signed, max value should not be used.
	//index type differenciation for input clarity, but it does not make much sense to have different types for them.
	using id_type = unsigned int;	//index in the boundaries, if unsigned, the max value is reserved.
	using index_type = id_type;		//index in the underlying container, if unsigned, the max value is reserved.
	using pos_type = id_type;		//relative position in the filtration of all simplices currently represented in the matrix. If unsigned, the max value is reserved.

	static const bool is_z2 = is_z2_only;
	static const Column_types column_type = col_type;

	static const Column_indexation_types column_indexation_type = Column_indexation_types::CONTAINER;

	static const bool is_separated_by_dimension = false;					//not implemented yet
	static const bool is_parallelizable = false;							//not implemented yet
	static const bool is_double_linked = true;								//not implemented yet, it depends of the column type for now (all double linked except for UNORDERED_SET). usefull?

	static const bool has_matrix_maximal_dimension_access = true;			//ignored and put to false if matrix is not specialised, as the notion of dimension makes no sense for free columns. Also ignored but set to true for base `has_column_pairings`.
	static const bool has_row_access = false;
	static const bool has_intrusive_rows = column_type != Column_types::SET && column_type != Column_types::UNORDERED_SET;	//ignored if has_row_access = false
	static const bool has_removable_rows = false;							//ignored if has_row_access = false
	static const bool has_column_pairings = false;
	static const bool has_vine_update = false;
	static const bool can_retrieve_representative_cycles = false;
	static const bool has_map_column_container = false;
	static const bool has_removable_columns = false;
	static const bool is_of_boundary_type = true;							//ignored if not at least one specialised method is enabled: has_column_pairings, has_vine_update, can_retrieve_representative_cycles
	static const bool has_column_compression = false;						//can be enabled only if no specialised method is enabled: has_column_pairings, has_vine_update, can_retrieve_representative_cycles, has_map_column_container
	static const bool has_column_and_row_swaps = false;						//ignored if has_vine_update or can_retrieve_representative_cycles is true.
};

template<Column_types column_type = Column_types::INTRUSIVE_LIST>
struct Zigzag_options : Default_options<column_type, true>{
	static const bool has_row_access = true;
	static const bool has_column_pairings = false;
	static const bool has_vine_update = true;
	static const bool is_of_boundary_type = false;
	static const bool has_map_column_container = true;
	static const bool has_removable_rows = true;
};

template<Column_types col_type = Column_types::INTRUSIVE_SET>
struct Representative_cycles_options : Default_options<col_type, true>{
	static const bool has_column_pairings = true;
	static const bool can_retrieve_representative_cycles = true;
};

template<Column_types column_type = Column_types::INTRUSIVE_SET>
struct Multi_persistence_options : Default_options<column_type, true>{
	static const bool has_column_pairings = true;
	static const bool has_vine_update = true;
};

template<Column_types column_type = Column_types::INTRUSIVE_LIST, bool is_z2_only = true, class Field_type = Zp_field_operators<> >
struct Cohomology_persistence_options : Default_options<column_type, is_z2_only, Field_type>{
	static const bool has_row_access = true;
	static const bool has_column_compression = true;
	static const bool has_removable_rows = true;
};

} //namespace persistence_matrix
} //namespace Gudhi

#endif // OPTIONS_INCLUDED


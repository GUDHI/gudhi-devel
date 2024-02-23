/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_MATRIX_TESTS_OPTIONS_H
#define PM_MATRIX_TESTS_OPTIONS_H

#include <gudhi/persistence_matrix_options.h>
#include <gudhi/matrix.h>
#include "gudhi/Fields/Zp_field_operators.h"

using Gudhi::persistence_matrix::Column_types;
using Gudhi::persistence_matrix::Matrix;
using Gudhi::persistence_matrix::Column_indexation_types;

using Zp = Gudhi::persistence_matrix::Zp_field_operators<> ;

template<bool is_z2_only, Column_types col_type, bool rem_col, bool swaps>
struct Base_options{
	using field_coeff_operators = Zp;
	using id_type = unsigned int;
	using index_type = id_type;
	using pos_type = id_type;
	using dimension_type = int;	//needs to be signed.

	static const bool is_z2 = is_z2_only;
	static const Column_types column_type = col_type;
	static const Column_indexation_types column_indexation_type = Column_indexation_types::CONTAINER;

	static const bool has_matrix_maximal_dimension_access = false;
	static const bool has_column_pairings = false;
	static const bool has_vine_update = false;
	static const bool can_retrieve_representative_cycles = false;
	static const bool is_of_boundary_type = true;
	static const bool has_column_compression = false;

	static const bool has_row_access = false;
	static const bool has_intrusive_rows = false;
	static const bool has_removable_rows = false;

	static const bool has_removable_columns = true;

	static const bool has_map_column_container = rem_col;
	static const bool has_column_and_row_swaps = swaps;
};

template<bool is_z2_only, Column_types col_type, bool rem_row, bool intr_row, bool rem_col, bool swaps>
struct Base_options_with_row_access{
	using field_coeff_operators = Zp;
	using id_type = unsigned int;
	using index_type = id_type;
	using pos_type = id_type;
	using dimension_type = int;	//needs to be signed.

	static const bool is_z2 = is_z2_only;
	static const Column_types column_type = col_type;
	static const Column_indexation_types column_indexation_type = Column_indexation_types::CONTAINER;

	static const bool has_matrix_maximal_dimension_access = false;
	static const bool has_column_pairings = false;
	static const bool has_vine_update = false;
	static const bool can_retrieve_representative_cycles = false;
	static const bool is_of_boundary_type = true;
	static const bool has_column_compression = false;

	static const bool has_row_access = true;
	static const bool has_intrusive_rows = intr_row;
	static const bool has_removable_rows = rem_row;

	static const bool has_removable_columns = true;

	static const bool has_map_column_container = rem_col;
	static const bool has_column_and_row_swaps = swaps;
};

template<bool is_z2_only, Column_types col_type>
struct Column_compression_options{
	using field_coeff_operators = Zp;
	using id_type = unsigned int;
	using index_type = id_type;
	using pos_type = id_type;
	using dimension_type = int;	//needs to be signed.

	static const bool is_z2 = is_z2_only;
	static const Column_types column_type = col_type;
	static const Column_indexation_types column_indexation_type = Column_indexation_types::CONTAINER;

	static const bool has_matrix_maximal_dimension_access = false;
	static const bool has_column_pairings = false;
	static const bool has_vine_update = false;
	static const bool can_retrieve_representative_cycles = false;
	static const bool is_of_boundary_type = true;
	static const bool has_map_column_container = false;
	static const bool has_column_and_row_swaps = false;
	static const bool has_removable_columns = false;

	static const bool has_column_compression = true;

	static const bool has_row_access = false;
	static const bool has_intrusive_rows = false;
	static const bool has_removable_rows = false;

};

template<bool is_z2_only, Column_types col_type, bool rem_row, bool intr_row>
struct Column_compression_options_with_row_access{
	using field_coeff_operators = Zp;
	using id_type = unsigned int;
	using index_type = id_type;
	using pos_type = id_type;
	using dimension_type = int;	//needs to be signed.

	static const bool is_z2 = is_z2_only;
	static const Column_types column_type = col_type;
	static const Column_indexation_types column_indexation_type = Column_indexation_types::CONTAINER;

	static const bool has_matrix_maximal_dimension_access = false;
	static const bool has_column_pairings = false;
	static const bool has_vine_update = false;
	static const bool can_retrieve_representative_cycles = false;
	static const bool is_of_boundary_type = true;
	static const bool has_map_column_container = false;
	static const bool has_column_and_row_swaps = false;
	static const bool has_removable_columns = false;

	static const bool has_column_compression = true;

	static const bool has_row_access = true;
	static const bool has_intrusive_rows = intr_row;
	static const bool has_removable_rows = rem_row;
};

template<bool is_z2_only, Column_types col_type, bool rem_col, bool swaps, bool pos_idx>
struct Boundary_options{
	using field_coeff_operators = Zp;
	using id_type = unsigned int;
	using index_type = id_type;
	using pos_type = id_type;
	using dimension_type = int;	//needs to be signed.

	static const bool is_z2 = is_z2_only;
	static const Column_types column_type = col_type;
	static const Column_indexation_types column_indexation_type = pos_idx ? Column_indexation_types::CONTAINER : Column_indexation_types::IDENTIFIER;

	static const bool has_matrix_maximal_dimension_access = true;
	static const bool has_column_pairings = true;
	static const bool has_vine_update = false;
	static const bool can_retrieve_representative_cycles = false;
	static const bool is_of_boundary_type = true;
	static const bool has_column_compression = false;

	static const bool has_row_access = false;
	static const bool has_intrusive_rows = false;
	static const bool has_removable_rows = false;

	static const bool has_removable_columns = true;

	static const bool has_map_column_container = rem_col;
	static const bool has_column_and_row_swaps = swaps;
};

template<bool is_z2_only, Column_types col_type, bool rem_row, bool intr_row, bool rem_col, bool swaps, bool pos_idx>
struct Boundary_options_with_row_access{
	using field_coeff_operators = Zp;
	using id_type = unsigned int;
	using index_type = id_type;
	using pos_type = id_type;
	using dimension_type = int;	//needs to be signed.

	static const bool is_z2 = is_z2_only;
	static const Column_types column_type = col_type;
	static const Column_indexation_types column_indexation_type = pos_idx ? Column_indexation_types::CONTAINER : Column_indexation_types::IDENTIFIER;

	static const bool has_matrix_maximal_dimension_access = true;
	static const bool has_column_pairings = true;
	static const bool has_vine_update = false;
	static const bool can_retrieve_representative_cycles = false;
	static const bool is_of_boundary_type = true;
	static const bool has_column_compression = false;

	static const bool has_row_access = true;
	static const bool has_intrusive_rows = intr_row;
	static const bool has_removable_rows = rem_row;

	static const bool has_removable_columns = true;

	static const bool has_map_column_container = rem_col;
	static const bool has_column_and_row_swaps = swaps;
};

template<Column_types col_type, bool rep, bool rem_col, bool pos_idx, bool dim, bool barcode>
struct RU_vine_options{
	using field_coeff_operators = Zp;
	using id_type = unsigned int;
	using index_type = id_type;
	using pos_type = id_type;
	using dimension_type = int;	//needs to be signed.

	static const bool is_z2 = true;
	static const Column_types column_type = col_type;
	static const Column_indexation_types column_indexation_type = pos_idx ? Column_indexation_types::CONTAINER : Column_indexation_types::IDENTIFIER;

	static const bool is_of_boundary_type = true;
	static const bool has_column_compression = false;
	static const bool has_column_and_row_swaps = false;

	//at least one of these two has to be true for the matrix to be RU
	static const bool has_vine_update = true;
	static const bool can_retrieve_representative_cycles = rep;

	static const bool has_row_access = false;
	static const bool has_intrusive_rows = false;
	static const bool has_removable_rows = false;

	static const bool has_removable_columns = true;

	static const bool has_map_column_container = rem_col;
	static const bool has_matrix_maximal_dimension_access = dim;
	static const bool has_column_pairings = barcode;
};

template<bool is_z2_only, Column_types col_type, bool rem_col, bool pos_idx, bool dim, bool barcode>
struct RU_rep_options{
	using field_coeff_operators = Zp;
	using id_type = unsigned int;
	using index_type = id_type;
	using pos_type = id_type;
	using dimension_type = int;	//needs to be signed.

	static const bool is_z2 = is_z2_only;
	static const Column_types column_type = col_type;
	static const Column_indexation_types column_indexation_type = pos_idx ? Column_indexation_types::CONTAINER : Column_indexation_types::IDENTIFIER;

	static const bool is_of_boundary_type = true;
	static const bool has_column_compression = false;
	static const bool has_column_and_row_swaps = false;

	//at least one of these two has to be true for the matrix to be RU
	static const bool has_vine_update = false;	//true-true combination already tested by RU_vine_options
	static const bool can_retrieve_representative_cycles = true;

	static const bool has_row_access = false;
	static const bool has_intrusive_rows = false;
	static const bool has_removable_rows = false;

	static const bool has_removable_columns = true;

	static const bool has_map_column_container = rem_col;
	static const bool has_matrix_maximal_dimension_access = dim;
	static const bool has_column_pairings = barcode;
};

template<Column_types col_type, bool rep, bool rem_row, bool intr_row, bool rem_col, bool pos_idx, bool dim, bool barcode>
struct RU_vine_options_with_row_access{
	using field_coeff_operators = Zp;
	using id_type = unsigned int;
	using index_type = id_type;
	using pos_type = id_type;
	using dimension_type = int;	//needs to be signed.

	static const bool is_z2 = true;
	static const Column_types column_type = col_type;
	static const Column_indexation_types column_indexation_type = pos_idx ? Column_indexation_types::CONTAINER : Column_indexation_types::IDENTIFIER;

	static const bool is_of_boundary_type = true;
	static const bool has_column_compression = false;
	static const bool has_column_and_row_swaps = false;

	//at least one of these two has to be true for the matrix to be RU
	static const bool has_vine_update = true;
	static const bool can_retrieve_representative_cycles = rep;

	static const bool has_row_access = true;
	static const bool has_intrusive_rows = intr_row;
	static const bool has_removable_rows = rem_row;

	static const bool has_removable_columns = true;

	static const bool has_map_column_container = rem_col;
	static const bool has_matrix_maximal_dimension_access = dim;
	static const bool has_column_pairings = barcode;
};

template<bool is_z2_only, Column_types col_type, bool rem_row, bool intr_row, bool rem_col, bool pos_idx, bool dim, bool barcode>
struct RU_rep_options_with_row_access{
	using field_coeff_operators = Zp;
	using id_type = unsigned int;
	using index_type = id_type;
	using pos_type = id_type;
	using dimension_type = int;	//needs to be signed.

	static const bool is_z2 = is_z2_only;
	static const Column_types column_type = col_type;
	static const Column_indexation_types column_indexation_type = pos_idx ? Column_indexation_types::CONTAINER : Column_indexation_types::IDENTIFIER;

	static const bool is_of_boundary_type = true;
	static const bool has_column_compression = false;
	static const bool has_column_and_row_swaps = false;

	//at least one of these two has to be true for the matrix to be RU
	static const bool has_vine_update = false;	//true-true combination already tested by RU_vine_options
	static const bool can_retrieve_representative_cycles = true;

	static const bool has_row_access = true;
	static const bool has_intrusive_rows = intr_row;
	static const bool has_removable_rows = rem_row;

	static const bool has_removable_columns = true;

	static const bool has_map_column_container = rem_col;
	static const bool has_matrix_maximal_dimension_access = dim;
	static const bool has_column_pairings = barcode;
};

template<bool is_z2_only, Column_types col_type, bool rem_col, bool pos_idx, bool dim>
struct Chain_barcode_options{
	using field_coeff_operators = Zp;
	using id_type = unsigned int;
	using index_type = id_type;
	using pos_type = id_type;
	using dimension_type = int;	//needs to be signed.

	static const bool is_z2 = is_z2_only;
	static const Column_types column_type = col_type;
	static const Column_indexation_types column_indexation_type = pos_idx ? Column_indexation_types::POSITION : Column_indexation_types::CONTAINER;

	static const bool is_of_boundary_type = false;
	static const bool has_column_compression = false;
	static const bool has_column_and_row_swaps = false;

	//at least one of these three has to be true for the matrix to be chain
	static const bool has_column_pairings = true;
	static const bool has_vine_update = false;	//true combinations already tested by Chain_vine_options
	static const bool can_retrieve_representative_cycles = false;	//true combinations already tested by Chain_rep_options

	static const bool has_row_access = false;
	static const bool has_intrusive_rows = false;
	static const bool has_removable_rows = false;

	static const bool has_removable_columns = true;

	static const bool has_map_column_container = rem_col;
	static const bool has_matrix_maximal_dimension_access = dim;
};

template<Column_types col_type, bool rep, bool barcode, bool rem_col, bool pos_idx, bool dim>
struct Chain_vine_options{
	using field_coeff_operators = Zp;
	using id_type = unsigned int;
	using index_type = id_type;
	using pos_type = id_type;
	using dimension_type = int;	//needs to be signed.

	static const bool is_z2 = true;
	static const Column_types column_type = col_type;
	static const Column_indexation_types column_indexation_type = pos_idx ? Column_indexation_types::POSITION : Column_indexation_types::CONTAINER;

	static const bool is_of_boundary_type = false;
	static const bool has_column_compression = false;
	static const bool has_column_and_row_swaps = false;

	//at least one of these three has to be true for the matrix to be chain
	static const bool has_column_pairings = barcode;
	static const bool has_vine_update = true;
	static const bool can_retrieve_representative_cycles = rep;

	static const bool has_row_access = false;
	static const bool has_intrusive_rows = false;
	static const bool has_removable_rows = false;

	static const bool has_removable_columns = rem_col;

	static const bool has_map_column_container = rem_col;
	static const bool has_matrix_maximal_dimension_access = dim;
};

template<bool is_z2_only, Column_types col_type, bool barcode, bool rem_col, bool pos_idx, bool dim>
struct Chain_rep_options{
	using field_coeff_operators = Zp;
	using id_type = unsigned int;
	using index_type = id_type;
	using pos_type = id_type;
	using dimension_type = int;	//needs to be signed.

	static const bool is_z2 = is_z2_only;
	static const Column_types column_type = col_type;
	static const Column_indexation_types column_indexation_type = pos_idx ? Column_indexation_types::POSITION : Column_indexation_types::CONTAINER;

	static const bool is_of_boundary_type = false;
	static const bool has_column_compression = false;
	static const bool has_column_and_row_swaps = false;

	//at least one of these three has to be true for the matrix to be chain
	static const bool has_column_pairings = barcode;
	static const bool has_vine_update = false;	//true combinations already tested by Chain_vine_options
	static const bool can_retrieve_representative_cycles = true;

	static const bool has_row_access = false;
	static const bool has_intrusive_rows = false;
	static const bool has_removable_rows = false;

	static const bool has_removable_columns = true;

	static const bool has_map_column_container = rem_col;
	static const bool has_matrix_maximal_dimension_access = dim;
};

template<bool is_z2_only, Column_types col_type, bool rem_row, bool intr_row, bool rem_col, bool pos_idx, bool dim>
struct Chain_barcode_options_with_row_access{
	using field_coeff_operators = Zp;
	using id_type = unsigned int;
	using index_type = id_type;
	using pos_type = id_type;
	using dimension_type = int;	//needs to be signed.

	static const bool is_z2 = is_z2_only;
	static const Column_types column_type = col_type;
	static const Column_indexation_types column_indexation_type = pos_idx ? Column_indexation_types::POSITION : Column_indexation_types::CONTAINER;

	static const bool is_of_boundary_type = false;
	static const bool has_column_compression = false;
	static const bool has_column_and_row_swaps = false;

	//at least one of these three has to be true for the matrix to be chain
	static const bool has_column_pairings = true;
	static const bool has_vine_update = false;	//true combinations already tested by Chain_vine_options
	static const bool can_retrieve_representative_cycles = false;	//true combinations already tested by Chain_rep_options

	static const bool has_row_access = true;
	static const bool has_intrusive_rows = intr_row;
	static const bool has_removable_rows = rem_row;

	static const bool has_removable_columns = true;

	static const bool has_map_column_container = rem_col;
	static const bool has_matrix_maximal_dimension_access = dim;
};

template<Column_types col_type, bool rep, bool barcode, bool rem_row, bool intr_row, bool rem_col, bool pos_idx, bool dim>
struct Chain_vine_options_with_row_access{
	using field_coeff_operators = Zp;
	using id_type = unsigned int;
	using index_type = id_type;
	using pos_type = id_type;
	using dimension_type = int;	//needs to be signed.

	static const bool is_z2 = true;
	static const Column_types column_type = col_type;
	static const Column_indexation_types column_indexation_type = pos_idx ? Column_indexation_types::POSITION : Column_indexation_types::CONTAINER;

	static const bool is_of_boundary_type = false;
	static const bool has_column_compression = false;
	static const bool has_column_and_row_swaps = false;

	//at least one of these three has to be true for the matrix to be chain
	static const bool has_column_pairings = barcode;
	static const bool has_vine_update = true;
	static const bool can_retrieve_representative_cycles = rep;

	static const bool has_row_access = true;
	static const bool has_intrusive_rows = intr_row;
	static const bool has_removable_rows = rem_row;

	static const bool has_removable_columns = rem_col;

	static const bool has_map_column_container = rem_col;
	static const bool has_matrix_maximal_dimension_access = dim;
};

template<bool is_z2_only, Column_types col_type, bool barcode, bool rem_row, bool intr_row, bool rem_col, bool pos_idx, bool dim>
struct Chain_rep_options_with_row_access{
	using field_coeff_operators = Zp;
	using id_type = unsigned int;
	using index_type = id_type;
	using pos_type = id_type;
	using dimension_type = int;	//needs to be signed.

	static const bool is_z2 = is_z2_only;
	static const Column_types column_type = col_type;
	static const Column_indexation_types column_indexation_type = pos_idx ? Column_indexation_types::POSITION : Column_indexation_types::CONTAINER;

	static const bool is_of_boundary_type = false;
	static const bool has_column_compression = false;
	static const bool has_column_and_row_swaps = false;

	//at least one of these three has to be true for the matrix to be chain
	static const bool has_column_pairings = barcode;
	static const bool has_vine_update = false;	//true combinations already tested by Chain_vine_options
	static const bool can_retrieve_representative_cycles = true;

	static const bool has_row_access = true;
	static const bool has_intrusive_rows = intr_row;
	static const bool has_removable_rows = rem_row;

	static const bool has_removable_columns = true;

	static const bool has_map_column_container = rem_col;
	static const bool has_matrix_maximal_dimension_access = dim;
};

#endif // PM_MATRIX_TESTS_OPTIONS_H

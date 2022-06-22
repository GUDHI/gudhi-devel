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

#include "Z2_field.h"

namespace Gudhi {
namespace persistence_matrix {

enum Column_type {
	LIST,
	SET,
	HEAP,
	VECTOR,
	UNORDERED_SET,
	BOOST_LIST,
	BOOST_SET
};

template<class Field_type = Z2_field_element, Column_type col_type = Column_type::SET, bool separated_by_dimension = false, bool parallelizable = false>
struct Default_options{
	using field_coeff_type = Field_type;
	static const Column_type column_type = col_type;

	static const bool is_separated_by_dimension = separated_by_dimension;
	static const bool is_parallelizable = parallelizable;

	static const bool has_row = false;
	static const bool has_column_pairings = false;
	static const bool has_general_vine_update = false;
	static const bool has_z_eq_1_vine_update = false;
	static const bool can_retrieve_representative_cycles = false;
	static const bool has_column_compression = false;
	static const bool is_double_linked = true;
	static const bool is_of_boundary_type = true;
	static const bool is_of_chain_type = false;
	static const bool has_removable_columns = false;
};

template<Column_type col_type = Column_type::BOOST_SET, bool separated_by_dimension = false, bool parallelizable = false>
struct Zigzag_options : Default_options<Z2_field_element, col_type, separated_by_dimension, parallelizable>{
	static const bool has_row = true;
	static const bool has_column_pairings = true;
	static const bool has_z_eq_1_vine_update = true;
	static const bool is_of_boundary_type = false;
	static const bool is_of_chain_type = true;
	static const bool has_removable_columns = true;
};

template<class Field_type = Z2_field_element, Column_type col_type = Column_type::SET, bool separated_by_dimension = false, bool parallelizable = false>
struct Representative_cycles_options : Default_options<Field_type, col_type, separated_by_dimension, parallelizable>{
	static const bool can_retrieve_representative_cycles = true;
};

template<Column_type col_type = Column_type::SET, bool separated_by_dimension = false, bool parallelizable = false>
struct Multi_persistence_options : Default_options<Z2_field_element, col_type, separated_by_dimension, parallelizable>{
	static const bool has_column_pairings = true;
	static const bool has_general_vine_update = true;
};

template<class Field_type = Z2_field_element>
struct Cohomology_persistence_options : Default_options<Field_type, Column_type::BOOST_LIST>{
	static const bool has_row = true;
	static const bool has_column_compression = true;
	static const bool is_of_boundary_type = false;
	static const bool is_of_chain_type = true;
	static const bool has_removable_columns = true;
};

} //namespace persistence_matrix
} //namespace Gudhi

#endif // OPTIONS_INCLUDED


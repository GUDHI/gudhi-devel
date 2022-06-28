/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef Z2_CELL_SET_COLUMN_H
#define Z2_CELL_SET_COLUMN_H

#include "z2_reduced_cell_column_with_row.h"

namespace Gudhi {
namespace persistence_matrix {

class Z2_reduced_cell_set_column_with_row : Z2_reduced_cell_column_with_row<Column_types::SET>
{
public:
	Z2_reduced_cell_set_column_with_row();
	template<class Chain_type>
	Z2_reduced_cell_set_column_with_row(index chainIndex, Chain_type& chain);
	template<class Chain_type>
	Z2_reduced_cell_set_column_with_row(index chainIndex, Chain_type& chain, index pairedColumnIndex);
	Z2_reduced_cell_set_column_with_row(const Z2_reduced_cell_set_column_with_row& other);
};

inline Z2_reduced_cell_set_column_with_row::Z2_reduced_cell_set_column_with_row()
	: Z2_reduced_cell_column_with_row<Column_types::SET>()
{}

template<class Chain_type>
inline Z2_reduced_cell_set_column_with_row::Z2_reduced_cell_set_column_with_row(
		index chainIndex, Chain_type& chain)
	: Z2_reduced_cell_column_with_row<Column_types::SET>(chainIndex, chain)
{}

template<class Chain_type>
inline Z2_reduced_cell_set_column_with_row::Z2_reduced_cell_set_column_with_row(
		index chainIndex, Chain_type& chain, index pairedColumnIndex)
	: Z2_reduced_cell_column_with_row<Column_types::SET>(chainIndex, chain, pairedColumnIndex)
{}

inline Z2_reduced_cell_set_column_with_row::Z2_reduced_cell_set_column_with_row(
		const Z2_reduced_cell_set_column_with_row& other)
	: Z2_reduced_cell_column_with_row<Column_types::SET>(other)
{}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // Z2_CELL_LIST_COLUMN_H

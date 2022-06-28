/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CELL_SET_COLUMN_H
#define CELL_SET_COLUMN_H

#include "reduced_cell_column_with_row.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Field_element_type = Zp_field_element<11> >
class Reduced_cell_set_column_with_row : Reduced_cell_column_with_row<Column_types::SET,Field_element_type>
{
public:
	Reduced_cell_set_column_with_row();
	template<class Chain_type>
	Reduced_cell_set_column_with_row(index chainIndex, Chain_type& rowIndices, Chain_type& values);
	template<class Chain_type>
	Reduced_cell_set_column_with_row(index chainIndex, Chain_type& rowIndices, Chain_type& values, index pairedColumnIndex);
	Reduced_cell_set_column_with_row(const Reduced_cell_set_column_with_row& other);
};

template<class Field_element_type>
inline Reduced_cell_set_column_with_row<Field_element_type>::Reduced_cell_set_column_with_row()
	: Reduced_cell_column_with_row<Column_types::SET,Field_element_type>()
{}

template<class Field_element_type>
template<class Chain_type>
inline Reduced_cell_set_column_with_row<Field_element_type>::Reduced_cell_set_column_with_row(
		index chainIndex, Chain_type& rowIndices, Chain_type &values)
	: Reduced_cell_column_with_row<Column_types::SET,Field_element_type>(chainIndex, rowIndices, values)
{}

template<class Field_element_type>
template<class Chain_type>
inline Reduced_cell_set_column_with_row<Field_element_type>::Reduced_cell_set_column_with_row(
		index chainIndex, Chain_type& rowIndices, Chain_type &values, index pairedColumnIndex)
	: Reduced_cell_column_with_row<Column_types::SET,Field_element_type>(chainIndex, rowIndices, values, pairedColumnIndex)
{}

template<class Field_element_type>
inline Reduced_cell_set_column_with_row<Field_element_type>::Reduced_cell_set_column_with_row(
		const Reduced_cell_set_column_with_row& other)
	: Reduced_cell_column_with_row<Column_types::SET,Field_element_type>(other)
{}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // CELL_SET_COLUMN_H

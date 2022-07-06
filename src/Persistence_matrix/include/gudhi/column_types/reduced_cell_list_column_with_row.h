/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CELL_LIST_COLUMN_H
#define CELL_LIST_COLUMN_H

#include "reduced_cell_column_with_row.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Field_element_type = Zp_field_element<11> >
class Reduced_cell_list_column_with_row : Reduced_cell_column_with_row<Column_types::LIST, Field_element_type>
{
public:
	Reduced_cell_list_column_with_row();
	template<class Chain_type>
	Reduced_cell_list_column_with_row(index chainIndex, Chain_type& chain);
	template<class Chain_type>
	Reduced_cell_list_column_with_row(index chainIndex, Chain_type& chain, index pairedColumnIndex);
	Reduced_cell_list_column_with_row(const Reduced_cell_list_column_with_row& other);
};

template<class Field_element_type>
inline Reduced_cell_list_column_with_row<Field_element_type>::Reduced_cell_list_column_with_row()
	: Reduced_cell_column_with_row<Column_types::LIST,Field_element_type>()
{}

template<class Field_element_type>
template<class Chain_type>
inline Reduced_cell_list_column_with_row<Field_element_type>::Reduced_cell_list_column_with_row(
		index chainIndex, Chain_type& chain)
	: Reduced_cell_column_with_row<Column_types::LIST,Field_element_type>(chainIndex, chain)
{}

template<class Field_element_type>
template<class Chain_type>
inline Reduced_cell_list_column_with_row<Field_element_type>::Reduced_cell_list_column_with_row(
		index chainIndex, Chain_type& chain, index pairedColumnIndex)
	: Reduced_cell_column_with_row<Column_types::LIST,Field_element_type>(chainIndex, chain, pairedColumnIndex)
{}

template<class Field_element_type>
inline Reduced_cell_list_column_with_row<Field_element_type>::Reduced_cell_list_column_with_row(
		const Reduced_cell_list_column_with_row& other)
	: Reduced_cell_column_with_row<Column_types::LIST,Field_element_type>(other)
{}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // CELL_LIST_COLUMN_H

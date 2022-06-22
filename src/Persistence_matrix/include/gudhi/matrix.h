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

#include <vector>

#include "options.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Options = Default_options<> >
class Matrix
{
public:
//	using Options::Field_type;
//	using Options::column_type;
//	using Options::is_double_linked;

//	using Options::is_of_boundary_type;
//	using Options::is_of_chain_type;

//	using Options::is_separated_by_dimension;	//later
//	using Options::is_parallelizable;			//later
//	using Options::has_column_compression;		//later
//	using Options::has_row;
//	using Options::has_removable_columns;

//	using Options::has_general_vine_update;
//	using Options::has_z_eq_1_vine_update;
//	using Options::can_retrieve_representative_cycles;
//	using Options::has_column_pairings;








//	using index = unsigned int;
//	using boundary_type = std::vector<index>;
//	using boundary_matrix = std::vector<boundary_type>;

//	Matrix();
//	Matrix(boundary_matrix& boundaries);	//simplex indices have to start at 0 and be consecutifs
//	Matrix(int numberOfColumns);
//	Matrix(Matrix &matrixToCopy);
//	Matrix(Matrix&& other) noexcept;

//	template<class Chain_type>
//	void insert_boundary(index simplexIndex, Chain_type& boundary);
//	template<class Chain_type>
//	void insert_boundary(index simplexIndex, Chain_type& boundary, std::vector<index>& chains_in_F);
//	dimension_type get_dimension(index simplexIndex) const;
//	unsigned int get_number_of_simplices() const;

//	index vine_swap_with_z_eq_1_case(index columnIndex1, index columnIndex2);	//returns index which was not modified, ie new i+1
//	index vine_swap(index columnIndex1, index columnIndex2);					//returns index which was not modified, ie new i+1
//	const barcode_type& get_current_barcode();
//	void erase_last(index simplexIndex);

//	void add_to(index sourceColumnIndex, index targetColumnIndex);

//	index get_column_with_pivot(index simplexIndex);
//	void get_row(index rowIndex, std::vector<index>& row);
//	index get_pivot(index columnIndex);

//	bool is_paired(index columnIndex);
//	index get_paired_column(index columnIndex);

//	void print_matrices();  //for debug

//	Matrix& operator=(Matrix other);
//	template<class Friend_options>
//	friend void swap(Matrix& matrix1, Matrix& matrix2);

private:

};

} //namespace persistence_matrix
} //namespace Gudhi

#endif // MASTER_MATRIX_H

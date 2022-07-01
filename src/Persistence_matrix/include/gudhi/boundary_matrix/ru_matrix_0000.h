/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef RU_MATRIX_0000_H
#define RU_MATRIX_0000_H

#include "../utilities.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class RU_matrix;

template<class Master_matrix>
class RU_matrix : Master_matrix::RU_pairing_option, Master_matrix::RU_vine_swap_option, Master_matrix::RU_representative_cycles_option{
public:
	using Column_type = typename Master_matrix::Column_type;
	using Row_type = void;

	RU_matrix()
		: Master_matrix::RU_pairing_option(),
		  Master_matrix::RU_vine_swap_option(reducedMatrixR_, mirrorMatrixU_),
		  Master_matrix::RU_representative_cycles_option(reducedMatrixR_, mirrorMatrixU_){};
	RU_matrix(std::vector<boundary_type>& orderedBoundaries);
	RU_matrix(unsigned int numberOfColumns);
	RU_matrix(RU_matrix& matrixToCopy);
	RU_matrix(RU_matrix&& other) noexcept;

	void insert_boundary(boundary_type& boundary);
	Column_type& get_column(index columnIndex);
	Row_type get_row(index rowIndex);
	void erase_last(index simplexIndex);

	dimension_type get_max_dimension() const;
	unsigned int get_number_of_columns() const;

	dimension_type get_column_dimension(index columnIndex) const;

	void add_to(index sourceColumnIndex, index targetColumnIndex);

	void zero_cell(index columnIndex, index rowIndex);
	void zero_column(index columnIndex);
	bool is_zero_cell(index columnIndex, index rowIndex) const;
	bool is_zero_column(index columnIndex);

	index get_column_with_pivot(index simplexIndex);
	index get_pivot(index columnIndex);

	RU_matrix& operator=(RU_matrix other);
	template<class Friend_master_matrix>
	friend void swap(RU_matrix<Friend_master_matrix>& matrix1,
					 RU_matrix<Friend_master_matrix>& matrix2);

	void print();  //for debug

private:
	typename Master_matrix::Base_matrix reducedMatrixR_;
	typename Master_matrix::Base_matrix mirrorMatrixU_;
};

} //namespace persistence_matrix
} //namespace Gudhi

#endif // RU_MATRIX_0000_H

/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022-23 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_BASE_MATRIX_RA_H
#define PM_BASE_MATRIX_RA_H

#include <utility>	//std::move

namespace Gudhi {
namespace persistence_matrix {

struct Dummy_matrix_row_access{
	Dummy_matrix_row_access(){};
	Dummy_matrix_row_access([[maybe_unused]] unsigned int numberOfColumns){};

	friend void swap([[maybe_unused]] Dummy_matrix_row_access& d1, [[maybe_unused]] Dummy_matrix_row_access& d2){}
};

template<typename Row_type, typename Row_container_type, bool has_removable_rows, typename id_index>
class Matrix_row_access
{
public:
	Matrix_row_access() : rows_(new Row_container_type()){};
	Matrix_row_access(unsigned int numberOfColumns) : rows_(new Row_container_type())
	{
		if constexpr (!has_removable_rows){
			rows_->resize(numberOfColumns);
		}
	}
	Matrix_row_access(const Matrix_row_access& toCopy) : rows_(new Row_container_type())	//as the matrix is rebuild, the rows should not be copied.
	{
		if constexpr (!has_removable_rows){
			rows_->resize(toCopy.rows_->size());
		}
	}
	Matrix_row_access(Matrix_row_access&& other) noexcept : rows_(std::exchange(other.rows_, nullptr)){}
	~Matrix_row_access(){
		delete rows_;
	}

	//get_row(rowIndex) --> simplex ID (=/= columnIndex)
	Row_type& get_row(id_index rowIndex){
		if constexpr (has_removable_rows) {
			return rows_->at(rowIndex);
		} else {
			return rows_->operator[](rowIndex);
		}
	}
	const Row_type& get_row(id_index rowIndex) const {
		if constexpr (has_removable_rows) {
			return rows_->at(rowIndex);
		} else {
			return rows_->operator[](rowIndex);
		}
	}
	//assumes that row is empty.
	void erase_row(id_index rowIndex){
		static_assert(has_removable_rows, "'erase_row' is not implemented for the chosen options.");

		auto it = rows_->find(rowIndex);
		if (it != rows_->end() && it->second.empty()) rows_->erase(it);
	}

	Matrix_row_access& operator=(const Matrix_row_access& other){
		if constexpr (has_removable_rows)
			rows_->reserve(other.rows_->size());
		else 
			rows_->resize(other.rows_->size());
		return *this;
	}
	friend void swap(Matrix_row_access& matrix1, Matrix_row_access& matrix2){
		std::swap(matrix1.rows_, matrix2.rows_);
	}

protected:
	//A pointer to faciliate column swaps when two matrices are swapped.
	Row_container_type *rows_;	//has to be destroyed after matrix_, therefore has to be inherited.
};

} //namespace persistence_matrix
} //namespace Gudhi

#endif // PM_BASE_MATRIX_RA_H

/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef BASE_SWAP_H
#define BASE_SWAP_H

#include <utility>
#include <vector>

#include "../utilities.h"  //type definitions

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class Base_swap
{
public:
	void swap_columns(index columnIndex1, index columnIndex2);
	void swap_rows(index rowIndex1, index rowIndex2);
	void swap_at_indices(index index1, index index2);

	Base_swap& operator=(Base_swap other);
	template<class Friend_master_matrix>
	friend void swap(Base_swap<Friend_master_matrix>& base1,
					 Base_swap<Friend_master_matrix>& base2);

	using matrix_type = typename Master_matrix::column_container_type;

	Base_swap(matrix_type &matrix);
	Base_swap(matrix_type &matrix, unsigned int numberOfColumns);
	Base_swap(const Base_swap& matrixToCopy);
	Base_swap(Base_swap&& other) noexcept;

protected:
	using index_dictionnary_type = typename Master_matrix::template dictionnary_type<unsigned int>;
	using row_dictionnary_type = typename Master_matrix::template dictionnary_type<index>;

	index_dictionnary_type indexToRow_;
	row_dictionnary_type rowToIndex_;
	bool rowSwapped_;
	static constexpr bool isActive_ = true;

	void _orderRows();

private:
	matrix_type &matrix_;
};

template<class Master_matrix>
inline Base_swap<Master_matrix>::Base_swap(matrix_type &matrix)
	: rowSwapped_(false), matrix_(matrix)
{}

template<class Master_matrix>
inline Base_swap<Master_matrix>::Base_swap(matrix_type &matrix, unsigned int numberOfColumns)
	: indexToRow_(numberOfColumns),
	  rowToIndex_(numberOfColumns),
	  rowSwapped_(false),
	  matrix_(matrix)
{}

template<class Master_matrix>
inline Base_swap<Master_matrix>::Base_swap(const Base_swap<Master_matrix>& matrixToCopy)
	: indexToRow_(matrixToCopy.indexToRow_),
	  rowToIndex_(matrixToCopy.rowToIndex_),
	  rowSwapped_(matrixToCopy.rowSwapped_),
	  matrix_(matrixToCopy.matrix_)
{}

template<class Master_matrix>
inline Base_swap<Master_matrix>::Base_swap(Base_swap<Master_matrix> &&other) noexcept
	: indexToRow_(std::move(other.indexToRow_)),
	  rowToIndex_(std::move(other.rowToIndex_)),
	  rowSwapped_(std::exchange(other.rowSwapped_, 0)),
	  matrix_(other.matrix_)
{}

template<class Master_matrix>
inline void Base_swap<Master_matrix>::swap_columns(index columnIndex1, index columnIndex2)
{
	std::swap(matrix_.at(columnIndex1), matrix_.at(columnIndex2));
}

template<class Master_matrix>
inline void Base_swap<Master_matrix>::swap_rows(index rowIndex1, index rowIndex2)
{
	rowSwapped_ = true;
	std::swap(rowToIndex_.at(indexToRow_.at(rowIndex1)), rowToIndex_.at(indexToRow_.at(rowIndex2)));
	std::swap(indexToRow_.at(rowIndex1), indexToRow_.at(rowIndex2));
}

template<class Master_matrix>
inline void Base_swap<Master_matrix>::swap_at_indices(index index1, index index2)
{
	swap_columns(index1, index2);
	swap_rows(index1, index2);
}

template<class Master_matrix>
inline Base_swap<Master_matrix> &Base_swap<Master_matrix>::operator=(Base_swap<Master_matrix> other)
{
	std::swap(matrix_, other.matrix_);
	std::swap(indexToRow_, other.indexToRow_);
	std::swap(rowToIndex_, other.rowToIndex_);
	std::swap(rowSwapped_, other.rowSwapped_);
	return *this;
}

template<class Master_matrix>
inline void Base_swap<Master_matrix>::_orderRows()
{
	for (unsigned int i = 0; i < matrix_.size(); i++){
		matrix_.at(i).reorder(rowToIndex_);
	}
	for (unsigned int i = 0; i < matrix_.size(); i++){
		indexToRow_.at(i) = i;
		rowToIndex_.at(i) = i;
	}
	rowSwapped_ = false;
}

template<class Friend_master_matrix>
inline void swap(Base_swap<Friend_master_matrix>& base1, Base_swap<Friend_master_matrix>& base2)
{
	std::swap(base1.matrix_, base2.matrix_);
	std::swap(base1.indexToRow_, base2.indexToRow_);
	std::swap(base1.rowToIndex_, base2.rowToIndex_);
	std::swap(base1.rowSwapped_, base2.rowSwapped_);
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // BASE_SWAP_H

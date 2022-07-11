/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef Chain_VINE_SWAP_H
#define Chain_VINE_SWAP_H

#include <utility>
#include <set>

#include "../utilities.h"
#include "chain_pairing.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class Chain_vine_swap : public Chain_pairing<Master_matrix>
{
public:
	index vine_swap_with_z_eq_1_case(index columnIndex1, index columnIndex2);	//returns index which was not modified, ie new i+1
	index vine_swap(index columnIndex1, index columnIndex2);					//returns index which was not modified, ie new i+1

	Chain_vine_swap& operator=(Chain_vine_swap other);
	template<class Friend_master_matrix>
	friend void swap(Chain_vine_swap<Friend_master_matrix>& swap1,
					 Chain_vine_swap<Friend_master_matrix>& swap2);

protected:
	using matrix_type = typename Master_matrix::column_container_type;
	using dictionnary_type = typename Master_matrix::template dictionnary_type<index>;

	Chain_vine_swap(matrix_type& matrix);
	Chain_vine_swap(Chain_vine_swap &matrixToCopy);
	Chain_vine_swap(Chain_vine_swap&& other) noexcept;

	dictionnary_type pivotToPosition_;
	static constexpr bool isActive_ = true;

private:
	using Chain_pairing<Master_matrix>::barcode_;
	using Chain_pairing<Master_matrix>::indexToBar_;
	using Column_type = typename Master_matrix::Column_type;

	matrix_type& matrix_;

	void _add_to(const typename Column_type::Column& column, std::set<index>& set);
	bool _is_negative_in_pair(index simplexIndex);

	void _positive_transpose(index simplexIndex1, index simplexIndex2);
	void _negative_transpose(index simplexIndex1, index simplexIndex2);
	void _positive_negative_transpose(index simplexIndex1, index simplexIndex2);
	void _negative_positive_transpose(index simplexIndex1, index simplexIndex2);
	index _positive_vine_swap(index columnIndex1, index columnIndex2);
	index _positive_negative_vine_swap(index columnIndex1, index columnIndex2);
	index _negative_positive_vine_swap(index columnIndex1, index columnIndex2);
	index _negative_vine_swap(index columnIndex1, index columnIndex2);

	//TODO: update use of matrix and indexToBar
};

template<class Master_matrix>
inline Chain_vine_swap<Master_matrix>::Chain_vine_swap(matrix_type &matrix)
	: matrix_(matrix), pivotToPosition_(matrix.size()), Chain_pairing<Master_matrix>()
{}

template<class Master_matrix>
inline Chain_vine_swap<Master_matrix>::Chain_vine_swap(
		Chain_vine_swap &matrixToCopy)
	: matrix_(matrixToCopy.matrix_), pivotToPosition_(matrixToCopy.pivotToPosition_), Chain_pairing<Master_matrix>(matrixToCopy)
{}

template<class Master_matrix>
inline Chain_vine_swap<Master_matrix>::Chain_vine_swap(Chain_vine_swap<Master_matrix> &&other) noexcept
	: matrix_(std::move(other.reducedMatrixR_)),
	  pivotToPosition_(std::move(other.mirrorMatrixU_)),
	  Chain_pairing<Master_matrix>(std::move(other))
{}

template<class Master_matrix>
inline index Chain_vine_swap<Master_matrix>::vine_swap_with_z_eq_1_case(index columnIndex1, index columnIndex2)
{
	index pivot1 = pivotToPosition_.at(matrix_.get_pivot(columnIndex1));
	index pivot2 = pivotToPosition_.at(matrix_.get_pivot(columnIndex2));

	if (_is_negative_in_pair(pivot1) && _is_negative_in_pair(pivot2))
		return _negative_vine_swap(columnIndex1, columnIndex2);

	if (_is_negative_in_pair(pivot1))
		return _negative_positive_vine_swap(columnIndex1, columnIndex2);

	if (_is_negative_in_pair(pivot2))
		return _positive_negative_vine_swap(columnIndex1, columnIndex2);

	return _positive_vine_swap(columnIndex1, columnIndex2);
}

template<class Master_matrix>
inline index Chain_vine_swap<Master_matrix>::vine_swap(index columnIndex1, index columnIndex2)
{
	index pivot1 = matrix_.get_pivot(columnIndex1);
	index pivot2 = matrix_.get_pivot(columnIndex2);
	index rowIndex1 = pivotToPosition_.at(pivot1);
	index rowIndex2 = pivotToPosition_.at(pivot2);

	if (_is_negative_in_pair(rowIndex1) && _is_negative_in_pair(rowIndex2)){
		if (!matrix_.contains(columnIndex2, matrix_.get_lowest_simplex_index(columnIndex1))){
			_negative_transpose(rowIndex1, rowIndex2);
			matrix_.swap_independent_rows(pivot1, pivot2);
			return columnIndex1;
		}
		return _negative_vine_swap(columnIndex1, columnIndex2);
	}

	if (_is_negative_in_pair(rowIndex1)){
		if (!matrix_.contains(columnIndex2, matrix_.get_lowest_simplex_index(columnIndex1))){
			_negative_positive_transpose(rowIndex1, rowIndex2);
			matrix_.swap_independent_rows(pivot1, pivot2);
			return columnIndex1;
		}
		return _negative_positive_vine_swap(columnIndex1, columnIndex2);
	}

	if (_is_negative_in_pair(rowIndex2)){
		if (!matrix_.contains(columnIndex2, matrix_.get_lowest_simplex_index(columnIndex1))){
			_positive_negative_transpose(rowIndex1, rowIndex2);
			matrix_.swap_independent_rows(pivot1, pivot2);
			return columnIndex1;
		}
		return _positive_negative_vine_swap(columnIndex1, columnIndex2);
	}

	if (!matrix_.contains(columnIndex2, matrix_.get_lowest_simplex_index(columnIndex1))){
		_positive_transpose(rowIndex1, rowIndex2);
		matrix_.swap_independent_rows(pivot1, pivot2);
		return columnIndex1;
	}
	return _positive_vine_swap(columnIndex1, columnIndex2);
}

template<class Master_matrix>
inline Chain_vine_swap<Master_matrix> &Chain_vine_swap<Master_matrix>::operator=(
		Chain_vine_swap<Master_matrix> other)
{
	std::swap(matrix_, other.matrix_);
	std::swap(pivotToPosition_, other.pivotToPosition_);
	Chain_pairing<Master_matrix>::operator=(other);
	return *this;
}

template<class Master_matrix>
inline void Chain_vine_swap<Master_matrix>::_add_to(
		const typename Column_type::Column& column, std::set<index>& set)
{
	std::pair<std::set<index>::iterator,bool> res_insert;
	for (const typename Column_type::Cell &cell : column) {
		res_insert = set.insert(cell.get_row_index());
		if (!res_insert.second) {
			set.erase(res_insert.first);
		}
	}
}

template<class Master_matrix>
inline bool Chain_vine_swap<Master_matrix>::_is_negative_in_pair(index simplexIndex)
{
	return indexToBar_.at(simplexIndex)->death == static_cast<int>(simplexIndex);
}

template<class Master_matrix>
inline void Chain_vine_swap<Master_matrix>::_positive_transpose(index simplexIndex1, index simplexIndex2)
{
	indexToBar_.at(simplexIndex1)->birth = simplexIndex2;
	indexToBar_.at(simplexIndex2)->birth = simplexIndex1;
	std::swap(indexToBar_.at(simplexIndex1), indexToBar_.at(simplexIndex2));
}

template<class Master_matrix>
inline void Chain_vine_swap<Master_matrix>::_negative_transpose(index simplexIndex1, index simplexIndex2)
{
	indexToBar_.at(simplexIndex1)->death = simplexIndex2;
	indexToBar_.at(simplexIndex2)->death = simplexIndex1;
	std::swap(indexToBar_.at(simplexIndex1), indexToBar_.at(simplexIndex2));
}

template<class Master_matrix>
inline void Chain_vine_swap<Master_matrix>::_positive_negative_transpose(index simplexIndex1, index simplexIndex2)
{
	indexToBar_.at(simplexIndex1)->birth = simplexIndex2;
	indexToBar_.at(simplexIndex2)->death = simplexIndex1;
	std::swap(indexToBar_.at(simplexIndex1), indexToBar_.at(simplexIndex2));
}

template<class Master_matrix>
inline void Chain_vine_swap<Master_matrix>::_negative_positive_transpose(index simplexIndex1, index simplexIndex2)
{
	indexToBar_.at(simplexIndex1)->death = simplexIndex2;
	indexToBar_.at(simplexIndex2)->birth = simplexIndex1;
	std::swap(indexToBar_.at(simplexIndex1), indexToBar_.at(simplexIndex2));
}

template<class Master_matrix>
inline index Chain_vine_swap<Master_matrix>::_positive_vine_swap(index columnIndex1, index columnIndex2)
{
	index pivot1 = pivotToPosition_.at(matrix_.get_pivot(columnIndex1));
	index pivot2 = pivotToPosition_.at(matrix_.get_pivot(columnIndex2));

	std::swap(pivotToPosition_.at(matrix_.get_pivot(columnIndex1)), pivotToPosition_.at(matrix_.get_pivot(columnIndex2)));

	if ((indexToBar_.at(pivot1)->death < indexToBar_.at(pivot2)->death && indexToBar_.at(pivot1)->death != -1)
			|| (indexToBar_.at(pivot1)->death != -1 && indexToBar_.at(pivot2)->death == -1)
			|| (indexToBar_.at(pivot1)->birth < indexToBar_.at(pivot2)->birth && indexToBar_.at(pivot1)->death == -1 && indexToBar_.at(pivot2)->death == -1))
	{
		if (indexToBar_.at(pivot2)->death != -1)
			matrix_.add_to(matrix_.get_paired_column(columnIndex1), matrix_.get_paired_column(columnIndex2));
		matrix_.add_to(columnIndex1, columnIndex2);

		_positive_transpose(pivot1, pivot2);

		return columnIndex1;
	}

	if (indexToBar_.at(pivot1)->death != -1 && indexToBar_.at(pivot2)->death != -1)
		matrix_.add_to(matrix_.get_paired_column(columnIndex2), matrix_.get_paired_column(columnIndex1));
	matrix_.add_to(columnIndex2, columnIndex1);

	return columnIndex2;
}

template<class Master_matrix>
inline index Chain_vine_swap<Master_matrix>::_positive_negative_vine_swap(index columnIndex1, index columnIndex2)
{
	matrix_.add_to(columnIndex1, columnIndex2);
	_positive_negative_transpose(pivotToPosition_.at(matrix_.get_pivot(columnIndex1)), pivotToPosition_.at(matrix_.get_pivot(columnIndex2)));
	std::swap(pivotToPosition_.at(matrix_.get_pivot(columnIndex1)), pivotToPosition_.at(matrix_.get_pivot(columnIndex2)));

	return columnIndex1;
}

template<class Master_matrix>
inline index Chain_vine_swap<Master_matrix>::_negative_positive_vine_swap(index columnIndex1, index columnIndex2)
{
	matrix_.add_to(columnIndex2, columnIndex1);
	std::swap(pivotToPosition_.at(matrix_.get_pivot(columnIndex1)), pivotToPosition_.at(matrix_.get_pivot(columnIndex2)));

	return columnIndex2;
}

template<class Master_matrix>
inline index Chain_vine_swap<Master_matrix>::_negative_vine_swap(index columnIndex1, index columnIndex2)
{
	index pivot1 = pivotToPosition_.at(matrix_.get_pivot(columnIndex1));
	index pivot2 = pivotToPosition_.at(matrix_.get_pivot(columnIndex2));

	index pairedIndex1 = matrix_.get_paired_column(columnIndex1);
	index pairedIndex2 = matrix_.get_paired_column(columnIndex2);

	std::swap(pivotToPosition_.at(matrix_.get_pivot(columnIndex1)), pivotToPosition_.at(matrix_.get_pivot(columnIndex2)));

	if (indexToBar_.at(pivot1)->birth < indexToBar_.at(pivot2)->birth)
	{
		matrix_.add_to(pairedIndex1, pairedIndex2);
		matrix_.add_to(columnIndex1, columnIndex2);
		_negative_transpose(pivot1, pivot2);

		return columnIndex1;
	}

	matrix_.add_to(pairedIndex2, pairedIndex1);
	matrix_.add_to(columnIndex2, columnIndex1);

	return columnIndex2;
}

template<class Friend_master_matrix>
inline void swap(Chain_vine_swap<Friend_master_matrix>& swap1,
				 Chain_vine_swap<Friend_master_matrix>& swap2)
{
	std::swap(swap1.reducedMatrixR_, swap2.reducedMatrixR_);
	std::swap(swap1.mirrorMatrixU_, swap2.mirrorMatrixU_);
	std::swap(static_cast<Chain_pairing<Friend_master_matrix> >(swap1),
			  static_cast<Chain_pairing<Friend_master_matrix> >(swap2));
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // Chain_VINE_SWAP_H

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

#include "../utilities/utilities.h"
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

	using matrix_type = typename Master_matrix::column_container_type;
	using dictionnary_type = typename Master_matrix::template dictionnary_type<index>;

	Chain_vine_swap(matrix_type& matrix);
	Chain_vine_swap(const Chain_vine_swap &matrixToCopy);
	Chain_vine_swap(Chain_vine_swap&& other) noexcept;

protected:
	dictionnary_type pivotToPosition_;
	static constexpr bool isActive_ = true;

private:
	using CP = Chain_pairing<Master_matrix>;
	using Column_type = typename Master_matrix::Column_type;

	matrix_type* matrix_;

	void _add_to(const typename Column_type::Column_type& column, std::set<index>& set);
	bool _is_negative_in_pair(index simplexIndex);

	void _positive_transpose(index simplexIndex1, index simplexIndex2);
	void _negative_transpose(index simplexIndex1, index simplexIndex2);
	void _positive_negative_transpose(index simplexIndex1, index simplexIndex2);
	void _negative_positive_transpose(index simplexIndex1, index simplexIndex2);
	index _positive_vine_swap(index columnIndex1, index columnIndex2);
	index _positive_negative_vine_swap(index columnIndex1, index columnIndex2);
	index _negative_positive_vine_swap(index columnIndex1, index columnIndex2);
	index _negative_vine_swap(index columnIndex1, index columnIndex2);

	int& _death(index simplexIndex);
	int& _birth(index simplexIndex);

	//TODO: update use of matrix and indexToBar
};

template<class Master_matrix>
inline Chain_vine_swap<Master_matrix>::Chain_vine_swap(matrix_type &matrix)
	: Chain_pairing<Master_matrix>(), matrix_(&matrix)
{}

template<class Master_matrix>
inline Chain_vine_swap<Master_matrix>::Chain_vine_swap(
		const Chain_vine_swap &matrixToCopy)
	: Chain_pairing<Master_matrix>(matrixToCopy), pivotToPosition_(matrixToCopy.pivotToPosition_), matrix_(matrixToCopy.matrix_)
{}

template<class Master_matrix>
inline Chain_vine_swap<Master_matrix>::Chain_vine_swap(Chain_vine_swap<Master_matrix> &&other) noexcept
	: Chain_pairing<Master_matrix>(std::move(other)),
	  pivotToPosition_(std::move(other.pivotToPosition_)),
	  matrix_(std::exchange(other.matrix_, nullptr))
{}

template<class Master_matrix>
inline index Chain_vine_swap<Master_matrix>::vine_swap_with_z_eq_1_case(index columnIndex1, index columnIndex2)
{
	index pivot1 = pivotToPosition_.at(matrix_->at(columnIndex1).get_pivot());
	index pivot2 = pivotToPosition_.at(matrix_->at(columnIndex2).get_pivot());

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
	index pivot1 = matrix_->at(columnIndex1).get_pivot();
	index pivot2 = matrix_->at(columnIndex2).get_pivot();
	index rowIndex1 = pivotToPosition_.at(pivot1);
	index rowIndex2 = pivotToPosition_.at(pivot2);

	if (_is_negative_in_pair(rowIndex1) && _is_negative_in_pair(rowIndex2)){
		if (!matrix_->at(columnIndex2).is_non_zero(matrix_->at(columnIndex1).get_lowest_simplex_index())){
			_negative_transpose(rowIndex1, rowIndex2);
			matrix_->at(pivot1).swap_independent_rows(pivot2);
			return columnIndex1;
		}
		return _negative_vine_swap(columnIndex1, columnIndex2);
	}

	if (_is_negative_in_pair(rowIndex1)){
		if (!matrix_->at(columnIndex2).is_non_zero(matrix_->at(columnIndex1).get_lowest_simplex_index())){
			_negative_positive_transpose(rowIndex1, rowIndex2);
			matrix_->at(pivot1).swap_independent_rows(pivot2);
			return columnIndex1;
		}
		return _negative_positive_vine_swap(columnIndex1, columnIndex2);
	}

	if (_is_negative_in_pair(rowIndex2)){
		if (!matrix_->at(columnIndex2).is_non_zero(matrix_->at(columnIndex1).get_lowest_simplex_index())){
			_positive_negative_transpose(rowIndex1, rowIndex2);
			matrix_->at(pivot1).swap_independent_rows(pivot2);
			return columnIndex1;
		}
		return _positive_negative_vine_swap(columnIndex1, columnIndex2);
	}

	if (!matrix_->at(columnIndex2).is_non_zero(matrix_->at(columnIndex1).get_lowest_simplex_index())){
		_positive_transpose(rowIndex1, rowIndex2);
		matrix_->at(pivot1).swap_independent_rows(pivot2);
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
		const typename Column_type::Column_type& column, std::set<index>& set)
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
	return _death(simplexIndex) == static_cast<int>(simplexIndex);
}

template<class Master_matrix>
inline void Chain_vine_swap<Master_matrix>::_positive_transpose(index simplexIndex1, index simplexIndex2)
{
	_birth(simplexIndex1) = simplexIndex2;
	_birth(simplexIndex2) = simplexIndex1;
	std::swap(CP::indexToBar_.at(simplexIndex1), CP::indexToBar_.at(simplexIndex2));
}

template<class Master_matrix>
inline void Chain_vine_swap<Master_matrix>::_negative_transpose(index simplexIndex1, index simplexIndex2)
{
	_death(simplexIndex1) = simplexIndex2;
	_death(simplexIndex2) = simplexIndex1;
	std::swap(CP::indexToBar_.at(simplexIndex1), CP::indexToBar_.at(simplexIndex2));
}

template<class Master_matrix>
inline void Chain_vine_swap<Master_matrix>::_positive_negative_transpose(index simplexIndex1, index simplexIndex2)
{
	_birth(simplexIndex1) = simplexIndex2;
	_death(simplexIndex2) = simplexIndex1;
	std::swap(CP::indexToBar_.at(simplexIndex1), CP::indexToBar_.at(simplexIndex2));
}

template<class Master_matrix>
inline void Chain_vine_swap<Master_matrix>::_negative_positive_transpose(index simplexIndex1, index simplexIndex2)
{
	_death(simplexIndex1) = simplexIndex2;
	_birth(simplexIndex2) = simplexIndex1;
	std::swap(CP::indexToBar_.at(simplexIndex1), CP::indexToBar_.at(simplexIndex2));
}

template<class Master_matrix>
inline index Chain_vine_swap<Master_matrix>::_positive_vine_swap(index columnIndex1, index columnIndex2)
{
	index pivot1 = pivotToPosition_.at(matrix_->at(columnIndex1).get_pivot());
	index pivot2 = pivotToPosition_.at(matrix_->at(columnIndex2).get_pivot());

	std::swap(pivotToPosition_.at(matrix_->at(columnIndex1).get_pivot()), pivotToPosition_.at(matrix_->at(columnIndex2).get_pivot()));

	if ((_death(pivot1) < _death(pivot2) && _death(pivot1) != -1)
			|| (_death(pivot1) != -1 && _death(pivot2) == -1)
			|| (_birth(pivot1) < _birth(pivot2) && _death(pivot1) == -1 && _death(pivot2) == -1))
	{
		if (_death(pivot2) != -1)
			matrix_->at(matrix_->at(columnIndex2).get_paired_chain_index()) += matrix_->at(matrix_->at(columnIndex1).get_paired_chain_index());
		matrix_->at(columnIndex2) += matrix_->at(columnIndex1);

		_positive_transpose(pivot1, pivot2);

		return columnIndex1;
	}

	if (_death(pivot1) != -1 && _death(pivot2) != -1)
		matrix_->at(matrix_->at(columnIndex1).get_paired_chain_index()) += matrix_->at(matrix_->at(columnIndex2).get_paired_chain_index());
	matrix_->at(columnIndex1) += matrix_->at(columnIndex2);

	return columnIndex2;
}

template<class Master_matrix>
inline index Chain_vine_swap<Master_matrix>::_positive_negative_vine_swap(index columnIndex1, index columnIndex2)
{
	matrix_->at(columnIndex2) += matrix_->at(columnIndex1);
	_positive_negative_transpose(pivotToPosition_.at(matrix_->at(columnIndex1).get_pivot()), pivotToPosition_.at(matrix_->at(columnIndex2).get_pivot()));
	std::swap(pivotToPosition_.at(matrix_->at(columnIndex1).get_pivot()), pivotToPosition_.at(matrix_->at(columnIndex2).get_pivot()));

	return columnIndex1;
}

template<class Master_matrix>
inline index Chain_vine_swap<Master_matrix>::_negative_positive_vine_swap(index columnIndex1, index columnIndex2)
{
	matrix_->at(columnIndex1) += matrix_->at(columnIndex2);
	std::swap(pivotToPosition_.at(matrix_->at(columnIndex1).get_pivot()), pivotToPosition_.at(matrix_->at(columnIndex2).get_pivot()));

	return columnIndex2;
}

template<class Master_matrix>
inline index Chain_vine_swap<Master_matrix>::_negative_vine_swap(index columnIndex1, index columnIndex2)
{
	index pivot1 = pivotToPosition_.at(matrix_->at(columnIndex1).get_pivot());
	index pivot2 = pivotToPosition_.at(matrix_->at(columnIndex2).get_pivot());

	index pairedIndex1 = matrix_->at(columnIndex1).get_paired_chain_index();
	index pairedIndex2 = matrix_->at(columnIndex2).get_paired_chain_index();

	std::swap(pivotToPosition_.at(matrix_->at(columnIndex1).get_pivot()), pivotToPosition_.at(matrix_->at(columnIndex2).get_pivot()));

	if (_birth(pivot1) < _birth(pivot2))
	{
		matrix_->at(pairedIndex2) += matrix_->at(pairedIndex1);
		matrix_->at(columnIndex2) += matrix_->at(columnIndex1);
		_negative_transpose(pivot1, pivot2);

		return columnIndex1;
	}

	matrix_->at(pairedIndex1) += matrix_->at(pairedIndex2);
	matrix_->at(columnIndex1) += matrix_->at(columnIndex2);

	return columnIndex2;
}

template<class Master_matrix>
inline int& Chain_vine_swap<Master_matrix>::_death(index simplexIndex)
{
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		return CP::indexToBar_.at(simplexIndex)->death;
	} else {
		return CP::barcode_.at(CP::indexToBar_.at(simplexIndex)).death;
	}
}

template<class Master_matrix>
inline int& Chain_vine_swap<Master_matrix>::_birth(index simplexIndex)
{
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		return CP::indexToBar_.at(simplexIndex)->birth;
	} else {
		return CP::barcode_.at(CP::indexToBar_.at(simplexIndex)).birth;
	}
}

template<class Friend_master_matrix>
inline void swap(Chain_vine_swap<Friend_master_matrix>& swap1,
				 Chain_vine_swap<Friend_master_matrix>& swap2)
{
	std::swap(swap1.matrix_, swap2.matrix_);
	std::swap(swap1.pivotToPosition_, swap2.pivotToPosition_);
	std::swap(static_cast<Chain_pairing<Friend_master_matrix> >(swap1),
			  static_cast<Chain_pairing<Friend_master_matrix> >(swap2));
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // Chain_VINE_SWAP_H

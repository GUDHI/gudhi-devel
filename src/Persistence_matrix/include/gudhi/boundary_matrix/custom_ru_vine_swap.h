/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CUSTOM_RU_VINE_SWAP_H
#define CUSTOM_RU_VINE_SWAP_H

#include <utility>
#include <functional>

#include "../utilities/utilities.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class Custom_RU_vine_swap
{
public:
	using Base_matrix = typename Master_matrix::Boundary_matrix_type;
	using dictionnary_type = typename Master_matrix::template dictionnary_type<int>;

	Custom_RU_vine_swap(Base_matrix &matrixR, Base_matrix &matrixU, 
						dictionnary_type &pivotToColumn);
	Custom_RU_vine_swap(const Custom_RU_vine_swap &matrixToCopy);
	Custom_RU_vine_swap(Custom_RU_vine_swap&& other) noexcept;

	bool vine_swap_with_z_eq_1_case(index index);	//returns true if barcode was changed
	bool vine_swap(index index);					//returns true if barcode was changed

	Custom_RU_vine_swap& operator=(Custom_RU_vine_swap other);
	friend void swap(Custom_RU_vine_swap& swap1, Custom_RU_vine_swap& swap2){
		return;
	}

protected:
	Base_matrix *reducedMatrixR_;
	Base_matrix *mirrorMatrixU_;
	dictionnary_type *pivotToColumnIndex_;
	static constexpr bool isActive_ = true;

private:
	bool _is_paired(index index);
	int _get_pair(index index);

	void _swap_at_index(index index);
	void _add_to(index sourceIndex, index targetIndex);
	void _positive_transpose(index index);
	void _negative_transpose(index index);
	void _positive_negative_transpose(index index);
	void _negative_positive_transpose(index index);
	bool _positive_vine_swap(index index);
	bool _negative_vine_swap(index index);
	bool _positive_negative_vine_swap(index index);
	bool _negative_positive_vine_swap(index index);
};

template <class Master_matrix>
inline Custom_RU_vine_swap<Master_matrix>::Custom_RU_vine_swap(
		Base_matrix &matrixR, Base_matrix &matrixU,
		dictionnary_type &pivotToColumn)
	: reducedMatrixR_(&matrixR),
	  mirrorMatrixU_(&matrixU),
	  pivotToColumnIndex_(&pivotToColumn)
{}

template <class Master_matrix>
inline Custom_RU_vine_swap<Master_matrix>::Custom_RU_vine_swap(const Custom_RU_vine_swap &matrixToCopy)
	: reducedMatrixR_(matrixToCopy.reducedMatrixR_),
	  mirrorMatrixU_(matrixToCopy.mirrorMatrixU_),
	  pivotToColumnIndex_(matrixToCopy.pivotToColumnIndex_)
{}

template <class Master_matrix>
inline Custom_RU_vine_swap<Master_matrix>::Custom_RU_vine_swap(Custom_RU_vine_swap<Master_matrix> &&other) noexcept
	: reducedMatrixR_(other.reducedMatrixR_),
	  mirrorMatrixU_(other.mirrorMatrixU_),
	  pivotToColumnIndex_(other.pivotToColumnIndex_)
{}

template<class Master_matrix>
inline bool Custom_RU_vine_swap<Master_matrix>::vine_swap_with_z_eq_1_case(index index)
{
	assert(index < reducedMatrixR_->get_number_of_columns() - 1 && "Index to swap out of bound.");

	bool iIsPositive = reducedMatrixR_->is_zero_column(index);
	bool iiIsPositive = reducedMatrixR_->is_zero_column(index + 1);

	if (iIsPositive && iiIsPositive) {
		mirrorMatrixU_->zero_cell(index + 1, index);
		return _positive_vine_swap(index);
	} else if (!iIsPositive && !iiIsPositive)
		return _negative_vine_swap(index);
	else if (iIsPositive && !iiIsPositive)
		return _positive_negative_vine_swap(index);
	else
		return _negative_positive_vine_swap(index);
}

template<class Master_matrix>
inline bool Custom_RU_vine_swap<Master_matrix>::vine_swap(index index)
{
	assert(index < reducedMatrixR_->get_number_of_columns() - 1 && "Index to swap out of bound.");

	bool iIsPositive = reducedMatrixR_->is_zero_column(index);
	bool iiIsPositive = reducedMatrixR_->is_zero_column(index + 1);

	if (iIsPositive && iiIsPositive) {
		if (reducedMatrixR_->get_column_dimension(index) != reducedMatrixR_->get_column_dimension(index + 1)){
			_swap_at_index(index);
			_positive_transpose(index);
			return true;
		}
		if (!mirrorMatrixU_->is_zero_cell(index + 1, index)){
			mirrorMatrixU_->zero_cell(index + 1, index);
		}
		return _positive_vine_swap(index);
	} else if (!iIsPositive && !iiIsPositive) {
		if (reducedMatrixR_->get_column_dimension(index) != reducedMatrixR_->get_column_dimension(index + 1) || mirrorMatrixU_->is_zero_cell(index + 1, index)){
			_swap_at_index(index);
			_negative_transpose(index);
			return true;
		}
		return _negative_vine_swap(index);
	} else if (iIsPositive && !iiIsPositive) {
		if (reducedMatrixR_->get_column_dimension(index) != reducedMatrixR_->get_column_dimension(index + 1) || mirrorMatrixU_->is_zero_cell(index + 1, index)){
			_swap_at_index(index);
			_positive_negative_transpose(index);
			return true;
		}
		return _positive_negative_vine_swap(index);
	} else {
		if (reducedMatrixR_->get_column_dimension(index) != reducedMatrixR_->get_column_dimension(index + 1) || mirrorMatrixU_->is_zero_cell(index + 1, index)){
			_swap_at_index(index);
			_negative_positive_transpose(index);
			return true;
		}
		return _negative_positive_vine_swap(index);
	}
}

template<class Master_matrix>
inline Custom_RU_vine_swap<Master_matrix> &Custom_RU_vine_swap<Master_matrix>::operator=(
		Custom_RU_vine_swap<Master_matrix> other)
{
	return *this;
}

template<class Master_matrix>
inline bool Custom_RU_vine_swap<Master_matrix>::_is_paired(index index){
	if (!reducedMatrixR_->is_zero_column(index)) return true;

	if constexpr (Master_matrix::Option_list::has_removable_columns){
		if (pivotToColumnIndex_->find(index) == pivotToColumnIndex_->end()) return false;
	} else {
		if (pivotToColumnIndex_->operator[](index) == -1) return false;
	}

	return true;
}

template <class Master_matrix>
inline int Custom_RU_vine_swap<Master_matrix>::_get_pair(index index) {
	if (!reducedMatrixR_->is_zero_column(index)) 
		return reducedMatrixR_->get_column(index).get_pivot();

	if constexpr (Master_matrix::Option_list::has_removable_columns) {
		auto it = pivotToColumnIndex_->find(index);
		if (it == pivotToColumnIndex_->end()) return -1;
		return *it;
	} else {
		return pivotToColumnIndex_->operator[](index);
	}
}

template<class Master_matrix>
inline void Custom_RU_vine_swap<Master_matrix>::_swap_at_index(index index){
	reducedMatrixR_->swap_at_indices(index, index + 1);
	mirrorMatrixU_->swap_at_indices(index + 1, index);
}

template<class Master_matrix>
inline void Custom_RU_vine_swap<Master_matrix>::_add_to(index sourceIndex, index targetIndex)
{
	reducedMatrixR_->add_to(sourceIndex, targetIndex);
	mirrorMatrixU_->add_to(sourceIndex, targetIndex);
}

template<class Master_matrix>
inline void Custom_RU_vine_swap<Master_matrix>::_positive_transpose(index index)
{
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		if (_is_paired(index) && _is_paired(index + 1)){
			std::swap(pivotToColumnIndex_->at(index), pivotToColumnIndex_->at(index + 1));
		} else if (_is_paired(index)){
			pivotToColumnIndex_->emplace(index + 1, pivotToColumnIndex_->at(index));
			pivotToColumnIndex_->erase(index);
		} else if (_is_paired(index + 1)){
			pivotToColumnIndex_->emplace(index, pivotToColumnIndex_->at(index + 1));
			pivotToColumnIndex_->erase(index + 1);
		}
	} else {
		std::swap(pivotToColumnIndex_->operator[](index), pivotToColumnIndex_->operator[](index + 1));
	}
}

template<class Master_matrix>
inline void Custom_RU_vine_swap<Master_matrix>::_negative_transpose(index index)
{
	auto pivot1 = reducedMatrixR_->get_column(index).get_pivot();
	auto pivot2 = reducedMatrixR_->get_column(index + 1).get_pivot();

	std::swap(pivotToColumnIndex_->at(pivot1), pivotToColumnIndex_->at(pivot2));
}

template<class Master_matrix>
inline void Custom_RU_vine_swap<Master_matrix>::_positive_negative_transpose(index index)
{
	pivotToColumnIndex_->at(reducedMatrixR_->get_column(index + 1).get_pivot()) = index;
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		if (_is_paired(index)){
			pivotToColumnIndex_->emplace(index + 1, pivotToColumnIndex_->at(index));
			pivotToColumnIndex_->erase(index);
		}
	} else {
		pivotToColumnIndex_->operator[](index + 1) = pivotToColumnIndex_->operator[](index);
		pivotToColumnIndex_->operator[](index) = -1;
	}
}

template<class Master_matrix>
inline void Custom_RU_vine_swap<Master_matrix>::_negative_positive_transpose(index index)
{
	pivotToColumnIndex_->at(reducedMatrixR_->get_column(index).get_pivot()) = index + 1;
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		if (_is_paired(index + 1)){
			pivotToColumnIndex_->emplace(index, pivotToColumnIndex_->at(index + 1));
			pivotToColumnIndex_->erase(index + 1);
		}
	} else {
		pivotToColumnIndex_->operator[](index) = pivotToColumnIndex_->operator[](index + 1);
		pivotToColumnIndex_->operator[](index + 1) = -1;
	}
}

template<class Master_matrix>
inline bool Custom_RU_vine_swap<Master_matrix>::_positive_vine_swap(index index)
{
	int iDeath = _get_pair(index);
	int iiDeath = _get_pair(index + 1);

	if (iDeath != -1 && iiDeath != -1 && !(reducedMatrixR_->is_zero_cell(iiDeath, index)))
	{
		if (iDeath < iiDeath) {
			_swap_at_index(index);
			_add_to(iDeath, iiDeath);
			_positive_transpose(index);
			return true;
		}

		_swap_at_index(index);
		_add_to(iiDeath, iDeath);
		return false;
	}

	_swap_at_index(index);

	if (iDeath != -1 || iiDeath == -1 || reducedMatrixR_->is_zero_cell(iiDeath, index + 1)) {
		_positive_transpose(index);
		return true;
	}

	return false;
}

template<class Master_matrix>
inline bool Custom_RU_vine_swap<Master_matrix>::_negative_vine_swap(index index)
{
	auto pivot1 = reducedMatrixR_->get_column(index).get_pivot();
	auto pivot2 = reducedMatrixR_->get_column(index + 1).get_pivot();

	_add_to(index, index + 1);
	_swap_at_index(index);

	if (pivot1 < pivot2)
	{
		_negative_transpose(index);
		return true;
	}

	_add_to(index, index + 1);

	return false;
}

template<class Master_matrix>
inline bool Custom_RU_vine_swap<Master_matrix>::_positive_negative_vine_swap(index index)
{
	mirrorMatrixU_->zero_cell(index + 1, index);

	_swap_at_index(index);
	_positive_negative_transpose(index);

	return true;
}

template<class Master_matrix>
inline bool Custom_RU_vine_swap<Master_matrix>::_negative_positive_vine_swap(index index)
{
	_add_to(index, index + 1);	//useless for R?
	_swap_at_index(index);		//if additions not made for R, do not swap R columns
	_add_to(index, index + 1);	//useless for R?

	return false;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // CUSTOM_RU_VINE_SWAP_H

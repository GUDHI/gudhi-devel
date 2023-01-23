/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef RU_VINE_SWAP_H
#define RU_VINE_SWAP_H

#include <utility>

#include "../utilities/utilities.h"
#include "ru_pairing.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class RU_vine_swap : public RU_pairing<Master_matrix>
{
public:
	using Base_matrix = typename Master_matrix::Boundary_matrix_type;
	using dictionnary_type = typename Master_matrix::template dictionnary_type<int>;

	RU_vine_swap(Base_matrix &matrixR, Base_matrix &matrixU, dictionnary_type &pivotToColumn);
	RU_vine_swap(const RU_vine_swap &matrixToCopy);
	RU_vine_swap(RU_vine_swap&& other) noexcept;

	bool vine_swap_with_z_eq_1_case(index index);	//returns true if barcode was changed
	bool vine_swap(index index);					//returns true if barcode was changed

	RU_vine_swap& operator=(RU_vine_swap other);
	friend void swap(RU_vine_swap& swap1, RU_vine_swap& swap2){
//		std::swap(swap1.reducedMatrixR_, swap2.reducedMatrixR_);
//		std::swap(swap1.mirrorMatrixU_, swap2.mirrorMatrixU_);
//		std::swap(swap1.pivotToColumnIndex_, swap2.pivotToColumnIndex_);
		swap(static_cast<RU_pairing<Master_matrix>&>(swap1),
			 static_cast<RU_pairing<Master_matrix>&>(swap2));
	}

protected:
	Base_matrix *reducedMatrixR_;
	Base_matrix *mirrorMatrixU_;
	dictionnary_type *pivotToColumnIndex_;
	static constexpr bool isActive_ = true;

private:
	using RUP = RU_pairing<Master_matrix>;

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

	int& _death(index simplexIndex);
	int& _birth(index simplexIndex);
	int _death(index simplexIndex) const;
	int _birth(index simplexIndex) const;
};

template<class Master_matrix>
inline RU_vine_swap<Master_matrix>::RU_vine_swap(Base_matrix &matrixR, Base_matrix &matrixU, dictionnary_type &pivotToColumn)
	: RU_pairing<Master_matrix>(), reducedMatrixR_(&matrixR), mirrorMatrixU_(&matrixU), pivotToColumnIndex_(&pivotToColumn)
{}

template<class Master_matrix>
inline RU_vine_swap<Master_matrix>::RU_vine_swap(
		const RU_vine_swap &matrixToCopy)
	: RU_pairing<Master_matrix>(matrixToCopy),
	  reducedMatrixR_(matrixToCopy.reducedMatrixR_),
	  mirrorMatrixU_(matrixToCopy.mirrorMatrixU_),
	  pivotToColumnIndex_(matrixToCopy.pivotToColumnIndex_)
{}

template<class Master_matrix>
inline RU_vine_swap<Master_matrix>::RU_vine_swap(RU_vine_swap<Master_matrix> &&other) noexcept
	: RU_pairing<Master_matrix>(std::move(other)),
	  reducedMatrixR_(other.reducedMatrixR_),
	  mirrorMatrixU_(other.mirrorMatrixU_),
	  pivotToColumnIndex_(other.pivotToColumnIndex_)
{}

template<class Master_matrix>
inline bool RU_vine_swap<Master_matrix>::vine_swap_with_z_eq_1_case(index index)
{
	assert(index < reducedMatrixR_->get_number_of_columns() - 1 && "Index to swap out of bound.");

	bool iIsPositive = (_birth(index) == static_cast<int>(index));
	bool iiIsPositive = (_birth(index + 1) == static_cast<int>(index) + 1);

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
inline bool RU_vine_swap<Master_matrix>::vine_swap(index index)
{
	assert(index < reducedMatrixR_->get_number_of_columns() - 1 && "Index to swap out of bound.");

	bool iIsPositive = (_birth(index) == static_cast<int>(index));
	bool iiIsPositive = (_birth(index + 1) == static_cast<int>(index) + 1);

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
inline RU_vine_swap<Master_matrix> &RU_vine_swap<Master_matrix>::operator=(
		RU_vine_swap<Master_matrix> other)
{
//	std::swap(reducedMatrixR_, other.reducedMatrixR_);
//	std::swap(mirrorMatrixU_, other.mirrorMatrixU_);
//	std::swap(pivotToColumnIndex_, other.pivotToColumnIndex_);
	RU_pairing<Master_matrix>::operator=(other);
	return *this;
}

template<class Master_matrix>
inline void RU_vine_swap<Master_matrix>::_swap_at_index(index index){
	reducedMatrixR_->swap_at_indices(index, index + 1);
	mirrorMatrixU_->swap_at_indices(index + 1, index);
}

template<class Master_matrix>
inline void RU_vine_swap<Master_matrix>::_add_to(index sourceIndex, index targetIndex)
{
	reducedMatrixR_->add_to(sourceIndex, targetIndex);
	mirrorMatrixU_->add_to(sourceIndex, targetIndex);
}

template<class Master_matrix>
inline void RU_vine_swap<Master_matrix>::_positive_transpose(index index)
{
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		if (_death(index) != -1 && _death(index + 1) != -1){
			std::swap(pivotToColumnIndex_->at(index), pivotToColumnIndex_->at(index + 1));
		} else if (_death(index) != -1){
			pivotToColumnIndex_->emplace(index + 1, pivotToColumnIndex_->at(index));
			pivotToColumnIndex_->erase(index);
		} else if (_death(index + 1) != -1){
			pivotToColumnIndex_->emplace(index, pivotToColumnIndex_->at(index + 1));
			pivotToColumnIndex_->erase(index + 1);
		}
	} else {
		std::swap(pivotToColumnIndex_->operator[](index), pivotToColumnIndex_->operator[](index + 1));
	}

	_birth(index) = index + 1;
	_birth(index + 1) = index;
	std::swap(RUP::indexToBar_.at(index), RUP::indexToBar_.at(index + 1));
}

template<class Master_matrix>
inline void RU_vine_swap<Master_matrix>::_negative_transpose(index index)
{
	_death(index) = index + 1;
	_death(index + 1) = index;
	std::swap(RUP::indexToBar_.at(index), RUP::indexToBar_.at(index + 1));
	std::swap(pivotToColumnIndex_->at(_birth(index)), pivotToColumnIndex_->at(_birth(index + 1)));
}

template<class Master_matrix>
inline void RU_vine_swap<Master_matrix>::_positive_negative_transpose(index index)
{
	pivotToColumnIndex_->at(_birth(index + 1)) = index;
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		if (_death(index) != -1){
			pivotToColumnIndex_->emplace(index + 1, pivotToColumnIndex_->at(index));
			pivotToColumnIndex_->erase(index);
		}
	} else {
		pivotToColumnIndex_->operator[](index + 1) = pivotToColumnIndex_->operator[](index);
		pivotToColumnIndex_->operator[](index) = -1;
	}

	_birth(index) = index + 1;
	_death(index + 1) = index;
	std::swap(RUP::indexToBar_.at(index), RUP::indexToBar_.at(index + 1));
}

template<class Master_matrix>
inline void RU_vine_swap<Master_matrix>::_negative_positive_transpose(index index)
{
	pivotToColumnIndex_->at(_birth(index)) = index + 1;
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		if (_death(index + 1) != -1){
			pivotToColumnIndex_->emplace(index, pivotToColumnIndex_->at(index + 1));
			pivotToColumnIndex_->erase(index + 1);
		}
	} else {
		pivotToColumnIndex_->operator[](index) = pivotToColumnIndex_->operator[](index + 1);
		pivotToColumnIndex_->operator[](index + 1) = -1;
	}

	_death(index) = index + 1;
	_birth(index + 1) = index;
	std::swap(RUP::indexToBar_.at(index), RUP::indexToBar_.at(index + 1));
}

template<class Master_matrix>
inline bool RU_vine_swap<Master_matrix>::_positive_vine_swap(index index)
{
	int iDeath = _death(index);
	int iiDeath = _death(index + 1);

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
inline bool RU_vine_swap<Master_matrix>::_negative_vine_swap(index index)
{
	_add_to(index, index + 1);
	_swap_at_index(index);

	if (_birth(index) < _birth(index + 1))
	{
		_negative_transpose(index);
		return true;
	}

	_add_to(index, index + 1);

	return false;
}

template<class Master_matrix>
inline bool RU_vine_swap<Master_matrix>::_positive_negative_vine_swap(index index)
{
	mirrorMatrixU_->zero_cell(index + 1, index);

	_swap_at_index(index);
	_positive_negative_transpose(index);

	return true;
}

template<class Master_matrix>
inline bool RU_vine_swap<Master_matrix>::_negative_positive_vine_swap(index index)
{
	_add_to(index, index + 1);	//useless for R?
	_swap_at_index(index);		//if additions not made for R, do not swap R columns
	_add_to(index, index + 1);	//useless for R?

	return false;
}

template<class Master_matrix>
inline int& RU_vine_swap<Master_matrix>::_death(index simplexIndex)
{
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		return RUP::indexToBar_.at(simplexIndex)->death;
	} else {
		return RUP::barcode_.at(RUP::indexToBar_.at(simplexIndex)).death;
	}
}

template<class Master_matrix>
inline int& RU_vine_swap<Master_matrix>::_birth(index simplexIndex)
{
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		return RUP::indexToBar_.at(simplexIndex)->birth;
	} else {
		return RUP::barcode_.at(RUP::indexToBar_.at(simplexIndex)).birth;
	}
}

template<class Master_matrix>
inline int RU_vine_swap<Master_matrix>::_death(index simplexIndex) const
{
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		return RUP::indexToBar_.at(simplexIndex)->death;
	} else {
		return RUP::barcode_.at(RUP::indexToBar_.at(simplexIndex)).death;
	}
}

template<class Master_matrix>
inline int RU_vine_swap<Master_matrix>::_birth(index simplexIndex) const
{
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		return RUP::indexToBar_.at(simplexIndex)->birth;
	} else {
		return RUP::barcode_.at(RUP::indexToBar_.at(simplexIndex)).birth;
	}
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // RU_VINE_SWAP_H

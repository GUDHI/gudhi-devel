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

#include "../utilities.h"
#include "ru_pairing.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class RU_vine_swap : public RU_pairing<Master_matrix>
{
public:
	void vine_swap_with_z_eq_1_case(index index);
	void vine_swap(index index);

	RU_vine_swap& operator=(RU_vine_swap other);
	template<class Friend_master_matrix>
	friend void swap(RU_vine_swap<Friend_master_matrix>& swap1,
					 RU_vine_swap<Friend_master_matrix>& swap2);

protected:
	using Base_matrix = typename Master_matrix::Base_matrix;

	RU_vine_swap(Base_matrix &matrixR, Base_matrix &matrixU);
	RU_vine_swap(const RU_vine_swap &matrixToCopy);
	RU_vine_swap(RU_vine_swap&& other) noexcept;

	static constexpr bool isActive_ = true;

private:
	using RUP = RU_pairing<Master_matrix>;

	Base_matrix &reducedMatrixR_;
	Base_matrix &mirrorMatrixU_;

	void _swap_at_index(index index);
	void _add_to(index sourceIndex, index targetIndex);
	void _positive_transpose(index index);
	void _negative_transpose(index index);
	void _positive_negative_transpose(index index);
	void _negative_positive_transpose(index index);
	void _positive_vine_swap(index index);
	void _negative_vine_swap(index index);
	void _positive_negative_vine_swap(index index);
	void _negative_positive_vine_swap(index index);

	int& _death(index simplexIndex);
	int& _birth(index simplexIndex);
};

template<class Master_matrix>
inline RU_vine_swap<Master_matrix>::RU_vine_swap(Base_matrix &matrixR, Base_matrix &matrixU)
	: reducedMatrixR_(matrixR), mirrorMatrixU_(matrixU), RU_pairing<Master_matrix>()
{}

template<class Master_matrix>
inline RU_vine_swap<Master_matrix>::RU_vine_swap(
		const RU_vine_swap &matrixToCopy)
	: reducedMatrixR_(matrixToCopy.reducedMatrixR_),
	  mirrorMatrixU_(matrixToCopy.mirrorMatrixU_),
	  RU_pairing<Master_matrix>(matrixToCopy)
{}

template<class Master_matrix>
inline RU_vine_swap<Master_matrix>::RU_vine_swap(RU_vine_swap<Master_matrix> &&other) noexcept
	: reducedMatrixR_(std::move(other.reducedMatrixR_)),
	  mirrorMatrixU_(std::move(other.mirrorMatrixU_)),
	  RU_pairing<Master_matrix>(std::move(other))
{}

template<class Master_matrix>
inline void RU_vine_swap<Master_matrix>::vine_swap_with_z_eq_1_case(index index)
{
	if (index >= reducedMatrixR_.get_number_of_columns() - 1) return;

	bool iIsPositive = (_birth(index) == static_cast<int>(index));
	bool iiIsPositive = (_birth(index + 1) == static_cast<int>(index) + 1);

	if (iIsPositive && iiIsPositive) {
		mirrorMatrixU_.zero_cell(index + 1, index);
		_positive_vine_swap(index);
	} else if (!iIsPositive && !iiIsPositive)
		_negative_vine_swap(index);
	else if (iIsPositive && !iiIsPositive)
		_positive_negative_vine_swap(index);
	else
		_negative_positive_vine_swap(index);
}

template<class Master_matrix>
inline void RU_vine_swap<Master_matrix>::vine_swap(index index)
{
	if (index >= reducedMatrixR_.get_number_of_columns() - 1) return;

	bool iIsPositive = (_birth(index) == static_cast<int>(index));
	bool iiIsPositive = (_birth(index + 1) == static_cast<int>(index) + 1);

	if (iIsPositive && iiIsPositive) {
		if (reducedMatrixR_.get_column_dimension(index) != reducedMatrixR_.get_column_dimension(index + 1)){
			_swap_at_index(index);
			_positive_transpose(index);
			return;
		}
		if (!mirrorMatrixU_.is_zero_cell(index + 1, index)){
			mirrorMatrixU_.zero_cell(index + 1, index);
		}
		_positive_vine_swap(index);
	} else if (!iIsPositive && !iiIsPositive) {
		if (reducedMatrixR_.get_column_dimension(index) != reducedMatrixR_.get_column_dimension(index + 1) || mirrorMatrixU_.is_zero_cell(index + 1, index)){
			_swap_at_index(index);
			_negative_transpose(index);
			return;
		}
		_negative_vine_swap(index);
	} else if (iIsPositive && !iiIsPositive) {
		if (reducedMatrixR_.get_column_dimension(index) != reducedMatrixR_.get_column_dimension(index + 1) || mirrorMatrixU_.is_zero_cell(index + 1, index)){
			_swap_at_index(index);
			_positive_negative_transpose(index);
			return;
		}
		_positive_negative_vine_swap(index);
	} else {
		if (reducedMatrixR_.get_column_dimension(index) != reducedMatrixR_.get_column_dimension(index + 1) || mirrorMatrixU_.is_zero_cell(index + 1, index)){
			_swap_at_index(index);
			_negative_positive_transpose(index);
			return;
		}
		_negative_positive_vine_swap(index);
	}
}

template<class Master_matrix>
inline RU_vine_swap<Master_matrix> &RU_vine_swap<Master_matrix>::operator=(
		RU_vine_swap<Master_matrix> other)
{
	std::swap(reducedMatrixR_, other.reducedMatrixR_);
	std::swap(mirrorMatrixU_, other.mirrorMatrixU_);
	RU_pairing<Master_matrix>::operator=(other);
	return *this;
}

template<class Master_matrix>
inline void RU_vine_swap<Master_matrix>::_swap_at_index(index index){
	reducedMatrixR_.swap_at_indices(index, index + 1);
	mirrorMatrixU_.swap_at_indices(index + 1, index);
}

template<class Master_matrix>
inline void RU_vine_swap<Master_matrix>::_add_to(index sourceIndex, index targetIndex)
{
	reducedMatrixR_.add_to(sourceIndex, targetIndex);
	mirrorMatrixU_.add_to(sourceIndex, targetIndex);
}

template<class Master_matrix>
inline void RU_vine_swap<Master_matrix>::_positive_transpose(index index)
{
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
}

template<class Master_matrix>
inline void RU_vine_swap<Master_matrix>::_positive_negative_transpose(index index)
{
	_birth(index) = index + 1;
	_death(index + 1) = index;
	std::swap(RUP::indexToBar_.at(index), RUP::indexToBar_.at(index + 1));
}

template<class Master_matrix>
inline void RU_vine_swap<Master_matrix>::_negative_positive_transpose(index index)
{
	_death(index) = index + 1;
	_birth(index + 1) = index;
	std::swap(RUP::indexToBar_.at(index), RUP::indexToBar_.at(index + 1));
}

template<class Master_matrix>
inline void RU_vine_swap<Master_matrix>::_positive_vine_swap(index index)
{
	int iDeath = _death(index);
	int iiDeath = _death(index + 1);

	if (iDeath != -1 && iiDeath != -1 && !(reducedMatrixR_.is_zero_cell(iiDeath, index)))
	{
		if (iDeath < iiDeath) {
			_swap_at_index(index);
			_add_to(iDeath, iiDeath);
			_positive_transpose(index);
			return;
		}

		_swap_at_index(index);
		_add_to(iiDeath, iDeath);
		return;
	}

	_swap_at_index(index);

	if (iDeath != -1 || iiDeath == -1 || reducedMatrixR_.is_zero_cell(iiDeath, index + 1)) {
		_positive_transpose(index);
	}
}

template<class Master_matrix>
inline void RU_vine_swap<Master_matrix>::_negative_vine_swap(index index)
{
	_add_to(index, index + 1);
	_swap_at_index(index);

	if (_birth(index) < _birth(index + 1))
	{
		_negative_transpose(index);
		return;
	}

	_add_to(index, index + 1);
}

template<class Master_matrix>
inline void RU_vine_swap<Master_matrix>::_positive_negative_vine_swap(index index)
{
	mirrorMatrixU_.zero_cell(index + 1, index);

	_swap_at_index(index);
	_positive_negative_transpose(index);
}

template<class Master_matrix>
inline void RU_vine_swap<Master_matrix>::_negative_positive_vine_swap(index index)
{
	_add_to(index, index + 1);
	_swap_at_index(index);
	_add_to(index, index + 1);
}

template<class Master_matrix>
inline int& RU_vine_swap<Master_matrix>::_death(index simplexIndex)
{
	if constexpr (Master_matrix::Options::has_removable_columns){
		return RUP::indexToBar_.at(simplexIndex)->death;
	} else {
		return RUP::barcode_.at(RUP::indexToBar_.at(simplexIndex)).death;
	}
}

template<class Master_matrix>
inline int& RU_vine_swap<Master_matrix>::_birth(index simplexIndex)
{
	if constexpr (Master_matrix::Options::has_removable_columns){
		return RUP::indexToBar_.at(simplexIndex)->birth;
	} else {
		return RUP::barcode_.at(RUP::indexToBar_.at(simplexIndex)).birth;
	}
}

template<class Friend_master_matrix>
inline void swap(RU_vine_swap<Friend_master_matrix>& swap1,
				 RU_vine_swap<Friend_master_matrix>& swap2)
{
	std::swap(swap1.reducedMatrixR_, swap2.reducedMatrixR_);
	std::swap(swap1.mirrorMatrixU_, swap2.mirrorMatrixU_);
	std::swap(static_cast<RU_pairing<Friend_master_matrix> >(swap1),
			  static_cast<RU_pairing<Friend_master_matrix> >(swap2));
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // RU_VINE_SWAP_H

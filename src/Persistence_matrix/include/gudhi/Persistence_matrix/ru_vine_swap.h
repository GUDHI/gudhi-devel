/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022-23 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_RU_VINE_SWAP_H
#define PM_RU_VINE_SWAP_H

#include <utility>	//std::move
#include <type_traits>	//std::conditional
#include <cassert>
#include <vector>

#include "ru_pairing.h"

namespace Gudhi {
namespace persistence_matrix {

struct Dummy_ru_vine_swap{
	friend void swap([[maybe_unused]] Dummy_ru_vine_swap& d1, [[maybe_unused]] Dummy_ru_vine_swap& d2){}

	Dummy_ru_vine_swap(){}
	Dummy_ru_vine_swap([[maybe_unused]] const Dummy_ru_vine_swap& matrixToCopy){}
	Dummy_ru_vine_swap([[maybe_unused]] Dummy_ru_vine_swap&& other) noexcept{}
};

struct Dummy_ru_vine_pairing{
	friend void swap([[maybe_unused]] Dummy_ru_vine_pairing& d1, [[maybe_unused]] Dummy_ru_vine_pairing& d2){}

	Dummy_ru_vine_pairing(){}
	Dummy_ru_vine_pairing([[maybe_unused]] const Dummy_ru_vine_pairing& matrixToCopy){}
	Dummy_ru_vine_pairing([[maybe_unused]] Dummy_ru_vine_pairing&& other) noexcept{}
};

template<class Master_matrix>
class RU_vine_swap : public std::conditional<
								Master_matrix::Option_list::has_column_pairings,
								RU_pairing<Master_matrix>,
								Dummy_ru_vine_pairing
							>::type
{
public:
	using Boundary_matrix = typename Master_matrix::Boundary_matrix_type;
	using Base_matrix = typename Master_matrix::Base_matrix_type;
	// using dictionnary_type = typename Master_matrix::template dictionnary_type<int>;
	using index = typename Master_matrix::index;
	using pos_index = typename Master_matrix::pos_index;
	using id_index = typename Master_matrix::id_index;

	RU_vine_swap();
	RU_vine_swap(const RU_vine_swap &matrixToCopy);
	RU_vine_swap(RU_vine_swap&& other) noexcept;

	bool vine_swap_with_z_eq_1_case(pos_index index);	//returns true if barcode was changed
	bool vine_swap(pos_index index);					//returns true if barcode was changed

	RU_vine_swap& operator=(RU_vine_swap other);
	friend void swap(RU_vine_swap& swap1, RU_vine_swap& swap2){
		if constexpr (Master_matrix::Option_list::has_column_pairings){
			swap(static_cast<RU_pairing<Master_matrix>&>(swap1), static_cast<RU_pairing<Master_matrix>&>(swap2));
		}
		swap1.positionToRowIdx_.swap(swap2.positionToRowIdx_);
	}

protected:
	//only usefull when simplex id does not corresponds to position, so feels kinda useless most of the time...
	//TODO: as it takes up some non trivial memory, see if this should not be optional
	//or only remember the positions with a difference. but then a map is needed, ie find instead of []. 
	std::vector<id_index> positionToRowIdx_;

private:
	using RUP = typename std::conditional<
								Master_matrix::Option_list::has_column_pairings,
								RU_pairing<Master_matrix>,
								Dummy_ru_vine_pairing
							>::type;
	using ru_matrix = typename Master_matrix::RU_matrix_type;

	bool _is_paired(index index);
	// int _get_pair(index index);

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

	pos_index& _death(pos_index simplexIndex);
	pos_index& _birth(pos_index simplexIndex);
	pos_index _get_death(index simplexIndex);
	pos_index _get_birth(index simplexIndex);

	constexpr ru_matrix* _matrix() { return static_cast<ru_matrix*>(this); }
	constexpr const ru_matrix* _matrix() const { return static_cast<const ru_matrix*>(this); }
};

template<class Master_matrix>
inline RU_vine_swap<Master_matrix>::RU_vine_swap()
	: RUP()
{}

template<class Master_matrix>
inline RU_vine_swap<Master_matrix>::RU_vine_swap(
		const RU_vine_swap &matrixToCopy)
	: RUP(static_cast<const RUP&>(matrixToCopy)), positionToRowIdx_(matrixToCopy.positionToRowIdx_)
{}

template<class Master_matrix>
inline RU_vine_swap<Master_matrix>::RU_vine_swap(RU_vine_swap<Master_matrix> &&other) noexcept
	: RUP(std::move(static_cast<RUP&>(other))), positionToRowIdx_(std::move(other.positionToRowIdx_))
{}

template<class Master_matrix>
inline bool RU_vine_swap<Master_matrix>::vine_swap_with_z_eq_1_case(pos_index index)
{
	assert(index < _matrix()->reducedMatrixR_.get_number_of_columns() - 1 && "Index to swap out of bound.");

	bool iIsPositive = _matrix()->reducedMatrixR_.is_zero_column(index);
	bool iiIsPositive = _matrix()->reducedMatrixR_.is_zero_column(index + 1);

	if (iIsPositive && iiIsPositive) {
		_matrix()->mirrorMatrixU_.zero_cell(index, positionToRowIdx_[index + 1]);
		return _positive_vine_swap(index);
	} else if (!iIsPositive && !iiIsPositive)
		return _negative_vine_swap(index);
	else if (iIsPositive && !iiIsPositive)
		return _positive_negative_vine_swap(index);
	else
		return _negative_positive_vine_swap(index);
}

template<class Master_matrix>
inline bool RU_vine_swap<Master_matrix>::vine_swap(pos_index index)
{
	assert(index < _matrix()->reducedMatrixR_.get_number_of_columns() - 1 && "Index to swap out of bound.");

	bool iIsPositive = _matrix()->reducedMatrixR_.is_zero_column(index);
	bool iiIsPositive = _matrix()->reducedMatrixR_.is_zero_column(index + 1);

	if (iIsPositive && iiIsPositive) {
		if (_matrix()->reducedMatrixR_.get_column_dimension(index) != _matrix()->reducedMatrixR_.get_column_dimension(index + 1)){
			_positive_transpose(index);
			_swap_at_index(index);
			return true;
		}
		if (!_matrix()->mirrorMatrixU_.is_zero_cell(index, positionToRowIdx_[index + 1])){
			_matrix()->mirrorMatrixU_.zero_cell(index, positionToRowIdx_[index + 1]);
		}
		return _positive_vine_swap(index);
	} else if (!iIsPositive && !iiIsPositive) {
		if (_matrix()->reducedMatrixR_.get_column_dimension(index) != _matrix()->reducedMatrixR_.get_column_dimension(index + 1) || _matrix()->mirrorMatrixU_.is_zero_cell(index, positionToRowIdx_[index + 1])){
			_negative_transpose(index);
			_swap_at_index(index);
			return true;
		}
		return _negative_vine_swap(index);
	} else if (iIsPositive && !iiIsPositive) {
		if (_matrix()->reducedMatrixR_.get_column_dimension(index) != _matrix()->reducedMatrixR_.get_column_dimension(index + 1) || _matrix()->mirrorMatrixU_.is_zero_cell(index, positionToRowIdx_[index + 1])){
			_positive_negative_transpose(index);
			_swap_at_index(index);
			return true;
		}
		return _positive_negative_vine_swap(index);
	} else {
		if (_matrix()->reducedMatrixR_.get_column_dimension(index) != _matrix()->reducedMatrixR_.get_column_dimension(index + 1) || _matrix()->mirrorMatrixU_.is_zero_cell(index, positionToRowIdx_[index + 1])){
			_negative_positive_transpose(index);
			_swap_at_index(index);
			return true;
		}
		return _negative_positive_vine_swap(index);
	}
}

template<class Master_matrix>
inline RU_vine_swap<Master_matrix> &RU_vine_swap<Master_matrix>::operator=(
		RU_vine_swap<Master_matrix> other)
{
	RUP::operator=(other);
	positionToRowIdx_.swap(other.positionToRowIdx_);
	return *this;
}

template<class Master_matrix>
inline bool RU_vine_swap<Master_matrix>::_is_paired(index index){
	if constexpr (Master_matrix::Option_list::has_column_pairings){
		return _get_death(index) != -1;
	} else {
		if (!_matrix()->reducedMatrixR_.is_zero_column(index)) return true;

		if constexpr (Master_matrix::Option_list::has_map_column_container){
			if (_matrix()->pivotToColumnIndex_.find(index) == _matrix()->pivotToColumnIndex_.end()) return false;
		} else {
			if (_matrix()->pivotToColumnIndex_.operator[](index) == -1) return false;
		}

		return true;
	}
}

// template <class Master_matrix>
// inline int RU_vine_swap<Master_matrix>::_get_pair(index index) {
// 	if (!_matrix()->reducedMatrixR_.is_zero_column(index)) 
// 		return _matrix()->reducedMatrixR_.get_column(index).get_pivot();

// 	if constexpr (Master_matrix::Option_list::has_map_column_container) {
// 		auto it = _matrix()->pivotToColumnIndex_.find(index);
// 		if (it == _matrix()->pivotToColumnIndex_.end()) return -1;
// 		return *it;
// 	} else {
// 		return _matrix()->pivotToColumnIndex_.operator[](index);
// 	}
// }

template<class Master_matrix>
inline void RU_vine_swap<Master_matrix>::_swap_at_index(index index){
	_matrix()->reducedMatrixR_.swap_columns(index, index + 1);
	_matrix()->reducedMatrixR_.swap_rows(positionToRowIdx_[index], positionToRowIdx_[index + 1]);
	_matrix()->mirrorMatrixU_.swap_columns(index, index + 1);
	_matrix()->mirrorMatrixU_.swap_rows(positionToRowIdx_[index], positionToRowIdx_[index + 1]);
}

template<class Master_matrix>
inline void RU_vine_swap<Master_matrix>::_add_to(index sourceIndex, index targetIndex)
{
	_matrix()->reducedMatrixR_.add_to(sourceIndex, targetIndex);
	_matrix()->mirrorMatrixU_.add_to(targetIndex, sourceIndex);
}

template<class Master_matrix>
inline void RU_vine_swap<Master_matrix>::_positive_transpose(index index)
{
	if constexpr (Master_matrix::Option_list::has_map_column_container){
		if (_is_paired(index) && _is_paired(index + 1)){
			std::swap(_matrix()->pivotToColumnIndex_.at(index), _matrix()->pivotToColumnIndex_.at(index + 1));
		} else if (_is_paired(index)){
			_matrix()->pivotToColumnIndex_.emplace(index + 1, _matrix()->pivotToColumnIndex_.at(index));
			_matrix()->pivotToColumnIndex_.erase(index);
		} else if (_is_paired(index + 1)){
			_matrix()->pivotToColumnIndex_.emplace(index, _matrix()->pivotToColumnIndex_.at(index + 1));
			_matrix()->pivotToColumnIndex_.erase(index + 1);
		}
	} else {
		std::swap(_matrix()->pivotToColumnIndex_.operator[](index), _matrix()->pivotToColumnIndex_.operator[](index + 1));
	}

	if constexpr (Master_matrix::Option_list::has_column_pairings){
		_birth(index) = index + 1;
		_birth(index + 1) = index;
		std::swap(RUP::indexToBar_.at(index), RUP::indexToBar_.at(index + 1));
	}
}

template<class Master_matrix>
inline void RU_vine_swap<Master_matrix>::_negative_transpose(index index)
{
	if constexpr (Master_matrix::Option_list::has_column_pairings){
		_death(index) = index + 1;
		_death(index + 1) = index;
		std::swap(RUP::indexToBar_.at(index), RUP::indexToBar_.at(index + 1));
	}
	std::swap(_matrix()->pivotToColumnIndex_.at(_get_birth(index)), _matrix()->pivotToColumnIndex_.at(_get_birth(index + 1)));
}

template<class Master_matrix>
inline void RU_vine_swap<Master_matrix>::_positive_negative_transpose(index index)
{
	_matrix()->pivotToColumnIndex_.at(_get_birth(index + 1)) = index;
	if constexpr (Master_matrix::Option_list::has_map_column_container){
		if (_is_paired(index)){
			_matrix()->pivotToColumnIndex_.emplace(index + 1, _matrix()->pivotToColumnIndex_.at(index));
			_matrix()->pivotToColumnIndex_.erase(index);
		}
	} else {
		_matrix()->pivotToColumnIndex_.operator[](index + 1) = _matrix()->pivotToColumnIndex_.operator[](index);
		_matrix()->pivotToColumnIndex_.operator[](index) = -1;
	}

	if constexpr (Master_matrix::Option_list::has_column_pairings){
		_birth(index) = index + 1;
		_death(index + 1) = index;
		std::swap(RUP::indexToBar_.at(index), RUP::indexToBar_.at(index + 1));
	}
}

template<class Master_matrix>
inline void RU_vine_swap<Master_matrix>::_negative_positive_transpose(index index)
{
	_matrix()->pivotToColumnIndex_.at(_get_birth(index)) = index + 1;
	if constexpr (Master_matrix::Option_list::has_map_column_container){
		if (_is_paired(index + 1)){
			_matrix()->pivotToColumnIndex_.emplace(index, _matrix()->pivotToColumnIndex_.at(index + 1));
			_matrix()->pivotToColumnIndex_.erase(index + 1);
		}
	} else {
		_matrix()->pivotToColumnIndex_.operator[](index) = _matrix()->pivotToColumnIndex_.operator[](index + 1);
		_matrix()->pivotToColumnIndex_.operator[](index + 1) = -1;
	}

	if constexpr (Master_matrix::Option_list::has_column_pairings){
		_death(index) = index + 1;
		_birth(index + 1) = index;
		std::swap(RUP::indexToBar_.at(index), RUP::indexToBar_.at(index + 1));
	}
}

template<class Master_matrix>
inline bool RU_vine_swap<Master_matrix>::_positive_vine_swap(index index)
{
	const pos_index iDeath = _get_death(index);
	const pos_index iiDeath = _get_death(index + 1);

	if (iDeath != static_cast<pos_index>(-1) && iiDeath != static_cast<pos_index>(-1) && !(_matrix()->reducedMatrixR_.is_zero_cell(iiDeath, positionToRowIdx_[index])))
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

	if (iDeath != static_cast<pos_index>(-1) || iiDeath == static_cast<pos_index>(-1) || _matrix()->reducedMatrixR_.is_zero_cell(iiDeath, positionToRowIdx_[index + 1])) {
		_positive_transpose(index);
		return true;
	}
	return false;
}

template<class Master_matrix>
inline bool RU_vine_swap<Master_matrix>::_negative_vine_swap(index index)
{
	const pos_index iBirth = _get_birth(index);
	const pos_index iiBirth = _get_birth(index + 1);

	_add_to(index, index + 1);
	_swap_at_index(index);

	if (iBirth < iiBirth)
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
	_matrix()->mirrorMatrixU_.zero_cell(index, positionToRowIdx_[index + 1]);

	_swap_at_index(index);
	_positive_negative_transpose(index);

	return true;
}

template<class Master_matrix>
inline bool RU_vine_swap<Master_matrix>::_negative_positive_vine_swap(index index)
{
	_add_to(index, index + 1);	//useless for R?
	_swap_at_index(index);		//if additions not made for R, do not swap R columns, just rows
	_add_to(index, index + 1);	//useless for R?

	return false;
}

template<class Master_matrix>
inline typename RU_vine_swap<Master_matrix>::pos_index& RU_vine_swap<Master_matrix>::_death(pos_index simplexIndex)
{
	static_assert(Master_matrix::Option_list::has_column_pairings, "Pairing necessary to modify death value.");

	if constexpr (Master_matrix::Option_list::has_removable_columns){
		return RUP::indexToBar_.at(simplexIndex)->death;
	} else {
		return RUP::barcode_.at(RUP::indexToBar_.at(simplexIndex)).death;
	}
}

template<class Master_matrix>
inline typename RU_vine_swap<Master_matrix>::pos_index& RU_vine_swap<Master_matrix>::_birth(pos_index simplexIndex)
{
	static_assert(Master_matrix::Option_list::has_column_pairings, "Pairing necessary to modify birth value.");

	if constexpr (Master_matrix::Option_list::has_removable_columns){
		return RUP::indexToBar_.at(simplexIndex)->birth;
	} else {
		return RUP::barcode_.at(RUP::indexToBar_.at(simplexIndex)).birth;
	}
}

template<class Master_matrix>
inline typename RU_vine_swap<Master_matrix>::pos_index RU_vine_swap<Master_matrix>::_get_death(index simplexIndex)
{
	if constexpr (Master_matrix::Option_list::has_column_pairings){
		if constexpr (Master_matrix::Option_list::has_removable_columns){
			return RUP::indexToBar_.at(simplexIndex)->death;
		} else {
			return RUP::barcode_.at(RUP::indexToBar_.at(simplexIndex)).death;
		}
	} else {
		if (!_matrix()->reducedMatrixR_.is_zero_column(simplexIndex)) 
			return _matrix()->reducedMatrixR_.get_column(simplexIndex).get_pivot();

		if constexpr (Master_matrix::Option_list::has_map_column_container) {
			auto it = _matrix()->pivotToColumnIndex_.find(simplexIndex);
			if (it == _matrix()->pivotToColumnIndex_.end()) return -1;
			return it->second;
		} else {
			return _matrix()->pivotToColumnIndex_.operator[](simplexIndex);
		}
	}
}

template<class Master_matrix>
inline typename RU_vine_swap<Master_matrix>::pos_index RU_vine_swap<Master_matrix>::_get_birth(index negativeSimplexIndex)
{
	if constexpr (Master_matrix::Option_list::has_column_pairings){
		if constexpr (Master_matrix::Option_list::has_removable_columns){
			return RUP::indexToBar_.at(negativeSimplexIndex)->birth;
		} else {
			return RUP::barcode_.at(RUP::indexToBar_.at(negativeSimplexIndex)).birth;
		}
	} else {
		return _matrix()->reducedMatrixR_.get_pivot(negativeSimplexIndex);
	}
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // PM_RU_VINE_SWAP_H

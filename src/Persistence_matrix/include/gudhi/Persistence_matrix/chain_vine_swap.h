/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022-23 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_CHAIN_VINE_SWAP_H
#define PM_CHAIN_VINE_SWAP_H

#include <utility>	//std::swap & std::move
#include <cmath>	//std::abs
#include <cassert>
#include <set>

#include <iostream> 	//for debug

#include "chain_pairing.h"

namespace Gudhi {
namespace persistence_matrix {

static constexpr bool _no_G_death_comparator(unsigned int columnIndex1, unsigned int columnIndex2){
	return false;
}

struct Dummy_chain_vine_swap{
	friend void swap([[maybe_unused]] Dummy_chain_vine_swap& d1, [[maybe_unused]] Dummy_chain_vine_swap& d2){}

	Dummy_chain_vine_swap(){}
	template<typename BirthComparatorFunction, typename DeathComparatorFunction>
	Dummy_chain_vine_swap([[maybe_unused]] BirthComparatorFunction&& birthComparator, [[maybe_unused]] DeathComparatorFunction&& deathComparator){}
	Dummy_chain_vine_swap([[maybe_unused]] const Dummy_chain_vine_swap &matrixToCopy){}
	Dummy_chain_vine_swap([[maybe_unused]] Dummy_chain_vine_swap&& other) noexcept{}
};

struct Dummy_chain_vine_pairing{
	friend void swap([[maybe_unused]] Dummy_chain_vine_pairing& d1, [[maybe_unused]] Dummy_chain_vine_pairing& d2){}

	Dummy_chain_vine_pairing(){}
};

template<typename Master_matrix>
class Chain_barcode_swap : public Chain_pairing<Master_matrix>
{
public:
	using index = typename Master_matrix::index;
	using CP = Chain_pairing<Master_matrix>;

	Chain_barcode_swap(){};
	Chain_barcode_swap(const Chain_barcode_swap &toCopy) 
		: CP(static_cast<const CP&>(toCopy)), 
		  pivotToPosition_(toCopy.pivotToPosition_) {};
	Chain_barcode_swap(Chain_barcode_swap &&other) 
		: CP(std::move(static_cast<CP&>(other))), 
		  pivotToPosition_(std::move(other.pivotToPosition_)) {};

protected:
	using dictionnary_type = typename Master_matrix::template dictionnary_type<index>;

	dictionnary_type pivotToPosition_;	//necessary to keep track of the barcode changes
	
	void swap_positions(index pivot1, index pivot2){
		if constexpr (Master_matrix::Option_list::has_removable_columns){
			std::swap(pivotToPosition_.at(pivot1), pivotToPosition_.at(pivot2));
		} else {
			std::swap(pivotToPosition_[pivot1], pivotToPosition_[pivot2]);
		}
	}

	bool is_negative_in_pair(index pivot) const
	{
		index pos = _get_pivot_position(pivot);
		return death(pivot) == static_cast<int>(pos);
	}

	void positive_transpose(index pivot1, index pivot2)
	{
		index simplexIndex1 = _get_pivot_position(pivot1);
		index simplexIndex2 = _get_pivot_position(pivot2);

		_birth(simplexIndex1) = simplexIndex2;
		_birth(simplexIndex2) = simplexIndex1;
		std::swap(CP::indexToBar_.at(simplexIndex1), CP::indexToBar_.at(simplexIndex2));
	}

	void negative_transpose(index pivot1, index pivot2)
	{
		index simplexIndex1 = _get_pivot_position(pivot1);
		index simplexIndex2 = _get_pivot_position(pivot2);

		_death(simplexIndex1) = simplexIndex2;
		_death(simplexIndex2) = simplexIndex1;
		std::swap(CP::indexToBar_.at(simplexIndex1), CP::indexToBar_.at(simplexIndex2));
	}

	void positive_negative_transpose(index pivot1, index pivot2)
	{
		index simplexIndex1 = _get_pivot_position(pivot1);
		index simplexIndex2 = _get_pivot_position(pivot2);

		_birth(simplexIndex1) = simplexIndex2;
		_death(simplexIndex2) = simplexIndex1;
		std::swap(CP::indexToBar_.at(simplexIndex1), CP::indexToBar_.at(simplexIndex2));
	}

	void negative_positive_transpose(index pivot1, index pivot2)
	{
		index simplexIndex1 = _get_pivot_position(pivot1);
		index simplexIndex2 = _get_pivot_position(pivot2);

		_death(simplexIndex1) = simplexIndex2;
		_birth(simplexIndex2) = simplexIndex1;
		std::swap(CP::indexToBar_.at(simplexIndex1), CP::indexToBar_.at(simplexIndex2));
	}

	bool are_adjacent(index pivot1, index pivot2) const{
		return std::abs(static_cast<int>(CP::_get_pivot_position(pivot2)) 
						- static_cast<int>(CP::_get_pivot_position(pivot1))) == 1;
	}

	Chain_barcode_swap& operator=(Chain_barcode_swap other){
		Chain_pairing<Master_matrix>::operator=(other);
		pivotToPosition_.swap(other.pivotToPosition_);
	}
	friend void swap(Chain_barcode_swap& swap1, Chain_barcode_swap& swap2){
		swap(static_cast<Chain_pairing<Master_matrix>&>(swap1), 
			 static_cast<Chain_pairing<Master_matrix>&>(swap2));
		swap1.pivotToPosition_.swap(swap2.pivotToPosition_);
	}

	int death(index pivot) const
	{
		index simplexIndex = _get_pivot_position(pivot);

		if constexpr (Master_matrix::Option_list::has_removable_columns){
			return CP::indexToBar_.at(simplexIndex)->death;
		} else {
			return CP::barcode_.at(CP::indexToBar_.at(simplexIndex)).death;
		}
	}

	int birth(index pivot) const
	{
		index simplexIndex = _get_pivot_position(pivot);

		if constexpr (Master_matrix::Option_list::has_removable_columns){
			return CP::indexToBar_.at(simplexIndex)->birth;
		} else {
			return CP::barcode_.at(CP::indexToBar_.at(simplexIndex)).birth;
		}
	}

private:
	index _get_pivot_position(index pivot) const{
		if constexpr (Master_matrix::Option_list::has_removable_columns){
			return pivotToPosition_.at(pivot);		//quite often called, make public and pass position instead of pivot to avoid find() everytime?
		} else {
			return pivotToPosition_[pivot];
		}
	}

	int& _death(index simplexIndex)
	{
		if constexpr (Master_matrix::Option_list::has_removable_columns){
			return CP::indexToBar_.at(simplexIndex)->death;
		} else {
			return CP::barcode_.at(CP::indexToBar_.at(simplexIndex)).death;
		}
	}

	int& _birth(index simplexIndex)
	{
		if constexpr (Master_matrix::Option_list::has_removable_columns){
			return CP::indexToBar_.at(simplexIndex)->birth;
		} else {
			return CP::barcode_.at(CP::indexToBar_.at(simplexIndex)).birth;
		}
	}
};

template<class Master_matrix>
class Chain_vine_swap : public std::conditional<
									Master_matrix::Option_list::has_column_pairings,
									Chain_barcode_swap<Master_matrix>,
									Dummy_chain_vine_pairing
								>::type
{
public:
	using index = typename Master_matrix::index;
	using matrix_type = typename Master_matrix::column_container_type;
	using dictionnary_type = typename Master_matrix::template dictionnary_type<index>;
	using Column_type = typename Master_matrix::Column_type;
	typedef bool (*BirthCompFuncPointer)(index,index);
	typedef bool (*DeathCompFuncPointer)(index,index);

	Chain_vine_swap();
	template<typename BirthComparatorFunction, typename DeathComparatorFunction>
	Chain_vine_swap(BirthComparatorFunction&& birthComparator, 
					DeathComparatorFunction&& deathComparator = _no_G_death_comparator);
	Chain_vine_swap(const Chain_vine_swap &matrixToCopy);
	Chain_vine_swap(Chain_vine_swap&& other) noexcept;

	index vine_swap_with_z_eq_1_case(index columnIndex1, index columnIndex2);	//returns index which was not modified, ie new i+1
	index vine_swap(index columnIndex1, index columnIndex2);					//returns index which was not modified, ie new i+1

	Chain_vine_swap& operator=(Chain_vine_swap other);
	friend void swap(Chain_vine_swap& swap1, Chain_vine_swap& swap2){
		if constexpr (Master_matrix::Option_list::has_column_pairings){
			swap(static_cast<Chain_barcode_swap<Master_matrix>&>(swap1), 
				 static_cast<Chain_barcode_swap<Master_matrix>&>(swap2));
		}
		std::swap(swap1.birthComp_, swap2.birthComp_);
		std::swap(swap1.deathComp_, swap2.deathComp_);
	}

private:
	using CP = typename std::conditional<
							Master_matrix::Option_list::has_column_pairings,
							Chain_barcode_swap<Master_matrix>,
							Dummy_chain_vine_pairing
						>::type;
	using chain_matrix = typename Master_matrix::Chain_matrix_type;

	BirthCompFuncPointer birthComp_;	// for F x F & H x H
	DeathCompFuncPointer deathComp_;	// for G x G

	bool _is_negative_in_pair(index columnIndex);

	index _positive_vine_swap(index columnIndex1, index columnIndex2);
	index _positive_negative_vine_swap(index columnIndex1, index columnIndex2);
	index _negative_positive_vine_swap(index columnIndex1, index columnIndex2);
	index _negative_vine_swap(index columnIndex1, index columnIndex2);

	constexpr chain_matrix* _matrix() { return static_cast<chain_matrix*>(this); }
	constexpr const chain_matrix* _matrix() const { return static_cast<const chain_matrix*>(this); }
};

template<class Master_matrix>
inline Chain_vine_swap<Master_matrix>::Chain_vine_swap()
	: CP(), birthComp_(nullptr), deathComp_(nullptr)
{
	static_assert(Master_matrix::Option_list::has_column_pairings, 
				  "If barcode is not stored, at least a birth comparator has to be specified.");
}

template<class Master_matrix>
template<typename BirthComparatorFunction, typename DeathComparatorFunction>
inline Chain_vine_swap<Master_matrix>::Chain_vine_swap(
	BirthComparatorFunction&& birthComparator, 
	DeathComparatorFunction&& deathComparator)
	: CP(), birthComp_(&birthComparator), deathComp_(&deathComparator)
{}

template<class Master_matrix>
inline Chain_vine_swap<Master_matrix>::Chain_vine_swap(
		const Chain_vine_swap &matrixToCopy)
	: CP(static_cast<const CP&>(matrixToCopy)),
	  birthComp_(matrixToCopy.birthComp_),
	  deathComp_(matrixToCopy.deathComp_)
{}

template<class Master_matrix>
inline Chain_vine_swap<Master_matrix>::Chain_vine_swap(Chain_vine_swap<Master_matrix> &&other) noexcept
	: CP(std::move(static_cast<CP&>(other))),
	  birthComp_(std::move(other.birthComp_)),
	  deathComp_(std::move(other.deathComp_))
{}

template<class Master_matrix>
inline typename Chain_vine_swap<Master_matrix>::index Chain_vine_swap<Master_matrix>::vine_swap_with_z_eq_1_case(index columnIndex1, index columnIndex2)
{
	if constexpr (Master_matrix::Option_list::has_column_pairings){
		assert(CP::are_adjacent(_matrix()->get_pivot(columnIndex1), _matrix()->get_pivot(columnIndex2))
			&& "Columns to be swaped need to be adjacent in the 'real' matrix.");
	}

	const bool col1IsNeg = _is_negative_in_pair(columnIndex1);
	const bool col2IsNeg = _is_negative_in_pair(columnIndex2);

	if (col1IsNeg && col2IsNeg)
		return _negative_vine_swap(columnIndex1, columnIndex2);

	if (col1IsNeg)
		return _negative_positive_vine_swap(columnIndex1, columnIndex2);

	if (col2IsNeg)
		return _positive_negative_vine_swap(columnIndex1, columnIndex2);

	return _positive_vine_swap(columnIndex1, columnIndex2);
}

template<class Master_matrix>
inline typename Chain_vine_swap<Master_matrix>::index Chain_vine_swap<Master_matrix>::vine_swap(index columnIndex1, index columnIndex2)
{
	if constexpr (Master_matrix::Option_list::has_column_pairings){
		assert(CP::are_adjacent(_matrix()->get_pivot(columnIndex1), _matrix()->get_pivot(columnIndex2)) && "Columns to be swaped need to be adjacent in the 'real' matrix.");
	}

	const bool col1IsNeg = _is_negative_in_pair(columnIndex1);
	const bool col2IsNeg = _is_negative_in_pair(columnIndex2);

	if (col1IsNeg && col2IsNeg){
		if (_matrix()->is_zero_cell(columnIndex2, _matrix()->get_pivot(columnIndex1))){
			if constexpr (Master_matrix::Option_list::has_column_pairings){
				index pivot1 = _matrix()->get_pivot(columnIndex1);
				index pivot2 = _matrix()->get_pivot(columnIndex2);

				CP::negative_transpose(pivot1, pivot2);
				CP::swap_positions(pivot1, pivot2);
			}
			return columnIndex1;
		}
		return _negative_vine_swap(columnIndex1, columnIndex2);
	}

	if (col1IsNeg){
		if (_matrix()->is_zero_cell(columnIndex2, _matrix()->get_pivot(columnIndex1))){
			if constexpr (Master_matrix::Option_list::has_column_pairings){
				index pivot1 = _matrix()->get_pivot(columnIndex1);
				index pivot2 = _matrix()->get_pivot(columnIndex2);

				CP::negative_positive_transpose(pivot1, pivot2);
				CP::swap_positions(pivot1, pivot2);
			}
			return columnIndex1;
		}
		return _negative_positive_vine_swap(columnIndex1, columnIndex2);
	}

	if (col2IsNeg){
		if (_matrix()->is_zero_cell(columnIndex2, _matrix()->get_pivot(columnIndex1))){
			if constexpr (Master_matrix::Option_list::has_column_pairings){
				index pivot1 = _matrix()->get_pivot(columnIndex1);
				index pivot2 = _matrix()->get_pivot(columnIndex2);

				CP::positive_negative_transpose(pivot1, pivot2);
				CP::swap_positions(pivot1, pivot2);
			}
			return columnIndex1;
		}
		return _positive_negative_vine_swap(columnIndex1, columnIndex2);
	}

	if (_matrix()->is_zero_cell(columnIndex2, _matrix()->get_pivot(columnIndex1))){
		if constexpr (Master_matrix::Option_list::has_column_pairings){
			index pivot1 = _matrix()->get_pivot(columnIndex1);
			index pivot2 = _matrix()->get_pivot(columnIndex2);

			CP::positive_transpose(pivot1, pivot2);
			CP::swap_positions(pivot1, pivot2);
		}
		return columnIndex1;
	}
	return _positive_vine_swap(columnIndex1, columnIndex2);
}

template<class Master_matrix>
inline Chain_vine_swap<Master_matrix> &Chain_vine_swap<Master_matrix>::operator=(
		Chain_vine_swap<Master_matrix> other)
{
	CP::operator=(other);
	std::swap(birthComp_, other.birthComp_);
	std::swap(deathComp_, other.deathComp_);
	return *this;
}

template<class Master_matrix>
inline bool Chain_vine_swap<Master_matrix>::_is_negative_in_pair(index columnIndex)
{
	if constexpr (Master_matrix::Option_list::has_column_pairings){
		return CP::is_negative_in_pair(_matrix()->get_pivot(columnIndex));
	} else {
		auto& col = _matrix()->get_column(columnIndex);
		if (!col.is_paired()) return false;
		return col.get_pivot() > _matrix()->get_pivot(col.get_paired_chain_index());
	}
}

template<class Master_matrix>
inline typename Chain_vine_swap<Master_matrix>::index Chain_vine_swap<Master_matrix>::_positive_vine_swap(index columnIndex1, index columnIndex2)
{
	auto& col1 = _matrix()->get_column(columnIndex1);
	auto& col2 = _matrix()->get_column(columnIndex2);

	if constexpr (Master_matrix::Option_list::has_column_pairings){
		CP::swap_positions(col1.get_pivot(), col2.get_pivot());
	}
//TODO: factorize the cases. But for debug it is much more easier to understand what is happening splitted like this
	if (!col1.is_paired()){		// F x *
		bool hasSmallerBirth;
		if constexpr (Master_matrix::Option_list::has_column_pairings){
			hasSmallerBirth = (CP::birth(col1.get_pivot()) < CP::birth(col2.get_pivot()));
		} else {
			hasSmallerBirth = birthComp_(columnIndex1, columnIndex2);
		}
		if (!col2.is_paired() && hasSmallerBirth){
			_matrix()->add_to(columnIndex1, columnIndex2);
			if constexpr (Master_matrix::Option_list::has_column_pairings){
				CP::positive_transpose(col1.get_pivot(), col2.get_pivot());
			}
			return columnIndex1;
		}
		_matrix()->add_to(columnIndex2, columnIndex1);

		return columnIndex2;
	}

	if (!col2.is_paired()){		// G x F
		static_cast<chain_matrix*>(this)->add_to(columnIndex1, columnIndex2);
		if constexpr (Master_matrix::Option_list::has_column_pairings){
			CP::positive_transpose(col1.get_pivot(), col2.get_pivot());
		}
		return columnIndex1;
	}

	bool hasSmallerDeath;
	if constexpr (Master_matrix::Option_list::has_column_pairings){
		hasSmallerDeath = (CP::death(col1.get_pivot()) < CP::death(col2.get_pivot()));
	} else {
		hasSmallerDeath = deathComp_(columnIndex1, columnIndex2);
	}

	// G x G
	if (hasSmallerDeath) // == if (matrix_->at(pairedIndex1).get_pivot() < matrix_->at(pairedIndex2).get_pivot()) ???
	{
		_matrix()->add_to(col1.get_paired_chain_index(), col2.get_paired_chain_index());
		_matrix()->add_to(columnIndex1, columnIndex2);
		if constexpr (Master_matrix::Option_list::has_column_pairings){
			CP::positive_transpose(col1.get_pivot(), col2.get_pivot());
		}
		return columnIndex1;
	}

	_matrix()->add_to(col2.get_paired_chain_index(), col1.get_paired_chain_index());
	_matrix()->add_to(columnIndex2, columnIndex1);

	return columnIndex2;
}

template<class Master_matrix>
inline typename Chain_vine_swap<Master_matrix>::index Chain_vine_swap<Master_matrix>::_positive_negative_vine_swap(index columnIndex1, index columnIndex2)
{
	_matrix()->add_to(columnIndex1, columnIndex2);

	if constexpr (Master_matrix::Option_list::has_column_pairings){
		index pivot1 = _matrix()->get_pivot(columnIndex1);
		index pivot2 = _matrix()->get_pivot(columnIndex2);

		CP::positive_negative_transpose(pivot1, pivot2);
		CP::swap_positions(pivot1, pivot2);
	}

	return columnIndex1;
}

template<class Master_matrix>
inline typename Chain_vine_swap<Master_matrix>::index Chain_vine_swap<Master_matrix>::_negative_positive_vine_swap(index columnIndex1, index columnIndex2)
{
	_matrix()->add_to(columnIndex2, columnIndex1);

	if constexpr (Master_matrix::Option_list::has_column_pairings){
		CP::swap_positions(_matrix()->get_pivot(columnIndex1), _matrix()->get_pivot(columnIndex2));
	}

	return columnIndex2;
}

template<class Master_matrix>
inline typename Chain_vine_swap<Master_matrix>::index Chain_vine_swap<Master_matrix>::_negative_vine_swap(index columnIndex1, index columnIndex2)
{
	auto& col1 = _matrix()->get_column(columnIndex1);
	auto& col2 = _matrix()->get_column(columnIndex2);

	index pairedIndex1 = col1.get_paired_chain_index();
	index pairedIndex2 = col2.get_paired_chain_index();

	bool hasSmallerBirth;
	if constexpr (Master_matrix::Option_list::has_column_pairings){
		hasSmallerBirth = (CP::birth(col1.get_pivot()) < CP::birth(col2.get_pivot()));
	} else {
		hasSmallerBirth = birthComp_(columnIndex1, columnIndex2);

		//for debug, to remove
		if (hasSmallerBirth != (_matrix()->get_pivot(pairedIndex1) < _matrix()->get_pivot(pairedIndex2)))
			std::cout << "!!!!!!!!!!!!!!!!!! not equal\n";
	}

	if constexpr (Master_matrix::Option_list::has_column_pairings){
		CP::swap_positions(col1.get_pivot(), col2.get_pivot());
	}

	if (hasSmallerBirth)	//== matrix_->at(pairedIndex1).get_pivot() < matrix_->at(pairedIndex2).get_pivot() ?
	{
		_matrix()->add_to(pairedIndex1, pairedIndex2);
		_matrix()->add_to(columnIndex1, columnIndex2);

		if constexpr (Master_matrix::Option_list::has_column_pairings){
			CP::negative_transpose(col1.get_pivot(), col2.get_pivot());
		}

		return columnIndex1;
	}

	_matrix()->add_to(pairedIndex2, pairedIndex1);
	_matrix()->add_to(columnIndex2, columnIndex1);

	return columnIndex2;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // PM_CHAIN_VINE_SWAP_H

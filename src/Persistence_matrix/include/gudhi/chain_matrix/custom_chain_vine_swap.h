/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CUSTOM_CHAIN_VINE_SWAP_H
#define CUSTOM_CHAIN_VINE_SWAP_H

#include <functional>
#include <utility>
#include <set>

#include "../utilities/utilities.h"

namespace Gudhi {
namespace persistence_matrix {

static constexpr bool _no_G_death_comparator(index columnIndex1, index columnIndex2){
	return false;
}

template<class Master_matrix>
class Custom_chain_vine_swap
{
public:
	using matrix_type = typename Master_matrix::column_container_type;
	using dictionnary_type = typename Master_matrix::template dictionnary_type<index>;

	Custom_chain_vine_swap(matrix_type& matrix, 
						   std::function<bool(index,index)> birthComparator, 
						   std::function<bool(index,index)> deathComparator = _no_G_death_comparator);
	Custom_chain_vine_swap(const Custom_chain_vine_swap &matrixToCopy);
	Custom_chain_vine_swap(Custom_chain_vine_swap&& other) noexcept;

	index vine_swap_with_z_eq_1_case(index columnIndex1, index columnIndex2);	//returns index which was not modified, ie new i+1
	index vine_swap(index columnIndex1, index columnIndex2);					//returns index which was not modified, ie new i+1

	Custom_chain_vine_swap& operator=(Custom_chain_vine_swap other);
	friend void swap(Custom_chain_vine_swap& swap1, Custom_chain_vine_swap& swap2){
		std::swap(swap1.birthComp_, swap2.birthComp_);
	}

protected:
	matrix_type* matrix_;
	static constexpr bool isActive_ = true;

private:
	using Column_type = typename Master_matrix::Column_type;

	void _add_to(const typename Column_type::Column_type& column, std::set<index>& set);
	bool _is_negative_in_pair(index columnIndex) const;
	bool _is_paired(index columnIndex) const;

	index _positive_vine_swap(index columnIndex1, index columnIndex2);
	index _positive_negative_vine_swap(index columnIndex1, index columnIndex2);
	index _negative_positive_vine_swap(index columnIndex1, index columnIndex2);
	index _negative_vine_swap(index columnIndex1, index columnIndex2);

	std::function<bool(index,index)> birthComp_;	// for F x F
	std::function<bool(index,index)> deathComp_;	// for G x G
};

template<class Master_matrix>
inline Custom_chain_vine_swap<Master_matrix>::Custom_chain_vine_swap(
		matrix_type& matrix, 
		std::function<bool(index,index)> birthComparator, 
		std::function<bool(index,index)> deathComparator)
	: matrix_(&matrix), birthComp_(std::move(birthComparator)), deathComp_(std::move(deathComparator))
{}

template<class Master_matrix>
inline Custom_chain_vine_swap<Master_matrix>::Custom_chain_vine_swap(
		const Custom_chain_vine_swap &matrixToCopy)
	: matrix_(matrixToCopy.matrix_),
	  birthComp_(matrixToCopy.birthComp_),
	  deathComp_(matrixToCopy.deathComp_)
{}

template<class Master_matrix>
inline Custom_chain_vine_swap<Master_matrix>::Custom_chain_vine_swap(Custom_chain_vine_swap<Master_matrix> &&other) noexcept
	: matrix_(other.matrix_),
	  birthComp_(std::move(other.birthComp_)),
	  deathComp_(std::move(other.deathComp_))
{}

template<class Master_matrix>
inline index Custom_chain_vine_swap<Master_matrix>::vine_swap_with_z_eq_1_case(index columnIndex1, index columnIndex2)
{
	if (_is_negative_in_pair(columnIndex1) && _is_negative_in_pair(columnIndex2))
		return _negative_vine_swap(columnIndex1, columnIndex2);

	if (_is_negative_in_pair(columnIndex1))
		return _negative_positive_vine_swap(columnIndex1, columnIndex2);

	if (_is_negative_in_pair(columnIndex2))
		return _positive_negative_vine_swap(columnIndex1, columnIndex2);

	return _positive_vine_swap(columnIndex1, columnIndex2);
}

template<class Master_matrix>
inline index Custom_chain_vine_swap<Master_matrix>::vine_swap(index columnIndex1, index columnIndex2)
{
	if (_is_negative_in_pair(columnIndex1) && _is_negative_in_pair(columnIndex2)){
//		if (!matrix_->at(columnIndex2).is_non_zero(matrix_->at(columnIndex1).get_lowest_simplex_index())){
		if (!matrix_->at(columnIndex2).is_non_zero(matrix_->at(columnIndex1).get_pivot())){
			return columnIndex1;
		}
		return _negative_vine_swap(columnIndex1, columnIndex2);
	}

	if (_is_negative_in_pair(columnIndex1)){
		if (!matrix_->at(columnIndex2).is_non_zero(matrix_->at(columnIndex1).get_pivot())){
			return columnIndex1;
		}
		return _negative_positive_vine_swap(columnIndex1, columnIndex2);
	}

	if (_is_negative_in_pair(columnIndex2)){
		if (!matrix_->at(columnIndex2).is_non_zero(matrix_->at(columnIndex1).get_pivot())){
			return columnIndex1;
		}
		return _positive_negative_vine_swap(columnIndex1, columnIndex2);
	}

	if (!matrix_->at(columnIndex2).is_non_zero(matrix_->at(columnIndex1).get_pivot())){
		return columnIndex1;
	}
	return _positive_vine_swap(columnIndex1, columnIndex2);
}

template<class Master_matrix>
inline Custom_chain_vine_swap<Master_matrix> &Custom_chain_vine_swap<Master_matrix>::operator=(
		Custom_chain_vine_swap<Master_matrix> other)
{
	std::swap(birthComp_, other.birthComp_);
	std::swap(deathComp_, other.deathComp_);
	return *this;
}

template<class Master_matrix>
inline void Custom_chain_vine_swap<Master_matrix>::_add_to(
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
inline bool Custom_chain_vine_swap<Master_matrix>::_is_negative_in_pair(index columnIndex) const
{
	const auto& col = matrix_->at(columnIndex);
	if (!col.is_paired()) return false;
	return col.get_pivot() > matrix_->at(col.get_paired_chain_index()).get_pivot();
}

template<class Master_matrix>
inline bool Custom_chain_vine_swap<Master_matrix>::_is_paired(index columnIndex) const
{
	return matrix_->at(columnIndex).is_paired();
}

template<class Master_matrix>
inline index Custom_chain_vine_swap<Master_matrix>::_positive_vine_swap(index columnIndex1, index columnIndex2)
{
	if (!_is_paired(columnIndex1)){		// F x *
		if (!_is_paired(columnIndex2) && birthComp_(columnIndex1, columnIndex2)){
			matrix_->at(columnIndex2) += matrix_->at(columnIndex1);
			return columnIndex1;
		}
		matrix_->at(columnIndex1) += matrix_->at(columnIndex2);
		return columnIndex2;
	}

	if (!_is_paired(columnIndex2)){		// G x F
		matrix_->at(columnIndex2) += matrix_->at(columnIndex1);
		return columnIndex1;
	}

	// G x G
	if (deathComp_(columnIndex1, columnIndex2)) // == if (matrix_->at(pairedIndex1).get_pivot() < matrix_->at(pairedIndex2).get_pivot()) ???
	{
		matrix_->at(matrix_->at(columnIndex2).get_paired_chain_index()) += matrix_->at(matrix_->at(columnIndex1).get_paired_chain_index());
		matrix_->at(columnIndex2) += matrix_->at(columnIndex1);

		return columnIndex1;
	}

	matrix_->at(matrix_->at(columnIndex1).get_paired_chain_index()) += matrix_->at(matrix_->at(columnIndex2).get_paired_chain_index());
	matrix_->at(columnIndex1) += matrix_->at(columnIndex2);

	return columnIndex2;
}

template<class Master_matrix>
inline index Custom_chain_vine_swap<Master_matrix>::_positive_negative_vine_swap(index columnIndex1, index columnIndex2)
{
	matrix_->at(columnIndex2) += matrix_->at(columnIndex1);

	return columnIndex1;
}

template<class Master_matrix>
inline index Custom_chain_vine_swap<Master_matrix>::_negative_positive_vine_swap(index columnIndex1, index columnIndex2)
{
	matrix_->at(columnIndex1) += matrix_->at(columnIndex2);

	return columnIndex2;
}

template<class Master_matrix>
inline index Custom_chain_vine_swap<Master_matrix>::_negative_vine_swap(index columnIndex1, index columnIndex2)
{
	index pairedIndex1 = matrix_->at(columnIndex1).get_paired_chain_index();
	index pairedIndex2 = matrix_->at(columnIndex2).get_paired_chain_index();

	if (matrix_->at(pairedIndex1).get_pivot() < matrix_->at(pairedIndex2).get_pivot())
	{
		matrix_->at(pairedIndex2) += matrix_->at(pairedIndex1);
		matrix_->at(columnIndex2) += matrix_->at(columnIndex1);

		return columnIndex1;
	}

	matrix_->at(pairedIndex1) += matrix_->at(pairedIndex2);
	matrix_->at(columnIndex1) += matrix_->at(columnIndex2);

	return columnIndex2;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // CUSTOM_CHAIN_VINE_SWAP_H

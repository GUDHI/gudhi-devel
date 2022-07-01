/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef Chain_REP_CYCLES_H
#define Chain_REP_CYCLES_H

#include <utility>
#include <vector>

#include "../utilities.h"  //type definitions

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class Chain_representative_cycles
{
public:
	using cycle_type = std::vector<index>;

	void update_representative_cycles();

	const std::vector<cycle_type>& get_representative_cycles();
	const cycle_type& get_representative_cycle(Bar& bar);

	Chain_representative_cycles& operator=(Chain_representative_cycles other);
	template<class Friend_master_matrix>
	friend void swap(Chain_representative_cycles<Friend_master_matrix>& base1,
					 Chain_representative_cycles<Friend_master_matrix>& base2);

protected:
	using matrix_type = typename Master_matrix::column_container_type;
	using dictionnary_type = typename Master_matrix::template dictionnary_type<index>;
	using column_type = typename Master_matrix::Column_type;

	Chain_representative_cycles(matrix_type& matrix, dictionnary_type& pivotToPosition);
	Chain_representative_cycles(Chain_representative_cycles& matrixToCopy);
	Chain_representative_cycles(Chain_representative_cycles&& other) noexcept;

	static constexpr bool isActive_ = true;

private:
	matrix_type& matrix_;
	dictionnary_type& pivotToPosition_;
	std::vector<cycle_type> representativeCycles_;
	std::vector<int> birthToCycle_;
};

template<class Master_matrix>
inline Chain_representative_cycles<Master_matrix>::Chain_representative_cycles(matrix_type &matrix, dictionnary_type &pivotToPosition)
	: matrix_(matrix), pivotToPosition_(pivotToPosition)
{}

template<class Master_matrix>
inline Chain_representative_cycles<Master_matrix>::Chain_representative_cycles(Chain_representative_cycles<Master_matrix>& matrixToCopy)
	: matrix_(matrixToCopy.matrix_),
	  pivotToPosition_(matrixToCopy.pivotToPosition_),
	  representativeCycles_(matrixToCopy.representativeCycles_),
	  birthToCycle_(matrixToCopy.birthToCycle_)
{}

template<class Master_matrix>
inline Chain_representative_cycles<Master_matrix>::Chain_representative_cycles(Chain_representative_cycles<Master_matrix>&& other) noexcept
	: matrix_(std::move(other.matrix_)),
	  pivotToPosition_(std::move(other.pivotToPosition_)),
	  representativeCycles_(std::move(other.representativeCycles_)),
	  birthToCycle_(std::move(other.birthToCycle_))
{}

template<class Master_matrix>
inline void Chain_representative_cycles<Master_matrix>::update_representative_cycles()
{
	birthToCycle_.clear();
	birthToCycle_.resize(matrix_.size(), -1);
	representativeCycles_.clear();

	for (index i = 0; i < matrix_.size(); i++){
		column_type &col = matrix_.at(pivotToPosition_.at(i));
		if (!col.is_paired() || i < col.get_paired_chain_index()){
			cycle_type cycle;
			for (typename column_type::Cell c : col){
				cycle.push_back(c.get_row_index());
			}
			representativeCycles_.push_back(cycle);
			birthToCycle_.at(i) = representativeCycles_.size() - 1;
		}
	}
}

template<class Master_matrix>
inline const std::vector<typename Chain_representative_cycles<Master_matrix>::cycle_type> &
Chain_representative_cycles<Master_matrix>::get_representative_cycles()
{
	if (representativeCycles_.empty()) update_representative_cycles();
	return representativeCycles_;
}

template<class Master_matrix>
inline const typename Chain_representative_cycles<Master_matrix>::cycle_type &
Chain_representative_cycles<Master_matrix>::get_representative_cycle(Bar &bar)
{
	if (representativeCycles_.empty()) update_representative_cycles();
	return representativeCycles_.at(birthToCycle_.at(bar.birth));
}

template<class Master_matrix>
inline Chain_representative_cycles<Master_matrix> &Chain_representative_cycles<Master_matrix>::operator=(Chain_representative_cycles<Master_matrix> other)
{
	std::swap(matrix_, other.matrix_);
	std::swap(pivotToPosition_, other.pivotToPosition_);
	std::swap(representativeCycles_, other.representativeCycles_);
	std::swap(birthToCycle_, other.birthToCycle_);
	return *this;
}

template<class Friend_master_matrix>
inline void swap(Chain_representative_cycles<Friend_master_matrix>& base1, Chain_representative_cycles<Friend_master_matrix>& base2)
{
	std::swap(base1.matrix_, base2.matrix_);
	std::swap(base1.pivotToPosition_, base2.pivotToPosition_);
	std::swap(base1.representativeCycles_, base2.representativeCycles_);
	std::swap(base1.birthToCycle_, base2.birthToCycle_);
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // Chain_REP_CYCLES_H

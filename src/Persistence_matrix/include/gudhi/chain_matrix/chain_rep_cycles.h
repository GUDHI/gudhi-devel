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

#include "../utilities/utilities.h"  //type definitions

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class Chain_representative_cycles
{
public:
	using cycle_type = std::vector<index>;
	using matrix_type = typename Master_matrix::column_container_type;
	using dictionnary_type = typename Master_matrix::template dictionnary_type<index>;

	Chain_representative_cycles(matrix_type& matrix, dictionnary_type& pivotToPosition);
	Chain_representative_cycles(const Chain_representative_cycles& matrixToCopy);
	Chain_representative_cycles(Chain_representative_cycles&& other) noexcept;

	void update_representative_cycles();

	const std::vector<cycle_type>& get_representative_cycles();
	const cycle_type& get_representative_cycle(const Bar& bar);

	Chain_representative_cycles& operator=(Chain_representative_cycles other);
	friend void swap(Chain_representative_cycles& base1,
					 Chain_representative_cycles& base2){
		std::swap(base1.matrix_, base2.matrix_);
		std::swap(base1.pivotToPosition_, base2.pivotToPosition_);
		base1.representativeCycles_.swap(base2.representativeCycles_);
		base1.birthToCycle_.swap(base2.birthToCycle_);
	}

protected:
	using column_type = typename Master_matrix::Column_type;

	static constexpr bool isActive_ = true;

private:
	matrix_type* matrix_;
	dictionnary_type* pivotToPosition_;
	std::vector<cycle_type> representativeCycles_;
	std::vector<int> birthToCycle_;
};

template<class Master_matrix>
inline Chain_representative_cycles<Master_matrix>::Chain_representative_cycles(matrix_type &matrix, dictionnary_type &pivotToPosition)
	: matrix_(&matrix), pivotToPosition_(&pivotToPosition)
{}

template<class Master_matrix>
inline Chain_representative_cycles<Master_matrix>::Chain_representative_cycles(const Chain_representative_cycles<Master_matrix>& matrixToCopy)
	: matrix_(matrixToCopy.matrix_),
	  pivotToPosition_(matrixToCopy.pivotToPosition_),
	  representativeCycles_(matrixToCopy.representativeCycles_),
	  birthToCycle_(matrixToCopy.birthToCycle_)
{}

template<class Master_matrix>
inline Chain_representative_cycles<Master_matrix>::Chain_representative_cycles(Chain_representative_cycles<Master_matrix>&& other) noexcept
	: matrix_(std::exchange(other.matrix_, nullptr)),
	  pivotToPosition_(std::exchange(other.pivotToPosition_, nullptr)),
	  representativeCycles_(std::move(other.representativeCycles_)),
	  birthToCycle_(std::move(other.birthToCycle_))
{}

template<class Master_matrix>
inline void Chain_representative_cycles<Master_matrix>::update_representative_cycles()
{
	birthToCycle_.clear();
	birthToCycle_.resize(matrix_->size(), -1);
	representativeCycles_.clear();

	for (index i = 0; i < matrix_->size(); i++){
		column_type &col = matrix_->at(pivotToPosition_->at(i));
		if (!col.is_paired() || i < col.get_paired_chain_index()){
			cycle_type cycle;
			for (typename column_type::Cell c : col){
				cycle.push_back(c.get_row_index());
			}
			representativeCycles_.push_back(cycle);
			birthToCycle_[i] = representativeCycles_.size() - 1;
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
Chain_representative_cycles<Master_matrix>::get_representative_cycle(const Bar &bar)
{
	if (representativeCycles_.empty()) update_representative_cycles();
	return representativeCycles_[birthToCycle_[bar.birth]];
}

template<class Master_matrix>
inline Chain_representative_cycles<Master_matrix> &Chain_representative_cycles<Master_matrix>::operator=(Chain_representative_cycles<Master_matrix> other)
{
	std::swap(matrix_, other.matrix_);
	std::swap(pivotToPosition_, other.pivotToPosition_);
	representativeCycles_.swap(other.representativeCycles_);
	birthToCycle_.swap(other.birthToCycle_);
	return *this;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // Chain_REP_CYCLES_H

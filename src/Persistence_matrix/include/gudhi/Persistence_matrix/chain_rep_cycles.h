/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022-23 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_CHAIN_REP_CYCLES_H
#define PM_CHAIN_REP_CYCLES_H

#include <utility>		//std::move
#include <algorithm>	//std::sort
#include <vector>

namespace Gudhi {
namespace persistence_matrix {

struct Dummy_chain_representative_cycles{
	friend void swap([[maybe_unused]] Dummy_chain_representative_cycles& d1, [[maybe_unused]] Dummy_chain_representative_cycles& d2){}

	Dummy_chain_representative_cycles(){}
	Dummy_chain_representative_cycles([[maybe_unused]] const Dummy_chain_representative_cycles& matrixToCopy){}
	Dummy_chain_representative_cycles([[maybe_unused]] Dummy_chain_representative_cycles&& other) noexcept{}
};

template<class Master_matrix>
class Chain_representative_cycles
{
public:
	using index = typename Master_matrix::index;
	using Bar = typename Master_matrix::Bar;
	using cycle_type = std::vector<index>;	//TODO: add coefficients
	using matrix_type = typename Master_matrix::column_container_type;
	using dictionnary_type = typename Master_matrix::template dictionnary_type<index>;

	Chain_representative_cycles();
	Chain_representative_cycles(const Chain_representative_cycles& matrixToCopy);
	Chain_representative_cycles(Chain_representative_cycles&& other) noexcept;

	void update_representative_cycles();

	const std::vector<cycle_type>& get_representative_cycles();
	const cycle_type& get_representative_cycle(const Bar& bar);

	Chain_representative_cycles& operator=(Chain_representative_cycles other);
	friend void swap(Chain_representative_cycles& base1,
					 Chain_representative_cycles& base2){
		base1.representativeCycles_.swap(base2.representativeCycles_);
		base1.birthToCycle_.swap(base2.birthToCycle_);
	}

private:
	using chain_matrix = typename Master_matrix::Chain_matrix_type;

	std::vector<cycle_type> representativeCycles_;
	std::vector<int> birthToCycle_;

	constexpr chain_matrix* _matrix() { return static_cast<chain_matrix*>(this); }
	constexpr const chain_matrix* _matrix() const { return static_cast<const chain_matrix*>(this); }
};

template<class Master_matrix>
inline Chain_representative_cycles<Master_matrix>::Chain_representative_cycles()
{}

template<class Master_matrix>
inline Chain_representative_cycles<Master_matrix>::Chain_representative_cycles(const Chain_representative_cycles<Master_matrix>& matrixToCopy)
	: representativeCycles_(matrixToCopy.representativeCycles_),
	  birthToCycle_(matrixToCopy.birthToCycle_)
{}

template<class Master_matrix>
inline Chain_representative_cycles<Master_matrix>::Chain_representative_cycles(Chain_representative_cycles<Master_matrix>&& other) noexcept
	: representativeCycles_(std::move(other.representativeCycles_)),
	  birthToCycle_(std::move(other.birthToCycle_))
{}

template<class Master_matrix>
inline void Chain_representative_cycles<Master_matrix>::update_representative_cycles()
{
	birthToCycle_.clear();
	birthToCycle_.resize(_matrix()->get_number_of_columns(), -1);
	representativeCycles_.clear();

	for (index i = 0; i < _matrix()->get_number_of_columns(); i++){
		auto &col = _matrix()->get_column(_matrix()->get_column_with_pivot(i));
		if (!col.is_paired() || static_cast<int>(i) < col.get_paired_chain_index()){
			cycle_type cycle;
			for (auto& c : col){
				cycle.push_back(c.get_row_index());
			}
			if constexpr (std::is_same_v<typename Master_matrix::Column_type, typename Master_matrix::Heap_column_type> 
							|| std::is_same_v<typename Master_matrix::Column_type, typename Master_matrix::Unordered_set_column_type>)
				std::sort(cycle.begin(), cycle.end());
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
	representativeCycles_.swap(other.representativeCycles_);
	birthToCycle_.swap(other.birthToCycle_);
	return *this;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // PM_CHAIN_REP_CYCLES_H

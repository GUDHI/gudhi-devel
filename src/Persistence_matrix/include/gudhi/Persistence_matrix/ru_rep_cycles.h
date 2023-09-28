/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022-23 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_RU_REP_CYCLES_H
#define PM_RU_REP_CYCLES_H

#include <utility>		//std::move
#include <algorithm>	//std::sort
#include <vector>

namespace Gudhi {
namespace persistence_matrix {

struct Dummy_ru_representative_cycles{
	friend void swap([[maybe_unused]] Dummy_ru_representative_cycles& d1, [[maybe_unused]] Dummy_ru_representative_cycles& d2){}

	Dummy_ru_representative_cycles(){}
	Dummy_ru_representative_cycles([[maybe_unused]] const Dummy_ru_representative_cycles& matrixToCopy){}
	Dummy_ru_representative_cycles([[maybe_unused]] Dummy_ru_representative_cycles&& other) noexcept{}
};

template<class Master_matrix>
class RU_representative_cycles
{
public:
	using index = typename Master_matrix::index;
	using Bar = typename Master_matrix::Bar;
	using Boundary_matrix = typename Master_matrix::Boundary_matrix_type;
	using Base_matrix = typename Master_matrix::Base_matrix_type;
	using cycle_type = std::vector<index>;	//TODO: add coefficients

	RU_representative_cycles();
	RU_representative_cycles(const RU_representative_cycles& matrixToCopy);
	RU_representative_cycles(RU_representative_cycles&& other) noexcept;

	void update_representative_cycles();

	const std::vector<cycle_type>& get_representative_cycles();
	const cycle_type& get_representative_cycle(const Bar& bar);

	RU_representative_cycles& operator=(RU_representative_cycles other);
	friend void swap(RU_representative_cycles& base1, RU_representative_cycles& base2){
		base1.representativeCycles_.swap(base2.representativeCycles_);
		base1.birthToCycle_.swap(base2.birthToCycle_);
	}

private:
	using ru_matrix = typename Master_matrix::RU_matrix_type;

	std::vector<cycle_type> representativeCycles_;
	std::vector<int> birthToCycle_;

	constexpr ru_matrix* _matrix() { return static_cast<ru_matrix*>(this); }
	constexpr const ru_matrix* _matrix() const { return static_cast<const ru_matrix*>(this); }
};

template<class Master_matrix>
inline RU_representative_cycles<Master_matrix>::RU_representative_cycles()
{}

template<class Master_matrix>
inline RU_representative_cycles<Master_matrix>::RU_representative_cycles(const RU_representative_cycles<Master_matrix>& matrixToCopy)
	: representativeCycles_(matrixToCopy.representativeCycles_),
	  birthToCycle_(matrixToCopy.birthToCycle_)
{}

template<class Master_matrix>
inline RU_representative_cycles<Master_matrix>::RU_representative_cycles(RU_representative_cycles<Master_matrix> &&other) noexcept
	: representativeCycles_(std::move(other.representativeCycles_)),
	  birthToCycle_(std::move(other.birthToCycle_))
{}

template<class Master_matrix>
inline void RU_representative_cycles<Master_matrix>::update_representative_cycles()
{
	birthToCycle_.clear();
	birthToCycle_.resize(_matrix()->reducedMatrixR_.get_number_of_columns(), -1);
	for (unsigned int i = 0; i < _matrix()->reducedMatrixR_.get_number_of_columns(); i++){
		if (_matrix()->reducedMatrixR_.is_zero_column(i)){
			representativeCycles_.push_back(cycle_type());
			for (const auto& cell : _matrix()->mirrorMatrixU_.get_column(i)){
				representativeCycles_.back().push_back(cell.get_row_index());
			}
			if constexpr (std::is_same_v<typename Master_matrix::Column_type, typename Master_matrix::Heap_column_type> 
							|| std::is_same_v<typename Master_matrix::Column_type, typename Master_matrix::Unordered_set_column_type>)
				std::sort(representativeCycles_.back().begin(), representativeCycles_.back().end());
			birthToCycle_[i] = representativeCycles_.size() - 1;
		}
	}
}

template<class Master_matrix>
inline const std::vector<typename RU_representative_cycles<Master_matrix>::cycle_type> &
RU_representative_cycles<Master_matrix>::get_representative_cycles()
{
	if (representativeCycles_.empty()) update_representative_cycles();
	return representativeCycles_;
}

template<class Master_matrix>
inline const typename RU_representative_cycles<Master_matrix>::cycle_type &
RU_representative_cycles<Master_matrix>::get_representative_cycle(const Bar &bar)
{
	if (representativeCycles_.empty()) update_representative_cycles();
	return representativeCycles_[birthToCycle_[bar.birth]];
}

template<class Master_matrix>
inline RU_representative_cycles<Master_matrix> &RU_representative_cycles<Master_matrix>::operator=(RU_representative_cycles<Master_matrix> other)
{
	representativeCycles_.swap(other.representativeCycles_);
	birthToCycle_.swap(other.birthToCycle_);
	return *this;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // PM_RU_REP_CYCLES_H

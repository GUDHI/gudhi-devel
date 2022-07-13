/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef RU_REP_CYCLES_H
#define RU_REP_CYCLES_H

#include <utility>
#include <vector>

#include "../utilities.h"  //type definitions

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class RU_representative_cycles
{
public:
	using cycle_type = std::vector<index>;

	void update_representative_cycles();

	const std::vector<cycle_type>& get_representative_cycles();
	const cycle_type& get_representative_cycle(Bar& bar);

	RU_representative_cycles& operator=(RU_representative_cycles other);
	template<class Friend_master_matrix>
	friend void swap(RU_representative_cycles<Friend_master_matrix>& base1,
					 RU_representative_cycles<Friend_master_matrix>& base2);

protected:
	using Base_matrix = typename Master_matrix::Base_matrix;

	RU_representative_cycles(Base_matrix &matrixR, Base_matrix &matrixU);
	RU_representative_cycles(const RU_representative_cycles& matrixToCopy);
	RU_representative_cycles(RU_representative_cycles&& other) noexcept;

	static constexpr bool isActive_ = true;

private:
	Base_matrix &reducedMatrixR_;
	Base_matrix &mirrorMatrixU_;
	std::vector<cycle_type> representativeCycles_;
	std::vector<int> birthToCycle_;
};

template<class Master_matrix>
inline RU_representative_cycles<Master_matrix>::RU_representative_cycles(Base_matrix &matrixR, Base_matrix &matrixU)
	: reducedMatrixR_(matrixR), mirrorMatrixU_(matrixU)
{}

template<class Master_matrix>
inline RU_representative_cycles<Master_matrix>::RU_representative_cycles(const RU_representative_cycles<Master_matrix>& matrixToCopy)
	: reducedMatrixR_(matrixToCopy.reducedMatrixR_),
	  mirrorMatrixU_(matrixToCopy.mirrorMatrixU_),
	  representativeCycles_(matrixToCopy.representativeCycles_),
	  birthToCycle_(matrixToCopy.birthToCycle_)
{}

template<class Master_matrix>
inline RU_representative_cycles<Master_matrix>::RU_representative_cycles(RU_representative_cycles<Master_matrix> &&other) noexcept
	: reducedMatrixR_(std::move(other.reducedMatrixR_)),
	  mirrorMatrixU_(std::move(other.mirrorMatrixU_)),
	  representativeCycles_(std::move(other.representativeCycles_)),
	  birthToCycle_(std::move(other.birthToCycle_))
{}

template<class Master_matrix>
inline void RU_representative_cycles<Master_matrix>::update_representative_cycles()
{
	birthToCycle_.clear();
	birthToCycle_.resize(reducedMatrixR_.get_number_of_columns(), -1);
	for (unsigned int i = 0; i < reducedMatrixR_.get_number_of_columns(); i++){
		if (reducedMatrixR_.is_zero_column(i)){
			representativeCycles_.push_back(cycle_type());
			auto column = mirrorMatrixU_.get_column(i);
			for (auto cell : column){
				representativeCycles_.back().push_back(cell.get_row_index());
			}
			birthToCycle_.at(i) = representativeCycles_.size() - 1;
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
RU_representative_cycles<Master_matrix>::get_representative_cycle(Bar &bar)
{
	if (representativeCycles_.empty()) update_representative_cycles();
	return representativeCycles_.at(birthToCycle_.at(bar.birth));
}

template<class Master_matrix>
inline RU_representative_cycles<Master_matrix> &RU_representative_cycles<Master_matrix>::operator=(RU_representative_cycles<Master_matrix> other)
{
	std::swap(reducedMatrixR_, other.reducedMatrixR_);
	std::swap(mirrorMatrixU_, other.mirrorMatrixU_);
	std::swap(representativeCycles_, other.representativeCycles_);
	std::swap(birthToCycle_, other.birthToCycle_);
	return *this;
}

template<class Friend_master_matrix>
inline void swap(RU_representative_cycles<Friend_master_matrix>& base1, RU_representative_cycles<Friend_master_matrix>& base2)
{
	std::swap(base1.reducedMatrixR_, base2.reducedMatrixR_);
	std::swap(base1.mirrorMatrixU_, base2.mirrorMatrixU_);
	std::swap(base1.representativeCycles_, base2.representativeCycles_);
	std::swap(base1.birthToCycle_, base2.birthToCycle_);
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // RU_REP_CYCLES_H

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
#include <algorithm>

#include "../utilities/utilities.h"  //type definitions
#include "../options.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class RU_representative_cycles
{
public:
	using cycle_type = std::vector<index>;	//TODO: add coefficients
	using Base_matrix = typename Master_matrix::Boundary_matrix_type;

	RU_representative_cycles(Base_matrix &matrixR, Base_matrix &matrixU);
	RU_representative_cycles(const RU_representative_cycles& matrixToCopy);
	RU_representative_cycles(RU_representative_cycles&& other) noexcept;

	void update_representative_cycles();

	const std::vector<cycle_type>& get_representative_cycles();
	const cycle_type& get_representative_cycle(const Bar& bar);

	RU_representative_cycles& operator=(RU_representative_cycles other);
	friend void swap(RU_representative_cycles& base1, RU_representative_cycles& base2){
//		std::swap(base1.reducedMatrixR_, base2.reducedMatrixR_);
//		std::swap(base1.mirrorMatrixU_, base2.mirrorMatrixU_);
		base1.representativeCycles_.swap(base2.representativeCycles_);
		base1.birthToCycle_.swap(base2.birthToCycle_);
	}

protected:
	Base_matrix *reducedMatrixR_;
	Base_matrix *mirrorMatrixU_;
	static constexpr bool isActive_ = true;

private:
	std::vector<cycle_type> representativeCycles_;
	std::vector<int> birthToCycle_;
};

template<class Master_matrix>
inline RU_representative_cycles<Master_matrix>::RU_representative_cycles(Base_matrix &matrixR, Base_matrix &matrixU)
	: reducedMatrixR_(&matrixR), mirrorMatrixU_(&matrixU)
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
	: reducedMatrixR_(other.reducedMatrixR_),
	  mirrorMatrixU_(other.mirrorMatrixU_),
	  representativeCycles_(std::move(other.representativeCycles_)),
	  birthToCycle_(std::move(other.birthToCycle_))
{}

template<class Master_matrix>
inline void RU_representative_cycles<Master_matrix>::update_representative_cycles()
{
	birthToCycle_.clear();
	birthToCycle_.resize(reducedMatrixR_->get_number_of_columns(), -1);
	for (unsigned int i = 0; i < reducedMatrixR_->get_number_of_columns(); i++){
		if (reducedMatrixR_->is_zero_column(i)){
			representativeCycles_.push_back(cycle_type());
			for (auto cell : mirrorMatrixU_->get_column(i)){
				representativeCycles_.back().push_back(cell.get_row_index());
				if constexpr (Master_matrix::Option_list::column_type == Column_types::HEAP || Master_matrix::Option_list::column_type == Column_types::UNORDERED_SET)
						std::sort(representativeCycles_.back().begin(), representativeCycles_.back().end());
			}
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
//	std::swap(reducedMatrixR_, other.reducedMatrixR_);
//	std::swap(mirrorMatrixU_, other.mirrorMatrixU_);
	representativeCycles_.swap(other.representativeCycles_);
	birthToCycle_.swap(other.birthToCycle_);
	return *this;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // RU_REP_CYCLES_H

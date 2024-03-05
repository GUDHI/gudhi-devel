/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022-23 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_RU_MATRIX_H
#define PM_RU_MATRIX_H

#include <vector>
#include <utility>	//std::swap, std::move & std::exchange
#include <iostream>	//print() only

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class RU_matrix
		: public Master_matrix::RU_pairing_option,
		  public Master_matrix::RU_vine_swap_option,
		  public Master_matrix::RU_representative_cycles_option
{
public:
	using Field_operators = typename Master_matrix::Field_operators;
	using Field_element_type = typename Master_matrix::element_type;
	using Column_type = typename Master_matrix::Column_type;
	using Row_type = typename Master_matrix::Row_type;
	using Cell_constructor = typename Master_matrix::Cell_constructor;
	using boundary_type = typename Master_matrix::boundary_type;
	using index = typename Master_matrix::index;
	using id_index = typename Master_matrix::id_index;
	using pos_index = typename Master_matrix::pos_index;
	using dimension_type = typename Master_matrix::dimension_type;

	RU_matrix(Field_operators* operators, Cell_constructor* cellConstructor);
	template<class Boundary_type = boundary_type>
	RU_matrix(const std::vector<Boundary_type>& orderedBoundaries, Field_operators* operators, Cell_constructor* cellConstructor);
	RU_matrix(unsigned int numberOfColumns, Field_operators* operators, Cell_constructor* cellConstructor);
	RU_matrix(const RU_matrix& matrixToCopy, Field_operators* operators = nullptr, Cell_constructor* cellConstructor = nullptr);
	RU_matrix(RU_matrix&& other) noexcept;

	template<class Boundary_type = boundary_type>
	void insert_boundary(const Boundary_type& boundary, dimension_type dim = -1);
	template<class Boundary_type = boundary_type>
	void insert_boundary(id_index simplexIndex, const Boundary_type& boundary, dimension_type dim = -1);
	Column_type& get_column(index columnIndex, bool inR = true);
	//get_row(rowIndex) --> simplex ID (=/= columnIndex)
	Row_type& get_row(index rowIndex, bool inR = true);
	void erase_row(index rowIndex);	//only erase row in R, as U will never have an empty row
	void remove_maximal_face(index columnIndex);
	void remove_last();

	dimension_type get_max_dimension() const;
	index get_number_of_columns() const;

	dimension_type get_column_dimension(index columnIndex) const;

	//avoid calling with specialized options or make it such that it makes sense for persistence
	//=================================================================
	void add_to(index sourceColumnIndex, index targetColumnIndex);
	void multiply_target_and_add_to(index sourceColumnIndex, const Field_element_type& coefficient, index targetColumnIndex);	//do not call with vine updates
	void multiply_source_and_add_to(const Field_element_type& coefficient, index sourceColumnIndex, index targetColumnIndex);	//do not call with vine updates

	void zero_cell(index columnIndex, index rowIndex, bool inR = true);
	void zero_column(index columnIndex, bool inR = true);
	//=================================================================
	bool is_zero_cell(index columnIndex, index rowIndex, bool inR = true) const;
	bool is_zero_column(index columnIndex, bool inR = true);

	index get_column_with_pivot(index simplexIndex) const;	//assumes that pivot exists
	index get_pivot(index columnIndex);

	void reset(Field_operators* operators, Cell_constructor* cellConstructor){
		reducedMatrixR_.reset(operators, cellConstructor);
		mirrorMatrixU_.reset(operators, cellConstructor);
		pivotToColumnIndex_.clear();
		nextEventIndex_ = 0;
		operators_ = operators;
	}

	// void set_operators(Field_operators* operators){ 
	// 	operators_ = operators;
	// 	reducedMatrixR_.set_operators(operators);
	// 	mirrorMatrixU_.set_operators(operators);
	// }

	RU_matrix& operator=(const RU_matrix& other);
	friend void swap(RU_matrix& matrix1, RU_matrix& matrix2){
		swap(static_cast<typename Master_matrix::RU_pairing_option&>(matrix1),
			 static_cast<typename Master_matrix::RU_pairing_option&>(matrix2));
		swap(static_cast<typename Master_matrix::RU_vine_swap_option&>(matrix1),
			 static_cast<typename Master_matrix::RU_vine_swap_option&>(matrix2));
		swap(static_cast<typename Master_matrix::RU_representative_cycles_option&>(matrix1),
			 static_cast<typename Master_matrix::RU_representative_cycles_option&>(matrix2));
		swap(matrix1.reducedMatrixR_, matrix2.reducedMatrixR_);
		swap(matrix1.mirrorMatrixU_, matrix2.mirrorMatrixU_);
		matrix1.pivotToColumnIndex_.swap(matrix2.pivotToColumnIndex_);
		std::swap(matrix1.nextEventIndex_, matrix2.nextEventIndex_);
		std::swap(matrix1.operators_, matrix2.operators_);
	}

	void print();  //for debug

private:
	using swap_opt = typename Master_matrix::RU_vine_swap_option;
	using pair_opt = typename Master_matrix::RU_pairing_option;
	using rep_opt = typename Master_matrix::RU_representative_cycles_option;
	using dictionnary_type = typename Master_matrix::template dictionnary_type<index>;
	using barcode_type = typename Master_matrix::barcode_type;
	using bar_dictionnary_type = typename Master_matrix::bar_dictionnary_type;
	using r_matrix_type = typename Master_matrix::Boundary_matrix_type;
	using u_matrix_type = typename Master_matrix::Base_matrix_type;

	friend rep_opt;		//direct access to the two matrices
	friend swap_opt;	//direct access to the two matrices

	r_matrix_type reducedMatrixR_;
	u_matrix_type mirrorMatrixU_;	//make U not accessible by default and add option to enable access? Inaccessible, it needs less options and we could avoid some ifs.
	dictionnary_type pivotToColumnIndex_;
	pos_index nextEventIndex_;
	Field_operators* operators_;

	void _insert_boundary(index currentIndex);
	void _initialize_U();
	void _reduce();
	void _reduce_last_column(index lastIndex);
	void _update_barcode(pos_index birth, pos_index death);
	void _add_bar(dimension_type dim, pos_index birth);

	constexpr barcode_type& _barcode();
	constexpr bar_dictionnary_type& _indexToBar();
};

template<class Master_matrix>
inline RU_matrix<Master_matrix>::RU_matrix(Field_operators* operators, Cell_constructor* cellConstructor)
	: pair_opt(),
	  swap_opt(),
	  rep_opt(),
	  reducedMatrixR_(operators, cellConstructor),
	  mirrorMatrixU_(operators, cellConstructor),
	  nextEventIndex_(0),
	  operators_(operators)
{}

template<class Master_matrix>
template<class Boundary_type>
inline RU_matrix<Master_matrix>::RU_matrix(const std::vector<Boundary_type> &orderedBoundaries, Field_operators* operators, Cell_constructor* cellConstructor)
	: pair_opt(),
	  swap_opt(),
	  rep_opt(),
	  reducedMatrixR_(orderedBoundaries, operators, cellConstructor),
	  mirrorMatrixU_(orderedBoundaries.size(), operators, cellConstructor),
	  nextEventIndex_(orderedBoundaries.size()),
	  operators_(operators)
{
	if constexpr (Master_matrix::Option_list::has_map_column_container){
		pivotToColumnIndex_.reserve(orderedBoundaries.size());
	} else {
		pivotToColumnIndex_.resize(orderedBoundaries.size(), -1);
	}

	_initialize_U();
	_reduce();
}

template<class Master_matrix>
inline RU_matrix<Master_matrix>::RU_matrix(unsigned int numberOfColumns, Field_operators* operators, Cell_constructor* cellConstructor)
	: pair_opt(),
	  swap_opt(),
	  rep_opt(),
	  reducedMatrixR_(numberOfColumns, operators, cellConstructor),
	  mirrorMatrixU_(numberOfColumns, operators, cellConstructor),
	  nextEventIndex_(0),
	  operators_(operators)
{
	if constexpr (Master_matrix::Option_list::has_map_column_container){
		pivotToColumnIndex_.reserve(numberOfColumns);
	} else {
		pivotToColumnIndex_.resize(numberOfColumns, -1);
	}
	if constexpr (Master_matrix::Option_list::has_column_pairings){
		_indexToBar().reserve(numberOfColumns);
	}
	if constexpr (Master_matrix::Option_list::has_vine_update){
		swap_opt::positionToRowIdx_.reserve(numberOfColumns);
	}
}

template<class Master_matrix>
inline RU_matrix<Master_matrix>::RU_matrix(const RU_matrix &matrixToCopy, Field_operators* operators, Cell_constructor* cellConstructor)
	: pair_opt(static_cast<const pair_opt&>(matrixToCopy)),
	  swap_opt(static_cast<const swap_opt&>(matrixToCopy)),
	  rep_opt(static_cast<const rep_opt&>(matrixToCopy)),
	  reducedMatrixR_(matrixToCopy.reducedMatrixR_, operators, cellConstructor),
	  mirrorMatrixU_(matrixToCopy.mirrorMatrixU_, operators, cellConstructor),
	  pivotToColumnIndex_(matrixToCopy.pivotToColumnIndex_),
	  nextEventIndex_(matrixToCopy.nextEventIndex_),
	  operators_(operators == nullptr ? matrixToCopy.operators_ : operators)
{}

template<class Master_matrix>
inline RU_matrix<Master_matrix>::RU_matrix(RU_matrix &&other) noexcept
	: pair_opt(std::move(static_cast<pair_opt&>(other))),
	  swap_opt(std::move(static_cast<swap_opt&>(other))),
	  rep_opt(std::move(static_cast<rep_opt&>(other))),
	  reducedMatrixR_(std::move(other.reducedMatrixR_)),
	  mirrorMatrixU_(std::move(other.mirrorMatrixU_)),
	  pivotToColumnIndex_(std::move(other.pivotToColumnIndex_)),
	  nextEventIndex_(std::exchange(other.nextEventIndex_, 0)),
	  operators_(std::exchange(other.operators_, nullptr))
{}

template<class Master_matrix>
template<class Boundary_type>
inline void RU_matrix<Master_matrix>::insert_boundary(const Boundary_type &boundary, dimension_type dim)
{
	if constexpr (Master_matrix::Option_list::has_vine_update){
		auto id = reducedMatrixR_.insert_boundary(boundary, dim);
		swap_opt::positionToRowIdx_.push_back(id);
		_insert_boundary(id);
	} else {
		_insert_boundary(reducedMatrixR_.insert_boundary(boundary, dim));
	}
}

template<class Master_matrix>
template<class Boundary_type>
inline void RU_matrix<Master_matrix>::insert_boundary(id_index simplexIndex, const Boundary_type& boundary, dimension_type dim){
	if constexpr (Master_matrix::Option_list::has_vine_update){
		swap_opt::positionToRowIdx_.push_back(simplexIndex);
	}
	_insert_boundary(reducedMatrixR_.insert_boundary(simplexIndex, boundary, dim));
}

template<class Master_matrix>
inline typename RU_matrix<Master_matrix>::Column_type &
RU_matrix<Master_matrix>::get_column(index columnIndex, bool inR)
{
	if (inR){
		return reducedMatrixR_.get_column(columnIndex);
	}
	return mirrorMatrixU_.get_column(columnIndex);
}

template<class Master_matrix>
inline typename RU_matrix<Master_matrix>::Row_type&
RU_matrix<Master_matrix>::get_row(index rowIndex, bool inR)
{
	static_assert(Master_matrix::Option_list::has_row_access,
			"'get_row' is not implemented for the chosen options.");

	if (inR){
		return reducedMatrixR_.get_row(rowIndex);
	}
	return mirrorMatrixU_.get_row(rowIndex);
}

template<class Master_matrix>
inline void RU_matrix<Master_matrix>::erase_row(index rowIndex)
{
	reducedMatrixR_.erase_row(rowIndex);
}

template<class Master_matrix>
inline void RU_matrix<Master_matrix>::remove_maximal_face(index columnIndex)
{
	static_assert(Master_matrix::Option_list::has_removable_columns && Master_matrix::Option_list::has_vine_update,
			"'remove_maximal_face' is not implemented for the chosen options.");

	//TODO: is there an easy test to verify maximality even without row access?

	for (index curr = columnIndex; curr < nextEventIndex_ - 1; ++curr) {
		swap_opt::vine_swap(curr);
	}

	remove_last();
}

template<class Master_matrix>
inline void RU_matrix<Master_matrix>::remove_last()
{
	static_assert(Master_matrix::Option_list::has_removable_columns,
			"'remove_last' is not implemented for the chosen options.");

	--nextEventIndex_;

	//assumes PosIdx == MatIdx for boundary matrices.
	if constexpr (Master_matrix::Option_list::has_column_pairings){
		if constexpr (Master_matrix::hasFixedBarcode){
			auto& bar = _barcode()[_indexToBar()[nextEventIndex_]];
			if (bar.death == -1) {	//birth
				_barcode().pop_back();	//sorted by birth and nextEventIndex_ has to be the heighest one
			} else {				//death
				bar.death = -1;
			};
			_indexToBar().pop_back();
		} else {	//birth order eventually shuffled by vine updates. No sort possible to keep the matchings.
			auto it = _indexToBar().find(nextEventIndex_);
			typename barcode_type::iterator bar = it->second;

			if (bar->death == -1) _barcode().erase(bar);
			else bar->death = -1;

			_indexToBar().erase(it);
		}
	}

	mirrorMatrixU_.remove_last();
	if constexpr (Master_matrix::Option_list::has_map_column_container){
		pivotToColumnIndex_.erase(reducedMatrixR_.remove_last());
	} else {
		id_index lastPivot = reducedMatrixR_.remove_last();
		if (lastPivot != -1) pivotToColumnIndex_[lastPivot] = -1;
	}

	if constexpr (Master_matrix::Option_list::has_vine_update){
		swap_opt::positionToRowIdx_.pop_back();
	}
}

template<class Master_matrix>
inline typename RU_matrix<Master_matrix>::dimension_type RU_matrix<Master_matrix>::get_max_dimension() const
{
	return reducedMatrixR_.get_max_dimension();
}

template<class Master_matrix>
inline typename RU_matrix<Master_matrix>::index RU_matrix<Master_matrix>::get_number_of_columns() const
{
	return reducedMatrixR_.get_number_of_columns();
}

template<class Master_matrix>
inline typename RU_matrix<Master_matrix>::dimension_type RU_matrix<Master_matrix>::get_column_dimension(index columnIndex) const
{
	return reducedMatrixR_.get_column_dimension(columnIndex);
}

template<class Master_matrix>
inline void RU_matrix<Master_matrix>::add_to(index sourceColumnIndex, index targetColumnIndex)
{
	reducedMatrixR_.add_to(sourceColumnIndex, targetColumnIndex);
	if constexpr (Master_matrix::Option_list::has_vine_update)
		mirrorMatrixU_.add_to(targetColumnIndex, sourceColumnIndex);
	else
		mirrorMatrixU_.add_to(sourceColumnIndex, targetColumnIndex);
}

template<class Master_matrix>
inline void RU_matrix<Master_matrix>::multiply_target_and_add_to(index sourceColumnIndex, const Field_element_type& coefficient, index targetColumnIndex)
{
	reducedMatrixR_.multiply_target_and_add_to(sourceColumnIndex, coefficient, targetColumnIndex);
	mirrorMatrixU_.multiply_target_and_add_to(sourceColumnIndex, coefficient, targetColumnIndex);
}

template<class Master_matrix>
inline void RU_matrix<Master_matrix>::multiply_source_and_add_to(const Field_element_type& coefficient, index sourceColumnIndex, index targetColumnIndex)
{
	reducedMatrixR_.multiply_source_and_add_to(coefficient, sourceColumnIndex, targetColumnIndex);
	mirrorMatrixU_.multiply_source_and_add_to(coefficient, sourceColumnIndex, targetColumnIndex);
}

template<class Master_matrix>
inline void RU_matrix<Master_matrix>::zero_cell(index columnIndex, index rowIndex, bool inR)
{
	if (inR){
		return reducedMatrixR_.zero_cell(columnIndex, rowIndex);
	}
	return mirrorMatrixU_.zero_cell(columnIndex, rowIndex);
}

template<class Master_matrix>
inline void RU_matrix<Master_matrix>::zero_column(index columnIndex, bool inR)
{
	if (inR){
		return reducedMatrixR_.zero_column(columnIndex);
	}
	return mirrorMatrixU_.zero_column(columnIndex);
}

template<class Master_matrix>
inline bool RU_matrix<Master_matrix>::is_zero_cell(index columnIndex, index rowIndex, bool inR) const
{
	if (inR){
		return reducedMatrixR_.is_zero_cell(columnIndex, rowIndex);
	}
	return mirrorMatrixU_.is_zero_cell(columnIndex, rowIndex);
}

template<class Master_matrix>
inline bool RU_matrix<Master_matrix>::is_zero_column(index columnIndex, bool inR)
{
	if (inR){
		return reducedMatrixR_.is_zero_column(columnIndex);
	}
	return mirrorMatrixU_.is_zero_column(columnIndex);
}

template<class Master_matrix>
inline typename RU_matrix<Master_matrix>::index RU_matrix<Master_matrix>::get_column_with_pivot(index simplexIndex) const
{
	if constexpr (Master_matrix::Option_list::has_map_column_container){
		return pivotToColumnIndex_.at(simplexIndex);
	} else {
		return pivotToColumnIndex_[simplexIndex];
	}
}

template<class Master_matrix>
inline typename RU_matrix<Master_matrix>::index RU_matrix<Master_matrix>::get_pivot(index columnIndex)
{
	return reducedMatrixR_.get_column(columnIndex).get_pivot();
}

template<class Master_matrix>
inline RU_matrix<Master_matrix> &RU_matrix<Master_matrix>::operator=(const RU_matrix& other)
{
	swap_opt::operator=(other);
	pair_opt::operator=(other);
	rep_opt::operator=(other);
	reducedMatrixR_ = other.reducedMatrixR_;
	mirrorMatrixU_ = other.mirrorMatrixU_;
	pivotToColumnIndex_ = other.pivotToColumnIndex_;
	nextEventIndex_ = other.nextEventIndex_;
	operators_ = other.operators_;
	return *this;
}

template<class Master_matrix>
inline void RU_matrix<Master_matrix>::print()
{
	std::cout << "R_matrix:\n";
	reducedMatrixR_.print();
	std::cout << "U_matrix:\n";
	mirrorMatrixU_.print();
}

template<class Master_matrix>
inline void RU_matrix<Master_matrix>::_insert_boundary(index currentIndex)
{
	if constexpr (Master_matrix::Option_list::is_z2) {
		mirrorMatrixU_.insert_column({currentIndex});
	} else {
		mirrorMatrixU_.insert_column({{currentIndex, 1}});
	}

	if constexpr (!Master_matrix::Option_list::has_map_column_container){
		while (pivotToColumnIndex_.size() <= currentIndex)
			pivotToColumnIndex_.resize(pivotToColumnIndex_.size()*2, -1);
	}
	
	_reduce_last_column(currentIndex);
	++nextEventIndex_;
}

template<class Master_matrix>
inline void RU_matrix<Master_matrix>::_initialize_U()
{
	typename std::conditional<Master_matrix::Option_list::is_z2, index, std::pair<index,Field_element_type> >::type id;
	if constexpr (!Master_matrix::Option_list::is_z2) id.second = 1;

	for (id_index i = 0; i < reducedMatrixR_.get_number_of_columns(); i++){
		if constexpr (Master_matrix::Option_list::is_z2) id = i;
		else id.first = i;
		mirrorMatrixU_.insert_column({id});
	}
}

template<class Master_matrix>
inline void RU_matrix<Master_matrix>::_reduce()
{
	auto get_column_with_pivot_ = [&](id_index pivot)->index{
		if (pivot == -1) return -1;
		if constexpr (Master_matrix::Option_list::has_map_column_container){
			auto it = pivotToColumnIndex_.find(pivot);
			if (it == pivotToColumnIndex_.end()) return -1;
			else return it->second;
		} else {
			return pivotToColumnIndex_[pivot];
		}
	};

	if constexpr (Master_matrix::Option_list::has_column_pairings){
		_indexToBar().reserve(reducedMatrixR_.get_number_of_columns());
	}
	if constexpr (Master_matrix::Option_list::has_vine_update){
		swap_opt::positionToRowIdx_.reserve(reducedMatrixR_.get_number_of_columns());
	}

	for (index i = 0; i < reducedMatrixR_.get_number_of_columns(); i++){
		if constexpr (Master_matrix::Option_list::has_vine_update){
			swap_opt::positionToRowIdx_.push_back(i);
		}
		if (!(reducedMatrixR_.is_zero_column(i)))
		{
			Column_type &curr = reducedMatrixR_.get_column(i);
			id_index pivot = curr.get_pivot();
			index currIndex = get_column_with_pivot_(pivot);

			while (pivot != -1 && currIndex != -1){
				if constexpr (Master_matrix::Option_list::is_z2){
					curr += reducedMatrixR_.get_column(currIndex);
					if constexpr (Master_matrix::Option_list::has_vine_update)
						mirrorMatrixU_.get_column(currIndex) += mirrorMatrixU_.get_column(i);
					else
						mirrorMatrixU_.get_column(i) += mirrorMatrixU_.get_column(currIndex);
				} else {
					Column_type &toadd = reducedMatrixR_.get_column(currIndex);
					Field_element_type coef = curr.get_pivot_value();
					coef = operators_->get_inverse(coef);
					coef = operators_->multiply(coef, operators_->get_characteristic() - toadd.get_pivot_value());

					curr.multiply_and_add(coef, toadd);
					mirrorMatrixU_.multiply_target_and_add_to(currIndex, coef, i);
				}

				pivot = curr.get_pivot();
				currIndex = get_column_with_pivot_(pivot);
			}

			if (pivot != -1){
				if constexpr (Master_matrix::Option_list::has_map_column_container){
					pivotToColumnIndex_.try_emplace(pivot, i);
				} else {
					pivotToColumnIndex_[pivot] = i;
				}
				_update_barcode(pivot, i);
			} else {
				_add_bar(get_column_dimension(i), i);
			}
		} else {
			_add_bar(get_column_dimension(i), i);
		}
	}
}

template<class Master_matrix>
inline void RU_matrix<Master_matrix>::_reduce_last_column(index lastIndex)
{
	auto get_column_with_pivot_ = [&](id_index pivot)->index{
		if (pivot == static_cast<id_index>(-1)) return -1;
		if constexpr (Master_matrix::Option_list::has_map_column_container){
			auto it = pivotToColumnIndex_.find(pivot);
			if (it == pivotToColumnIndex_.end()) return -1;
			else return it->second;
		} else {
			return pivotToColumnIndex_[pivot];
		}
	};

	Column_type &curr = reducedMatrixR_.get_column(lastIndex);
	if (curr.is_empty()) {
		_add_bar(get_column_dimension(lastIndex), nextEventIndex_);
		return;
	}

	id_index pivot = curr.get_pivot();
	index currIndex = get_column_with_pivot_(pivot);

	while (pivot != static_cast<id_index>(-1) && currIndex != static_cast<index>(-1)){
		if constexpr (Master_matrix::Option_list::is_z2){
			curr += reducedMatrixR_.get_column(currIndex);
			if constexpr (Master_matrix::Option_list::has_vine_update)
				mirrorMatrixU_.get_column(currIndex) += mirrorMatrixU_.get_column(lastIndex);
			else
				mirrorMatrixU_.get_column(lastIndex) += mirrorMatrixU_.get_column(currIndex);
		} else {
			Column_type &toadd = reducedMatrixR_.get_column(currIndex);
			Field_element_type coef = curr.get_pivot_value();
			coef = operators_->get_inverse(coef);
			coef = operators_->multiply(coef, operators_->get_characteristic() - toadd.get_pivot_value());

			curr.multiply_and_add(coef, toadd);
			mirrorMatrixU_.get_column(lastIndex).multiply_and_add(coef, mirrorMatrixU_.get_column(currIndex));
		}

		pivot = curr.get_pivot();
		currIndex = get_column_with_pivot_(pivot);
	}

	if (pivot != static_cast<id_index>(-1)){
		if constexpr (Master_matrix::Option_list::has_map_column_container){
			pivotToColumnIndex_.try_emplace(pivot, lastIndex);
		} else {
			pivotToColumnIndex_[pivot] = lastIndex;
		}
		_update_barcode(pivot, nextEventIndex_);
	} else {
		_add_bar(get_column_dimension(lastIndex), nextEventIndex_);
	}
}

template<class Master_matrix>
inline void RU_matrix<Master_matrix>::_update_barcode(pos_index birth, pos_index death)
{
	if constexpr (Master_matrix::Option_list::has_column_pairings){
		if constexpr (Master_matrix::hasFixedBarcode || !Master_matrix::Option_list::has_removable_columns){
			_barcode()[_indexToBar()[birth]].death = death;
			_indexToBar().push_back(_indexToBar()[birth]);
		} else {
			auto& barIt = _indexToBar().at(birth);
			barIt->death = death;
			_indexToBar().try_emplace(death, barIt);	//list so iterators are stable
		}
	}
}

template<class Master_matrix>
inline void RU_matrix<Master_matrix>::_add_bar(dimension_type dim, pos_index birth)
{
	if constexpr (Master_matrix::Option_list::has_column_pairings){
		_barcode().emplace_back(dim, birth, -1);
		if constexpr (Master_matrix::hasFixedBarcode || !Master_matrix::Option_list::has_removable_columns){
			_indexToBar().push_back(_barcode().size() - 1);
		} else {
			_indexToBar().try_emplace(birth, --_barcode().end());
		}
	}
}

template<class Master_matrix>
inline constexpr typename RU_matrix<Master_matrix>::barcode_type &RU_matrix<Master_matrix>::_barcode()
{
	if constexpr (Master_matrix::Option_list::has_vine_update)
		return swap_opt::template RU_pairing<Master_matrix>::barcode_;
	else
		return pair_opt::barcode_;
}

template<class Master_matrix>
inline constexpr typename RU_matrix<Master_matrix>::bar_dictionnary_type &RU_matrix<Master_matrix>::_indexToBar()
{
	if constexpr (Master_matrix::Option_list::has_vine_update)
		return swap_opt::template RU_pairing<Master_matrix>::indexToBar_;
	else
		return pair_opt::indexToBar_;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // PM_RU_MATRIX_H

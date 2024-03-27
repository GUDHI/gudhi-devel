/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022-23 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_POS_TO_ID_TRANSLATION_H
#define PM_POS_TO_ID_TRANSLATION_H

#include <vector>
#include <utility>		//std::swap, std::move & std::exchange
#include <algorithm>	//std::transform

namespace Gudhi {
namespace persistence_matrix {

template<class Matrix_type, class Master_matrix_type>
class Position_to_index_overlay
{
public:
	using index = typename Master_matrix_type::index;
	using id_index = typename Master_matrix_type::id_index;
	using pos_index = typename Master_matrix_type::pos_index;
	using dimension_type = typename Master_matrix_type::dimension_type;
	using Field_operators = typename Master_matrix_type::Field_operators;
	using Field_element_type = typename Master_matrix_type::element_type;
	using boundary_type = typename Master_matrix_type::boundary_type;
	using Column_type = typename Master_matrix_type::Column_type;
	using Row_type = typename Master_matrix_type::Row_type;
	using bar_type = typename Master_matrix_type::Bar;
	using barcode_type = typename Master_matrix_type::barcode_type;
	using cycle_type = typename Master_matrix_type::cycle_type;
	using cell_rep_type = typename Master_matrix_type::cell_rep_type;
	using Cell_constructor = typename Master_matrix_type::Cell_constructor;

	Position_to_index_overlay(Field_operators* operators, Cell_constructor* cellConstructor);
	template<class Boundary_type = boundary_type>
	Position_to_index_overlay(const std::vector<Boundary_type>& orderedBoundaries, Field_operators* operators, Cell_constructor* cellConstructor);
	Position_to_index_overlay(unsigned int numberOfColumns, Field_operators* operators, Cell_constructor* cellConstructor);
	//chain
	template<typename EventComparatorFunction>
	Position_to_index_overlay(
		Field_operators* operators, Cell_constructor* cellConstructor,
		EventComparatorFunction&& birthComparator, 
		EventComparatorFunction&& deathComparator);
	//chain
	template<typename EventComparatorFunction, class Boundary_type>
	Position_to_index_overlay(
		const std::vector<Boundary_type>& orderedBoundaries, 
		Field_operators* operators, Cell_constructor* cellConstructor,
		EventComparatorFunction&& birthComparator, 
		EventComparatorFunction&& deathComparator);
	//chain
	template<typename EventComparatorFunction>
	Position_to_index_overlay(
		unsigned int numberOfColumns, 
		Field_operators* operators, Cell_constructor* cellConstructor,
		EventComparatorFunction&& birthComparator, 
		EventComparatorFunction&& deathComparator);
	Position_to_index_overlay(const Position_to_index_overlay& matrixToCopy, Field_operators* operators = nullptr, Cell_constructor* cellConstructor = nullptr);
	Position_to_index_overlay(Position_to_index_overlay&& other) noexcept;

	//chain: new simplex = new ID even if the same simplex was already inserted and then removed, ie., an ID cannot come back.
	template<class Boundary_type = boundary_type>
	std::vector<cell_rep_type> insert_boundary(const Boundary_type& boundary, dimension_type dim = -1);
	template<class Boundary_type = boundary_type>
	std::vector<cell_rep_type> insert_boundary(id_index faceIndex, const Boundary_type& boundary, dimension_type dim = -1);
	//boundary
	//ru: inR = true forced
	//chain
	Column_type& get_column(pos_index position);
	//boundary
	//ru: inR = true forced
	//chain
	const Column_type& get_column(pos_index position) const;
	//get_row(rowIndex) --> simplex ID (=/= columnIndex)
	//boundary
	//ru: inR = true forced
	//chain
	Row_type& get_row(id_index rowIndex);
	//boundary
	//ru: inR = true forced
	//chain
	const Row_type& get_row(id_index rowIndex) const;
	//boundary: indirect
	//ru: indirect
	//chain: indirect
	void erase_row(id_index rowIndex);
	//boundary: update barcode if already computed, does not verify if it really was maximal
	//ru
	//chain
	void remove_maximal_face(pos_index position);
	void remove_last();

	//boundary: indirect
	//ru
	//chain: indirect
	dimension_type get_max_dimension() const;
	//boundary
	//ru
	//chain
	index get_number_of_columns() const;
	//boundary
	//ru
	//chain
	dimension_type get_column_dimension(pos_index position) const;

	//boundary: avoid calling with pairing option or make it such that it makes sense for persistence
	//ru: avoid calling with specialized options or make it such that it makes sense for persistence
	//chain: avoid calling with specialized options or make it such that it makes sense for persistence
	void add_to(pos_index sourcePosition, pos_index targetPosition);
	//boundary: avoid calling with pairing option or make it such that it makes sense for persistence
	//ru: avoid calling with specialized options or make it such that it makes sense for persistence
	//chain: avoid calling with specialized options or make it such that it makes sense for persistence
	void multiply_target_and_add_to(pos_index sourcePosition, const Field_element_type& coefficient, pos_index targetPosition);
	//boundary: avoid calling with pairing option or make it such that it makes sense for persistence
	//ru: avoid calling with specialized options or make it such that it makes sense for persistence
	//chain: avoid calling with specialized options or make it such that it makes sense for persistence
	void multiply_source_and_add_to(const Field_element_type& coefficient, pos_index sourcePosition, pos_index targetPosition);

	//boundary
	//ru: inR = true forced
	//chain
	bool is_zero_cell(pos_index position, id_index rowIndex) const;
	//boundary
	//ru: inR = true forced
	//chain: just for sanity checks as a valid chain matrix never has an empty column.
	bool is_zero_column(pos_index position);

	//ru: assumes that pivot exists
	//chain
	pos_index get_column_with_pivot(id_index faceIndex) const;	//assumes that pivot exists
	//boundary
	//ru
	//chain
	id_index get_pivot(pos_index position);

	void reset(Field_operators* operators, Cell_constructor* cellConstructor){
		matrix_.reset(operators, cellConstructor);
		positionToIndex_.clear();
		nextPosition_ = 0;
		nextIndex_ = 0;
	}

	// void set_operators(Field_operators* operators){ 
	// 	matrix_.set_operators(operators);
	// }

	Position_to_index_overlay& operator=(const Position_to_index_overlay& other);
	friend void swap(Position_to_index_overlay& matrix1,
					 Position_to_index_overlay& matrix2){
		swap(matrix1.matrix_, matrix2.matrix_);
		matrix1.positionToIndex_.swap(matrix2.positionToIndex_);
		std::swap(matrix1.nextPosition_, matrix2.nextPosition_);
		std::swap(matrix1.nextIndex_, matrix2.nextIndex_);
	}

	void print();  //for debug

	//access to optionnal methods

	//boundary
	//chain
	//ru
	const barcode_type& get_current_barcode() const;
	//chain
	//ru
	void update_representative_cycles();
	//chain
	//ru
	const std::vector<cycle_type>& get_representative_cycles();
	//chain
	//ru
	const cycle_type& get_representative_cycle(const bar_type& bar);
	//ru: returns true if barcode was changed
	bool vine_swap_with_z_eq_1_case(pos_index position);					//by column position with ordered columns
	//ru: returns true if barcode was changed
	bool vine_swap(pos_index position);										//by column position with ordered columns

private:
	Matrix_type matrix_;
	std::vector<index> positionToIndex_;
	pos_index nextPosition_;
	index nextIndex_;
};

template<class Matrix_type, class Master_matrix_type>
inline Position_to_index_overlay<Matrix_type,Master_matrix_type>::Position_to_index_overlay(Field_operators* operators, Cell_constructor* cellConstructor)
	: matrix_(operators, cellConstructor),
	  nextPosition_(0),
	  nextIndex_(0)
{}

template<class Matrix_type, class Master_matrix_type>
template<class Boundary_type>
inline Position_to_index_overlay<Matrix_type,Master_matrix_type>::Position_to_index_overlay(
		const std::vector<Boundary_type>& orderedBoundaries, Field_operators* operators, Cell_constructor* cellConstructor)
	: matrix_(orderedBoundaries, operators, cellConstructor),
	  positionToIndex_(orderedBoundaries.size()),
	  nextPosition_(orderedBoundaries.size()),
	  nextIndex_(orderedBoundaries.size())
{
	for (index i = 0; i < orderedBoundaries.size(); i++){
		positionToIndex_[i] = i;
	}
}

template<class Matrix_type, class Master_matrix_type>
inline Position_to_index_overlay<Matrix_type,Master_matrix_type>::Position_to_index_overlay(
		unsigned int numberOfColumns, Field_operators* operators, Cell_constructor* cellConstructor)
	: matrix_(numberOfColumns, operators, cellConstructor), positionToIndex_(numberOfColumns), nextPosition_(0), nextIndex_(0)
{}

template<class Matrix_type, class Master_matrix_type>
template<typename EventComparatorFunction>
inline Position_to_index_overlay<Matrix_type,Master_matrix_type>::Position_to_index_overlay(
		Field_operators* operators, Cell_constructor* cellConstructor,
		EventComparatorFunction&& birthComparator, 
		EventComparatorFunction&& deathComparator) 
	: matrix_(operators, cellConstructor, birthComparator, deathComparator), 
	  nextPosition_(0), 
	  nextIndex_(0)
{}

template<class Matrix_type, class Master_matrix_type>
template<typename EventComparatorFunction, class Boundary_type>
inline Position_to_index_overlay<Matrix_type,Master_matrix_type>::Position_to_index_overlay(
		const std::vector<Boundary_type>& orderedBoundaries, 
		Field_operators* operators, Cell_constructor* cellConstructor,
		EventComparatorFunction&& birthComparator, 
		EventComparatorFunction&& deathComparator)
	: matrix_(orderedBoundaries, operators, cellConstructor, birthComparator, deathComparator),
	  positionToIndex_(orderedBoundaries.size()),
	  nextPosition_(orderedBoundaries.size()),
	  nextIndex_(orderedBoundaries.size())
{
	for (index i = 0; i < orderedBoundaries.size(); i++){
		positionToIndex_[i] = i;
	}
}

template<class Matrix_type, class Master_matrix_type>
template<typename EventComparatorFunction>
inline Position_to_index_overlay<Matrix_type,Master_matrix_type>::Position_to_index_overlay(
		unsigned int numberOfColumns, 
		Field_operators* operators, Cell_constructor* cellConstructor,
		EventComparatorFunction&& birthComparator, 
		EventComparatorFunction&& deathComparator)
	: matrix_(numberOfColumns, operators, cellConstructor, birthComparator, deathComparator), 
	  positionToIndex_(numberOfColumns), 
	  nextPosition_(0), 
	  nextIndex_(0)
{}

template<class Matrix_type, class Master_matrix_type>
inline Position_to_index_overlay<Matrix_type,Master_matrix_type>::Position_to_index_overlay(
		const Position_to_index_overlay &matrixToCopy, Field_operators* operators, Cell_constructor* cellConstructor)
	: matrix_(matrixToCopy.matrix_, operators, cellConstructor),
	  positionToIndex_(matrixToCopy.positionToIndex_),
	  nextPosition_(matrixToCopy.nextPosition_),
	  nextIndex_(matrixToCopy.nextIndex_)
{}

template<class Matrix_type, class Master_matrix_type>
inline Position_to_index_overlay<Matrix_type,Master_matrix_type>::Position_to_index_overlay(
		Position_to_index_overlay &&other) noexcept
	: matrix_(std::move(other.matrix_)),
	  positionToIndex_(std::move(other.positionToIndex_)),
	  nextPosition_(std::exchange(other.nextPosition_, 0)),
	  nextIndex_(std::exchange(other.nextIndex_, 0))
{}

template<class Matrix_type, class Master_matrix_type>
template<class Boundary_type>
inline std::vector<typename Position_to_index_overlay<Matrix_type,Master_matrix_type>::cell_rep_type> 
Position_to_index_overlay<Matrix_type,Master_matrix_type>::insert_boundary(const Boundary_type &boundary, dimension_type dim)
{
	if (positionToIndex_.size() <= nextPosition_) {
		positionToIndex_.resize(nextPosition_ * 2 + 1);
	}

	positionToIndex_[nextPosition_++] = nextIndex_++;

	return matrix_.insert_boundary(boundary, dim);
}

template<class Matrix_type, class Master_matrix_type>
template<class Boundary_type>
inline std::vector<typename Position_to_index_overlay<Matrix_type,Master_matrix_type>::cell_rep_type> 
Position_to_index_overlay<Matrix_type,Master_matrix_type>::insert_boundary(id_index faceIndex, const Boundary_type& boundary, dimension_type dim){
	if (positionToIndex_.size() <= nextPosition_) {
		positionToIndex_.resize(nextPosition_ * 2 + 1);
	}

	positionToIndex_[nextPosition_++] = nextIndex_++;

	return matrix_.insert_boundary(faceIndex, boundary, dim);
}

template<class Matrix_type, class Master_matrix_type>
inline typename Position_to_index_overlay<Matrix_type,Master_matrix_type>::Column_type &
Position_to_index_overlay<Matrix_type,Master_matrix_type>::get_column(pos_index position)
{
	return matrix_.get_column(positionToIndex_[position]);
}

template<class Matrix_type, class Master_matrix_type>
inline const typename Position_to_index_overlay<Matrix_type,Master_matrix_type>::Column_type &
Position_to_index_overlay<Matrix_type,Master_matrix_type>::get_column(pos_index position) const
{
	return matrix_.get_column(positionToIndex_[position]);
}

template<class Matrix_type, class Master_matrix_type>
inline typename Position_to_index_overlay<Matrix_type,Master_matrix_type>::Row_type &
Position_to_index_overlay<Matrix_type,Master_matrix_type>::get_row(id_index rowIndex)
{
	return matrix_.get_row(rowIndex);
}

template<class Matrix_type, class Master_matrix_type>
inline const typename Position_to_index_overlay<Matrix_type,Master_matrix_type>::Row_type &
Position_to_index_overlay<Matrix_type,Master_matrix_type>::get_row(id_index rowIndex) const
{
	return matrix_.get_row(rowIndex);
}

template<class Matrix_type, class Master_matrix_type>
inline void Position_to_index_overlay<Matrix_type,Master_matrix_type>::erase_row(id_index rowIndex)
{
	return matrix_.erase_row(rowIndex);
}

template<class Matrix_type, class Master_matrix_type>
inline void Position_to_index_overlay<Matrix_type,Master_matrix_type>::remove_maximal_face(pos_index position)
{
	--nextPosition_;

	id_index pivot = matrix_.get_pivot(positionToIndex_[position]);
	std::vector<index> columnsToSwap(nextPosition_ - position);

	if (nextPosition_ != position){
		positionToIndex_[position] = positionToIndex_[position + 1];
		for (pos_index p = position + 1; p < nextPosition_; ++p) {
			columnsToSwap[p - position - 1] = positionToIndex_[p];
			positionToIndex_[p] = positionToIndex_[p + 1];
		}
		columnsToSwap.back() = positionToIndex_[nextPosition_];
	}

	matrix_.remove_maximal_face(pivot, columnsToSwap);
}

template<class Matrix_type, class Master_matrix_type>
inline void Position_to_index_overlay<Matrix_type,Master_matrix_type>::remove_last()
{
	--nextPosition_;
	if constexpr (Master_matrix_type::Option_list::has_vine_update){
		std::vector<index> columnsToSwap;
		matrix_.remove_maximal_face(matrix_.get_pivot(positionToIndex_[nextPosition_]), columnsToSwap);
	} else {
		matrix_.remove_last();	//linear with vine updates, so it is better to use remove_maximal_face
	}
}

template<class Matrix_type, class Master_matrix_type>
inline typename Position_to_index_overlay<Matrix_type,Master_matrix_type>::dimension_type 
Position_to_index_overlay<Matrix_type,Master_matrix_type>::get_max_dimension() const
{
	return matrix_.get_max_dimension();
}

template<class Matrix_type, class Master_matrix_type>
inline typename Position_to_index_overlay<Matrix_type,Master_matrix_type>::index Position_to_index_overlay<Matrix_type,Master_matrix_type>::get_number_of_columns() const
{
	return matrix_.get_number_of_columns();
}

template<class Matrix_type, class Master_matrix_type>
inline typename Position_to_index_overlay<Matrix_type,Master_matrix_type>::dimension_type 
Position_to_index_overlay<Matrix_type,Master_matrix_type>::get_column_dimension(pos_index position) const
{
	return matrix_.get_column_dimension(positionToIndex_[position]);
}

template<class Matrix_type, class Master_matrix_type>
inline void Position_to_index_overlay<Matrix_type,Master_matrix_type>::add_to(pos_index sourcePosition, pos_index targetPosition)
{
	return matrix_.add_to(positionToIndex_[sourcePosition], positionToIndex_[targetPosition]);
}

template<class Matrix_type, class Master_matrix_type>
inline void Position_to_index_overlay<Matrix_type,Master_matrix_type>::multiply_target_and_add_to(pos_index sourcePosition, const Field_element_type& coefficient, pos_index targetPosition)
{
	return matrix_.multiply_target_and_add_to(positionToIndex_[sourcePosition], coefficient, positionToIndex_[targetPosition]);
}

template<class Matrix_type, class Master_matrix_type>
inline void Position_to_index_overlay<Matrix_type,Master_matrix_type>::multiply_source_and_add_to(const Field_element_type& coefficient, pos_index sourcePosition, pos_index targetPosition)
{
	return matrix_.multiply_source_and_add_to(coefficient, positionToIndex_[sourcePosition], positionToIndex_[targetPosition]);
}

template<class Matrix_type, class Master_matrix_type>
inline bool Position_to_index_overlay<Matrix_type,Master_matrix_type>::is_zero_cell(pos_index position, id_index rowIndex) const
{
	return matrix_.is_zero_cell(positionToIndex_[position], rowIndex);
}

template<class Matrix_type, class Master_matrix_type>
inline bool Position_to_index_overlay<Matrix_type,Master_matrix_type>::is_zero_column(pos_index position)
{
	return matrix_.is_zero_column(positionToIndex_[position]);
}

template<class Matrix_type, class Master_matrix_type>
inline typename Position_to_index_overlay<Matrix_type,Master_matrix_type>::pos_index 
Position_to_index_overlay<Matrix_type,Master_matrix_type>::get_column_with_pivot(id_index faceIndex) const
{
	index id = matrix_.get_column_with_pivot(faceIndex);
	pos_index i = 0;
	while (positionToIndex_[i] != id) ++i;
	return i;
}

template<class Matrix_type, class Master_matrix_type>
inline typename Position_to_index_overlay<Matrix_type,Master_matrix_type>::id_index Position_to_index_overlay<Matrix_type,Master_matrix_type>::get_pivot(pos_index position)
{
	return matrix_.get_pivot(positionToIndex_[position]);
}

template<class Matrix_type, class Master_matrix_type>
inline Position_to_index_overlay<Matrix_type,Master_matrix_type>&
Position_to_index_overlay<Matrix_type,Master_matrix_type>::operator=(const Position_to_index_overlay& other)
{
	matrix_ = other.matrix_;
	positionToIndex_ = other.positionToIndex_;
	nextPosition_ = other.nextPosition_;
	nextIndex_ = other.nextIndex_;

	return *this;
}

template<class Matrix_type, class Master_matrix_type>
inline void Position_to_index_overlay<Matrix_type,Master_matrix_type>::print()
{
	return matrix_.print();
}

template<class Matrix_type, class Master_matrix_type>
inline const typename Position_to_index_overlay<Matrix_type,Master_matrix_type>::barcode_type&
Position_to_index_overlay<Matrix_type,Master_matrix_type>::get_current_barcode() const
{
	return matrix_.get_current_barcode();
}

template<class Matrix_type, class Master_matrix_type>
inline void Position_to_index_overlay<Matrix_type,Master_matrix_type>::update_representative_cycles()
{
	matrix_.update_representative_cycles();
}

template<class Matrix_type, class Master_matrix_type>
inline const std::vector<typename Position_to_index_overlay<Matrix_type,Master_matrix_type>::cycle_type>&
Position_to_index_overlay<Matrix_type,Master_matrix_type>::get_representative_cycles()
{
	return matrix_.get_representative_cycles();
}

template<class Matrix_type, class Master_matrix_type>
inline const typename Position_to_index_overlay<Matrix_type,Master_matrix_type>::cycle_type&
Position_to_index_overlay<Matrix_type,Master_matrix_type>::get_representative_cycle(const bar_type& bar)
{
	return matrix_.get_representative_cycle(bar);
}

template<class Matrix_type, class Master_matrix_type>
inline bool Position_to_index_overlay<Matrix_type,Master_matrix_type>::vine_swap_with_z_eq_1_case(pos_index position)
{
	index next = matrix_.vine_swap_with_z_eq_1_case(positionToIndex_[position], positionToIndex_[position + 1]);
	if (next == positionToIndex_[position]){
		std::swap(positionToIndex_[position],
				  positionToIndex_[position + 1]);
		return true;
	}

	return false;
}

template<class Matrix_type, class Master_matrix_type>
inline bool Position_to_index_overlay<Matrix_type, Master_matrix_type>::vine_swap(pos_index position)
{
	index next = matrix_.vine_swap(positionToIndex_[position], positionToIndex_[position + 1]);
	if (next == positionToIndex_[position]){
		std::swap(positionToIndex_[position],
				  positionToIndex_[position + 1]);
		return true;
	}

	return false;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // PM_POS_TO_ID_TRANSLATION_H

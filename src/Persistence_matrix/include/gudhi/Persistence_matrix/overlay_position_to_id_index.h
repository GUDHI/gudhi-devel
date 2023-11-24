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
#include <utility>	//std::swap, std::move & std::exchange

namespace Gudhi {
namespace persistence_matrix {

template<class Matrix_type, class Master_matrix_type>
class Position_to_id_indexation_overlay
{
public:
	using index = typename Master_matrix_type::index;
	using dimension_type = typename Master_matrix_type::dimension_type;
	using Field_element_type = typename Master_matrix_type::Field_type;
	using boundary_type = typename Master_matrix_type::boundary_type;
	using Column_type = typename Master_matrix_type::Column_type;
	using Row_type = typename Master_matrix_type::Row_type;
	using bar_type = typename Master_matrix_type::Bar;
	using barcode_type = typename Master_matrix_type::barcode_type;
	using cycle_type = typename Master_matrix_type::cycle_type;
	using cell_rep_type = typename Master_matrix_type::cell_rep_type;

	Position_to_id_indexation_overlay();
	template<class Boundary_type = boundary_type>
	Position_to_id_indexation_overlay(const std::vector<Boundary_type>& orderedBoundaries);
	Position_to_id_indexation_overlay(unsigned int numberOfColumns);
	//chain
	template<typename BirthComparatorFunction, typename DeathComparatorFunction>
	Position_to_id_indexation_overlay(
		BirthComparatorFunction&& birthComparator, 
		DeathComparatorFunction&& deathComparator);
	//chain
	template<typename BirthComparatorFunction, typename DeathComparatorFunction, class Boundary_type>
	Position_to_id_indexation_overlay(
		const std::vector<Boundary_type>& orderedBoundaries,
		BirthComparatorFunction&& birthComparator, 
		DeathComparatorFunction&& deathComparator);
	//chain
	template<typename BirthComparatorFunction, typename DeathComparatorFunction>
	Position_to_id_indexation_overlay(
		unsigned int numberOfColumns,
		BirthComparatorFunction&& birthComparator, 
		DeathComparatorFunction&& deathComparator);
	Position_to_id_indexation_overlay(const Position_to_id_indexation_overlay& matrixToCopy);
	Position_to_id_indexation_overlay(Position_to_id_indexation_overlay&& other) noexcept;

	//chain: new simplex = new ID even if the same simplex was already inserted and then removed, ie., an ID cannot come back.
	template<class Boundary_type = boundary_type>
	std::vector<cell_rep_type> insert_boundary(const Boundary_type& boundary, dimension_type dim = -1);
	//boundary
	//ru: inR = true forced
	//chain
	Column_type& get_column(index columnIndex);
	//boundary
	//ru: inR = true forced
	//chain
	const Column_type& get_column(index columnIndex) const;
	//get_row(rowIndex) --> simplex ID (=/= columnIndex)
	//boundary
	//ru: inR = true forced
	//chain
	Row_type& get_row(index rowIndex);
	//boundary
	//ru: inR = true forced
	//chain
	const Row_type& get_row(index rowIndex) const;
	//boundary: indirect
	//ru: indirect
	//chain: indirect
	void erase_row(index rowIndex);
	//boundary: update barcode if already computed, does not verify if it really was maximal
	//ru
	//chain
	void remove_maximal_simplex(index columnIndex);

	//boundary: indirect
	//ru
	//chain: indirect
	dimension_type get_max_dimension() const;
	//boundary
	//ru
	//chain
	unsigned int get_number_of_columns() const;
	//boundary
	//ru
	//chain
	dimension_type get_column_dimension(index columnIndex) const;

	//boundary: avoid calling with pairing option or make it such that it makes sense for persistence
	//ru: avoid calling with specialized options or make it such that it makes sense for persistence
	//chain: avoid calling with specialized options or make it such that it makes sense for persistence
	void add_to(index sourceColumnIndex, index targetColumnIndex);
	//boundary: avoid calling with pairing option or make it such that it makes sense for persistence
	//ru: avoid calling with specialized options or make it such that it makes sense for persistence
	//chain: avoid calling with specialized options or make it such that it makes sense for persistence
	void add_to(index sourceColumnIndex, const Field_element_type& coefficient, index targetColumnIndex);
	//boundary: avoid calling with pairing option or make it such that it makes sense for persistence
	//ru: avoid calling with specialized options or make it such that it makes sense for persistence
	//chain: avoid calling with specialized options or make it such that it makes sense for persistence
	void add_to(const Field_element_type& coefficient, index sourceColumnIndex, index targetColumnIndex);

	//boundary
	//ru: inR = true forced
	//chain
	bool is_zero_cell(index columnIndex, index rowIndex) const;
	//boundary
	//ru: inR = true forced
	//chain: just for sanity checks as a valid chain matrix never has an empty column.
	bool is_zero_column(index columnIndex);

	//ru: assumes that pivot exists
	//chain
	index get_column_with_pivot(index simplexIndex) const;	//assumes that pivot exists
	//boundary
	//ru
	//chain
	int get_pivot(index columnIndex);

	Position_to_id_indexation_overlay& operator=(Position_to_id_indexation_overlay other);
	friend void swap(Position_to_id_indexation_overlay& matrix1,
					 Position_to_id_indexation_overlay& matrix2){
		swap(matrix1.matrix_, matrix2.matrix_);
		matrix1.columnPositionToID_.swap(matrix2.columnPositionToID_);
		std::swap(matrix1.nextIndex_, matrix2.nextIndex_);
		std::swap(matrix1.nextID_, matrix2.nextID_);
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
	bool vine_swap_with_z_eq_1_case(index position);					//by column position with ordered columns
	//ru: returns true if barcode was changed
	bool vine_swap(index position);										//by column position with ordered columns

private:
	Matrix_type matrix_;
	std::vector<index> columnPositionToID_;
	index nextIndex_;
	index nextID_;
};

template<class Matrix_type, class Master_matrix_type>
inline Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::Position_to_id_indexation_overlay()
	: nextIndex_(0),
	  nextID_(0)
{}

template<class Matrix_type, class Master_matrix_type>
template<class Boundary_type>
inline Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::Position_to_id_indexation_overlay(
		const std::vector<Boundary_type>& orderedBoundaries)
	: matrix_(orderedBoundaries),
	  columnPositionToID_(orderedBoundaries.size()),
	  nextIndex_(orderedBoundaries.size()),
	  nextID_(orderedBoundaries.size())
{
	for (unsigned int i = 0; i < orderedBoundaries.size(); i++){
		columnPositionToID_[i] = i;
	}
}

template<class Matrix_type, class Master_matrix_type>
inline Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::Position_to_id_indexation_overlay(
		unsigned int numberOfColumns)
	: matrix_(numberOfColumns), columnPositionToID_(numberOfColumns), nextIndex_(0), nextID_(0)
{}

template<class Matrix_type, class Master_matrix_type>
template<typename BirthComparatorFunction, typename DeathComparatorFunction>
inline Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::Position_to_id_indexation_overlay(
		BirthComparatorFunction&& birthComparator, 
		DeathComparatorFunction&& deathComparator) 
	: matrix_(birthComparator, deathComparator), 
	  nextIndex_(0), 
	  nextID_(0)
{}

template<class Matrix_type, class Master_matrix_type>
template<typename BirthComparatorFunction, typename DeathComparatorFunction, class Boundary_type>
inline Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::Position_to_id_indexation_overlay(
		const std::vector<Boundary_type>& orderedBoundaries,
		BirthComparatorFunction&& birthComparator, 
		DeathComparatorFunction&& deathComparator)
	: matrix_(orderedBoundaries, birthComparator, deathComparator),
	  columnPositionToID_(orderedBoundaries.size()),
	  nextIndex_(orderedBoundaries.size()),
	  nextID_(orderedBoundaries.size())
{
	for (unsigned int i = 0; i < orderedBoundaries.size(); i++){
		columnPositionToID_[i] = i;
	}
}

template<class Matrix_type, class Master_matrix_type>
template<typename BirthComparatorFunction, typename DeathComparatorFunction>
inline Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::Position_to_id_indexation_overlay(
		unsigned int numberOfColumns,
		BirthComparatorFunction&& birthComparator, 
		DeathComparatorFunction&& deathComparator)
	: matrix_(numberOfColumns, birthComparator, deathComparator), 
	  columnPositionToID_(numberOfColumns), 
	  nextIndex_(0), 
	  nextID_(0)
{}

template<class Matrix_type, class Master_matrix_type>
inline Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::Position_to_id_indexation_overlay(
		const Position_to_id_indexation_overlay &matrixToCopy)
	: matrix_(matrixToCopy.matrix_),
	  columnPositionToID_(matrixToCopy.columnPositionToID_),
	  nextIndex_(matrixToCopy.nextIndex_),
	  nextID_(matrixToCopy.nextID_)
{}

template<class Matrix_type, class Master_matrix_type>
inline Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::Position_to_id_indexation_overlay(
		Position_to_id_indexation_overlay &&other) noexcept
	: matrix_(std::move(other.matrix_)),
	  columnPositionToID_(std::move(other.columnPositionToID_)),
	  nextIndex_(std::exchange(other.nextIndex_, 0)),
	  nextID_(std::exchange(other.nextID_, 0))
{}

template<class Matrix_type, class Master_matrix_type>
template<class Boundary_type>
inline std::vector<typename Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::cell_rep_type> 
Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::insert_boundary(const Boundary_type &boundary, dimension_type dim)
{
	if (columnPositionToID_.size() <= nextIndex_) {
		columnPositionToID_.resize(nextIndex_ * 2 + 1);
	}

	columnPositionToID_[nextIndex_++] = nextID_++;

	return matrix_.insert_boundary(boundary, dim);
}

template<class Matrix_type, class Master_matrix_type>
inline typename Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::Column_type &
Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::get_column(index columnIndex)
{
	return matrix_.get_column(columnPositionToID_[columnIndex]);
}

template<class Matrix_type, class Master_matrix_type>
inline const typename Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::Column_type &
Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::get_column(index columnIndex) const
{
	return matrix_.get_column(columnPositionToID_[columnIndex]);
}

template<class Matrix_type, class Master_matrix_type>
inline typename Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::Row_type &
Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::get_row(index rowIndex)
{
	return matrix_.get_row(matrix_.get_pivot(columnPositionToID_[rowIndex]));		//TODO: ???
}

template<class Matrix_type, class Master_matrix_type>
inline const typename Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::Row_type &
Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::get_row(index rowIndex) const
{
	return matrix_.get_row(matrix_.get_pivot(columnPositionToID_[rowIndex]));		//TODO: ???
}

template<class Matrix_type, class Master_matrix_type>
inline void Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::erase_row(index rowIndex)
{
	return matrix_.erase_row(matrix_.get_pivot(columnPositionToID_[rowIndex]));		//TODO: ???
}

template<class Matrix_type, class Master_matrix_type>
inline void Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::remove_maximal_simplex(index columnIndex)
{
	matrix_.remove_maximal_simplex(matrix_.get_pivot(columnPositionToID_[columnIndex]));
	--nextIndex_;
	for (unsigned int i = columnIndex; i < nextIndex_; ++i){
		columnPositionToID_[i] = columnPositionToID_[i + 1];
	}
	// if (columnPositionToID_.size() > nextIndex_ * 2)
	// 	columnPositionToID_.resize(nextIndex_ + 1);
}

template<class Matrix_type, class Master_matrix_type>
inline typename Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::dimension_type 
Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::get_max_dimension() const
{
	return matrix_.get_max_dimension();
}

template<class Matrix_type, class Master_matrix_type>
inline unsigned int Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::get_number_of_columns() const
{
	return matrix_.get_number_of_columns();
}

template<class Matrix_type, class Master_matrix_type>
inline typename Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::dimension_type 
Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::get_column_dimension(index columnIndex) const
{
	return matrix_.get_column_dimension(columnPositionToID_[columnIndex]);
}

template<class Matrix_type, class Master_matrix_type>
inline void Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::add_to(index sourceColumnIndex, index targetColumnIndex)
{
	return matrix_.add_to(columnPositionToID_[sourceColumnIndex], columnPositionToID_[targetColumnIndex]);
}

template<class Matrix_type, class Master_matrix_type>
inline void Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::add_to(index sourceColumnIndex, const Field_element_type& coefficient, index targetColumnIndex)
{
	return matrix_.add_to(columnPositionToID_[sourceColumnIndex], coefficient, columnPositionToID_[targetColumnIndex]);
}

template<class Matrix_type, class Master_matrix_type>
inline void Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::add_to(const Field_element_type& coefficient, index sourceColumnIndex, index targetColumnIndex)
{
	return matrix_.add_to(coefficient, columnPositionToID_[sourceColumnIndex], columnPositionToID_[targetColumnIndex]);
}

template<class Matrix_type, class Master_matrix_type>
inline bool Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::is_zero_cell(index columnIndex, index rowIndex) const
{
	return matrix_.is_zero_cell(columnPositionToID_[columnIndex], columnPositionToID_[rowIndex]);
}

template<class Matrix_type, class Master_matrix_type>
inline bool Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::is_zero_column(index columnIndex)
{
	return matrix_.is_zero_column(columnPositionToID_[columnIndex]);
}

template<class Matrix_type, class Master_matrix_type>
inline typename Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::index 
Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::get_column_with_pivot(index simplexIndex) const
{
	index id = matrix_.get_column_with_pivot(simplexIndex);
	unsigned int i = 0;
	while (columnPositionToID_[i] != id) ++i;
	return i;
}

template<class Matrix_type, class Master_matrix_type>
inline int Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::get_pivot(index columnIndex)
{
	return matrix_.get_pivot(columnPositionToID_[columnIndex]);
}

template<class Matrix_type, class Master_matrix_type>
inline Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>&
Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::operator=(Position_to_id_indexation_overlay other)
{
	swap(matrix_, other.matrix_);
	columnPositionToID_.swap(other.columnPositionToID_);
	std::swap(nextIndex_, other.nextIndex_);
	std::swap(nextID_, other.nextID_);

	return *this;
}

template<class Matrix_type, class Master_matrix_type>
inline void Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::print()
{
	return matrix_.print();
}

template<class Matrix_type, class Master_matrix_type>
inline const typename Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::barcode_type&
Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::get_current_barcode() const
{
	return matrix_.get_current_barcode();
}

template<class Matrix_type, class Master_matrix_type>
inline void Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::update_representative_cycles()
{
	matrix_.update_representative_cycles();
}

template<class Matrix_type, class Master_matrix_type>
inline const std::vector<typename Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::cycle_type>&
Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::get_representative_cycles()
{
	return matrix_.get_representative_cycles();
}

template<class Matrix_type, class Master_matrix_type>
inline const typename Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::cycle_type&
Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::get_representative_cycle(const bar_type& bar)
{
	return matrix_.get_representative_cycle(bar);
}

template<class Matrix_type, class Master_matrix_type>
inline bool Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::vine_swap_with_z_eq_1_case(index position)
{
	index next = matrix_.vine_swap_with_z_eq_1_case(columnPositionToID_[position], columnPositionToID_[position + 1]);
	if (next == columnPositionToID_[position]){
		std::swap(columnPositionToID_[position],
				  columnPositionToID_[position + 1]);
		return true;
	}

	return false;
}

template<class Matrix_type, class Master_matrix_type>
inline bool Position_to_id_indexation_overlay<Matrix_type, Master_matrix_type>::vine_swap(index position)
{
	index next = matrix_.vine_swap(columnPositionToID_[position], columnPositionToID_[position + 1]);
	if (next == columnPositionToID_[position]){
		std::swap(columnPositionToID_[position],
				  columnPositionToID_[position + 1]);
		return true;
	}

	return false;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // PM_POS_TO_ID_TRANSLATION_H

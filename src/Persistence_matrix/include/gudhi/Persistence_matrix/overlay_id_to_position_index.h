/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022-23 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_ID_TO_POS_TRANSLATION_H
#define PM_ID_TO_POS_TRANSLATION_H

#include <vector>
#include <cassert>
#include <utility>	//std::swap, std::move & std::exchange
#include <cmath>	//std::abs

namespace Gudhi {
namespace persistence_matrix {

template<class Matrix_type, class Master_matrix_type>
class Id_to_position_indexation_overlay
{
public:
	using index = typename Master_matrix_type::index;
	using dimension_type = typename Master_matrix_type::dimension_type;
	using Field_element_type = typename Master_matrix_type::Field_type;
	using boundary_type = typename Master_matrix_type::boundary_type;
	using Column_type = typename Master_matrix_type::Column_type;	//cast to parent for compression?
	using Row_type = typename Master_matrix_type::Row_type;
	using bar_type = typename Master_matrix_type::Bar;
	using barcode_type = typename Master_matrix_type::barcode_type;
	using cycle_type = typename Master_matrix_type::cycle_type;

	Id_to_position_indexation_overlay();
	template<class Boundary_type = boundary_type>
	Id_to_position_indexation_overlay(const std::vector<Boundary_type>& boundaries);
	Id_to_position_indexation_overlay(unsigned int numberOfColumns);
	Id_to_position_indexation_overlay(const Id_to_position_indexation_overlay& matrixToCopy);
	Id_to_position_indexation_overlay(Id_to_position_indexation_overlay&& other) noexcept;

	//boundary: does not update barcode as it needs reduction
	//ru
	template<class Boundary_type = boundary_type>
	void insert_boundary(const Boundary_type& boundary, dimension_type dim = -1);
	//chain: new simplex = new ID even if the same simplex was already inserted and then removed, ie., an ID cannot come back.
	template<class Boundary_type = boundary_type>
	void insert_boundary(index simplexIndex, const Boundary_type& boundary, dimension_type dim = -1);
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
	void remove_maximal_simplex(index simplexIndex);

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
	void add_to(index sourceColumnIndex, index targetColumnIndex);
	//boundary: avoid calling with pairing option or make it such that it makes sense for persistence
	//ru: avoid calling with specialized options or make it such that it makes sense for persistence
	void add_to(index sourceColumnIndex, const Field_element_type& coefficient, index targetColumnIndex);
	//boundary: avoid calling with pairing option or make it such that it makes sense for persistence
	//ru: avoid calling with specialized options or make it such that it makes sense for persistence
	void add_to(const Field_element_type& coefficient, index sourceColumnIndex, index targetColumnIndex);

	//boundary: avoid calling with pairing option or make it such that it makes sense for persistence
	//ru: inR = true forced, avoid calling with specialized options or make it such that it makes sense for persistence
	void zero_cell(index columnIndex, index rowIndex);
	//boundary: avoid calling with pairing option or make it such that it makes sense for persistence
	//ru: inR = true forced, avoid calling with specialized options or make it such that it makes sense for persistence
	void zero_column(index columnIndex);
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

	Id_to_position_indexation_overlay& operator=(Id_to_position_indexation_overlay other);
	friend void swap(Id_to_position_indexation_overlay& matrix1,
					 Id_to_position_indexation_overlay& matrix2){
		swap(matrix1.matrix_, matrix2.matrix_);
		matrix1.columnIDToPosition_.swap(matrix2.columnIDToPosition_);
		std::swap(matrix1.nextIndex_, matrix2.nextIndex_);
	}

	void print();  //for debug

	//access to optionnal methods

	//boundary
	//ru
	//chain
	const barcode_type& get_current_barcode();
	// //boundary
	// void swap_columns(index columnIndex1, index columnIndex2);	//to enable, row and column position have to be separated
	// //boundary
	// void swap_rows(index rowIndex1, index rowIndex2);	//to enable, row and column position have to be separated
	//boundary
	void swap_at_indices(index index1, index index2);
	//chain
	//ru
	void update_representative_cycles();
	//chain
	//ru
	const std::vector<cycle_type>& get_representative_cycles();
	//chain
	//ru
	const cycle_type& get_representative_cycle(const bar_type& bar);
	//chain: returns index which was not modified, ie new i+1
	index vine_swap_with_z_eq_1_case(index columnIndex1, index columnIndex2);	//by column id with potentielly unordered columns
	//chain: returns index which was not modified, ie new i+1
	index vine_swap(index columnIndex1, index columnIndex2);					//by column id with potentielly unordered columns;

private:
	using dictionnary_type = typename Master_matrix_type::template dictionnary_type<int>;

	Matrix_type matrix_;
	dictionnary_type columnIDToPosition_;
	index nextIndex_;
	
	int _column_id_to_position(index id) const;
};

template<class Matrix_type, class Master_matrix_type>
inline Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::Id_to_position_indexation_overlay()
{}

template<class Matrix_type, class Master_matrix_type>
template<class Boundary_type>
inline Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::Id_to_position_indexation_overlay(
		const std::vector<Boundary_type> &boundaries)
	: matrix_(boundaries),
	  columnIDToPosition_(boundaries.size()),
	  nextIndex_(boundaries.size())
{
	for (unsigned int i = 0; i < boundaries.size(); i++){
		columnIDToPosition_[i] = i;
	}
}

template<class Matrix_type, class Master_matrix_type>
inline Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::Id_to_position_indexation_overlay(
		unsigned int numberOfColumns)
	: matrix_(numberOfColumns), columnIDToPosition_(numberOfColumns), nextIndex_(0)//, lastID_(0)
{}

template<class Matrix_type, class Master_matrix_type>
inline Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::Id_to_position_indexation_overlay(
		const Id_to_position_indexation_overlay &matrixToCopy)
	: matrix_(matrixToCopy.matrix_),
	  columnIDToPosition_(matrixToCopy.columnIDToPosition_),
	  nextIndex_(matrixToCopy.nextIndex_)
{}

template<class Matrix_type, class Master_matrix_type>
inline Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::Id_to_position_indexation_overlay(
		Id_to_position_indexation_overlay &&other) noexcept
	: matrix_(std::move(other.matrix_)),
	  columnIDToPosition_(std::move(other.columnIDToPosition_)),
	  nextIndex_(std::exchange(other.nextIndex_, 0))
{}

template<class Matrix_type, class Master_matrix_type>
template<class Boundary_type>
inline void Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::insert_boundary(const Boundary_type& boundary, dimension_type dim)
{
	matrix_.insert_boundary(boundary, dim);
	if constexpr (Master_matrix_type::Option_list::has_removable_columns){
		columnIDToPosition_[nextIndex_] = columnIDToPosition_.size();
	} else {
		if (columnIDToPosition_.size() == nextIndex_) {
			columnIDToPosition_.push_back(nextIndex_);
		} else {
			columnIDToPosition_[nextIndex_] = nextIndex_;
		}
	}
	nextIndex_++;
}

template<class Matrix_type, class Master_matrix_type>
template<class Boundary_type>
inline void Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::insert_boundary(index simplexIndex, const Boundary_type& boundary, dimension_type dim)
{
	if constexpr (Master_matrix_type::Option_list::has_removable_columns){
		assert(columnIDToPosition_.find(simplexIndex) == columnIDToPosition_.end() && "Index for simplex already chosen!");
	} else {
		assert((columnIDToPosition_.size() <= simplexIndex || columnIDToPosition_[simplexIndex] == -1)  && "Index for simplex already chosen!");
	}
	matrix_.insert_boundary(boundary, dim);
	if constexpr (Master_matrix_type::Option_list::has_removable_columns){
		columnIDToPosition_[simplexIndex] = matrix_.get_number_of_columns();
	} else {
		if (columnIDToPosition_.size() <= simplexIndex) {
			columnIDToPosition_.resize(simplexIndex + 1, -1);
		}
		columnIDToPosition_[simplexIndex] = matrix_.get_number_of_columns();
	}
	if (nextIndex_ <= simplexIndex) nextIndex_ = simplexIndex + 1;
}

template<class Matrix_type, class Master_matrix_type>
inline typename Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::Column_type &
Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::get_column(index columnIndex)
{
	return matrix_.get_column(_column_id_to_position(columnIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline const typename Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::Column_type &
Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::get_column(index columnIndex) const
{
	return matrix_.get_column(_column_id_to_position(columnIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline typename Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::Row_type &
Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::get_row(index rowIndex)
{
	return matrix_.get_row(_column_id_to_position(rowIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline const typename Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::Row_type &
Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::get_row(index rowIndex) const
{
	return matrix_.get_row(_column_id_to_position(rowIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline void Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::erase_row(index rowIndex)
{
	return matrix_.erase_row(_column_id_to_position(rowIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline void Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::remove_maximal_simplex(index simplexIndex)
{
	//assumes that for boudary and ru, columnID == simplexIndex
	matrix_.remove_maximal_simplex(_column_id_to_position(simplexIndex));
	columnIDToPosition_.erase(simplexIndex);
}

template<class Matrix_type, class Master_matrix_type>
inline typename Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::dimension_type 
Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::get_max_dimension() const
{
	return matrix_.get_max_dimension();
}

template<class Matrix_type, class Master_matrix_type>
inline unsigned int Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::get_number_of_columns() const
{
	return matrix_.get_number_of_columns();
}

template<class Matrix_type, class Master_matrix_type>
inline typename Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::dimension_type 
Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::get_column_dimension(index columnIndex) const
{
	return matrix_.get_column_dimension(_column_id_to_position(columnIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline void Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::add_to(index sourceColumnIndex, index targetColumnIndex)
{
	return matrix_.add_to(_column_id_to_position(sourceColumnIndex), _column_id_to_position(targetColumnIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline void Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::add_to(index sourceColumnIndex, const Field_element_type& coefficient, index targetColumnIndex)
{
	return matrix_.add_to(_column_id_to_position(sourceColumnIndex), coefficient, _column_id_to_position(targetColumnIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline void Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::add_to(const Field_element_type& coefficient, index sourceColumnIndex, index targetColumnIndex)
{
	return matrix_.add_to(coefficient, _column_id_to_position(sourceColumnIndex), _column_id_to_position(targetColumnIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline void Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::zero_cell(index columnIndex, index rowIndex)
{
	return matrix_.zero_cell(_column_id_to_position(columnIndex), _column_id_to_position(rowIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline void Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::zero_column(index columnIndex)
{
	return matrix_.zero_column(_column_id_to_position(columnIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline bool Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::is_zero_cell(index columnIndex, index rowIndex) const
{
	return matrix_.is_zero_cell(_column_id_to_position(columnIndex), _column_id_to_position(rowIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline bool Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::is_zero_column(index columnIndex)
{
	return matrix_.is_zero_column(_column_id_to_position(columnIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline typename Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::index 
Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::get_column_with_pivot(index simplexIndex) const
{
	int pos = matrix_.get_column_with_pivot(simplexIndex);
	unsigned int i = 0;
	while (_column_id_to_position(i) != pos) ++i;
	return i;
}

template<class Matrix_type, class Master_matrix_type>
inline int Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::get_pivot(index columnIndex)
{
	return matrix_.get_pivot(_column_id_to_position(columnIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>&
Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::operator=(Id_to_position_indexation_overlay other)
{
	swap(matrix_, other.matrix_);
	columnIDToPosition_.swap(other.columnIDToPosition_);
	std::swap(nextIndex_, other.nextIndex_);

	return *this;
}

template<class Matrix_type, class Master_matrix_type>
inline void Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::print()
{
	return matrix_.print();
}

template<class Matrix_type, class Master_matrix_type>
inline const typename Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::barcode_type& 
Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::get_current_barcode()
{
	return matrix_.get_current_barcode();
}

template<class Matrix_type, class Master_matrix_type>
inline void Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::update_representative_cycles()
{
	matrix_.update_representative_cycles();
}

template<class Matrix_type, class Master_matrix_type>
inline const std::vector<typename Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::cycle_type>& 
Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::get_representative_cycles()
{
	return matrix_.get_representative_cycles();
}

template<class Matrix_type, class Master_matrix_type>
inline const typename Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::cycle_type& 
Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::get_representative_cycle(const bar_type& bar)
{
	return matrix_.get_representative_cycle(bar);
}

// template<class Matrix_type, class Master_matrix_type>
// inline void Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::swap_columns(index columnIndex1, index columnIndex2)
// {
// 	return matrix_.swap_columns(_column_id_to_position(columnIndex1), _column_id_to_position(columnIndex2));
// }

// template<class Matrix_type, class Master_matrix_type>
// inline void Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::swap_rows(index rowIndex1, index rowIndex2)
// {
// 	return matrix_.swap_rows(_column_id_to_position(rowIndex1), _column_id_to_position(rowIndex2));
// }

template<class Matrix_type, class Master_matrix_type>
inline void Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::swap_at_indices(index columnIndex1, index columnIndex2)
{
	std::swap(columnIDToPosition_.at(columnIndex1), columnIDToPosition_.at(columnIndex2));
	return matrix_.swap_at_indices(_column_id_to_position(columnIndex1), _column_id_to_position(columnIndex2));
}

template<class Matrix_type, class Master_matrix_type>
inline typename Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::index 
Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::vine_swap_with_z_eq_1_case(index columnIndex1, index columnIndex2)
{
	assert(std::abs(_column_id_to_position(columnIndex1) - _column_id_to_position(columnIndex2)) == 1 && "The columns to swap are not contiguous.");
	index first = columnIDToPosition_.at(columnIndex1) < columnIDToPosition_.at(columnIndex2) ? columnIDToPosition_.at(columnIndex1) : columnIDToPosition_.at(columnIndex2);

	bool change = matrix_.vine_swap_with_z_eq_1_case(first);

	std::swap(columnIDToPosition_.at(columnIndex1), columnIDToPosition_.at(columnIndex2));

	if (change){
		return columnIndex1;
	}

	return columnIndex2;
}

template<class Matrix_type, class Master_matrix_type>
inline typename Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::index 
Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::vine_swap(index columnIndex1, index columnIndex2)
{
	assert(std::abs(_column_id_to_position(columnIndex1) - _column_id_to_position(columnIndex2)) == 1 && "The columns to swap are not contiguous.");
	index first = columnIDToPosition_.at(columnIndex1) < columnIDToPosition_.at(columnIndex2) ? columnIDToPosition_.at(columnIndex1) : columnIDToPosition_.at(columnIndex2);

	bool change = matrix_.vine_swap(first);

	std::swap(columnIDToPosition_.at(columnIndex1), columnIDToPosition_.at(columnIndex2));

	if (change){
		return columnIndex1;
	}

	return columnIndex2;
}

template<class Matrix_type, class Master_matrix_type>
inline int Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::_column_id_to_position(index id) const{
	if constexpr (Master_matrix_type::Option_list::has_removable_columns){
		return columnIDToPosition_.at(id);
	} else {
		return columnIDToPosition_[id];
	}
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // PM_ID_TO_POS_TRANSLATION_H

#ifndef ID_TO_POS_TRANSLATION_H
#define ID_TO_POS_TRANSLATION_H

#include <vector>
#include <cassert>

#include "utilities.h"  //type definitions

namespace Gudhi {
namespace persistence_matrix {

template<class Matrix_type, class Master_matrix_type>
class Id_to_position_indexation_overlay
{
public:
	using Field_element_type = typename Master_matrix_type::Field_type;
	using boundary_type = typename Master_matrix_type::boundary_type;
	using Column_type = typename Master_matrix_type::Column_type;
	using Row_type = typename Master_matrix_type::Row_type;
	using barcode_type = typename Master_matrix_type::barcode_type;
	using cycle_type = typename Master_matrix_type::cycle_type;

	Id_to_position_indexation_overlay();
	template<class Boundary_type = boundary_type>
	Id_to_position_indexation_overlay(const std::vector<Boundary_type>& boundaries);	//simplex indices have to start at 0 and be consecutifs
	Id_to_position_indexation_overlay(int numberOfColumns);
	Id_to_position_indexation_overlay(const Id_to_position_indexation_overlay &matrixToCopy);
	Id_to_position_indexation_overlay(Id_to_position_indexation_overlay&& other) noexcept;

	template<class Boundary_type = boundary_type>
	void insert_boundary(const Boundary_type& boundary);
	Column_type& get_column(index columnIndex);
	const Column_type& get_column(index columnIndex) const;
	Row_type& get_row(index rowIndex);
	const Row_type& get_row(index rowIndex) const;
	void remove_maximal_simplex(index columnIndex);

	dimension_type get_max_dimension() const;
	unsigned int get_number_of_columns() const;

	dimension_type get_column_dimension(index columnIndex) const;

	void add_to(index sourceColumnIndex, index targetColumnIndex);
	void add_to(const Column_type& sourceColumn, index targetColumnIndex);
	void add_to(const Column_type& sourceColumn, const Field_element_type& coefficient, index targetColumnIndex);
	void add_to(const Field_element_type& coefficient, const Column_type& sourceColumn, index targetColumnIndex);

	void zero_cell(index columnIndex, index rowIndex);
	void zero_column(index columnIndex);
	bool is_zero_cell(index columnIndex, index rowIndex) const;
	bool is_zero_column(index columnIndex);

	index get_column_with_pivot(index simplexIndex) const;
	index get_pivot(index columnIndex);

	Id_to_position_indexation_overlay& operator=(Id_to_position_indexation_overlay other);
	friend void swap(Id_to_position_indexation_overlay& matrix1,
					 Id_to_position_indexation_overlay& matrix2){
		swap(matrix1.matrix_, matrix2.matrix_);
		matrix1.columnIDToPosition_.swap(matrix2.columnIDToPosition_);
		std::swap(matrix1.nextIndex_, matrix2.nextIndex_);
		std::swap(matrix1.lastID_, matrix2.lastID_);
	}

	void print();  //for debug

	//access to optionnal methods
	const barcode_type& get_current_barcode();
	void update_representative_cycles();
	const std::vector<cycle_type>& get_representative_cycles();
	const cycle_type& get_representative_cycle(const Bar& bar);
	index vine_swap_with_z_eq_1_case(index columnIndex1, index columnIndex2);	//by column id with potentielly unordered columns
	index vine_swap(index columnIndex1, index columnIndex2);					//by column id with potentielly unordered columns

private:
	using dictionnary_type = typename Master_matrix_type::template dictionnary_type<index>;

	Matrix_type matrix_;
	dictionnary_type columnIDToPosition_;
	index nextIndex_;
	index lastID_;
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
	  nextIndex_(boundaries.size()),
	  lastID_(boundaries.size() - 1)
{
	for (unsigned int i = 0; i < boundaries.size(); i++){
		columnIDToPosition_[i] = i;
	}
}

template<class Matrix_type, class Master_matrix_type>
inline Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::Id_to_position_indexation_overlay(
		int numberOfColumns)
	: matrix_(numberOfColumns), columnIDToPosition_(numberOfColumns), nextIndex_(0), lastID_(0)
{}

template<class Matrix_type, class Master_matrix_type>
inline Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::Id_to_position_indexation_overlay(
		const Id_to_position_indexation_overlay &matrixToCopy)
	: matrix_(matrixToCopy.matrix_),
	  columnIDToPosition_(matrixToCopy.columnIDToPosition_),
	  nextIndex_(matrixToCopy.nextIndex_),
	  lastID_(matrixToCopy.lastID_)
{}

template<class Matrix_type, class Master_matrix_type>
inline Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::Id_to_position_indexation_overlay(
		Id_to_position_indexation_overlay &&other) noexcept
	: matrix_(std::move(other.matrix_)),
	  columnIDToPosition_(std::move(other.columnIDToPosition_)),
	  nextIndex_(std::exchange(other.nextIndex_, 0)),
	  lastID_(std::exchange(other.lastID_, 0))
{}

template<class Matrix_type, class Master_matrix_type>
template<class Boundary_type>
inline void Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::insert_boundary(const Boundary_type &boundary)
{
	matrix_.insert_boundary(boundary);
	if constexpr (Master_matrix_type::Option_list::has_removable_columns){
		columnIDToPosition_[nextIndex_] = columnIDToPosition_.size();
	} else {
		if (columnIDToPosition_.size() == nextIndex_) {
			columnIDToPosition_.push_back(nextIndex_);
		} else {
			columnIDToPosition_[nextIndex_] = nextIndex_;
		}
	}
	lastID_ = nextIndex_++;
}

template<class Matrix_type, class Master_matrix_type>
inline typename Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::Column_type &
Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::get_column(index columnIndex)
{
	return matrix_.get_column(columnIDToPosition_.at(columnIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline const typename Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::Column_type &
Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::get_column(index columnIndex) const
{
	return matrix_.get_column(columnIDToPosition_.at(columnIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline typename Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::Row_type &
Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::get_row(index rowIndex)
{
	return matrix_.get_row(columnIDToPosition_.at(rowIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline const typename Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::Row_type &
Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::get_row(index rowIndex) const
{
	return matrix_.get_row(columnIDToPosition_.at(rowIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline void Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::remove_maximal_simplex(index columnIndex)
{
	matrix_.remove_maximal_simplex(columnIndex);
	columnIDToPosition_.erase(lastID_);
	if constexpr (Master_matrix_type::Option_list::has_removable_columns){
		auto it = columnIDToPosition_.begin();
		while (it != columnIDToPosition_.end() && it->second != columnIDToPosition_.size() - 1)
			++it;
		lastID_ = it->first;
	} else {
		for (lastID_ = 0; columnIDToPosition_[lastID_] != columnIDToPosition_.size() - 1; ++lastID_);
	}
}

template<class Matrix_type, class Master_matrix_type>
inline dimension_type Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::get_max_dimension() const
{
	return matrix_.get_max_dimension();
}

template<class Matrix_type, class Master_matrix_type>
inline unsigned int Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::get_number_of_columns() const
{
	return matrix_.get_number_of_columns();
}

template<class Matrix_type, class Master_matrix_type>
inline dimension_type Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::get_column_dimension(index columnIndex) const
{
	return matrix_.get_column_dimension(columnIDToPosition_.at(columnIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline void Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::add_to(index sourceColumnIndex, index targetColumnIndex)
{
	return matrix_.add_to(columnIDToPosition_.at(sourceColumnIndex), columnIDToPosition_.at(targetColumnIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline void Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::add_to(const Column_type& sourceColumn, index targetColumnIndex)
{
	return matrix_.add_to(sourceColumn, columnIDToPosition_.at(targetColumnIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline void Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::add_to(const Column_type& sourceColumn, const Field_element_type& coefficient, index targetColumnIndex)
{
	return matrix_.add_to(sourceColumn, coefficient, columnIDToPosition_.at(targetColumnIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline void Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::add_to(const Field_element_type& coefficient, const Column_type& sourceColumn, index targetColumnIndex)
{
	return matrix_.add_to(coefficient, sourceColumn, columnIDToPosition_.at(targetColumnIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline void Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::zero_cell(index columnIndex, index rowIndex)
{
	return matrix_.zero_cell(columnIDToPosition_.at(columnIndex), columnIDToPosition_.at(rowIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline void Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::zero_column(index columnIndex)
{
	return matrix_.zero_column(columnIDToPosition_.at(columnIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline bool Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::is_zero_cell(index columnIndex, index rowIndex) const
{
	return matrix_.is_zero_cell(columnIDToPosition_.at(columnIndex), columnIDToPosition_.at(rowIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline bool Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::is_zero_column(index columnIndex)
{
	return matrix_.is_zero_column(columnIDToPosition_.at(columnIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline index Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::get_column_with_pivot(index simplexIndex) const
{
	index pos = matrix_.get_column_with_pivot(simplexIndex);
	unsigned int i = 0;
	while (columnIDToPosition_.at(i) != pos) ++i;
	return i;
}

template<class Matrix_type, class Master_matrix_type>
inline index Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::get_pivot(index columnIndex)
{
	return matrix_.get_pivot(columnIDToPosition_.at(columnIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>&
Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::operator=(Id_to_position_indexation_overlay other)
{
	swap(matrix_, other.matrix_);
	columnIDToPosition_.swap(other.columnIDToPosition_);
	std::swap(nextIndex_, other.nextIndex_);
	std::swap(lastID_, other.lastID_);

	return *this;
}

template<class Matrix_type, class Master_matrix_type>
inline void Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::print()
{
	return matrix_.print();
}

template<class Matrix_type, class Master_matrix_type>
inline const typename Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::barcode_type& Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::get_current_barcode()
{
	return matrix_.get_current_barcode();
}

template<class Matrix_type, class Master_matrix_type>
inline void Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::update_representative_cycles()
{
	matrix_.update_representative_cycles();
}

template<class Matrix_type, class Master_matrix_type>
inline const std::vector<typename Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::cycle_type>& Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::get_representative_cycles()
{
	return matrix_.get_representative_cycles();
}

template<class Matrix_type, class Master_matrix_type>
inline const typename Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::cycle_type& Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::get_representative_cycle(const Bar& bar)
{
	return matrix_.get_representative_cycle(bar);
}

template<class Matrix_type, class Master_matrix_type>
inline index Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::vine_swap_with_z_eq_1_case(index columnIndex1, index columnIndex2)
{
	assert(std::abs(static_cast<int>(columnIDToPosition_.at(columnIndex1)) - static_cast<int>(columnIDToPosition_.at(columnIndex2))) == 1 && "The columns to swap are not contiguous.");
	index first = columnIDToPosition_.at(columnIndex1) < columnIDToPosition_.at(columnIndex2) ? columnIDToPosition_.at(columnIndex1) : columnIDToPosition_.at(columnIndex2);

	bool change = matrix_.vine_swap_with_z_eq_1_case(first);

	std::swap(columnIDToPosition_.at(columnIndex1), columnIDToPosition_.at(columnIndex2));

	if (change){
		return columnIndex1;
	}

	return columnIndex2;
}

template<class Matrix_type, class Master_matrix_type>
inline index Id_to_position_indexation_overlay<Matrix_type,Master_matrix_type>::vine_swap(index columnIndex1, index columnIndex2)
{
	assert(std::abs(static_cast<int>(columnIDToPosition_.at(columnIndex1)) - static_cast<int>(columnIDToPosition_.at(columnIndex2))) == 1 && "The columns to swap are not contiguous.");
	index first = columnIDToPosition_.at(columnIndex1) < columnIDToPosition_.at(columnIndex2) ? columnIDToPosition_.at(columnIndex1) : columnIDToPosition_.at(columnIndex2);

	bool change = matrix_.vine_swap(first);

	std::swap(columnIDToPosition_.at(columnIndex1), columnIDToPosition_.at(columnIndex2));

	if (change){
		return columnIndex1;
	}

	return columnIndex2;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // ID_TO_POS_TRANSLATION_H

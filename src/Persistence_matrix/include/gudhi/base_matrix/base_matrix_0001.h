/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef BASE_MATRIX_0001_H
#define BASE_MATRIX_0001_H

#include <iostream>
#include <vector>

#include "../utilities/utilities.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class Base_matrix_with_removals
		: public Master_matrix::Base_swap_option
{
public:
	using Field_element_type = typename Master_matrix::Field_type;
	using Column_type = typename Master_matrix::Column_type;
	using Row_type = void;

	Base_matrix_with_removals();
	template<class Container_type>
	Base_matrix_with_removals(const std::vector<Container_type>& columns);
	Base_matrix_with_removals(unsigned int numberOfColumns);
	Base_matrix_with_removals(const Base_matrix_with_removals& matrixToCopy);
	Base_matrix_with_removals(Base_matrix_with_removals&& other) noexcept;

	template<class Container_type>
	void insert_column(const Container_type& column, int columnIndex = -1);
	template<class Boundary_type>
	void insert_boundary(const Boundary_type& boundary);
	Column_type& get_column(index columnIndex);
	const Column_type& get_column(index columnIndex) const;
	Row_type get_row(index rowIndex) const;
	void erase_column(index columnIndex);
	void erase_row(index rowIndex);		// /!\ assumes row is empty

	unsigned int get_number_of_columns() const;

	void add_to(index sourceColumnIndex, index targetColumnIndex);
	template<class Cell_range>
	void add_to(const Cell_range& sourceColumn, index targetColumnIndex);
	template<class Cell_range>
	void add_to(const Cell_range& sourceColumn, const Field_element_type& coefficient, index targetColumnIndex);
	template<class Cell_range>
	void add_to(const Field_element_type& coefficient, const Cell_range& sourceColumn, index targetColumnIndex);

	void zero_cell(index columnIndex, index rowIndex);
	void zero_column(index columnIndex);
	bool is_zero_cell(index columnIndex, index rowIndex) const;
	bool is_zero_column(index columnIndex);

	Base_matrix_with_removals& operator=(Base_matrix_with_removals other);
	friend void swap(Base_matrix_with_removals& matrix1,
					 Base_matrix_with_removals& matrix2){
		swap(static_cast<typename Master_matrix::Base_swap_option&>(matrix1),
			 static_cast<typename Master_matrix::Base_swap_option&>(matrix2));
		matrix1.matrix_.swap(matrix2.matrix_);
		std::swap(matrix1.nextInsertIndex_, matrix2.nextInsertIndex_);
	}

	void print();  //for debug

private:
	using swap_opt = typename Master_matrix::Base_swap_option;

	typename Master_matrix::column_container_type matrix_;
	index nextInsertIndex_;
};

template<class Master_matrix>
inline Base_matrix_with_removals<Master_matrix>::Base_matrix_with_removals()
	: Master_matrix::Base_swap_option(matrix_),
	  nextInsertIndex_(0)
{}

template<class Master_matrix>
template<class Container_type>
inline Base_matrix_with_removals<Master_matrix>::Base_matrix_with_removals(const std::vector<Container_type> &columns)
	: Master_matrix::Base_swap_option(matrix_, columns.size()),
	  matrix_(columns.size()),
	  nextInsertIndex_(columns.size())
{
	for (unsigned int i = 0; i < columns.size(); i++){
		matrix_.try_emplace(i, Column_type(columns[i]));
	}
}

template<class Master_matrix>
inline Base_matrix_with_removals<Master_matrix>::Base_matrix_with_removals(unsigned int numberOfColumns)
	: Master_matrix::Base_swap_option(matrix_, numberOfColumns),
	  matrix_(numberOfColumns),
	  nextInsertIndex_(0)
{}

template<class Master_matrix>
inline Base_matrix_with_removals<Master_matrix>::Base_matrix_with_removals(const Base_matrix_with_removals &matrixToCopy)
	: Master_matrix::Base_swap_option(matrixToCopy),
	  matrix_(matrixToCopy.matrix_),
	  nextInsertIndex_(matrixToCopy.nextInsertIndex_)
{
	if constexpr (swap_opt::isActive_)
		swap_opt::matrix_ = &matrix_;
}

template<class Master_matrix>
inline Base_matrix_with_removals<Master_matrix>::Base_matrix_with_removals(Base_matrix_with_removals &&other) noexcept
	: Master_matrix::Base_swap_option(std::move(other)),
	  matrix_(std::move(other.matrix_)),
	  nextInsertIndex_(std::exchange(other.nextInsertIndex_, 0))
{
	if constexpr (swap_opt::isActive_)
		swap_opt::matrix_ = &matrix_;
}

template<class Master_matrix>
template<class Container_type>
inline void Base_matrix_with_removals<Master_matrix>::insert_column(const Container_type &column, int columnIndex)
{
	index id = columnIndex < 0 ? nextInsertIndex_ : columnIndex;
	assert(matrix_.find(id) == matrix_.end() && "Column already existing at given index.");
	if (columnIndex > static_cast<int>(nextInsertIndex_)) nextInsertIndex_ = columnIndex + 1;
	else if (columnIndex < 0) nextInsertIndex_++;

	if constexpr (swap_opt::isActive_){
		if (swap_opt::rowSwapped_) swap_opt::_orderRows();
		swap_opt::indexToRow_[id] = id;
		swap_opt::rowToIndex_[id] = id;
	}

	matrix_.try_emplace(id, column);
}

template<class Master_matrix>
template<class Boundary_type>
inline void Base_matrix_with_removals<Master_matrix>::insert_boundary(const Boundary_type &boundary)
{
	insert_column(boundary);
}

template<class Master_matrix>
inline typename Base_matrix_with_removals<Master_matrix>::Column_type &
Base_matrix_with_removals<Master_matrix>::get_column(index columnIndex)
{
	if constexpr (swap_opt::isActive_){
		if (swap_opt::rowSwapped_) swap_opt::_orderRows();
	}

	return matrix_.at(columnIndex);
}

template<class Master_matrix>
inline const typename Base_matrix_with_removals<Master_matrix>::Column_type &
Base_matrix_with_removals<Master_matrix>::get_column(index columnIndex) const
{
	if constexpr (swap_opt::isActive_){
		if (swap_opt::rowSwapped_) swap_opt::_orderRows();
	}

	return matrix_.at(columnIndex);
}

template<class Master_matrix>
inline typename Base_matrix_with_removals<Master_matrix>::Row_type
Base_matrix_with_removals<Master_matrix>::get_row(index rowIndex) const
{
	static_assert(Master_matrix::Option_list::has_row_access,
			"'get_row' is not implemented for the chosen options.");
}

template<class Master_matrix>
inline void Base_matrix_with_removals<Master_matrix>::erase_column(index columnIndex)
{
	if (columnIndex == nextInsertIndex_ - 1) --nextInsertIndex_;

	matrix_.erase(columnIndex);
}

template<class Master_matrix>
inline void Base_matrix_with_removals<Master_matrix>::erase_row(index rowIndex)
{
	if constexpr (swap_opt::isActive_){
		auto it = swap_opt::indexToRow_.find(rowIndex);
		swap_opt::rowToIndex_.erase(it->second);
		swap_opt::indexToRow_.erase(it);
	}
}

template<class Master_matrix>
inline unsigned int Base_matrix_with_removals<Master_matrix>::get_number_of_columns() const
{
//	return nextInsertIndex_;	//if erased columns are viewed as zero columns, otherwise use matrix size.
	return matrix_.size();
}

template<class Master_matrix>
inline void Base_matrix_with_removals<Master_matrix>::add_to(index sourceColumnIndex, index targetColumnIndex)
{
	matrix_.at(targetColumnIndex) += matrix_.at(sourceColumnIndex);
}

template<class Master_matrix>
template<class Cell_range>
inline void Base_matrix_with_removals<Master_matrix>::add_to(const Cell_range& sourceColumn, index targetColumnIndex)
{
	matrix_.at(targetColumnIndex) += sourceColumn;
}

template<class Master_matrix>
template<class Cell_range>
inline void Base_matrix_with_removals<Master_matrix>::add_to(const Cell_range& sourceColumn, const Field_element_type& coefficient, index targetColumnIndex)
{
	matrix_.at(targetColumnIndex).multiply_and_add(coefficient, sourceColumn);
}

template<class Master_matrix>
template<class Cell_range>
inline void Base_matrix_with_removals<Master_matrix>::add_to(const Field_element_type& coefficient, const Cell_range& sourceColumn, index targetColumnIndex)
{
	matrix_.at(targetColumnIndex).multiply_and_add(sourceColumn, coefficient);
}

template<class Master_matrix>
inline void Base_matrix_with_removals<Master_matrix>::zero_cell(index columnIndex, index rowIndex)
{
	static_assert(Master_matrix::Option_list::has_row_access,
			"'zero_cell' is not implemented for the chosen options.");
}

template<class Master_matrix>
inline void Base_matrix_with_removals<Master_matrix>::zero_column(index columnIndex)
{
	static_assert(Master_matrix::Option_list::has_row_access,
			"'zero_column' is not implemented for the chosen options.");
}

template<class Master_matrix>
inline bool Base_matrix_with_removals<Master_matrix>::is_zero_cell(index columnIndex, index rowIndex) const
{
	if constexpr (swap_opt::isActive_){
		return !(matrix_.at(columnIndex).is_non_zero(swap_opt::indexToRow_.at(rowIndex)));
	} else {
		return !(matrix_.at(columnIndex).is_non_zero(rowIndex));
	}
}

template<class Master_matrix>
inline bool Base_matrix_with_removals<Master_matrix>::is_zero_column(index columnIndex)
{
	return matrix_.at(columnIndex).is_empty();
}

template<class Master_matrix>
inline Base_matrix_with_removals<Master_matrix> &Base_matrix_with_removals<Master_matrix>::operator=(Base_matrix_with_removals other)
{
	swap_opt::operator=(other);
	matrix_.swap(other.matrix_);
	std::swap(nextInsertIndex_, other.nextInsertIndex_);
	return *this;
}

template<class Master_matrix>
inline void Base_matrix_with_removals<Master_matrix>::print()
{
	std::cout << "Base_matrix_with_removals:\n";
	for (unsigned int i = 0; i < nextInsertIndex_; ++i){
		const Column_type& col = matrix_.at(i);
		for (const auto e : col.get_content(nextInsertIndex_)){
			if (e == 0u) std::cout << "- ";
			else std::cout << e << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // BASE_MATRIX_0001_H

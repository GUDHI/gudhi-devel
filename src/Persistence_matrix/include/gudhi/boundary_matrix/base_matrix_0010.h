/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef BASE_MATRIX_0010_H
#define BASE_MATRIX_0010_H

#include <iostream>
#include <vector>

#include "../utilities.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class Base_matrix_with_row_access : Master_matrix::Base_swap_option, Master_matrix::Base_pairing_option{
public:
	using Column_type = typename Master_matrix::Column_type;
	using Row_type = typename Master_matrix::Column_type::Row_type;

	Base_matrix_with_row_access();
	Base_matrix_with_row_access(std::vector<boundary_type>& orderedBoundaries);
	Base_matrix_with_row_access(unsigned int numberOfColumns);
	Base_matrix_with_row_access(Base_matrix_with_row_access& matrixToCopy);
	Base_matrix_with_row_access(Base_matrix_with_row_access&& other) noexcept;

	void insert_boundary(boundary_type& boundary);
	Column_type& get_column(index columnIndex);
	Row_type& get_row(index rowIndex);
	void erase_last();

	dimension_type get_max_dimension() const;
	unsigned int get_number_of_columns() const;

	dimension_type get_column_dimension(index columnIndex) const;

	void add_to(index sourceColumnIndex, index targetColumnIndex);

	void zero_cell(index columnIndex, index rowIndex);
	void zero_column(index columnIndex);
	bool is_zero_cell(index columnIndex, index rowIndex) const;
	bool is_zero_column(index columnIndex);

	Base_matrix_with_row_access& operator=(Base_matrix_with_row_access other);
	template<class Friend_master_matrix>
	friend void swap(Base_matrix_with_row_access<Friend_master_matrix>& matrix1,
					 Base_matrix_with_row_access<Friend_master_matrix>& matrix2);

	void print();  //for debug

private:
	using swap_opt = typename Master_matrix::Base_swap_option;
	using pair_opt = typename Master_matrix::Base_pairing_option;
	using dictionnary_type = typename Master_matrix::template dictionnary_type<index>;

	typename Master_matrix::column_container_type matrix_;
	dictionnary_type pivotToColumnIndex_;
	dimension_type maxDim_;
	index nextInsertIndex_;

	void _reduce_boundary(boundary_type& boundary);
};

template<class Master_matrix>
inline Base_matrix_with_row_access<Master_matrix>::Base_matrix_with_row_access()
	: Master_matrix::Base_swap_option(matrix_, maxDim_),
	  Master_matrix::Base_pairing_option(),
	  maxDim_(-1),
	  nextInsertIndex_(0)
{}

template<class Master_matrix>
inline Base_matrix_with_row_access<Master_matrix>::Base_matrix_with_row_access(std::vector<boundary_type> &orderedBoundaries)
	: Master_matrix::Base_swap_option(matrix_, maxDim_),
	  Master_matrix::Base_pairing_option(),
	  matrix_(orderedBoundaries.size()),
	  maxDim_(0),
	  nextInsertIndex_(orderedBoundaries.size())
{
	if constexpr (swap_opt::isActive_){
		swap_opt::indexToRow_.resize(orderedBoundaries.size());
		swap_opt::rowToIndex_.resize(orderedBoundaries.size());
	}

	for (unsigned int i = 0; i < orderedBoundaries.size(); i++){
		boundary_type& b = orderedBoundaries.at(i);
		matrix_.at(i) = Column_type(b);
		if (maxDim_ < static_cast<int>(b.size()) - 1) maxDim_ = b.size() - 1;

		if constexpr (swap_opt::isActive_){
			swap_opt::indexToRow_.at(i) = i;
			swap_opt::rowToIndex_.at(i) = i;
		}
	}
	if (maxDim_ == -1) maxDim_ = 0;
}

template<class Master_matrix>
inline Base_matrix_with_row_access<Master_matrix>::Base_matrix_with_row_access(unsigned int numberOfColumns)
	: Master_matrix::Base_swap_option(matrix_, maxDim_),
	  Master_matrix::Base_pairing_option(),
	  matrix_(numberOfColumns),
	  maxDim_(0),
	  nextInsertIndex_(0)
{
	if constexpr (swap_opt::isActive_){
		swap_opt::indexToRow_.resize(numberOfColumns);
		swap_opt::rowToIndex_.resize(numberOfColumns);

		for (unsigned int i = 0; i < numberOfColumns; i++){
			swap_opt::indexToRow_.at(i) = i;
			swap_opt::rowToIndex_.at(i) = i;
		}
	}
}

template<class Master_matrix>
inline Base_matrix_with_row_access<Master_matrix>::Base_matrix_with_row_access(Base_matrix_with_row_access &matrixToCopy)
	: Master_matrix::Base_swap_option(matrixToCopy),
	  Master_matrix::Base_pairing_option(matrixToCopy),
	  matrix_(matrixToCopy.matrix_),
	  maxDim_(matrixToCopy.maxDim_),
	  nextInsertIndex_(matrixToCopy.nextInsertIndex_)
{}

template<class Master_matrix>
inline Base_matrix_with_row_access<Master_matrix>::Base_matrix_with_row_access(Base_matrix_with_row_access &&other) noexcept
	: Master_matrix::Base_swap_option(std::move(other)),
	  Master_matrix::Base_pairing_option(std::move(other)),
	  matrix_(std::move(other.matrix_)),
	  maxDim_(std::exchange(other.maxDim_, 0)),
	  nextInsertIndex_(std::move(other.nextInsertIndex_))
{}

template<class Master_matrix>
inline void Base_matrix_with_row_access<Master_matrix>::insert_boundary(boundary_type &boundary)
{
	if constexpr (swap_opt::isActive_){
		if (swap_opt::rowSwapped_) swap_opt::_orderRows();
	}

	unsigned int size = matrix_.size();
	if (size <= nextInsertIndex_) {
		if constexpr (swap_opt::isActive_){
			for (unsigned int i = size; i <= size * 2; i++){
				swap_opt::indexToRow_.push_back(i);
				swap_opt::rowToIndex_.push_back(i);
			}
		}
		matrix_.resize(size * 2);
	}

	matrix_.at(nextInsertIndex_++) = Column_type(boundary);
	if (maxDim_ < boundary.size() - 1) maxDim_ = boundary.size() - 1;
}

template<class Master_matrix>
inline typename Base_matrix_with_row_access<Master_matrix>::Column_type &Base_matrix_with_row_access<Master_matrix>::get_column(index columnIndex)
{
	if constexpr (swap_opt::isActive_){
		if (swap_opt::rowSwapped_) swap_opt::_orderRows();
	}

	return matrix_.at(columnIndex);
}

template<class Master_matrix>
inline typename Base_matrix_with_row_access<Master_matrix>::Row_type &Base_matrix_with_row_access<Master_matrix>::get_row(index rowIndex)
{
	static_assert(static_cast<int>(Master_matrix::Field_type::get_characteristic()) == -1,
			"'get_row' is not implemented for the chosen options.");
}

template<class Master_matrix>
inline void Base_matrix_with_row_access<Master_matrix>::erase_last(index simplexIndex)
{
	static_assert(static_cast<int>(Master_matrix::Field_type::get_characteristic()) == -1,
			"'erase_last' is not implemented for the chosen options.");
}

template<class Master_matrix>
inline dimension_type Base_matrix_with_row_access<Master_matrix>::get_max_dimension() const
{
	return maxDim_;
}

template<class Master_matrix>
inline unsigned int Base_matrix_with_row_access<Master_matrix>::get_number_of_columns() const
{
	return matrix_.size();
}

template<class Master_matrix>
inline dimension_type Base_matrix_with_row_access<Master_matrix>::get_column_dimension(index columnIndex) const
{
	return matrix_.at(columnIndex).get_dimension();
}

template<class Master_matrix>
inline void Base_matrix_with_row_access<Master_matrix>::add_to(index sourceColumnIndex, index targetColumnIndex)
{
	matrix_.at(targetColumnIndex) += matrix_.at(sourceColumnIndex);
}

template<class Master_matrix>
inline void Base_matrix_with_row_access<Master_matrix>::zero_cell(index columnIndex, index rowIndex)
{
	if constexpr (swap_opt::isActive_){
		matrix_.at(columnIndex).clear(swap_opt::indexToRow_.at(rowIndex));
	} else {
		matrix_.at(columnIndex).clear(rowIndex);
	}
}

template<class Master_matrix>
inline void Base_matrix_with_row_access<Master_matrix>::zero_column(index columnIndex)
{
	matrix_.at(columnIndex).clear();
}

template<class Master_matrix>
inline bool Base_matrix_with_row_access<Master_matrix>::is_zero_cell(index columnIndex, index rowIndex) const
{
	if constexpr (swap_opt::isActive_){
		return !(matrix_.at(columnIndex).contains(swap_opt::indexToRow_.at(rowIndex)));
	} else {
		return !(matrix_.at(columnIndex).contains(rowIndex));
	}
}

template<class Master_matrix>
inline bool Base_matrix_with_row_access<Master_matrix>::is_zero_column(index columnIndex)
{
	return matrix_.at(columnIndex).is_empty();
}

template<class Master_matrix>
inline Base_matrix_with_row_access<Master_matrix> &Base_matrix_with_row_access<Master_matrix>::operator=(Base_matrix_with_row_access other)
{
	swap_opt::operator=(other);
	pair_opt::operator=(other);
	std::swap(matrix_, other.matrix_);
	std::swap(maxDim_, other.maxDim_);
	std::swap(nextInsertIndex_, other.nextInsertIndex_);
	return *this;
}

template<class Master_matrix>
inline void Base_matrix_with_row_access<Master_matrix>::print()
{
	std::cout << "Base_matrix:\n";
	for (unsigned int i = 0; i < nextInsertIndex_; ++i){
		const Column_type& col = matrix_.at(i);
		for (auto e : col.get_content(nextInsertIndex_)){
			if (e == 0) std::cout << "-\n";
			else std::cout << e << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";
}

template<class Friend_master_matrix>
void swap(Base_matrix_with_row_access<Friend_master_matrix>& matrix1,
		  Base_matrix_with_row_access<Friend_master_matrix>& matrix2)
{
	std::swap(static_cast<typename Friend_master_matrix::Base_swap_option>(matrix1),
			  static_cast<typename Friend_master_matrix::Base_swap_option>(matrix2));
	std::swap(static_cast<typename Friend_master_matrix::Base_pairing_option>(matrix1),
			  static_cast<typename Friend_master_matrix::Base_pairing_option>(matrix2));
	std::swap(matrix1.matrix_, matrix2.matrix_);
	std::swap(matrix1.maxDim_, matrix2.maxDim_);
	std::swap(matrix1.nextInsertIndex_, matrix2.nextInsertIndex_);
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // BASE_MATRIX_0010_H

/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef BASE_MATRIX_0000_H
#define BASE_MATRIX_0000_H

#include <iostream>
#include <vector>

#include "../utilities/utilities.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class Base_matrix
		: public Master_matrix::Base_swap_option,
		  public Master_matrix::Base_pairing_option
{
public:
	using Column_type = typename Master_matrix::Column_type;
	using boundary_type = typename Master_matrix::boundary_type;
	using Row_type = void;

	Base_matrix();
	template<class Boundary_type = boundary_type>
	Base_matrix(const std::vector<Boundary_type>& orderedBoundaries);
	Base_matrix(unsigned int numberOfColumns);
	Base_matrix(const Base_matrix& matrixToCopy);
	Base_matrix(Base_matrix&& other) noexcept;
//	~Base_matrix(){std::cout << "base matrix destr: " << maxDim_ << ", " << matrix_.size() << "\n";};

	template<class Boundary_type = boundary_type>
	void insert_boundary(const Boundary_type& boundary);
	Column_type& get_column(index columnIndex);
	const Column_type& get_column(index columnIndex) const;
	Row_type get_row(index rowIndex) const;
	void erase_last();

	dimension_type get_max_dimension() const;
	unsigned int get_number_of_columns() const;

	dimension_type get_column_dimension(index columnIndex) const;

	void add_to(index sourceColumnIndex, index targetColumnIndex);

	void zero_cell(index columnIndex, index rowIndex);
	void zero_column(index columnIndex);
	bool is_zero_cell(index columnIndex, index rowIndex) const;
	bool is_zero_column(index columnIndex);

	index get_column_with_pivot(index simplexIndex) const;
	index get_pivot(index columnIndex);

	Base_matrix& operator=(Base_matrix other);
	friend void swap(Base_matrix& matrix1, Base_matrix& matrix2){
		swap(static_cast<typename Master_matrix::Base_swap_option&>(matrix1),
			 static_cast<typename Master_matrix::Base_swap_option&>(matrix2));
		swap(static_cast<typename Master_matrix::Base_pairing_option&>(matrix1),
			 static_cast<typename Master_matrix::Base_pairing_option&>(matrix2));
		matrix1.matrix_.swap(matrix2.matrix_);
		std::swap(matrix1.maxDim_, matrix2.maxDim_);
		std::swap(matrix1.nextInsertIndex_, matrix2.nextInsertIndex_);
	}

	void print();  //for debug

private:
	using swap_opt = typename Master_matrix::Base_swap_option;
	using pair_opt = typename Master_matrix::Base_pairing_option;

	typename Master_matrix::column_container_type matrix_;
	dimension_type maxDim_;
	index nextInsertIndex_;
};

template<class Master_matrix>
inline Base_matrix<Master_matrix>::Base_matrix()
	: Master_matrix::Base_swap_option(matrix_),
	  Master_matrix::Base_pairing_option(matrix_, maxDim_),
	  maxDim_(-1),
	  nextInsertIndex_(0)
{}

template<class Master_matrix>
template<class Boundary_type>
inline Base_matrix<Master_matrix>::Base_matrix(const std::vector<Boundary_type> &orderedBoundaries)
	: Master_matrix::Base_swap_option(matrix_, orderedBoundaries.size()),
	  Master_matrix::Base_pairing_option(matrix_, maxDim_),
	  matrix_(orderedBoundaries.size()),
	  maxDim_(0),
	  nextInsertIndex_(orderedBoundaries.size())
{
	for (unsigned int i = 0; i < orderedBoundaries.size(); i++){
		const Boundary_type& b = orderedBoundaries[i];
		matrix_[i] = Column_type(b);
		if (maxDim_ < matrix_[i].get_dimension()) maxDim_ = matrix_[i].get_dimension();
	}
}

template<class Master_matrix>
inline Base_matrix<Master_matrix>::Base_matrix(unsigned int numberOfColumns)
	: Master_matrix::Base_swap_option(matrix_, numberOfColumns),
	  Master_matrix::Base_pairing_option(matrix_, maxDim_),
	  matrix_(numberOfColumns),
	  maxDim_(-1),
	  nextInsertIndex_(0)
{}

template<class Master_matrix>
inline Base_matrix<Master_matrix>::Base_matrix(const Base_matrix &matrixToCopy)
	: Master_matrix::Base_swap_option(matrixToCopy),
	  Master_matrix::Base_pairing_option(matrixToCopy),
	  matrix_(matrixToCopy.matrix_),
	  maxDim_(matrixToCopy.maxDim_),
	  nextInsertIndex_(matrixToCopy.nextInsertIndex_)
{}

template<class Master_matrix>
inline Base_matrix<Master_matrix>::Base_matrix(Base_matrix &&other) noexcept
	: Master_matrix::Base_swap_option(std::move(other)),
	  Master_matrix::Base_pairing_option(std::move(other)),
	  matrix_(std::move(other.matrix_)),
	  maxDim_(std::exchange(other.maxDim_, -1)),
	  nextInsertIndex_(std::exchange(other.nextInsertIndex_, 0))
{}

template<class Master_matrix>
template<class Boundary_type>
inline void Base_matrix<Master_matrix>::insert_boundary(const Boundary_type &boundary)
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

	matrix_[nextInsertIndex_++] = Column_type(boundary);
	if (maxDim_ < boundary.size() - 1) maxDim_ = boundary.size() - 1;
}

template<class Master_matrix>
inline typename Base_matrix<Master_matrix>::Column_type &Base_matrix<Master_matrix>::get_column(index columnIndex)
{
	if constexpr (swap_opt::isActive_){
		if (swap_opt::rowSwapped_) swap_opt::_orderRows();
	}

	return matrix_[columnIndex];
}

template<class Master_matrix>
inline const typename Base_matrix<Master_matrix>::Column_type &Base_matrix<Master_matrix>::get_column(index columnIndex) const
{
	if constexpr (swap_opt::isActive_){
		if (swap_opt::rowSwapped_) swap_opt::_orderRows();
	}

	return matrix_[columnIndex];
}

template<class Master_matrix>
inline typename Base_matrix<Master_matrix>::Row_type Base_matrix<Master_matrix>::get_row(index rowIndex) const
{
	static_assert(static_cast<int>(Master_matrix::Field_type::get_characteristic()) == -1,
			"'get_row' is not implemented for the chosen options.");
}

template<class Master_matrix>
inline void Base_matrix<Master_matrix>::erase_last()
{
	static_assert(static_cast<int>(Master_matrix::Field_type::get_characteristic()) == -1,
			"'erase_last' is not implemented for the chosen options.");
}

template<class Master_matrix>
inline dimension_type Base_matrix<Master_matrix>::get_max_dimension() const
{
	return maxDim_;
}

template<class Master_matrix>
inline unsigned int Base_matrix<Master_matrix>::get_number_of_columns() const
{
	return nextInsertIndex_;
}

template<class Master_matrix>
inline dimension_type Base_matrix<Master_matrix>::get_column_dimension(index columnIndex) const
{
	return matrix_[columnIndex].get_dimension();
}

template<class Master_matrix>
inline void Base_matrix<Master_matrix>::add_to(index sourceColumnIndex, index targetColumnIndex)
{
	matrix_[targetColumnIndex] += matrix_[sourceColumnIndex];
}

template<class Master_matrix>
inline void Base_matrix<Master_matrix>::zero_cell(index columnIndex, index rowIndex)
{
	if constexpr (swap_opt::isActive_){
		matrix_[columnIndex].clear(swap_opt::indexToRow_[rowIndex]);
	} else {
		matrix_[columnIndex].clear(rowIndex);
	}
}

template<class Master_matrix>
inline void Base_matrix<Master_matrix>::zero_column(index columnIndex)
{
	matrix_[columnIndex].clear();
}

template<class Master_matrix>
inline bool Base_matrix<Master_matrix>::is_zero_cell(index columnIndex, index rowIndex) const
{
	if constexpr (swap_opt::isActive_){
		return !(matrix_[columnIndex].is_non_zero(swap_opt::indexToRow_[rowIndex]));
	} else {
		return !(matrix_[columnIndex].is_non_zero(rowIndex));
	}
}

template<class Master_matrix>
inline bool Base_matrix<Master_matrix>::is_zero_column(index columnIndex)
{
	return matrix_[columnIndex].is_empty();
}

template<class Master_matrix>
inline index Base_matrix<Master_matrix>::get_column_with_pivot(index simplexIndex) const
{
	static_assert(static_cast<int>(Master_matrix::Field_type::get_characteristic()) == -1,
			"'get_column_with_pivot' is not implemented for the chosen options.");
}

template<class Master_matrix>
inline index Base_matrix<Master_matrix>::get_pivot(index columnIndex)
{
	return matrix_[columnIndex].get_pivot();
}

template<class Master_matrix>
inline Base_matrix<Master_matrix> &Base_matrix<Master_matrix>::operator=(Base_matrix other)
{
	swap_opt::operator=(other);
	pair_opt::operator=(other);
	matrix_.swap(other.matrix_);
	std::swap(maxDim_, other.maxDim_);
	std::swap(nextInsertIndex_, other.nextInsertIndex_);
	return *this;
}

template<class Master_matrix>
inline void Base_matrix<Master_matrix>::print()
{
	std::cout << "Base_matrix:\n";
	for (unsigned int i = 0; i < nextInsertIndex_; ++i){
		Column_type& col = matrix_[i];
		for (auto e : col.get_content(nextInsertIndex_)){
			if (e == 0u) std::cout << "- ";
			else std::cout << e << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // BASE_MATRIX_0000_H

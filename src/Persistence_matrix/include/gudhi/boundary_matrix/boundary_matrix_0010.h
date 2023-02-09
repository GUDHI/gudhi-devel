/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef BOUNDARY_MATRIX_0010_H
#define BOUNDARY_MATRIX_0010_H

#include <iostream>
#include <vector>
#include <set>
#include <algorithm>

#include "../utilities/utilities.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class Boundary_matrix_with_row_access
		: public Master_matrix::Base_swap_option,
		  public Master_matrix::Base_pairing_option
{
public:
	using Field_element_type = typename Master_matrix::Field_type;
	using Column_type = typename Master_matrix::Column_type;
	using boundary_type = typename Master_matrix::boundary_type;
	using Row_type = typename Master_matrix::Row_type;

	Boundary_matrix_with_row_access();
	template<class Boundary_type = boundary_type>
	Boundary_matrix_with_row_access(const std::vector<Boundary_type>& orderedBoundaries);
	Boundary_matrix_with_row_access(unsigned int numberOfColumns);
	Boundary_matrix_with_row_access(const Boundary_matrix_with_row_access& matrixToCopy);
	Boundary_matrix_with_row_access(Boundary_matrix_with_row_access&& other) noexcept;

	template<class Boundary_type = boundary_type>
	void insert_boundary(const Boundary_type& boundary);
	Column_type& get_column(index columnIndex);
	const Column_type& get_column(index columnIndex) const;
	//get_row(rowIndex) --> simplex ID (=/= columnIndex)
	Row_type& get_row(index rowIndex);
	const Row_type& get_row(index rowIndex) const;
	void erase_last();

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

	Boundary_matrix_with_row_access& operator=(const Boundary_matrix_with_row_access& other);
	friend void swap(Boundary_matrix_with_row_access& matrix1, Boundary_matrix_with_row_access& matrix2){
		swap(static_cast<typename Master_matrix::Base_swap_option&>(matrix1),
			 static_cast<typename Master_matrix::Base_swap_option&>(matrix2));
		swap(static_cast<typename Master_matrix::Base_pairing_option&>(matrix1),
			 static_cast<typename Master_matrix::Base_pairing_option&>(matrix2));
		matrix1.rows_.swap(matrix2.rows_);
		matrix1.matrix_.swap(matrix2.matrix_);
		std::swap(matrix1.maxDim_, matrix2.maxDim_);
		std::swap(matrix1.nextInsertIndex_, matrix2.nextInsertIndex_);
		for (Column_type& col : matrix1.matrix_){
			col.set_rows(&matrix1.rows_);
		}
		for (Column_type& col : matrix2.matrix_){
			col.set_rows(&matrix2.rows_);
		}
	}

	void print();  //for debug

private:
	using swap_opt = typename Master_matrix::Base_swap_option;
	using pair_opt = typename Master_matrix::Base_pairing_option;
	using matrix_type = typename Master_matrix::column_container_type;
	using rows_type = typename Master_matrix::row_container_type;
	using cell_rep_type = typename std::conditional<
								Master_matrix::Option_list::is_z2,
								index,
								std::pair<index,typename Master_matrix::Field_type>
							>::type;

	rows_type rows_;	//has to be destroyed after matrix_
	matrix_type matrix_;
	dimension_type maxDim_;
	index nextInsertIndex_;
};

template<class Master_matrix>
inline Boundary_matrix_with_row_access<Master_matrix>::Boundary_matrix_with_row_access()
	: Master_matrix::Base_swap_option(matrix_),
	  Master_matrix::Base_pairing_option(matrix_, maxDim_),
	  maxDim_(-1),
	  nextInsertIndex_(0)
{}

template<class Master_matrix>
template<class Boundary_type>
inline Boundary_matrix_with_row_access<Master_matrix>::Boundary_matrix_with_row_access(const std::vector<Boundary_type> &orderedBoundaries)
	: Master_matrix::Base_swap_option(matrix_, orderedBoundaries.size()),
	  Master_matrix::Base_pairing_option(matrix_, maxDim_),
	  rows_(orderedBoundaries.size()),
	  maxDim_(-1),
	  nextInsertIndex_(orderedBoundaries.size())
{
	matrix_.reserve(orderedBoundaries.size());
	for (unsigned int i = 0; i < orderedBoundaries.size(); i++){
		const Boundary_type& b = orderedBoundaries[i];
		matrix_.emplace_back(i, b, rows_);
		if (maxDim_ < matrix_.back().get_dimension()) maxDim_ = matrix_.back().get_dimension();
	}
}

template<class Master_matrix>
inline Boundary_matrix_with_row_access<Master_matrix>::Boundary_matrix_with_row_access(unsigned int numberOfColumns)
	: Master_matrix::Base_swap_option(matrix_, numberOfColumns),
	  Master_matrix::Base_pairing_option(matrix_, maxDim_),
	  rows_(numberOfColumns),
	  maxDim_(-1),
	  nextInsertIndex_(0)
{
	matrix_.reserve(numberOfColumns);
}

template<class Master_matrix>
inline Boundary_matrix_with_row_access<Master_matrix>::Boundary_matrix_with_row_access(const Boundary_matrix_with_row_access &matrixToCopy)
	: Master_matrix::Base_swap_option(matrixToCopy),
	  Master_matrix::Base_pairing_option(matrixToCopy),
	  rows_(matrixToCopy.rows_.size()),
	  maxDim_(matrixToCopy.maxDim_),
	  nextInsertIndex_(matrixToCopy.nextInsertIndex_)
{
	if constexpr (swap_opt::isActive_)
		swap_opt::matrix_ = &matrix_;
	if constexpr (pair_opt::isActive_){
		pair_opt::matrix_ = &matrix_;
		pair_opt::maxDim_ = &maxDim_;
	}
	matrix_.reserve(matrixToCopy.matrix_.size());
	for (const Column_type& col : matrixToCopy.matrix_){
		std::vector<cell_rep_type> tmp(col.begin(), col.end());
		matrix_.emplace_back(col.get_column_index(), tmp, col.get_dimension(), rows_);
	}
}

template<class Master_matrix>
inline Boundary_matrix_with_row_access<Master_matrix>::Boundary_matrix_with_row_access(Boundary_matrix_with_row_access &&other) noexcept
	: Master_matrix::Base_swap_option(std::move(other)),
	  Master_matrix::Base_pairing_option(std::move(other)),
	  rows_(std::move(other.rows_)),
	  matrix_(std::move(other.matrix_)),
	  maxDim_(std::exchange(other.maxDim_, -1)),
	  nextInsertIndex_(std::exchange(other.nextInsertIndex_, 0))
{
	if constexpr (swap_opt::isActive_)
		swap_opt::matrix_ = &matrix_;
	if constexpr (pair_opt::isActive_){
		pair_opt::matrix_ = &matrix_;
		pair_opt::maxDim_ = &maxDim_;
	}
	for (Column_type& col : matrix_){
		col.set_rows(&rows_);
	}
}

template<class Master_matrix>
template<class Boundary_type>
inline void Boundary_matrix_with_row_access<Master_matrix>::insert_boundary(const Boundary_type &boundary)
{
	if constexpr (swap_opt::isActive_){
		if (swap_opt::rowSwapped_) swap_opt::_orderRows();
	}

	if constexpr (swap_opt::isActive_){
		if (swap_opt::indexToRow_.size() <= nextInsertIndex_){
			swap_opt::indexToRow_.push_back(nextInsertIndex_);
			swap_opt::rowToIndex_.push_back(nextInsertIndex_);
		}
	}

	if (rows_.size() <= nextInsertIndex_) rows_.resize(nextInsertIndex_ + 1);
	matrix_.emplace_back(nextInsertIndex_++, boundary, rows_);
	if (maxDim_ < matrix_.back().get_dimension()) maxDim_ = matrix_.back().get_dimension();
}

template<class Master_matrix>
inline typename Boundary_matrix_with_row_access<Master_matrix>::Column_type &Boundary_matrix_with_row_access<Master_matrix>::get_column(index columnIndex)
{
	if constexpr (swap_opt::isActive_){
		if (swap_opt::rowSwapped_) swap_opt::_orderRows();
	}

	return matrix_[columnIndex];
}

template<class Master_matrix>
inline const typename Boundary_matrix_with_row_access<Master_matrix>::Column_type &Boundary_matrix_with_row_access<Master_matrix>::get_column(index columnIndex) const
{
	if constexpr (swap_opt::isActive_){
		if (swap_opt::rowSwapped_) swap_opt::_orderRows();
	}

	return matrix_[columnIndex];
}

template<class Master_matrix>
inline typename Boundary_matrix_with_row_access<Master_matrix>::Row_type& Boundary_matrix_with_row_access<Master_matrix>::get_row(index rowIndex)
{
	if constexpr (swap_opt::isActive_){
		if (swap_opt::rowSwapped_) swap_opt::_orderRows();
	}
	return rows_[rowIndex];
}

template<class Master_matrix>
inline const typename Boundary_matrix_with_row_access<Master_matrix>::Row_type& Boundary_matrix_with_row_access<Master_matrix>::get_row(index rowIndex) const
{
	if constexpr (swap_opt::isActive_){
		if (swap_opt::rowSwapped_) swap_opt::_orderRows();
	}
	return rows_[rowIndex];
}

template<class Master_matrix>
inline void Boundary_matrix_with_row_access<Master_matrix>::erase_last()
{
	static_assert(Master_matrix::Option_list::has_row_access,
			"'erase_last' is not implemented for the chosen options.");
}

template<class Master_matrix>
inline dimension_type Boundary_matrix_with_row_access<Master_matrix>::get_max_dimension() const
{
	return maxDim_;
}

template<class Master_matrix>
inline unsigned int Boundary_matrix_with_row_access<Master_matrix>::get_number_of_columns() const
{
	return nextInsertIndex_;
}

template<class Master_matrix>
inline dimension_type Boundary_matrix_with_row_access<Master_matrix>::get_column_dimension(index columnIndex) const
{
	return matrix_[columnIndex].get_dimension();
}

template<class Master_matrix>
inline void Boundary_matrix_with_row_access<Master_matrix>::add_to(index sourceColumnIndex, index targetColumnIndex)
{
	matrix_[targetColumnIndex] += matrix_[sourceColumnIndex];
}

template<class Master_matrix>
inline void Boundary_matrix_with_row_access<Master_matrix>::add_to(const Column_type& sourceColumn, index targetColumnIndex)
{
	matrix_[targetColumnIndex] += sourceColumn;
}

template<class Master_matrix>
inline void Boundary_matrix_with_row_access<Master_matrix>::add_to(const Column_type& sourceColumn, const Field_element_type& coefficient, index targetColumnIndex)
{
	matrix_[targetColumnIndex].multiply_and_add(coefficient, sourceColumn);
}

template<class Master_matrix>
inline void Boundary_matrix_with_row_access<Master_matrix>::add_to(const Field_element_type& coefficient, const Column_type& sourceColumn, index targetColumnIndex)
{
	matrix_[targetColumnIndex].multiply_and_add(sourceColumn, coefficient);
}

template<class Master_matrix>
inline void Boundary_matrix_with_row_access<Master_matrix>::zero_cell(index columnIndex, index rowIndex)
{
	if constexpr (swap_opt::isActive_){
		matrix_[columnIndex].clear(swap_opt::indexToRow_[rowIndex]);
	} else {
		matrix_[columnIndex].clear(rowIndex);
	}
}

template<class Master_matrix>
inline void Boundary_matrix_with_row_access<Master_matrix>::zero_column(index columnIndex)
{
	matrix_[columnIndex].clear();
}

template<class Master_matrix>
inline bool Boundary_matrix_with_row_access<Master_matrix>::is_zero_cell(index columnIndex, index rowIndex) const
{
	if constexpr (swap_opt::isActive_){
		return !(matrix_[columnIndex].is_non_zero(swap_opt::indexToRow_[rowIndex]));
	} else {
		return !(matrix_[columnIndex].is_non_zero(rowIndex));
	}
}

template<class Master_matrix>
inline bool Boundary_matrix_with_row_access<Master_matrix>::is_zero_column(index columnIndex)
{
	return matrix_[columnIndex].is_empty();
}

template<class Master_matrix>
inline index Boundary_matrix_with_row_access<Master_matrix>::get_column_with_pivot(index simplexIndex) const
{
	static_assert(Master_matrix::Option_list::has_row_access,
			"'get_column_with_pivot' is not implemented for the chosen options.");
}

template<class Master_matrix>
inline index Boundary_matrix_with_row_access<Master_matrix>::get_pivot(index columnIndex)
{
	return matrix_[columnIndex].get_pivot();
}

template<class Master_matrix>
inline Boundary_matrix_with_row_access<Master_matrix> &
Boundary_matrix_with_row_access<Master_matrix>::operator=(const Boundary_matrix_with_row_access &other)
{
	swap_opt::operator=(other);
	pair_opt::operator=(other);
	rows_.resize(other.rows_.size());
	maxDim_ = other.maxDim_;
	nextInsertIndex_ = other.nextInsertIndex_;
	matrix_.reserve(other.matrix_.size());
	for (const Column_type& col : other.matrix_){
		std::vector<cell_rep_type> tmp(col.begin(), col.end());
		matrix_.emplace_back(col.get_column_index(), tmp, col.get_dimension(), rows_);
	}
	return *this;
}

template<class Master_matrix>
inline void Boundary_matrix_with_row_access<Master_matrix>::print()
{
	std::cout << "Base_matrix_with_row_access:\n";
	for (unsigned int i = 0; i < nextInsertIndex_; ++i){
		Column_type& col = matrix_[i];
		for (auto e : col.get_content(nextInsertIndex_)){
			if (e == 0u) std::cout << "- ";
			else std::cout << e << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";
	std::cout << "Row Matrix:\n";
	for (unsigned int i = 0; i < nextInsertIndex_; ++i){
		const Row_type& row = rows_[i];
		for (const auto &cell : row){
			std::cout << cell.get_column_index() << " ";
		}
		std::cout << "(" << i << ")\n";
	}
	std::cout << "\n";
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // BOUNDARY_MATRIX_0010_H

/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CHAIN_MATRIX_0010_H
#define CHAIN_MATRIX_0010_H

#include <iostream>
#include <set>

#include "../utilities/utilities.h"
#include "custom_chain_vine_swap.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class Chain_matrix_with_row_access
		: public Master_matrix::Chain_pairing_option,
		  public Master_matrix::Chain_vine_swap_option,
		  public Master_matrix::Chain_representative_cycles_option
{
public:
	using Field_element_type = typename Master_matrix::Field_type;
	using Column_type = typename Master_matrix::Column_type;
	using Row_type = typename Master_matrix::Row_type;
	using Cell = typename Master_matrix::Cell_type;
	using boundary_type = typename Master_matrix::boundary_type;
	using cell_rep_type = typename Master_matrix::cell_rep_type;

	Chain_matrix_with_row_access();
	template<class Boundary_type = boundary_type>
	Chain_matrix_with_row_access(const std::vector<Boundary_type>& orderedBoundaries);
	Chain_matrix_with_row_access(unsigned int numberOfColumns);
	Chain_matrix_with_row_access(
		std::function<bool(index,index)> birthComparator, 
		std::function<bool(index,index)> deathComparator = _no_G_death_comparator);
	template<class Boundary_type = boundary_type>
	Chain_matrix_with_row_access(
		const std::vector<Boundary_type>& orderedBoundaries,
		std::function<bool(index,index)> birthComparator, 
		std::function<bool(index,index)> deathComparator = _no_G_death_comparator);
	Chain_matrix_with_row_access(
		unsigned int numberOfColumns,
		std::function<bool(index,index)> birthComparator, 
		std::function<bool(index,index)> deathComparator = _no_G_death_comparator);
	Chain_matrix_with_row_access(const Chain_matrix_with_row_access& matrixToCopy);
	Chain_matrix_with_row_access(Chain_matrix_with_row_access&& other) noexcept;

	//new simplex = new ID even if the same simplex was already inserted and then removed, ie., an ID cannot come back.
	template<class Boundary_type = boundary_type>
	std::vector<cell_rep_type> insert_boundary(const Boundary_type& boundary);
	template<class Boundary_type = boundary_type>
	std::vector<cell_rep_type> insert_boundary(index simplexIndex, const Boundary_type& boundary);
	Column_type& get_column(index columnIndex);
	const Column_type& get_column(index columnIndex) const;
	Row_type& get_row(index rowIndex);
	const Row_type& get_row(index rowIndex) const;
	void remove_maximal_simplex(index simplexIndex);

	dimension_type get_max_dimension() const;
	unsigned int get_number_of_columns() const;

	dimension_type get_column_dimension(index columnIndex) const;

	void add_to(index sourceColumnIndex, index targetColumnIndex);
	void add_to(Column_type& sourceColumn, index targetColumnIndex);
	void add_to(Column_type& sourceColumn, const Field_element_type& coefficient, index targetColumnIndex);
	void add_to(const Field_element_type& coefficient, Column_type& sourceColumn, index targetColumnIndex);

	void zero_cell(index columnIndex, index rowIndex);
	void zero_column(index columnIndex);
	bool is_zero_cell(index columnIndex, index rowIndex) const;
	bool is_zero_column(index columnIndex) const;

	index get_column_with_pivot(index simplexIndex) const;
	index get_pivot(index columnIndex) const;

	Chain_matrix_with_row_access& operator=(const Chain_matrix_with_row_access& other);
	friend void swap(Chain_matrix_with_row_access& matrix1,
					 Chain_matrix_with_row_access& matrix2){
		swap(static_cast<typename Master_matrix::Chain_pairing_option&>(matrix1),
			 static_cast<typename Master_matrix::Chain_pairing_option&>(matrix2));
		swap(static_cast<typename Master_matrix::Chain_vine_swap_option&>(matrix1),
			 static_cast<typename Master_matrix::Chain_vine_swap_option&>(matrix2));
		swap(static_cast<typename Master_matrix::Chain_representative_cycles_option&>(matrix1),
			 static_cast<typename Master_matrix::Chain_representative_cycles_option&>(matrix2));
		matrix1.matrix_.swap(matrix2.matrix_);
		matrix1.rows_.swap(matrix2.rows_);
		matrix1.pivotToColumnIndex_.swap(matrix2.pivotToColumnIndex_);
		std::swap(matrix1.nextInsertIndex_, matrix2.nextInsertIndex_);
		std::swap(matrix1.maxDim_, matrix2.maxDim_);
		for (Column_type& col : matrix1.matrix_){
			col.set_rows(&matrix1.rows_);
		}
		for (Column_type& col : matrix2.matrix_){
			col.set_rows(&matrix2.rows_);
		}
	}

	void print() const;  //for debug

private:
	using swap_opt = typename Master_matrix::Chain_vine_swap_option;
	using pair_opt = typename Master_matrix::Chain_pairing_option;
	using rep_opt = typename Master_matrix::Chain_representative_cycles_option;
	using matrix_type = typename Master_matrix::column_container_type;
	using rows_type = typename Master_matrix::row_container_type;
	using dictionnary_type = typename Master_matrix::template dictionnary_type<index>;
	using barcode_type = typename Master_matrix::barcode_type;
	using bar_dictionnary_type = typename Master_matrix::bar_dictionnary_type;
	using tmp_column_type = typename std::conditional<
								Master_matrix::Option_list::is_z2,
								std::set<index>,
								std::set<std::pair<index,Field_element_type>,CellPairComparator<Field_element_type> >
							>::type;

	rows_type rows_;	//has to be destroyed after matrix_
	matrix_type matrix_;
	dictionnary_type pivotToColumnIndex_;
	index nextInsertIndex_;
	dimension_type maxDim_;

	template<class Boundary_type>
	std::vector<cell_rep_type> _reduce_boundary(index simplexIndex, const Boundary_type& boundary);
	void _reduce_by_G(tmp_column_type& column,
					  std::vector<cell_rep_type>& chainsInH,
					  index currentPivot);
	void _reduce_by_F(tmp_column_type& column,
					  std::vector<cell_rep_type>& chainsInF,
					  index currentPivot);
	void _build_from_H(index simplexIndex, 
					   tmp_column_type& column,
					   std::vector<cell_rep_type>& chainsInH);
	void _update_largest_death_in_F(const std::vector<cell_rep_type>& chainsInF);
	void _insert_chain(const tmp_column_type& column, dimension_type dimension);
	void _insert_chain(const tmp_column_type& column, dimension_type dimension, index pair);
	void _add_to(const Column_type& column, tmp_column_type& set, unsigned int coef = 1u);

	static constexpr bool _barcode_option_is_active();
	constexpr barcode_type& _barcode();
	constexpr bar_dictionnary_type& _indexToBar();
};

template<class Master_matrix>
inline Chain_matrix_with_row_access<Master_matrix>::Chain_matrix_with_row_access()
	: Master_matrix::Chain_pairing_option(),
	  Master_matrix::Chain_vine_swap_option(matrix_),
	  Master_matrix::Chain_representative_cycles_option(matrix_, pivotToColumnIndex_),
	  nextInsertIndex_(0),
	  maxDim_(-1)
{}

template<class Master_matrix>
template<class Boundary_type>
inline Chain_matrix_with_row_access<Master_matrix>::Chain_matrix_with_row_access(
		const std::vector<Boundary_type> &orderedBoundaries)
	: Master_matrix::Chain_pairing_option(),
	  Master_matrix::Chain_vine_swap_option(matrix_),
	  Master_matrix::Chain_representative_cycles_option(matrix_, pivotToColumnIndex_),
//	  rows_(orderedBoundaries.size()),
	  pivotToColumnIndex_(orderedBoundaries.size(), -1),
	  nextInsertIndex_(0),
	  maxDim_(-1)
{
	matrix_.reserve(orderedBoundaries.size());
	if constexpr (!Master_matrix::Option_list::has_removable_rows){
		rows_.resize(orderedBoundaries.size());
	}
	for (const Boundary_type &b : orderedBoundaries){
		insert_boundary(b);
	}
}

template<class Master_matrix>
inline Chain_matrix_with_row_access<Master_matrix>::Chain_matrix_with_row_access(unsigned int numberOfColumns)
	: Master_matrix::Chain_pairing_option(),
	  Master_matrix::Chain_vine_swap_option(matrix_),
	  Master_matrix::Chain_representative_cycles_option(matrix_, pivotToColumnIndex_),
//	  rows_(numberOfColumns),
	  pivotToColumnIndex_(numberOfColumns, -1),
	  nextInsertIndex_(0),
	  maxDim_(-1)
{
	matrix_.reserve(numberOfColumns);
	if constexpr (!Master_matrix::Option_list::has_removable_rows){
		rows_.resize(numberOfColumns);
	}
}

template<class Master_matrix>
inline Chain_matrix_with_row_access<Master_matrix>::Chain_matrix_with_row_access(
		std::function<bool(index,index)> birthComparator, 
		std::function<bool(index,index)> deathComparator)
	: Master_matrix::Chain_pairing_option(),
	  Master_matrix::Chain_vine_swap_option(matrix_, birthComparator, deathComparator),
	  Master_matrix::Chain_representative_cycles_option(matrix_, pivotToColumnIndex_),
	  nextInsertIndex_(0),
	  maxDim_(-1)
{}

template<class Master_matrix>
template<class Boundary_type>
inline Chain_matrix_with_row_access<Master_matrix>::Chain_matrix_with_row_access(
		const std::vector<Boundary_type> &orderedBoundaries,
		std::function<bool(index,index)> birthComparator, 
		std::function<bool(index,index)> deathComparator)
	: Master_matrix::Chain_pairing_option(),
	  Master_matrix::Chain_vine_swap_option(matrix_, birthComparator, deathComparator),
	  Master_matrix::Chain_representative_cycles_option(matrix_, pivotToColumnIndex_),
//	  rows_(orderedBoundaries.size()),
	  matrix_(orderedBoundaries.size()),
	  pivotToColumnIndex_(orderedBoundaries.size(), -1),
	  nextInsertIndex_(0),
	  maxDim_(-1)
{
	matrix_.reserve(orderedBoundaries.size());
	if constexpr (!Master_matrix::Option_list::has_removable_rows){
		rows_.resize(orderedBoundaries.size());
	}
	for (const Boundary_type &b : orderedBoundaries){
		insert_boundary(b);
	}
}

template<class Master_matrix>
inline Chain_matrix_with_row_access<Master_matrix>::Chain_matrix_with_row_access(
		unsigned int numberOfColumns,
		std::function<bool(index,index)> birthComparator, 
		std::function<bool(index,index)> deathComparator)
	: Master_matrix::Chain_pairing_option(),
	  Master_matrix::Chain_vine_swap_option(matrix_, birthComparator, deathComparator),
	  Master_matrix::Chain_representative_cycles_option(matrix_, pivotToColumnIndex_),
//	  rows_(numberOfColumns),
	  pivotToColumnIndex_(numberOfColumns, -1),
	  nextInsertIndex_(0),
	  maxDim_(-1)
{
	matrix_.reserve(numberOfColumns);
	if constexpr (!Master_matrix::Option_list::has_removable_rows){
		rows_.resize(numberOfColumns);
	}
}

template<class Master_matrix>
inline Chain_matrix_with_row_access<Master_matrix>::Chain_matrix_with_row_access(
		const Chain_matrix_with_row_access &matrixToCopy)
	: Master_matrix::Chain_pairing_option(matrixToCopy),
	  Master_matrix::Chain_vine_swap_option(matrixToCopy),
	  Master_matrix::Chain_representative_cycles_option(matrixToCopy),
//	  rows_(matrixToCopy.rows_.size()),
	  pivotToColumnIndex_(matrixToCopy.pivotToColumnIndex_),
	  nextInsertIndex_(matrixToCopy.nextInsertIndex_),
	  maxDim_(matrixToCopy.maxDim_)
{
	if constexpr (!Master_matrix::Option_list::has_removable_rows){
		rows_.resize(matrixToCopy.rows_.size());
	}
	if constexpr (rep_opt::isActive_){
		rep_opt::matrix_ = &matrix_;
		rep_opt::pivotToPosition_ = &pivotToColumnIndex_;
	}
	if constexpr (swap_opt::isActive_){
		swap_opt::matrix_ = &matrix_;
	}
	matrix_.reserve(matrixToCopy.matrix_.size());
	for (const Column_type& col : matrixToCopy.matrix_){
		std::vector<cell_rep_type> tmp(col.begin(), col.end());
		matrix_.emplace_back(col.get_column_index(), tmp, col.get_dimension(), rows_, pivotToColumnIndex_);
	}
}

template<class Master_matrix>
inline Chain_matrix_with_row_access<Master_matrix>::Chain_matrix_with_row_access(
		Chain_matrix_with_row_access &&other) noexcept
	: Master_matrix::Chain_pairing_option(std::move(other)),
	  Master_matrix::Chain_vine_swap_option(std::move(other)),
	  Master_matrix::Chain_representative_cycles_option(std::move(other)),
	  rows_(std::move(other.rows_)),
	  matrix_(std::move(other.matrix_)),
	  pivotToColumnIndex_(std::move(other.pivotToColumnIndex_)),
	  nextInsertIndex_(std::exchange(other.nextInsertIndex_, 0)),
	  maxDim_(std::exchange(other.maxDim_, -1))
{
	if constexpr (rep_opt::isActive_){
		rep_opt::matrix_ = &matrix_;
		rep_opt::pivotToPosition_ = &pivotToColumnIndex_;
	}
	if constexpr (swap_opt::isActive_){
		swap_opt::matrix_ = &matrix_;
	}
	for (Column_type& col : matrix_){
		col.set_rows(&rows_);
	}
}

template<class Master_matrix>
template<class Boundary_type>
inline std::vector<typename Master_matrix::cell_rep_type> Chain_matrix_with_row_access<Master_matrix>::insert_boundary(
		const Boundary_type &boundary)
{
	if (pivotToColumnIndex_.size() <= nextInsertIndex_){
		pivotToColumnIndex_.resize(nextInsertIndex_ * 2, -1);
	}
	if constexpr (swap_opt::isActive_ && _barcode_option_is_active()){
		if (swap_opt::pivotToPosition_.size() <= nextInsertIndex_)
			swap_opt::pivotToPosition_.resize(pivotToColumnIndex_.size());
		swap_opt::pivotToPosition_[nextInsertIndex_] = nextInsertIndex_;
	}
	int dim = boundary.size() - 1;
	if (maxDim_ < dim) maxDim_ = dim;

	return _reduce_boundary(nextInsertIndex_, boundary);
}

template<class Master_matrix>
template<class Boundary_type>
inline std::vector<typename Master_matrix::cell_rep_type> Chain_matrix_with_row_access<Master_matrix>::insert_boundary(
		index simplexIndex, const Boundary_type &boundary)
{
	if (pivotToColumnIndex_.size() <= simplexIndex){
		pivotToColumnIndex_.resize(simplexIndex * 2, -1);
	}
	if constexpr (swap_opt::isActive_ && _barcode_option_is_active()){
		if (swap_opt::pivotToPosition_.size() <= simplexIndex)
			swap_opt::pivotToPosition_.resize(pivotToColumnIndex_.size());
		swap_opt::pivotToPosition_[simplexIndex] = simplexIndex;
	}
	int dim = boundary.size() - 1;
	if (maxDim_ < dim) maxDim_ = dim;

	return _reduce_boundary(simplexIndex, boundary);
}

template<class Master_matrix>
inline typename Chain_matrix_with_row_access<Master_matrix>::Column_type &
Chain_matrix_with_row_access<Master_matrix>::get_column(index columnIndex)
{
	return matrix_[columnIndex];
}

template<class Master_matrix>
inline const typename Chain_matrix_with_row_access<Master_matrix>::Column_type &
Chain_matrix_with_row_access<Master_matrix>::get_column(index columnIndex) const
{
	return matrix_[columnIndex];
}

template<class Master_matrix>
inline typename Chain_matrix_with_row_access<Master_matrix>::Row_type&
Chain_matrix_with_row_access<Master_matrix>::get_row(index rowIndex)
{
	return rows_[rowIndex];
}

template<class Master_matrix>
inline const typename Chain_matrix_with_row_access<Master_matrix>::Row_type&
Chain_matrix_with_row_access<Master_matrix>::get_row(index rowIndex) const
{
	return rows_[rowIndex];
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access<Master_matrix>::remove_maximal_simplex([[maybe_unused]] index simplexIndex)
{
	static_assert(Master_matrix::Option_list::has_removable_columns,
			"'remove_maximal_simplex' is not implemented for the chosen options.");
}

template<class Master_matrix>
inline dimension_type Chain_matrix_with_row_access<Master_matrix>::get_max_dimension() const
{
	return maxDim_;
}

template<class Master_matrix>
inline unsigned int Chain_matrix_with_row_access<Master_matrix>::get_number_of_columns() const
{
	return nextInsertIndex_;
}

template<class Master_matrix>
inline dimension_type Chain_matrix_with_row_access<Master_matrix>::get_column_dimension(index columnIndex) const
{
	return matrix_[columnIndex].get_dimension();
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access<Master_matrix>::add_to(index sourceColumnIndex, index targetColumnIndex)
{
	matrix_[targetColumnIndex] += matrix_[sourceColumnIndex];
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access<Master_matrix>::add_to(Column_type& sourceColumn, index targetColumnIndex)
{
	matrix_[targetColumnIndex] += sourceColumn;
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access<Master_matrix>::add_to(Column_type& sourceColumn, const Field_element_type& coefficient, index targetColumnIndex)
{
	matrix_[targetColumnIndex].multiply_and_add(coefficient, sourceColumn);
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access<Master_matrix>::add_to(const Field_element_type& coefficient, Column_type& sourceColumn, index targetColumnIndex)
{
	matrix_[targetColumnIndex].multiply_and_add(sourceColumn, coefficient);
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access<Master_matrix>::zero_cell(index columnIndex, index rowIndex)
{
	static_assert(Master_matrix::Option_list::has_removable_columns,
			"'zero_cell' is not implemented for the chosen options.");
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access<Master_matrix>::zero_column(index columnIndex)
{
	static_assert(Master_matrix::Option_list::has_removable_columns,
			"'zero_column' is not implemented for the chosen options.");
}

template<class Master_matrix>
inline bool Chain_matrix_with_row_access<Master_matrix>::is_zero_cell(index columnIndex, index rowIndex) const
{
	return !matrix_[columnIndex].is_non_zero(rowIndex);
}

template<class Master_matrix>
inline bool Chain_matrix_with_row_access<Master_matrix>::is_zero_column(index columnIndex) const
{
	return matrix_[columnIndex].is_empty();
}

template<class Master_matrix>
inline index Chain_matrix_with_row_access<Master_matrix>::get_column_with_pivot(index simplexIndex) const
{
	return pivotToColumnIndex_[simplexIndex];
}

template<class Master_matrix>
inline index Chain_matrix_with_row_access<Master_matrix>::get_pivot(index columnIndex) const
{
	return matrix_[columnIndex].get_pivot();
}

template<class Master_matrix>
inline Chain_matrix_with_row_access<Master_matrix> &
Chain_matrix_with_row_access<Master_matrix>::operator=(const Chain_matrix_with_row_access &other)
{
	swap_opt::operator=(other);
	pair_opt::operator=(other);
	rep_opt::operator=(other);
	rows_.resize(other.rows_.size());
	pivotToColumnIndex_ = other.pivotToColumnIndex_;
	nextInsertIndex_ = other.nextInsertIndex_;
	maxDim_ = other.maxDim_;
	matrix_.reserve(other.matrix_.size());
	for (const Column_type& col : other.matrix_){
		std::vector<cell_rep_type> tmp(col.begin(), col.end());
		matrix_.emplace_back(col.get_column_index(), tmp, col.get_dimension(), rows_, pivotToColumnIndex_);
	}
	return *this;
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access<Master_matrix>::print() const
{
	std::cout << "Column Matrix:\n";
	for (unsigned int i = 0; i < pivotToColumnIndex_.size() && pivotToColumnIndex_[i] != -1; ++i){
		index pos = pivotToColumnIndex_[i];
		const Column_type& col = matrix_[pos];
		for (const auto &cell : col){
			std::cout << cell.get_row_index() << " ";
		}
		std::cout << "(" << i << ", " << pos << ")\n";
	}
	std::cout << "\n";
	std::cout << "Row Matrix:\n";
	for (unsigned int i = 0; i < pivotToColumnIndex_.size() && pivotToColumnIndex_[i] != -1; ++i){
		index pos = pivotToColumnIndex_[i];
		const Row_type& row = rows_[pos];
		for (const auto &cell : row){
			std::cout << cell.get_column_index() << " ";
		}
		std::cout << "(" << i << ", " << pos << ")\n";
	}
	std::cout << "\n";
}

template<class Master_matrix>
template<class Boundary_type>
inline std::vector<typename Master_matrix::cell_rep_type> Chain_matrix_with_row_access<Master_matrix>::_reduce_boundary(
		index simplexIndex, const Boundary_type& boundary)
{
	tmp_column_type column(boundary.begin(), boundary.end());
	int dim = boundary.begin() == boundary.end() ? 0 : boundary.size() - 1;
	std::vector<cell_rep_type> chainsInH; //for corresponding indices in H (paired columns)
	std::vector<cell_rep_type> chainsInF; //for corresponding indices in F (unpaired, essential columns)

	auto get_last = [&column](){
		if constexpr (Master_matrix::Option_list::is_z2)
			return *(column.rbegin());
		else
			return column.rbegin()->first;
	};

	if (boundary.begin() == boundary.end())
	{
		if constexpr (Master_matrix::Option_list::is_z2)
			column.insert(simplexIndex);
		else
			column.emplace(simplexIndex, 1);
		_insert_chain(column, dim);
		return chainsInF;
	}

	index currentPivot = pivotToColumnIndex_[get_last()];

	while (matrix_[currentPivot].is_paired())
	{
		_reduce_by_G(column, chainsInH, currentPivot);

		if (column.empty()) {
			//produce the sum of all col_h in chains_in_H
			_build_from_H(simplexIndex, column, chainsInH);
			//create a new cycle (in F) sigma - \sum col_h
			_insert_chain(column, dim);
			return chainsInF;
		}

		currentPivot = pivotToColumnIndex_[get_last()];
	}

	while (!column.empty())
	{
		currentPivot = pivotToColumnIndex_[get_last()];

		if (!matrix_[currentPivot].is_paired()) {
			//only fills currentEssentialCycleIndices if Z2 coefficients, so chainsInF remains empty
			_reduce_by_F(column, chainsInF, currentPivot);
		} else {
			_reduce_by_G(column, chainsInH, currentPivot);
		}
	}

	_update_largest_death_in_F(chainsInF);

	//Compute the new column zzsh + \sum col_h, for col_h in chains_in_H
	_build_from_H(simplexIndex, column, chainsInH);

	//Create and insert (\sum col_h) + sigma (in H, paired with chain_fp) in matrix_
	if constexpr (Master_matrix::Option_list::is_z2)
		_insert_chain(column, dim, chainsInF[0]);
	else
		_insert_chain(column, dim, chainsInF[0].first);

	return chainsInF;
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access<Master_matrix>::_reduce_by_G(
		tmp_column_type &column,
		std::vector<cell_rep_type>& chainsInH,
		index currentPivot)
{
	Column_type& col = matrix_[currentPivot];
	if constexpr (Master_matrix::Option_list::is_z2){
		_add_to(col, column);	//Reduce with the column col_g
		chainsInH.push_back(col.get_paired_chain_index());//keep the col_h with which col_g is paired
	} else {
		Field_element_type coef = col.get_pivot_value();
		coef = coef.get_inverse();
		coef *= (Master_matrix::Field_type::get_characteristic() - static_cast<unsigned int>(column.rbegin()->second));

		_add_to(col, column, coef);	//Reduce with the column col_g
		chainsInH.emplace_back(col.get_paired_chain_index(), coef);//keep the col_h with which col_g is paired
	}
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access<Master_matrix>::_reduce_by_F(
		tmp_column_type& column,
		std::vector<cell_rep_type>& chainsInF,
		index currentPivot)
{
	Column_type& col = matrix_[currentPivot];
	if constexpr (Master_matrix::Option_list::is_z2){
		_add_to(col, column);	//Reduce with the column col_g
		chainsInF.push_back(currentPivot);
	} else {
		Field_element_type coef = col.get_pivot_value();
		coef = coef.get_inverse();
		coef *= (Master_matrix::Field_type::get_characteristic() - static_cast<unsigned int>(column.rbegin()->second));

		_add_to(col, column, coef);	//Reduce with the column col_g
		chainsInF.emplace_back(currentPivot, Master_matrix::Field_type::get_characteristic() - static_cast<unsigned int>(coef));
	}
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access<Master_matrix>::_build_from_H(
		index simplexIndex, 
		tmp_column_type& column,
		std::vector<cell_rep_type>& chainsInH)
{
	if constexpr (Master_matrix::Option_list::is_z2){
		column.insert(simplexIndex);
		for (index idx_h : chainsInH) {
			_add_to(matrix_[idx_h], column);
		}
	} else {
		column.emplace(simplexIndex, 1);
		for (std::pair<index,Field_element_type>& idx_h : chainsInH) {
			_add_to(matrix_[idx_h.first], column, idx_h.second);
		}
	}
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access<Master_matrix>::_update_largest_death_in_F(
		const std::vector<cell_rep_type>& chainsInF)
{
	if constexpr (Master_matrix::Option_list::is_z2){
		index toUpdate = chainsInF[0];
		for (auto other_col_it = chainsInF.begin() + 1;
			other_col_it != chainsInF.end();
			 ++other_col_it)
		{
			add_to(*other_col_it, toUpdate);
		}
	} else {
		index toUpdate = chainsInF[0].first;
		matrix_[toUpdate] *= chainsInF[0].second;
		for (auto other_col_it = chainsInF.begin() + 1;
			other_col_it != chainsInF.end();
			 ++other_col_it)
		{
			add_to(other_col_it->second, matrix_[other_col_it->first], toUpdate);
		}
	}
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access<Master_matrix>::_insert_chain(
		const tmp_column_type &column, dimension_type dimension)
{
	matrix_.emplace_back(nextInsertIndex_, column, dimension, rows_, pivotToColumnIndex_);
	pivotToColumnIndex_[nextInsertIndex_] = nextInsertIndex_;

	if constexpr (_barcode_option_is_active()){
		_barcode().emplace_back(dimension, nextInsertIndex_, -1);
		_indexToBar().push_back(_barcode().size() - 1);
	}
	++nextInsertIndex_;
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access<Master_matrix>::_insert_chain(
		const tmp_column_type &column, dimension_type dimension, index pair)
{
	matrix_.emplace_back(nextInsertIndex_, column, dimension, rows_, pivotToColumnIndex_);
	matrix_[nextInsertIndex_].assign_paired_chain(pair);
	matrix_[pair].assign_paired_chain(nextInsertIndex_);
	pivotToColumnIndex_[nextInsertIndex_] = nextInsertIndex_;

	if constexpr (_barcode_option_is_active()){
		_barcode()[_indexToBar()[matrix_[pair].get_pivot()]].death = nextInsertIndex_;
		_indexToBar().push_back(_indexToBar()[matrix_[pair].get_pivot()]);
	}
	++nextInsertIndex_;
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access<Master_matrix>::_add_to(
		const Column_type &column,
		tmp_column_type &set,
		[[maybe_unused]] unsigned int coef)
{
	if constexpr (Master_matrix::Option_list::is_z2){
		std::pair<std::set<index>::iterator,bool> res_insert;
		for (const Cell &cell : column) {
			res_insert = set.insert(cell.get_row_index());
			if (!res_insert.second) {
				set.erase(res_insert.first);
			}
		}
	} else {
		for (const Cell &cell : column) {
			std::pair<index,Field_element_type> p(cell.get_row_index(), cell.get_element());
			auto res_it = set.find(p);

			if (res_it != set.end()){
				p.second *= coef;
				p.second += res_it->second;
				set.erase(res_it);
				if (p.second != Field_element_type::get_additive_identity()){
					set.insert(p);
				}
			} else {
				p.second *= coef;
				set.insert(p);
			}
		}
	}
}

template<class Master_matrix>
inline constexpr bool Chain_matrix_with_row_access<Master_matrix>::_barcode_option_is_active()
{
	return Master_matrix::Option_list::has_column_pairings;
}

template<class Master_matrix>
inline constexpr typename Chain_matrix_with_row_access<Master_matrix>::barcode_type &
Chain_matrix_with_row_access<Master_matrix>::_barcode()
{
	if constexpr (swap_opt::isActive_)
		return swap_opt::template Chain_pairing<Master_matrix>::barcode_;
	else
		return pair_opt::barcode_;
}

template<class Master_matrix>
inline constexpr typename Chain_matrix_with_row_access<Master_matrix>::bar_dictionnary_type &
Chain_matrix_with_row_access<Master_matrix>::_indexToBar()
{
	if constexpr (swap_opt::isActive_)
		return swap_opt::template Chain_pairing<Master_matrix>::indexToBar_;
	else
		return pair_opt::indexToBar_;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // CHAIN_MATRIX_0010_H

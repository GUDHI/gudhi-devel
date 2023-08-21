/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CHAIN_MATRIX_0011_H
#define CHAIN_MATRIX_0011_H

#include <cassert>
#include <iostream>
#include <set>
#include <unordered_map>

#include <gudhi/Simple_object_pool.h>

#include "../utilities/utilities.h"
#include "custom_chain_vine_swap.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class Chain_matrix_with_row_access_with_removals
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

	Chain_matrix_with_row_access_with_removals();
	template<class Boundary_type = boundary_type>
	Chain_matrix_with_row_access_with_removals(const std::vector<Boundary_type>& orderedBoundaries);
	Chain_matrix_with_row_access_with_removals(unsigned int numberOfColumns);
	Chain_matrix_with_row_access_with_removals(
		std::function<bool(index,index)> birthComparator, 
		std::function<bool(index,index)> deathComparator = _no_G_death_comparator);
	template<class Boundary_type = boundary_type>
	Chain_matrix_with_row_access_with_removals(
		const std::vector<Boundary_type>& orderedBoundaries,
		std::function<bool(index,index)> birthComparator, 
		std::function<bool(index,index)> deathComparator = _no_G_death_comparator);
	Chain_matrix_with_row_access_with_removals(
		unsigned int numberOfColumns,
		std::function<bool(index,index)> birthComparator, 
		std::function<bool(index,index)> deathComparator = _no_G_death_comparator);
	Chain_matrix_with_row_access_with_removals(const Chain_matrix_with_row_access_with_removals& matrixToCopy);
	Chain_matrix_with_row_access_with_removals(Chain_matrix_with_row_access_with_removals&& other) noexcept;

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

	Chain_matrix_with_row_access_with_removals& operator=(const Chain_matrix_with_row_access_with_removals& other);
	friend void swap(Chain_matrix_with_row_access_with_removals& matrix1,
					 Chain_matrix_with_row_access_with_removals& matrix2){
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
		matrix1.dimensions_.swap(matrix2.dimensions_);
		std::swap(matrix1.maxDim_, matrix2.maxDim_);
		for (auto& p : matrix1.matrix_){
			Column_type& col = p.second;
			col.set_rows(&matrix1.rows_);
		}
		for (auto& p : matrix2.matrix_){
			Column_type& col = p.second;
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
								std::set<std::pair<index,Field_element_type>,
										 CellPairComparator<Field_element_type> >
							>::type;

	rows_type rows_;	//has to be destroyed after matrix_
	matrix_type matrix_;
	dictionnary_type pivotToColumnIndex_;
	index nextInsertIndex_;
	std::vector<unsigned int> dimensions_;
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
inline Chain_matrix_with_row_access_with_removals<Master_matrix>::Chain_matrix_with_row_access_with_removals()
	: Master_matrix::Chain_pairing_option(),
	  Master_matrix::Chain_vine_swap_option(matrix_),
	  Master_matrix::Chain_representative_cycles_option(matrix_, pivotToColumnIndex_),
	  nextInsertIndex_(0),
	  maxDim_(-1)
{}

template<class Master_matrix>
template<class Boundary_type>
inline Chain_matrix_with_row_access_with_removals<Master_matrix>::Chain_matrix_with_row_access_with_removals(
		const std::vector<Boundary_type> &orderedBoundaries)
	: Master_matrix::Chain_pairing_option(),
	  Master_matrix::Chain_vine_swap_option(matrix_),
	  Master_matrix::Chain_representative_cycles_option(matrix_, pivotToColumnIndex_),
//	  rows_(orderedBoundaries.size()),
	  matrix_(orderedBoundaries.size()),
	  pivotToColumnIndex_(orderedBoundaries.size()),
	  nextInsertIndex_(0),
	  maxDim_(-1)
{
	if constexpr (!Master_matrix::Option_list::has_removable_rows){
		rows_.resize(orderedBoundaries.size());
	}
	for (const Boundary_type &b : orderedBoundaries){
		insert_boundary(b);
	}
}

template<class Master_matrix>
inline Chain_matrix_with_row_access_with_removals<Master_matrix>::Chain_matrix_with_row_access_with_removals(
		unsigned int numberOfColumns)
	: Master_matrix::Chain_pairing_option(),
	  Master_matrix::Chain_vine_swap_option(matrix_),
	  Master_matrix::Chain_representative_cycles_option(matrix_, pivotToColumnIndex_),
//	  rows_(numberOfColumns),
	  matrix_(numberOfColumns),
	  pivotToColumnIndex_(numberOfColumns),
	  nextInsertIndex_(0),
	  maxDim_(-1)
{
	if constexpr (!Master_matrix::Option_list::has_removable_rows){
		rows_.resize(numberOfColumns);
	}
}

template<class Master_matrix>
inline Chain_matrix_with_row_access_with_removals<Master_matrix>::Chain_matrix_with_row_access_with_removals(
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
inline Chain_matrix_with_row_access_with_removals<Master_matrix>::Chain_matrix_with_row_access_with_removals(
		const std::vector<Boundary_type> &orderedBoundaries,
		std::function<bool(index,index)> birthComparator, 
		std::function<bool(index,index)> deathComparator)
	: Master_matrix::Chain_pairing_option(),
	  Master_matrix::Chain_vine_swap_option(matrix_, birthComparator, deathComparator),
	  Master_matrix::Chain_representative_cycles_option(matrix_, pivotToColumnIndex_),
//	  rows_(orderedBoundaries.size()),
	  matrix_(orderedBoundaries.size()),
	  pivotToColumnIndex_(orderedBoundaries.size()),
	  nextInsertIndex_(0),
	  maxDim_(-1)
{
	if constexpr (!Master_matrix::Option_list::has_removable_rows){
		rows_.resize(orderedBoundaries.size());
	}
	for (const Boundary_type &b : orderedBoundaries){
		insert_boundary(b);
	}
}

template<class Master_matrix>
inline Chain_matrix_with_row_access_with_removals<Master_matrix>::Chain_matrix_with_row_access_with_removals(
		unsigned int numberOfColumns,
		std::function<bool(index,index)> birthComparator, 
		std::function<bool(index,index)> deathComparator)
	: Master_matrix::Chain_pairing_option(),
	  Master_matrix::Chain_vine_swap_option(matrix_, birthComparator, deathComparator),
	  Master_matrix::Chain_representative_cycles_option(matrix_, pivotToColumnIndex_),
//	  rows_(numberOfColumns),
	  matrix_(numberOfColumns),
	  pivotToColumnIndex_(numberOfColumns),
	  nextInsertIndex_(0),
	  maxDim_(-1)
{
	if constexpr (!Master_matrix::Option_list::has_removable_rows){
		rows_.resize(numberOfColumns);
	}
}

template<class Master_matrix>
inline Chain_matrix_with_row_access_with_removals<Master_matrix>::Chain_matrix_with_row_access_with_removals(
		const Chain_matrix_with_row_access_with_removals &matrixToCopy)
	: Master_matrix::Chain_pairing_option(matrixToCopy),
	  Master_matrix::Chain_vine_swap_option(matrixToCopy),
	  Master_matrix::Chain_representative_cycles_option(matrixToCopy),
//	  rows_(matrixToCopy.rows_.size()),
	  matrix_(matrixToCopy.matrix_.size()),
	  pivotToColumnIndex_(matrixToCopy.pivotToColumnIndex_),
	  nextInsertIndex_(matrixToCopy.nextInsertIndex_),
	  dimensions_(matrixToCopy.dimensions_),
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
	for (const auto& p : matrixToCopy.matrix_){
		const Column_type& col = p.second;
		std::vector<cell_rep_type> tmp(col.begin(), col.end());
		matrix_.try_emplace(p.first,
						Column_type(col.get_column_index(),
									tmp,
									col.get_dimension(),
									rows_,
									pivotToColumnIndex_));
	}
}

template<class Master_matrix>
inline Chain_matrix_with_row_access_with_removals<Master_matrix>::Chain_matrix_with_row_access_with_removals(
		Chain_matrix_with_row_access_with_removals &&other) noexcept
	: Master_matrix::Chain_pairing_option(std::move(other)),
	  Master_matrix::Chain_vine_swap_option(std::move(other)),
	  Master_matrix::Chain_representative_cycles_option(std::move(other)),
	  rows_(std::move(other.rows_)),
	  matrix_(std::move(other.matrix_)),
	  pivotToColumnIndex_(std::move(other.pivotToColumnIndex_)),
	  nextInsertIndex_(std::exchange(other.nextInsertIndex_, 0)),
	  dimensions_(std::move(other.dimensions_)),
	  maxDim_(std::exchange(other.maxDim_, -1))
{
	if constexpr (rep_opt::isActive_){
		rep_opt::matrix_ = &matrix_;
		rep_opt::pivotToPosition_ = &pivotToColumnIndex_;
	}
	if constexpr (swap_opt::isActive_){
		swap_opt::matrix_ = &matrix_;
	}
	for (auto& p : matrix_){
		Column_type& col = p.second;
		col.set_rows(&rows_);
	}
}

template<class Master_matrix>
template<class Boundary_type>
inline std::vector<typename Master_matrix::cell_rep_type> 
Chain_matrix_with_row_access_with_removals<Master_matrix>::insert_boundary(
		const Boundary_type &boundary)
{
	if constexpr (swap_opt::isActive_ && _barcode_option_is_active()){
		swap_opt::pivotToPosition_.try_emplace(nextInsertIndex_, nextInsertIndex_);
	}
	unsigned int dim = boundary.size() == 0 ? 0 : boundary.size() - 1;
	if (maxDim_ < static_cast<int>(dim)) maxDim_ = dim;
	if (dimensions_.size() <= dim) dimensions_.resize(dim + 1);
	++(dimensions_[dim]);

	return _reduce_boundary(nextInsertIndex_, boundary);
}

template<class Master_matrix>
template<class Boundary_type>
inline std::vector<typename Master_matrix::cell_rep_type> 
Chain_matrix_with_row_access_with_removals<Master_matrix>::insert_boundary(
		index simplexIndex, const Boundary_type &boundary)
{
	if constexpr (swap_opt::isActive_ && _barcode_option_is_active()){
		swap_opt::pivotToPosition_.try_emplace(simplexIndex, simplexIndex);
	}
	unsigned int dim = boundary.size() == 0 ? 0 : boundary.size() - 1;
	if (maxDim_ < static_cast<int>(dim)) maxDim_ = dim;
	if (dimensions_.size() <= dim) dimensions_.resize(dim + 1);
	++(dimensions_[dim]);

	return _reduce_boundary(simplexIndex, boundary);
}

template<class Master_matrix>
inline typename Chain_matrix_with_row_access_with_removals<Master_matrix>::Column_type&
Chain_matrix_with_row_access_with_removals<Master_matrix>::get_column(index columnIndex)
{
	return matrix_.at(columnIndex);
}

template<class Master_matrix>
inline const typename Chain_matrix_with_row_access_with_removals<Master_matrix>::Column_type&
Chain_matrix_with_row_access_with_removals<Master_matrix>::get_column(index columnIndex) const
{
	return matrix_.at(columnIndex);
}

template<class Master_matrix>
inline typename Chain_matrix_with_row_access_with_removals<Master_matrix>::Row_type&
Chain_matrix_with_row_access_with_removals<Master_matrix>::get_row(index rowIndex)
{
	// if constexpr (Master_matrix::Option_list::has_intrusive_rows){
	// 	rows_.at(rowIndex).sort([this](const Cell& c1, const Cell& c2) { 
	// 			return get_pivot(c1.get_column_index()) < 
	// 				   get_pivot(c2.get_column_index());
	// 		});
	// }
	return rows_.at(rowIndex);
}

template<class Master_matrix>
inline const typename Chain_matrix_with_row_access_with_removals<Master_matrix>::Row_type&
Chain_matrix_with_row_access_with_removals<Master_matrix>::get_row(index rowIndex) const
{
	// if constexpr (Master_matrix::Option_list::has_intrusive_rows) {
	// 	rows_.at(rowIndex).sort([this](const Cell& c1, const Cell& c2) { 
	// 			return get_pivot(c1.get_column_index()) < 
	// 				   get_pivot(c2.get_column_index());
	// 		});
	// }
	return rows_.at(rowIndex);
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access_with_removals<Master_matrix>::remove_maximal_simplex(index simplexIndex)
{
	// --nextInsertIndex_;
	index toErase = pivotToColumnIndex_.at(simplexIndex);

	int dim = matrix_.at(toErase).get_dimension();
	--(dimensions_[dim]);
	while (dimensions_.back() == 0)
		dimensions_.pop_back();
	maxDim_ = dimensions_.size() - 1;

	if constexpr (_barcode_option_is_active()){
		index timeStamp;
		if constexpr (swap_opt::isActive_) timeStamp = swap_opt::pivotToPosition_[simplexIndex];
		else timeStamp = simplexIndex;

		typename barcode_type::iterator bar = _indexToBar().at(timeStamp);

		if (bar->death == -1) _barcode().erase(bar);
		else bar->death = -1;

		_indexToBar().erase(timeStamp);
		if constexpr (swap_opt::isActive_) swap_opt::pivotToPosition_.erase(simplexIndex);
	}

	Column_type& c = matrix_.at(toErase);
	index pivot = c.get_pivot();

	assert(rows_[pivot].size() == 1 && "Column asked to be removed do not corresponds to a maximal simplex.");

	if (c.is_paired()) matrix_.at(c.get_paired_chain_index()).unassign_paired_chain();
	pivotToColumnIndex_.erase(simplexIndex);
	matrix_.erase(toErase);
	if constexpr (Master_matrix::Option_list::has_removable_rows)
		rows_.erase(pivot);		//c should have been the only one with an entry in rows_[pivot]
}

template<class Master_matrix>
inline dimension_type Chain_matrix_with_row_access_with_removals<Master_matrix>::get_max_dimension() const
{
	return maxDim_;
}

template<class Master_matrix>
inline unsigned int Chain_matrix_with_row_access_with_removals<Master_matrix>::get_number_of_columns() const
{
	return matrix_.size();
}

template<class Master_matrix>
inline dimension_type Chain_matrix_with_row_access_with_removals<Master_matrix>::get_column_dimension(
		index columnIndex) const
{
	return matrix_.at(columnIndex).get_dimension();
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access_with_removals<Master_matrix>::add_to(
		index sourceColumnIndex, index targetColumnIndex)
{
	matrix_.at(targetColumnIndex) += matrix_.at(sourceColumnIndex);
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access_with_removals<Master_matrix>::add_to(Column_type& sourceColumn, index targetColumnIndex)
{
	matrix_.at(targetColumnIndex) += sourceColumn;
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access_with_removals<Master_matrix>::add_to(Column_type& sourceColumn, const Field_element_type& coefficient, index targetColumnIndex)
{
	matrix_.at(targetColumnIndex).multiply_and_add(coefficient, sourceColumn);
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access_with_removals<Master_matrix>::add_to(const Field_element_type& coefficient, Column_type& sourceColumn, index targetColumnIndex)
{
	matrix_.at(targetColumnIndex).multiply_and_add(sourceColumn, coefficient);
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access_with_removals<Master_matrix>::zero_cell(
		index columnIndex, index rowIndex)
{
	static_assert(!Master_matrix::Option_list::has_removable_columns,
			"'zero_cell' is not implemented for the chosen options.");
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access_with_removals<Master_matrix>::zero_column(index columnIndex)
{
	static_assert(!Master_matrix::Option_list::has_removable_columns,
			"'zero_column' is not implemented for the chosen options.");
}

template<class Master_matrix>
inline bool Chain_matrix_with_row_access_with_removals<Master_matrix>::is_zero_cell(
		index columnIndex, index rowIndex) const
{
	return !matrix_.at(columnIndex).is_non_zero(rowIndex);
}

template<class Master_matrix>
inline bool Chain_matrix_with_row_access_with_removals<Master_matrix>::is_zero_column(index columnIndex) const
{
	return matrix_.at(columnIndex).is_empty();
}

template<class Master_matrix>
inline index Chain_matrix_with_row_access_with_removals<Master_matrix>::get_column_with_pivot(
		index simplexIndex) const
{
	return pivotToColumnIndex_.at(simplexIndex);
}

template<class Master_matrix>
inline index Chain_matrix_with_row_access_with_removals<Master_matrix>::get_pivot(index columnIndex) const
{
	return matrix_.at(columnIndex).get_pivot();
}

template<class Master_matrix>
inline Chain_matrix_with_row_access_with_removals<Master_matrix> &
Chain_matrix_with_row_access_with_removals<Master_matrix>::operator=(const Chain_matrix_with_row_access_with_removals &other)
{
	swap_opt::operator=(other);
	pair_opt::operator=(other);
	rep_opt::operator=(other);
	rows_.reserve(other.rows_.size());
	matrix_.reserve(other.matrix_.size());
	pivotToColumnIndex_ = other.pivotToColumnIndex_;
	nextInsertIndex_ = other.nextInsertIndex_;
	dimensions_ = other.dimensions_;
	maxDim_ = other.maxDim_;
	for (const auto& p : other.matrix_){
		const Column_type& col = p.second;
		std::vector<cell_rep_type> tmp(col.begin(), col.end());
		matrix_.try_emplace(p.first,
						Column_type(col.get_column_index(),
									tmp,
									col.get_dimension(),
									rows_,
									pivotToColumnIndex_));
	}
	return *this;
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access_with_removals<Master_matrix>::print() const
{
	std::cout << "Column Matrix:\n";
	for (const auto& p : pivotToColumnIndex_){
		const Column_type& col = matrix_.at(p.second);
		for (const auto &cell : col){
			std::cout << cell.get_row_index() << " ";
		}
		std::cout << "(" << p.first << ", " << p.second << ")\n";
	}
	std::cout << "\n";
	std::cout << "Row Matrix:\n";
	for (const auto& p : pivotToColumnIndex_){
		const Row_type& row = rows_.at(p.first);
		for (const auto &cell : row){
			std::cout << cell.get_column_index() << " ";
		}
		std::cout << "(" << p.first << ", " << p.second << ")\n";
	}
	std::cout << "\n";
}

template<class Master_matrix>
template<class Boundary_type>
inline std::vector<typename Master_matrix::cell_rep_type> 
Chain_matrix_with_row_access_with_removals<Master_matrix>::_reduce_boundary(
		index simplexIndex, const Boundary_type& boundary)
{
	tmp_column_type column(boundary.begin(), boundary.end());
	int dim = boundary.empty() ? 0 : boundary.size() - 1;
	std::vector<cell_rep_type> chainsInH; //for corresponding indices in H (paired columns)
	std::vector<cell_rep_type> chainsInF; //for corresponding indices in F (unpaired, essential columns)

	auto get_last = [&column](){
		if constexpr (Master_matrix::Option_list::is_z2)
			return *(column.rbegin());
		else
			return column.rbegin()->first;
	};

	if (boundary.empty())
	{
		if constexpr (Master_matrix::Option_list::is_z2)
			column.insert(simplexIndex);
		else
			column.emplace(simplexIndex, 1);
		_insert_chain(column, dim);
		return chainsInF;
	}

	index currentPivot = pivotToColumnIndex_.at(get_last());

	while (matrix_.at(currentPivot).is_paired())
	{
		_reduce_by_G(column, chainsInH, currentPivot);

		if (column.empty()) {
			//produce the sum of all col_h in chains_in_H
			_build_from_H(simplexIndex, column, chainsInH);
			//create a new cycle (in F) sigma - \sum col_h
			_insert_chain(column, dim);
			return chainsInF;
		}

		currentPivot = pivotToColumnIndex_.at(get_last());
	}

	while (!column.empty())
	{
		currentPivot = pivotToColumnIndex_.at(get_last());

		if (!matrix_.at(currentPivot).is_paired()) {
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
inline void Chain_matrix_with_row_access_with_removals<Master_matrix>::_reduce_by_G(
		tmp_column_type &column,
		std::vector<cell_rep_type>& chainsInH,
		index currentPivot)
{
	Column_type& col = matrix_.at(currentPivot);
	if constexpr (Master_matrix::Option_list::is_z2){
		_add_to(col, column);	//Reduce with the column col_g
		chainsInH.push_back(col.get_paired_chain_index());//keep the col_h with which col_g is paired
	} else {
		Field_element_type coef = col.get_pivot_value();
		coef = coef.get_inverse();
		coef *= (Master_matrix::Field_type::get_characteristic()
				 - static_cast<unsigned int>(column.rbegin()->second));

		_add_to(col, column, coef);	//Reduce with the column col_g
		chainsInH.emplace_back(col.get_paired_chain_index(), coef);//keep the col_h with which col_g is paired
	}
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access_with_removals<Master_matrix>::_reduce_by_F(
		tmp_column_type& column,
		std::vector<cell_rep_type>& chainsInF,
		index currentPivot)
{
	Column_type& col = matrix_.at(currentPivot);
	if constexpr (Master_matrix::Option_list::is_z2){
		_add_to(col, column);	//Reduce with the column col_g
		chainsInF.push_back(currentPivot);
	} else {
		Field_element_type coef = col.get_pivot_value();
		coef = coef.get_inverse();
		coef *= (Master_matrix::Field_type::get_characteristic()
				 - static_cast<unsigned int>(column.rbegin()->second));

		_add_to(col, column, coef);	//Reduce with the column col_g
		chainsInF.emplace_back(currentPivot,
							   Master_matrix::Field_type::get_characteristic() - static_cast<unsigned int>(coef));
	}
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access_with_removals<Master_matrix>::_build_from_H(
		index simplexIndex, 
		tmp_column_type& column,
		std::vector<cell_rep_type>& chainsInH)
{
	if constexpr (Master_matrix::Option_list::is_z2){
		column.insert(simplexIndex);
		for (index idx_h : chainsInH) {
			_add_to(matrix_.at(idx_h), column);
		}
	} else {
		column.emplace(simplexIndex, 1);
		for (std::pair<index,Field_element_type>& idx_h : chainsInH) {
			_add_to(matrix_.at(idx_h.first), column, idx_h.second);
		}
	}
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access_with_removals<Master_matrix>::_update_largest_death_in_F(
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
		matrix_.at(toUpdate) *= chainsInF[0].second;
		for (auto other_col_it = chainsInF.begin() + 1;
			other_col_it != chainsInF.end();
			 ++other_col_it)
		{
			add_to(other_col_it->second, matrix_.at(other_col_it->first), toUpdate);
		}
	}
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access_with_removals<Master_matrix>::_insert_chain(
		const tmp_column_type &column, dimension_type dimension)
{
	index pivot;
	if constexpr (Master_matrix::Option_list::is_z2){
		pivot = *(column.rbegin());
	} else {
		pivot = column.rbegin()->first;
	}
	pivotToColumnIndex_.try_emplace(pivot, nextInsertIndex_);

	matrix_.try_emplace(nextInsertIndex_,
					Column_type(nextInsertIndex_, column, dimension, rows_, pivotToColumnIndex_));

	if constexpr (_barcode_option_is_active()){
		_barcode().emplace_back(dimension, pivot, -1);
		_indexToBar().try_emplace(pivot, --_barcode().end());
	}
	++nextInsertIndex_;
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access_with_removals<Master_matrix>::_insert_chain(
		const tmp_column_type &column, dimension_type dimension, index pair)
{
	index pivot;
	if constexpr (Master_matrix::Option_list::is_z2){
		pivot = *(column.rbegin());
	} else {
		pivot = column.rbegin()->first;
	}
	pivotToColumnIndex_.try_emplace(pivot, nextInsertIndex_);

	matrix_.try_emplace(nextInsertIndex_,
					Column_type(nextInsertIndex_, column, dimension, rows_, pivotToColumnIndex_));
	matrix_.at(nextInsertIndex_).assign_paired_chain(pair);
	auto& p = matrix_.at(pair);
	p.assign_paired_chain(nextInsertIndex_);

	if constexpr (_barcode_option_is_active() && swap_opt::isActive_){
		auto barIt = _indexToBar().at(swap_opt::pivotToPosition_[p.get_pivot()]);
		barIt->death = pivot;
		_indexToBar().try_emplace(pivot, barIt);
	} else if constexpr (_barcode_option_is_active()){
		auto barIt = _indexToBar().at(p.get_pivot());
		barIt->death = pivot;
		_indexToBar().try_emplace(pivot, barIt);
	}
	++nextInsertIndex_;
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access_with_removals<Master_matrix>::_add_to(
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
inline constexpr bool Chain_matrix_with_row_access_with_removals<Master_matrix>::_barcode_option_is_active()
{
	return  Master_matrix::Option_list::has_column_pairings;
}

template<class Master_matrix>
inline constexpr typename Chain_matrix_with_row_access_with_removals<Master_matrix>::barcode_type &
Chain_matrix_with_row_access_with_removals<Master_matrix>::_barcode()
{
	if constexpr (swap_opt::isActive_)
		return swap_opt::barcode_;
	else
		return pair_opt::barcode_;
}

template<class Master_matrix>
inline constexpr typename Chain_matrix_with_row_access_with_removals<Master_matrix>::bar_dictionnary_type &
Chain_matrix_with_row_access_with_removals<Master_matrix>::_indexToBar()
{
	if constexpr (swap_opt::isActive_)
		return swap_opt::indexToBar_;
	else
		return pair_opt::indexToBar_;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // CHAIN_MATRIX_0011_H
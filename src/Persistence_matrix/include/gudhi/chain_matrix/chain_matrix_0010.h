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

#include "../utilities.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class Chain_matrix_with_row_access : public Master_matrix::Chain_pairing_option, Master_matrix::Chain_vine_swap_option, Master_matrix::Chain_representative_cycles_option{
public:
	using Column_type = typename Master_matrix::Column_type;
	using Row_type = typename Master_matrix::Column_type::Row_type;
	using Cell = typename Master_matrix::Column_type::Cell;
	using boundary_type = typename Master_matrix::boundary_type;

	Chain_matrix_with_row_access();
	template<class Boundary_type = boundary_type>
	Chain_matrix_with_row_access(std::vector<Boundary_type>& orderedBoundaries);
	Chain_matrix_with_row_access(unsigned int numberOfColumns);
	Chain_matrix_with_row_access(const Chain_matrix_with_row_access& matrixToCopy);
	Chain_matrix_with_row_access(Chain_matrix_with_row_access&& other) noexcept;

	template<class Boundary_type = boundary_type>
	void insert_boundary(Boundary_type& boundary);
	template<class Boundary_type = boundary_type>
	void insert_boundary(Boundary_type& boundary, std::vector<index>& currentEssentialCycleIndices);
	Column_type& get_column(index columnIndex);
	Row_type& get_row(index rowIndex);
	void erase_last();

	dimension_type get_max_dimension() const;
	unsigned int get_number_of_columns() const;

	dimension_type get_column_dimension(index columnIndex) const;

	void add_to(index sourceColumnIndex, index targetColumnIndex);

	void zero_cell(index columnIndex, index rowIndex);
	void zero_column(index columnIndex);
	bool is_zero_cell(index columnIndex, index rowIndex);
	bool is_zero_column(index columnIndex);

	index get_column_with_pivot(index simplexIndex);
	index get_pivot(index columnIndex);

	Chain_matrix_with_row_access& operator=(Chain_matrix_with_row_access other);
	template<class Friend_master_matrix>
	friend void swap(Chain_matrix_with_row_access<Friend_master_matrix>& matrix1,
					 Chain_matrix_with_row_access<Friend_master_matrix>& matrix2);

	void print();  //for debug

private:
	using swap_opt = typename Master_matrix::Chain_vine_swap_option;
	using pair_opt = typename Master_matrix::Chain_pairing_option;
	using rep_opt = typename Master_matrix::Chain_representative_cycles_option;
	using matrix_type = typename Master_matrix::column_container_type;
	using dictionnary_type = typename Master_matrix::template dictionnary_type<index>;
	using barcode_type = typename Master_matrix::barcode_type;
	using bar_dictionnary_type = typename Master_matrix::bar_dictionnary_type;
	using Field_element_type = typename Master_matrix::Field_type;

	matrix_type matrix_;
	dictionnary_type pivotToColumnIndex_;
	index nextInsertIndex_;
	dimension_type maxDim_;

	void _insert_chain(std::set<index>& column, dimension_type dimension);
	void _insert_chain(std::set<index>& column, dimension_type dimension, index pair);
	void _add_to(Column_type& column, std::set<index>& set);
	void _add_to(Column_type &column, std::set<std::pair<index,Field_element_type> >& set, Field_element_type& coef);

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
inline Chain_matrix_with_row_access<Master_matrix>::Chain_matrix_with_row_access(std::vector<Boundary_type> &orderedBoundaries)
	: Master_matrix::Chain_pairing_option(),
	  Master_matrix::Chain_vine_swap_option(matrix_),
	  Master_matrix::Chain_representative_cycles_option(matrix_, pivotToColumnIndex_),
	  matrix_(orderedBoundaries.size(), Column_type(matrix_, pivotToColumnIndex_)),
	  pivotToColumnIndex_(orderedBoundaries.size(), -1),
	  nextInsertIndex_(orderedBoundaries.size()),
	  maxDim_(-1)
{
	for (Boundary_type &b : orderedBoundaries){
		insert_boundary(b);
	}
}

template<class Master_matrix>
inline Chain_matrix_with_row_access<Master_matrix>::Chain_matrix_with_row_access(unsigned int numberOfColumns)
	: Master_matrix::Chain_pairing_option(),
	  Master_matrix::Chain_vine_swap_option(matrix_),
	  Master_matrix::Chain_representative_cycles_option(matrix_, pivotToColumnIndex_),
	  matrix_(numberOfColumns, Column_type(matrix_, pivotToColumnIndex_)),
	  pivotToColumnIndex_(numberOfColumns, -1),
	  nextInsertIndex_(0),
	  maxDim_(-1)
{}

template<class Master_matrix>
inline Chain_matrix_with_row_access<Master_matrix>::Chain_matrix_with_row_access(const Chain_matrix_with_row_access &matrixToCopy)
	: Master_matrix::Chain_pairing_option(matrixToCopy),
	  Master_matrix::Chain_vine_swap_option(matrixToCopy),
	  Master_matrix::Chain_representative_cycles_option(matrixToCopy),
	  matrix_(matrixToCopy.matrix_),
	  pivotToColumnIndex_(matrixToCopy.pivotToColumnIndex_),
	  nextInsertIndex_(matrixToCopy.nextInsertIndex_),
	  maxDim_(matrixToCopy.maxDim_)
{}

template<class Master_matrix>
inline Chain_matrix_with_row_access<Master_matrix>::Chain_matrix_with_row_access(Chain_matrix_with_row_access &&other) noexcept
	: Master_matrix::Chain_pairing_option(std::move(other)),
	  Master_matrix::Chain_vine_swap_option(std::move(other)),
	  Master_matrix::Chain_representative_cycles_option(std::move(other)),
	  matrix_(std::move(other.matrix_)),
	  pivotToColumnIndex_(std::move(other.pivotToColumnIndex_)),
	  nextInsertIndex_(std::move(other.nextInsertIndex_)),
	  maxDim_(std::exchange(other.maxDim_, 0))
{}

template<class Master_matrix>
template<class Boundary_type>
inline void Chain_matrix_with_row_access<Master_matrix>::insert_boundary(Boundary_type &boundary)
{
	std::vector<index> chains_in_F;
	insert_boundary(boundary, chains_in_F);
}

template<class Master_matrix>
template<class Boundary_type>
inline void Chain_matrix_with_row_access<Master_matrix>::insert_boundary(Boundary_type &boundary, std::vector<index>& currentEssentialCycleIndices)
{
	if (matrix_.size() <= nextInsertIndex_){
		matrix_.resize(nextInsertIndex_ * 2, Column_type(matrix_, pivotToColumnIndex_));
		pivotToColumnIndex_.resize(matrix_.size(), -1);
	}
	if constexpr (swap_opt::isActive_){
		if (swap_opt::pivotToPosition_.size() <= nextInsertIndex_) swap_opt::pivotToPosition_.resize(matrix_.size());
		swap_opt::pivotToPosition_.at(nextInsertIndex_) = nextInsertIndex_;
	}
	int dim = boundary.size() - 1;
	if (maxDim_ < dim) maxDim_ = dim;

	if constexpr (Master_matrix::Field_type::get_characteristic() == 2){
		std::set<index> column(boundary.begin(), boundary.end());

		if (boundary.empty())
		{
			column.insert(nextInsertIndex_);
			_insert_chain(column, dim);
			if constexpr (_barcode_option_is_active()){
				_barcode().emplace_back(0, nextInsertIndex_, -1);
				_indexToBar().push_back(_barcode().size() - 1);
			}
			++nextInsertIndex_;
			return;
		}

		index col_low = pivotToColumnIndex_.at(*(column.rbegin()));
		std::vector<index> chains_in_H; //for corresponding indices in H

		while (matrix_.at(col_low).is_paired())
		{
			chains_in_H.push_back(matrix_.at(col_low).get_paired_chain_index());//keep the col_h with which col_g is paired
			_add_to(matrix_.at(col_low), column);	//Reduce with the column col_g

			if (column.empty()) {
				//produce the sum of all col_h in chains_in_H
				column.insert(nextInsertIndex_);
				for (index idx_h : chains_in_H) {
					_add_to(matrix_.at(idx_h), column);
				}
				//create a new cycle (in F) sigma - \sum col_h
				_insert_chain(column, dim);
				if constexpr (_barcode_option_is_active()){
					_barcode().emplace_back(dim, nextInsertIndex_, -1);
					_indexToBar().push_back(_barcode().size() - 1);
				}
				++nextInsertIndex_;
				return;
			}

			col_low = pivotToColumnIndex_.at(*(column.rbegin()));
		}

		while (!column.empty())
		{
			col_low = pivotToColumnIndex_.at(*(column.rbegin()));

			if (!matrix_.at(col_low).is_paired()) {
				currentEssentialCycleIndices.push_back(col_low);
			} else {
				chains_in_H.push_back(matrix_.at(col_low).get_paired_chain_index());
			}

			_add_to(matrix_.at(col_low), column);
		}

		index chain_fp = currentEssentialCycleIndices.front(); //col_fp, with largest death <d index.

		for (std::vector<index>::iterator other_col_it = currentEssentialCycleIndices.begin() + 1;
			other_col_it != currentEssentialCycleIndices.end(); ++other_col_it)
		{
			add_to(*other_col_it, chain_fp);
		}

		//Compute the new column zzsh + \sum col_h, for col_h in chains_in_H
		column.insert(nextInsertIndex_);
		for(index idx_h : chains_in_H){
			_add_to(matrix_.at(idx_h), column);
		}

		//Create and insert (\sum col_h) + sigma (in H, paired with chain_fp) in matrix_
		_insert_chain(column, dim, chain_fp);
		if constexpr (_barcode_option_is_active()){
			_barcode().at(_indexToBar().at(matrix_.at(chain_fp).get_pivot())).death = nextInsertIndex_;
			_indexToBar().push_back(_indexToBar().at(matrix_.at(chain_fp).get_pivot()));
		}
		++nextInsertIndex_;
	} else {
		std::set<std::pair<index,Field_element_type> > column(boundary.begin(), boundary.end());

		if (boundary.empty())
		{
			column.emplace(nextInsertIndex_, 1);
			_insert_chain(column, dim);
			if constexpr (_barcode_option_is_active()){
				_barcode().emplace_back(0, nextInsertIndex_, -1);
				_indexToBar().push_back(_barcode().size() - 1);
			}
			++nextInsertIndex_;
			return;
		}

		index col_low = pivotToColumnIndex_.at(column.rbegin()->first);
		std::vector<std::pair<index,Field_element_type> > chains_in_H; //for corresponding indices in H
		std::vector<std::pair<index,Field_element_type> > chains_in_F; //for corresponding indices in F

		while (matrix_.at(col_low).is_paired())
		{
			Field_element_type coef = matrix_.at(col_low).get_pivot_value();
			coef = coef.get_inverse();
			coef *= (Master_matrix::Field_type::get_characteristic() - static_cast<unsigned int>(column.rbegin()->second));

			_add_to(matrix_.at(col_low), column, coef);	//Reduce with the column col_g
			chains_in_H.emplace_back(matrix_.at(col_low).get_paired_chain_index(), coef);//keep the col_h with which col_g is paired

			if (column.empty()) {
				//produce the sum of all col_h in chains_in_H
				column.emplace(nextInsertIndex_, 1);
				for (std::pair<index,Field_element_type>& idx_h : chains_in_H) {
					_add_to(matrix_.at(idx_h.first), column, idx_h.second);
				}
				//create a new cycle (in F) sigma - \sum col_h
				_insert_chain(column, dim);
				if constexpr (_barcode_option_is_active()){
					_barcode().emplace_back(dim, nextInsertIndex_, -1);
					_indexToBar().push_back(_barcode().size() - 1);
				}
				++nextInsertIndex_;
				return;
			}

			col_low = pivotToColumnIndex_.at(column.rbegin()->first);
		}

		while (!column.empty())
		{
			Field_element_type coef = matrix_.at(col_low).get_pivot_value();
			coef = coef.get_inverse();
			coef *= (Master_matrix::Field_type::get_characteristic() - static_cast<unsigned int>(column.rbegin()->second));

			col_low = pivotToColumnIndex_.at(column.rbegin()->first);

			if (!matrix_.at(col_low).is_paired()) {
				currentEssentialCycleIndices.push_back(col_low);
				chains_in_F.emplace_back(col_low, coef);
			} else {
				chains_in_H.emplace_back(matrix_.at(col_low).get_paired_chain_index(), coef);
			}

			_add_to(matrix_.at(col_low), column);
		}

		index chain_fp = currentEssentialCycleIndices.front(); //col_fp, with largest death <d index.

		for (auto other_col_it = chains_in_F.begin() + 1;
			other_col_it != chains_in_F.end(); ++other_col_it)
		{
			matrix_.at(chain_fp) *= other_col_it->second.get_inverse();
			add_to(other_col_it->first, chain_fp);
		}

		//Compute the new column zzsh + \sum col_h, for col_h in chains_in_H
		column.emplace(nextInsertIndex_, 1);
		for (std::pair<index,Field_element_type>& idx_h : chains_in_H) {
			_add_to(matrix_.at(idx_h.first), column, idx_h.second);
		}

		//Create and insert (\sum col_h) + sigma (in H, paired with chain_fp) in matrix_
		_insert_chain(column, dim, chain_fp);
		if constexpr (_barcode_option_is_active()){
			_barcode().at(_indexToBar().at(matrix_.at(chain_fp).get_pivot())).death = nextInsertIndex_;
			_indexToBar().push_back(_indexToBar().at(matrix_.at(chain_fp).get_pivot()));
		}
		++nextInsertIndex_;
	}
}

template<class Master_matrix>
inline typename Chain_matrix_with_row_access<Master_matrix>::Column_type &Chain_matrix_with_row_access<Master_matrix>::get_column(index columnIndex)
{
	return matrix_.at(columnIndex);
}

template<class Master_matrix>
inline typename Chain_matrix_with_row_access<Master_matrix>::Row_type& Chain_matrix_with_row_access<Master_matrix>::get_row(index rowIndex)
{
	return matrix_.at(pivotToColumnIndex_.at(rowIndex)).get_row();
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access<Master_matrix>::erase_last()
{
	static_assert(static_cast<int>(Master_matrix::Field_type::get_characteristic()) == -1,
			"'erase_last' is not implemented for the chosen options.");
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
	return matrix_.at(columnIndex).get_dimension();
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access<Master_matrix>::add_to(index sourceColumnIndex, index targetColumnIndex)
{
	matrix_.at(targetColumnIndex) += matrix_.at(sourceColumnIndex);
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access<Master_matrix>::zero_cell(index columnIndex, index rowIndex)
{
	static_assert(static_cast<int>(Master_matrix::Field_type::get_characteristic()) == -1,
			"'zero_cell' is not implemented for the chosen options.");
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access<Master_matrix>::zero_column(index columnIndex)
{
	static_assert(static_cast<int>(Master_matrix::Field_type::get_characteristic()) == -1,
			"'zero_column' is not implemented for the chosen options.");
}

template<class Master_matrix>
inline bool Chain_matrix_with_row_access<Master_matrix>::is_zero_cell(index columnIndex, index rowIndex)
{
	return matrix_.at(columnIndex).is_non_zero(rowIndex);
}

template<class Master_matrix>
inline bool Chain_matrix_with_row_access<Master_matrix>::is_zero_column(index columnIndex)
{
	return matrix_.at(columnIndex).is_empty();
}

template<class Master_matrix>
inline index Chain_matrix_with_row_access<Master_matrix>::get_column_with_pivot(index simplexIndex)
{
	return pivotToColumnIndex_.at(simplexIndex);
}

template<class Master_matrix>
inline index Chain_matrix_with_row_access<Master_matrix>::get_pivot(index columnIndex)
{
	return matrix_.at(columnIndex).get_pivot();
}

template<class Master_matrix>
inline Chain_matrix_with_row_access<Master_matrix> &Chain_matrix_with_row_access<Master_matrix>::operator=(Chain_matrix_with_row_access other)
{
	swap_opt::operator=(other);
	pair_opt::operator=(other);
	rep_opt::operator=(other);
	std::swap(matrix_, other.matrix_);
	std::swap(pivotToColumnIndex_, other.pivotToColumnIndex_);
	std::swap(nextInsertIndex_, other.nextInsertIndex_);
	std::swap(maxDim_, other.maxDim_);
	return *this;
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access<Master_matrix>::print()
{
	std::cout << "Column Matrix:\n";
	for (auto p : pivotToColumnIndex_){
		typename Column_type::Column_type& col = matrix_.at(p.second).get_column();
		for (auto it = col.begin(); it != col.end(); it++){
			auto &cell = *it;
			std::cout << cell.get_row_index() << " ";
		}
		std::cout << "(" << p.first << ")\n";
	}
	std::cout << "\n";
	std::cout << "Row Matrix:\n";
	for (auto p : pivotToColumnIndex_){
		Row_type& col = matrix_.at(p.second).get_row();
		for (auto it = col.begin(); it != col.end(); it++){
			auto &cell = *it;
			std::cout << cell.get_column_index() << " ";
		}
		std::cout << "(" << p.second << ")\n";
	}
	std::cout << "\n";
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access<Master_matrix>::_insert_chain(std::set<index> &column, dimension_type dimension)
{
	matrix_.at(nextInsertIndex_) = Column_type(nextInsertIndex_, column, dimension, matrix_, pivotToColumnIndex_);
	pivotToColumnIndex_.at(nextInsertIndex_) = *(column.rbegin());
	for (Cell& cell : matrix_.at(nextInsertIndex_).get_column()){
		matrix_.at(pivotToColumnIndex_.at(cell.get_row_index())).get_row().push_back(cell);
	}
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access<Master_matrix>::_insert_chain(std::set<index> &column, dimension_type dimension, index pair)
{
	matrix_.at(nextInsertIndex_) = Column_type(nextInsertIndex_, column, dimension, matrix_, pivotToColumnIndex_);
	matrix_.at(nextInsertIndex_).assign_paired_chain(pair);
	matrix_.at(pair).assign_paired_chain(nextInsertIndex_);
	pivotToColumnIndex_.at(nextInsertIndex_) = *(column.rbegin());
	for (Cell& cell : matrix_.at(nextInsertIndex_).get_column()){
		matrix_.at(pivotToColumnIndex_.at(cell.get_row_index())).get_row().push_back(cell);
	}
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access<Master_matrix>::_add_to(Column_type &column, std::set<index> &set)
{
	std::pair<std::set<index>::iterator,bool> res_insert;
	for (const Cell &cell : column.get_column()) {
		res_insert = set.insert(cell.get_row_index());
		if (!res_insert.second) {
			set.erase(res_insert.first);
		}
	}
}

template<class Master_matrix>
inline void Chain_matrix_with_row_access<Master_matrix>::_add_to(
		Column_type &column, std::set<std::pair<index, Field_element_type> > &set, Field_element_type &coef)
{
	std::set<std::pair<index,Field_element_type> > newSet;

	for (const Cell &cell : column.get_column()) {
		std::pair<index,Field_element_type> p(cell.get_row_index(), cell.get_element());
		auto res_it = set.find(p);

		if (res_it != set.end()){
			p.second *= coef;
			p.second += res_it->second;
			if (p.second != 0){
				newSet.insert(p);
			}
		} else {
			newSet.insert(p);
		}
	}

	set.swap(newSet);
}

template<class Master_matrix>
inline constexpr bool Chain_matrix_with_row_access<Master_matrix>::_barcode_option_is_active()
{
	return swap_opt::isActive_ || pair_opt::isActive_;
}

template<class Master_matrix>
inline constexpr typename Chain_matrix_with_row_access<Master_matrix>::barcode_type &Chain_matrix_with_row_access<Master_matrix>::_barcode()
{
	if constexpr (swap_opt::isActive_)
		return swap_opt::template Chain_pairing<Master_matrix>::barcode_;
	else
		return pair_opt::barcode_;
}

template<class Master_matrix>
inline constexpr typename Chain_matrix_with_row_access<Master_matrix>::bar_dictionnary_type &Chain_matrix_with_row_access<Master_matrix>::_indexToBar()
{
	if constexpr (swap_opt::isActive_)
		return swap_opt::template Chain_pairing<Master_matrix>::indexToBar_;
	else
		return pair_opt::indexToBar_;
}

template<class Friend_master_matrix>
void swap(Chain_matrix_with_row_access<Friend_master_matrix>& matrix1,
		  Chain_matrix_with_row_access<Friend_master_matrix>& matrix2)
{
	std::swap(static_cast<typename Friend_master_matrix::Chain_pairing_option>(matrix1),
			  static_cast<typename Friend_master_matrix::Chain_pairing_option>(matrix2));
	std::swap(static_cast<typename Friend_master_matrix::Chain_vine_swap_option>(matrix1),
			  static_cast<typename Friend_master_matrix::Chain_vine_swap_option>(matrix2));
	std::swap(static_cast<typename Friend_master_matrix::Chain_representative_cycles_option>(matrix1),
			  static_cast<typename Friend_master_matrix::Chain_representative_cycles_option>(matrix2));
	std::swap(matrix1.matrix_, matrix2.matrix_);
	std::swap(matrix1.pivotToColumnIndex_, matrix2.pivotToColumnIndex_);
	std::swap(matrix1.nextInsertIndex_, matrix2.nextInsertIndex_);
	std::swap(matrix1.maxDim_, matrix2.maxDim_);
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // CHAIN_MATRIX_0010_H

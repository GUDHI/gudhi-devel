/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CHAIN_MATRIX_0000_H
#define CHAIN_MATRIX_0000_H

#include "../utilities.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class Chain_matrix : Master_matrix::Chain_pairing_option, Master_matrix::Chain_vine_swap_option, Master_matrix::Chain_representative_cycles_option{
public:
	using Column_type = typename Master_matrix::Column_type;
	using Row_type = void;
	using boundary_type = typename Master_matrix::boundary_type;

	Chain_matrix();
	template<class Boundary_type = boundary_type>
	Chain_matrix(std::vector<Boundary_type>& orderedBoundaries);
	Chain_matrix(unsigned int numberOfColumns);
	Chain_matrix(Chain_matrix& matrixToCopy);
	Chain_matrix(Chain_matrix&& other) noexcept;

	template<class Boundary_type = boundary_type>
	void insert_boundary(Boundary_type& boundary);
	template<class Boundary_type = boundary_type>
	void insert_boundary(Boundary_type& boundary, std::vector<index>& currentEssentialCycleIndices);
	Column_type& get_column(index columnIndex);
	Row_type get_row(index rowIndex);
	void erase_last();

	dimension_type get_max_dimension() const;
	unsigned int get_number_of_columns() const;

	dimension_type get_column_dimension(index columnIndex) const;

	void add_to(index sourceColumnIndex, index targetColumnIndex);

	void zero_cell(index columnIndex, index rowIndex);
	void zero_column(index columnIndex);
	bool is_zero_cell(index columnIndex, index rowIndex, bool inR = true) const;
	bool is_zero_column(index columnIndex, bool inR = true);

	index get_column_with_pivot(index simplexIndex);
	index get_pivot(index columnIndex);

	Chain_matrix& operator=(Chain_matrix other);
	template<class Friend_master_matrix>
	friend void swap(Chain_matrix<Friend_master_matrix>& matrix1,
					 Chain_matrix<Friend_master_matrix>& matrix2);

	void print();  //for debug

private:
	using swap_opt = typename Master_matrix::Chain_vine_swap_option;
	using pair_opt = typename Master_matrix::Chain_pairing_option;
	using rep_opt = typename Master_matrix::Chain_representative_cycles_option;
	using matrix_type = typename Master_matrix::column_container_type;
	using dictionnary_type = typename Master_matrix::template dictionnary_type<index>;
	using barcode_type = typename Master_matrix::barcode_type;
	using bar_dictionnary_type = typename Master_matrix::bar_dictionnary_type;

	matrix_type matrix_;
	dictionnary_type pivotToColumnIndex_;
	index nextInsertIndex_;
	dimension_type maxDim_;

	static constexpr bool _barcode_option_is_active();
	constexpr barcode_type& _barcode();
	constexpr bar_dictionnary_type& _indexToBar();
};

template<class Master_matrix>
inline Chain_matrix<Master_matrix>::Chain_matrix()
	: Master_matrix::Chain_pairing_option(),
	  Master_matrix::Chain_vine_swap_option(matrix_, pivotToColumnIndex_),
	  Master_matrix::Chain_representative_cycles_option(matrix_, pivotToColumnIndex_),
	  nextInsertIndex_(0),
	  maxDim_(-1)
{}

template<class Master_matrix>
template<class Boundary_type>
inline Chain_matrix<Master_matrix>::Chain_matrix(std::vector<Boundary_type> &orderedBoundaries)
	: Master_matrix::Chain_pairing_option(),
	  Master_matrix::Chain_vine_swap_option(matrix_, pivotToColumnIndex_),
	  Master_matrix::Chain_representative_cycles_option(matrix_, pivotToColumnIndex_),
	  matrix_(orderedBoundaries.size()),
	  pivotToColumnIndex_(orderedBoundaries.size(), -1),
	  nextInsertIndex_(orderedBoundaries.size()),
	  maxDim_(-1)
{
	for (Boundary_type &b : orderedBoundaries){
		insert_boundary(b);
	}
}

template<class Master_matrix>
inline Chain_matrix<Master_matrix>::Chain_matrix(unsigned int numberOfColumns)
	: Master_matrix::Chain_pairing_option(),
	  Master_matrix::Chain_vine_swap_option(matrix_, pivotToColumnIndex_),
	  Master_matrix::Chain_representative_cycles_option(matrix_, pivotToColumnIndex_),
	  matrix_(numberOfColumns),
	  pivotToColumnIndex_(numberOfColumns, -1),
	  nextInsertIndex_(0),
	  maxDim_(-1)
{}

template<class Master_matrix>
inline Chain_matrix<Master_matrix>::Chain_matrix(Chain_matrix &matrixToCopy)
	: Master_matrix::Chain_pairing_option(matrixToCopy),
	  Master_matrix::Chain_vine_swap_option(matrixToCopy),
	  Master_matrix::Chain_representative_cycles_option(matrixToCopy),
	  matrix_(matrixToCopy.matrix_),
	  pivotToColumnIndex_(matrixToCopy.pivotToColumnIndex_),
	  nextInsertIndex_(matrixToCopy.nextInsertIndex_),
	  maxDim_(matrixToCopy.maxDim_)
{}

template<class Master_matrix>
inline Chain_matrix<Master_matrix>::Chain_matrix(Chain_matrix &&other) noexcept
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
inline void Chain_matrix<Master_matrix>::insert_boundary(Boundary_type &boundary)
{
	std::vector<index> chains_in_F;
	insert_boundary(boundary, chains_in_F);
}

template<class Master_matrix>
template<class Boundary_type>
inline void Chain_matrix<Master_matrix>::insert_boundary(Boundary_type &boundary, std::vector<index>& currentEssentialCycleIndices)
{
	if (matrix_.size() <= nextInsertIndex_) matrix_.resize(nextInsertIndex_ * 2);
	pivotToColumnIndex_.emplace(nextInsertIndex_, nextInsertIndex_);

	int dim = boundary.size() - 1;
	if (maxDim_ < dim) maxDim_ = dim;

	if (boundary.empty())
	{
		matrix_.at(nextInsertIndex_) = Column_type({nextInsertIndex_}, 0);
		if constexpr (_barcode_option_is_active()){
			_barcode().emplace_back(0, nextInsertIndex_, -1);
			_indexToBar().push_back(_barcode().size() - 1);
		}
		return;
	}

	Column_type column(boundary);
	//std::set<index> column(boundary.begin(), boundary.end());

	index col_low = pivotToColumnIndex_.at(column.get_pivot());
	std::vector<index> chains_in_H; //for corresponding indices in H

	while (matrix_.at(col_low).is_paired())
	{
		chains_in_H.push_back(matrix_.at(col_low).get_paired_column());//keep the col_h with which col_g is paired
		column += matrix_.at(col_low);	//Reduce with the column col_g

		if (column.is_empty()) {
			column = Column_type({nextInsertIndex_}, dim);
			//produce the sum of all col_h in chains_in_H
			for (index idx_h : chains_in_H) {
				column += matrix_.at(idx_h);
			}
			//create a new cycle (in F) sigma - \sum col_h
			matrix_.at(nextInsertIndex_) = column;
			if constexpr (_barcode_option_is_active()){
				_barcode().emplace_back(dim, nextInsertIndex_, -1);
				_indexToBar().push_back(_barcode().size() - 1);
			}
			return;
		}

		col_low = pivotToColumnIndex_.at(column.get_pivot());
	}

	while (!column.empty())
	{
		col_low = pivotToColumnIndex_.at(column.get_pivot());

		if (!matrix_.at(col_low).is_paired()) {
			currentEssentialCycleIndices.push_back(col_low);
		} else {
			chains_in_H.push_back(matrix_.at(col_low).get_paired_column());
		}

		column += matrix_.at(col_low);
	}

	index chain_fp = currentEssentialCycleIndices.front(); //col_fp, with largest death <d index.

	for (std::vector<index>::iterator other_col_it = currentEssentialCycleIndices.begin() + 1;
		other_col_it != currentEssentialCycleIndices.end(); ++other_col_it)
	{
		matrix_.at(chain_fp) += matrix_.at(*other_col_it);
	}

	//Compute the new column zzsh + \sum col_h, for col_h in chains_in_H
	column = Column_type({nextInsertIndex_}, dim);
	for(index idx_h : chains_in_H){
		column += matrix_.at(idx_h);
	}

	//Create and insert (\sum col_h) + sigma (in H, paired with chain_fp) in matrix_
	column.assign_paired_chain(chain_fp);
	matrix_.at(nextInsertIndex_) = column;
	if constexpr (_barcode_option_is_active()){
		_barcode().at(_indexToBar().at(matrix_.at(chain_fp).get_pivot())).death = nextInsertIndex_;
		_indexToBar().push_back(_indexToBar().at(matrix_.at(chain_fp).get_pivot()));
	}
	++nextInsertIndex_;
}

template<class Master_matrix>
inline typename Chain_matrix<Master_matrix>::Column_type &Chain_matrix<Master_matrix>::get_column(index columnIndex)
{
	return matrix_.at(columnIndex);
}

template<class Master_matrix>
inline typename Chain_matrix<Master_matrix>::Row_type Chain_matrix<Master_matrix>::get_row(index rowIndex)
{
	static_assert(static_cast<int>(Master_matrix::Field_type::get_characteristic()) == -1,
			"'get_row' is not implemented for the chosen options.");
}

template<class Master_matrix>
inline void Chain_matrix<Master_matrix>::erase_last()
{
	static_assert(static_cast<int>(Master_matrix::Field_type::get_characteristic()) == -1,
			"'erase_last' is not implemented for the chosen options.");
}

template<class Master_matrix>
inline dimension_type Chain_matrix<Master_matrix>::get_max_dimension() const
{
	return maxDim_;
}

template<class Master_matrix>
inline unsigned int Chain_matrix<Master_matrix>::get_number_of_columns() const
{
	return matrix_.size();
}

template<class Master_matrix>
inline dimension_type Chain_matrix<Master_matrix>::get_column_dimension(index columnIndex) const
{
	return matrix_.at(columnIndex).get_dimension();
}

template<class Master_matrix>
inline void Chain_matrix<Master_matrix>::add_to(index sourceColumnIndex, index targetColumnIndex)
{
	matrix_.at(targetColumnIndex) += matrix_.at(sourceColumnIndex);
}

template<class Master_matrix>
inline void Chain_matrix<Master_matrix>::zero_cell(index columnIndex, index rowIndex)
{
	static_assert(static_cast<int>(Master_matrix::Field_type::get_characteristic()) == -1,
			"'zero_cell' is not implemented for the chosen options.");
}

template<class Master_matrix>
inline void Chain_matrix<Master_matrix>::zero_column(index columnIndex)
{
	static_assert(static_cast<int>(Master_matrix::Field_type::get_characteristic()) == -1,
			"'zero_column' is not implemented for the chosen options.");
}

template<class Master_matrix>
inline bool Chain_matrix<Master_matrix>::is_zero_cell(index columnIndex, index rowIndex, bool inR) const
{
	if (inR){
		return reducedMatrixR_.is_zero_cell(columnIndex, rowIndex);
	}
	return mirrorMatrixU_.is_zero_cell(columnIndex, rowIndex);
}

template<class Master_matrix>
inline bool Chain_matrix<Master_matrix>::is_zero_column(index columnIndex, bool inR)
{
	if (inR){
		return reducedMatrixR_.is_zero_column(columnIndex);
	}
	return mirrorMatrixU_.is_zero_column(columnIndex);
}

template<class Master_matrix>
inline Chain_matrix<Master_matrix> &Chain_matrix<Master_matrix>::operator=(Chain_matrix other)
{
	swap_opt::operator=(other);
	pair_opt::operator=(other);
	rep_opt::operator=(other);
	std::swap(reducedMatrixR_, other.reducedMatrixR_);
	std::swap(mirrorMatrixU_, other.mirrorMatrixU_);
	std::swap(pivotToColumnIndex_, other.pivotToColumnIndex_);
	std::swap(nextInsertIndex_, other.nextInsertIndex_);
	return *this;
}

template<class Master_matrix>
inline void Chain_matrix<Master_matrix>::print()
{
	std::cout << "R_matrix:\n";
	for (unsigned int i = 0; i < nextInsertIndex_; ++i){
		const Column_type& col = reducedMatrixR_.get_column(i);
		for (auto e : col.get_content(nextInsertIndex_)){
			if (e == 0) std::cout << "-\n";
			else std::cout << e << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";
	std::cout << "U_matrix:\n";
	for (unsigned int i = 0; i < nextInsertIndex_; ++i){
		const Column_type& col = mirrorMatrixU_.get_column(i);
		for (auto e : col.get_content(nextInsertIndex_)){
			if (e == 0) std::cout << "-\n";
			else std::cout << e << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";
}

template<class Master_matrix>
inline void Chain_matrix<Master_matrix>::_initialize_U()
{
	boundary_type id(1);
	if constexpr (Master_matrix::Field_type::get_characteristic() != 2) id.at(0).second = 1;

	for (unsigned int i = 0; i < reducedMatrixR_.get_number_of_columns(); i++){
		if constexpr (Master_matrix::Field_type::get_characteristic() == 2) id.at(0) = i;
		else id.at(0).first = i;
		mirrorMatrixU_.insert_boundary(id);
	}
}

template<class Master_matrix>
inline void Chain_matrix<Master_matrix>::_reduce()
{
	if constexpr (_barcode_option_is_active()){
		_indexToBar().resize(reducedMatrixR_.get_number_of_columns(), -1);
	}

	for (unsigned int i = 0; i < reducedMatrixR_.get_number_of_columns(); i++){
		if (!(reducedMatrixR_.is_zero_column(i)))
		{
			Column_type &curr = reducedMatrixR_.get_column(i);
			int pivot = curr.get_pivot();

			while (pivot != -1 && pivotToColumnIndex_.at(pivot) != -1){
				curr += reducedMatrixR_.get_column(pivotToColumnIndex_.at(pivot));
				mirrorMatrixU_.get_column(i) += mirrorMatrixU_.get_column(pivotToColumnIndex_.at(pivot));
				pivot = curr.get_pivot();
			}

			if (pivot != -1){
				pivotToColumnIndex_.at(pivot) = i;
				if constexpr (_barcode_option_is_active()){
					_barcode().at(_indexToBar().at(pivot)).death = i;
					_indexToBar().at(i) = _indexToBar().at(pivot);
				}
			} else if constexpr (_barcode_option_is_active()){
				_barcode().emplace_back(get_column_dimension(i), i, -1);
				_indexToBar().at(i) = _barcode().size() - 1;
			}
		} else if constexpr (_barcode_option_is_active()){
			_barcode().emplace_back(0, i, -1);
			_indexToBar().at(i) = _barcode().size() - 1;
		}
	}
}

template<class Master_matrix>
inline void Chain_matrix<Master_matrix>::_reduce_last_column()
{
	Column_type &curr = reducedMatrixR_.get_column(nextInsertIndex_);
	if (curr.is_empty()) {
		if constexpr (_barcode_option_is_active()){
			_barcode().emplace_back(0, nextInsertIndex_, -1);
			_indexToBar().push_back(_barcode().size() - 1);
		}
		return;
	}

	int pivot = curr.get_pivot();

	while (pivot != -1 && pivotToColumnIndex_.at(pivot) != -1){
		curr += reducedMatrixR_.get_column(pivotToColumnIndex_.at(pivot));
		mirrorMatrixU_.get_column(nextInsertIndex_) += mirrorMatrixU_.get_column(pivotToColumnIndex_.at(pivot));
		pivot = curr.get_pivot();
	}

	if (pivot != -1){
		pivotToColumnIndex_.at(pivot) = nextInsertIndex_;
		if constexpr (_barcode_option_is_active()){
			_barcode().at(_indexToBar().at(pivot)).death = nextInsertIndex_;
			_indexToBar().push_back(_indexToBar().at(pivot));
		}
	} else if constexpr (_barcode_option_is_active()){
		_barcode().emplace_back(get_column_dimension(nextInsertIndex_), nextInsertIndex_, -1);
		_indexToBar().push_back(_barcode().size() - 1);
	}
}

template<class Master_matrix>
inline constexpr bool Chain_matrix<Master_matrix>::_barcode_option_is_active()
{
	return swap_opt::isActive_ || pair_opt::isActive_;
}

template<class Master_matrix>
inline constexpr typename Chain_matrix<Master_matrix>::barcode_type &Chain_matrix<Master_matrix>::_barcode()
{
	return swap_opt::isActive_ ? swap_opt::indexToBar_ : pair_opt::indexToBar_;
}

template<class Master_matrix>
inline constexpr typename Chain_matrix<Master_matrix>::bar_dictionnary_type &Chain_matrix<Master_matrix>::_indexToBar()
{
	return swap_opt::isActive_ ? swap_opt::barcode_ : pair_opt::barcode_;
}

template<class Friend_master_matrix>
void swap(Chain_matrix<Friend_master_matrix>& matrix1,
		  Chain_matrix<Friend_master_matrix>& matrix2)
{
	std::swap(static_cast<typename Friend_master_matrix::Chain_pairing_option>(matrix1),
			  static_cast<typename Friend_master_matrix::Chain_pairing_option>(matrix2));
	std::swap(static_cast<typename Friend_master_matrix::Chain_vine_swap_option>(matrix1),
			  static_cast<typename Friend_master_matrix::Chain_vine_swap_option>(matrix2));
	std::swap(static_cast<typename Friend_master_matrix::Chain_representative_cycles_option>(matrix1),
			  static_cast<typename Friend_master_matrix::Chain_representative_cycles_option>(matrix2));
	std::swap(matrix1.reducedMatrixR_, matrix2.reducedMatrixR_);
	std::swap(matrix1.mirrorMatrixU_, matrix2.mirrorMatrixU_);
	std::swap(matrix1.pivotToColumnIndex_, matrix2.pivotToColumnIndex_);
	std::swap(matrix1.nextInsertIndex_, matrix2.nextInsertIndex_);
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // CHAIN_MATRIX_0000_H

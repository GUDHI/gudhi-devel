/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022-23 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_CHAIN_MATRIX_H
#define PM_CHAIN_MATRIX_H

#include <algorithm>
#include <iostream>	//print() only
#include <set>
#include <vector>
#include <utility>	//std::swap, std::move & std::exchange

#include <gudhi/Persistence_matrix/overlay_ididx_to_matidx.h>

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class Chain_matrix
		: public Master_matrix::Matrix_dimension_option,
		  public Master_matrix::Chain_pairing_option,
		  public Master_matrix::Chain_vine_swap_option,
		  public Master_matrix::Chain_representative_cycles_option,
		  public Master_matrix::Matrix_row_access_option
{
public:
	using Field_operators = typename Master_matrix::Field_operators;
	using Field_element_type = typename Master_matrix::element_type;
	using Column_type = typename Master_matrix::Column_type;
	using Row_type = typename Master_matrix::Row_type;
	using Cell = typename Master_matrix::Cell_type;
	using Cell_constructor = typename Master_matrix::Cell_constructor;
	using boundary_type = typename Master_matrix::boundary_type;
	using cell_rep_type = typename Master_matrix::cell_rep_type;
	using index = typename Master_matrix::index;
	using id_index = typename Master_matrix::id_index;
	using pos_index = typename Master_matrix::pos_index;
	using dimension_type = typename Master_matrix::dimension_type;

	Chain_matrix(Field_operators* operators, Cell_constructor* cellConstructor);
	template<class Boundary_type = boundary_type>
	Chain_matrix(const std::vector<Boundary_type>& orderedBoundaries, Field_operators* operators, Cell_constructor* cellConstructor);
	Chain_matrix(unsigned int numberOfColumns, Field_operators* operators, Cell_constructor* cellConstructor);
	template<typename EventComparatorFunction>
	Chain_matrix(
		Field_operators* operators, Cell_constructor* cellConstructor,
		EventComparatorFunction&& birthComparator, 
		EventComparatorFunction&& deathComparator);
	template<typename EventComparatorFunction, class Boundary_type = boundary_type>
	Chain_matrix(
		const std::vector<Boundary_type>& orderedBoundaries, 
		Field_operators* operators, Cell_constructor* cellConstructor,
		EventComparatorFunction&& birthComparator, 
		EventComparatorFunction&& deathComparator);
	template<typename EventComparatorFunction>
	Chain_matrix(
		unsigned int numberOfColumns, 
		Field_operators* operators, Cell_constructor* cellConstructor,
		EventComparatorFunction&& birthComparator, 
		EventComparatorFunction&& deathComparator);
	Chain_matrix(const Chain_matrix& matrixToCopy, Field_operators* operators = nullptr, Cell_constructor* cellConstructor = nullptr);
	Chain_matrix(Chain_matrix&& other) noexcept;

	//new simplex = new ID even if the same simplex was already inserted and then removed, ie., an ID cannot come back.
	template<class Boundary_type = boundary_type>
	std::vector<cell_rep_type> insert_boundary(const Boundary_type& boundary, dimension_type dim = -1);
	template<class Boundary_type = boundary_type>
	std::vector<cell_rep_type> insert_boundary(id_index faceID, const Boundary_type& boundary, dimension_type dim = -1);
	Column_type& get_column(index columnIndex);
	const Column_type& get_column(index columnIndex) const;
	void remove_maximal_face(id_index faceID);
	void remove_maximal_face(id_index faceID, const std::vector<index>& columnsToSwap);
	void remove_last();

	index get_number_of_columns() const;

	dimension_type get_column_dimension(index columnIndex) const;

	//avoid calling with specialized options or make it such that it makes sense for persistence
	//=================================================================
	void add_to(index sourceColumnIndex, index targetColumnIndex);
	void multiply_target_and_add_to(index sourceColumnIndex, const Field_element_type& coefficient, index targetColumnIndex);
	void multiply_source_and_add_to(const Field_element_type& coefficient, index sourceColumnIndex, index targetColumnIndex);
	//=================================================================

	bool is_zero_cell(index columnIndex, id_index rowIndex) const;
	bool is_zero_column(index columnIndex);	//just for sanity checks as a valid chain matrix never has an empty column.

	index get_column_with_pivot(id_index faceID) const;
	id_index get_pivot(index columnIndex);

	void reset(Field_operators* operators, Cell_constructor* cellConstructor){
		matrix_.clear();
		pivotToColumnIndex_.clear();
		nextIndex_ = 0;
		operators_ = operators;
		cellPool_ = cellConstructor;
	}

	// void set_operators(Field_operators* operators){ 
	// 	operators_ = operators;
	// 	if constexpr (Master_matrix::Option_list::has_map_column_container){
	// 		for (auto& p : matrix_){
	// 			p.second.set_operators(operators);
	// 		}
	// 	} else {
	// 		for (auto& col : matrix_){
	// 			col.set_operators(operators);
	// 		}
	// 	}
	// }

	Chain_matrix& operator=(const Chain_matrix& other);
	friend void swap(Chain_matrix& matrix1,
					 Chain_matrix& matrix2){
		swap(static_cast<typename Master_matrix::Matrix_dimension_option&>(matrix1),
			 static_cast<typename Master_matrix::Matrix_dimension_option&>(matrix2));
		swap(static_cast<typename Master_matrix::Chain_pairing_option&>(matrix1),
			 static_cast<typename Master_matrix::Chain_pairing_option&>(matrix2));
		swap(static_cast<typename Master_matrix::Chain_vine_swap_option&>(matrix1),
			 static_cast<typename Master_matrix::Chain_vine_swap_option&>(matrix2));
		swap(static_cast<typename Master_matrix::Chain_representative_cycles_option&>(matrix1),
			 static_cast<typename Master_matrix::Chain_representative_cycles_option&>(matrix2));
		matrix1.matrix_.swap(matrix2.matrix_);
		matrix1.pivotToColumnIndex_.swap(matrix2.pivotToColumnIndex_);
		std::swap(matrix1.nextIndex_, matrix2.nextIndex_);
		std::swap(matrix1.operators_, matrix2.operators_);
		std::swap(matrix1.cellPool_, matrix2.cellPool_);

		if constexpr (Master_matrix::Option_list::has_row_access){
			swap(static_cast<typename Master_matrix::Matrix_row_access_option&>(matrix1),
				 static_cast<typename Master_matrix::Matrix_row_access_option&>(matrix2));
			// if constexpr (Master_matrix::Option_list::has_map_column_container){
			// 	for (auto& p : matrix1.matrix_){
			// 		p.second.set_rows(&matrix1.rows_);
			// 	}
			// 	for (auto& p : matrix2.matrix_){
			// 		p.second.set_rows(&matrix2.rows_);
			// 	}
			// } else {
			// 	for (auto& col : matrix1.matrix_){
			// 		col.set_rows(&matrix1.rows_);
			// 	}
			// 	for (auto& col : matrix2.matrix_){
			// 		col.set_rows(&matrix2.rows_);
			// 	}
			// }
		}
	}

	void print() const;  //for debug

	friend class Id_to_index_overlay<Chain_matrix<Master_matrix>, Master_matrix>;

private:
	using dim_opt = typename Master_matrix::Matrix_dimension_option;
	using swap_opt = typename Master_matrix::Chain_vine_swap_option;
	using pair_opt = typename Master_matrix::Chain_pairing_option;
	using rep_opt = typename Master_matrix::Chain_representative_cycles_option;
	using ra_opt = typename Master_matrix::Matrix_row_access_option;
	using matrix_type = typename Master_matrix::column_container_type;
	using dictionnary_type = typename Master_matrix::template dictionnary_type<index>;
	using barcode_type = typename Master_matrix::barcode_type;
	using bar_dictionnary_type = typename Master_matrix::bar_dictionnary_type;
	using tmp_column_type = typename std::conditional<
								Master_matrix::Option_list::is_z2,
								std::set<id_index>,
								std::set<std::pair<id_index,Field_element_type>,typename Master_matrix::CellPairComparator>
							>::type;

	matrix_type matrix_;
	dictionnary_type pivotToColumnIndex_;
	index nextIndex_;
	Field_operators* operators_;
	Cell_constructor* cellPool_;

	template<class Boundary_type>
	std::vector<cell_rep_type> _reduce_boundary(id_index faceID, const Boundary_type& boundary, dimension_type dim);
	void _reduce_by_G(tmp_column_type& column,
					  std::vector<cell_rep_type>& chainsInH,
					  index currentPivot);
	void _reduce_by_F(tmp_column_type& column,
					  std::vector<cell_rep_type>& chainsInF,
					  index currentPivot);
	void _build_from_H(id_index faceID, 
					   tmp_column_type& column,
					   std::vector<cell_rep_type>& chainsInH);
	void _update_largest_death_in_F(const std::vector<cell_rep_type>& chainsInF);
	void _insert_chain(const tmp_column_type& column, dimension_type dimension);
	void _insert_chain(const tmp_column_type& column, dimension_type dimension, index pair);
	void _add_to(const Column_type& column, tmp_column_type& set, unsigned int coef);
	template<typename F>
	void _add_to(Column_type& target, F&& addition);
	void _remove_last(index lastIndex);
	void _update_barcode(pos_index birth);
	void _add_bar(dimension_type dim);

	constexpr barcode_type& _barcode();
	constexpr bar_dictionnary_type& _indexToBar();
	constexpr pos_index& _nextPosition();
};

template<class Master_matrix>
inline Chain_matrix<Master_matrix>::Chain_matrix(Field_operators* operators, Cell_constructor* cellConstructor)
	: dim_opt(-1),
	  pair_opt(),
	  swap_opt(),
	  rep_opt(),
	  ra_opt(),
	  nextIndex_(0),
	  operators_(operators),
	  cellPool_(cellConstructor)
{}

template<class Master_matrix>
template<class Boundary_type>
inline Chain_matrix<Master_matrix>::Chain_matrix(const std::vector<Boundary_type> &orderedBoundaries, Field_operators* operators, Cell_constructor* cellConstructor)
	: dim_opt(-1),
	  pair_opt(),
	  swap_opt(),
	  rep_opt(),
	  ra_opt(orderedBoundaries.size()),
	  nextIndex_(0),
	  operators_(operators),
	  cellPool_(cellConstructor)
{
	matrix_.reserve(orderedBoundaries.size());
	if constexpr (Master_matrix::Option_list::has_map_column_container){
		pivotToColumnIndex_.reserve(orderedBoundaries.size());
	} else {
		pivotToColumnIndex_.resize(orderedBoundaries.size(), -1);
	}

	for (const Boundary_type &b : orderedBoundaries){
		insert_boundary(b);
	}
}

template<class Master_matrix>
inline Chain_matrix<Master_matrix>::Chain_matrix(unsigned int numberOfColumns, Field_operators* operators, Cell_constructor* cellConstructor)
	: dim_opt(-1),
	  pair_opt(),
	  swap_opt(),
	  rep_opt(),
	  ra_opt(numberOfColumns),
	  nextIndex_(0),
	  operators_(operators),
	  cellPool_(cellConstructor)
{
	matrix_.reserve(numberOfColumns);
	if constexpr (Master_matrix::Option_list::has_map_column_container){
		pivotToColumnIndex_.reserve(numberOfColumns);
	} else {
		pivotToColumnIndex_.resize(numberOfColumns, -1);
	}
}

template<class Master_matrix>
template<typename EventComparatorFunction>
inline Chain_matrix<Master_matrix>::Chain_matrix(
		Field_operators* operators, Cell_constructor* cellConstructor,
		EventComparatorFunction&& birthComparator, 
		EventComparatorFunction&& deathComparator)
	: dim_opt(-1),
	  pair_opt(),
	  swap_opt(birthComparator, deathComparator),
	  rep_opt(),
	  ra_opt(),
	  nextIndex_(0),
	  operators_(operators),
	  cellPool_(cellConstructor)
{}

template<class Master_matrix>
template<typename EventComparatorFunction, class Boundary_type>
inline Chain_matrix<Master_matrix>::Chain_matrix(
		const std::vector<Boundary_type> &orderedBoundaries, 
		Field_operators* operators, Cell_constructor* cellConstructor,
		EventComparatorFunction&& birthComparator, 
		EventComparatorFunction&& deathComparator)
	: dim_opt(-1),
	  pair_opt(),
	  swap_opt(birthComparator, deathComparator),
	  rep_opt(),
	  ra_opt(orderedBoundaries.size()),
	  nextIndex_(0),
	  operators_(operators),
	  cellPool_(cellConstructor)
{
	matrix_.reserve(orderedBoundaries.size());
	if constexpr (Master_matrix::Option_list::has_map_column_container){
		pivotToColumnIndex_.reserve(orderedBoundaries.size());
	} else {
		pivotToColumnIndex_.resize(orderedBoundaries.size(), -1);
	}
	for (const Boundary_type &b : orderedBoundaries){
		insert_boundary(b);
	}
}

template<class Master_matrix>
template<typename EventComparatorFunction>
inline Chain_matrix<Master_matrix>::Chain_matrix(
		unsigned int numberOfColumns, 
		Field_operators* operators, Cell_constructor* cellConstructor,
		EventComparatorFunction&& birthComparator, 
		EventComparatorFunction&& deathComparator)
	: dim_opt(-1),
	  pair_opt(),
	  swap_opt(birthComparator, deathComparator),
	  rep_opt(),
	  ra_opt(numberOfColumns),
	  nextIndex_(0),
	  operators_(operators),
	  cellPool_(cellConstructor)
{
	matrix_.reserve(numberOfColumns);
	if constexpr (Master_matrix::Option_list::has_map_column_container){
		pivotToColumnIndex_.reserve(numberOfColumns);
	} else {
		pivotToColumnIndex_.resize(numberOfColumns, -1);
	}
}

template<class Master_matrix>
inline Chain_matrix<Master_matrix>::Chain_matrix(const Chain_matrix &matrixToCopy, Field_operators* operators, Cell_constructor* cellConstructor)
	: dim_opt(static_cast<const dim_opt&>(matrixToCopy)),
	  pair_opt(static_cast<const pair_opt&>(matrixToCopy)),
	  swap_opt(static_cast<const swap_opt&>(matrixToCopy)),
	  rep_opt(static_cast<const rep_opt&>(matrixToCopy)),
	  ra_opt(static_cast<const ra_opt&>(matrixToCopy)),
	  pivotToColumnIndex_(matrixToCopy.pivotToColumnIndex_),
	  nextIndex_(matrixToCopy.nextIndex_),
	  operators_(operators == nullptr ? matrixToCopy.operators_ : operators),
	  cellPool_(cellConstructor == nullptr ? matrixToCopy.cellPool_ : cellConstructor)
{
	matrix_.reserve(matrixToCopy.matrix_.size());
	if constexpr (Master_matrix::Option_list::has_row_access){
		if constexpr (Master_matrix::Option_list::has_map_column_container){
			for (const auto& p : matrixToCopy.matrix_){
				const Column_type& col = p.second;
				matrix_.try_emplace(p.first, Column_type(col, col.get_column_index(), ra_opt::rows_, operators_, cellPool_));
			}
		} else {
			for (const auto& col : matrixToCopy.matrix_){
				matrix_.emplace_back(col, col.get_column_index(), ra_opt::rows_, operators_, cellPool_);
			}
		}
	} else {
		if constexpr (Master_matrix::Option_list::has_map_column_container){
			for (const auto& p : matrixToCopy.matrix_){
				const Column_type& col = p.second;
				matrix_.try_emplace(p.first, Column_type(col, operators_, cellPool_));
			}
		} else {
			for (const auto& col : matrixToCopy.matrix_){
				matrix_.emplace_back(col, operators_, cellPool_);
			}
		}
	}
}

template<class Master_matrix>
inline Chain_matrix<Master_matrix>::Chain_matrix(
		Chain_matrix &&other) noexcept
	: dim_opt(std::move(static_cast<dim_opt&>(other))),
	  pair_opt(std::move(static_cast<pair_opt&>(other))),
	  swap_opt(std::move(static_cast<swap_opt&>(other))),
	  rep_opt(std::move(static_cast<rep_opt&>(other))),
	  ra_opt(std::move(static_cast<ra_opt&>(other))),
	  matrix_(std::move(other.matrix_)),
	  pivotToColumnIndex_(std::move(other.pivotToColumnIndex_)),
	  nextIndex_(std::exchange(other.nextIndex_, 0)),
	  operators_(std::exchange(other.operators_, nullptr)),
	  cellPool_(std::exchange(other.cellPool_, nullptr))
{
	//TODO: not sur this is necessary, as the address of rows_ should not change from the move, no?
	// if constexpr (Master_matrix::Option_list::has_map_column_container){
	// 	for (auto& p : matrix_){
	// 		if constexpr (Master_matrix::Option_list::has_row_access){
	// 			p.second.set_rows(&this->rows_);
	// 		}
	// 	}
	// } else {
	// 	for (auto& col : matrix_){
	// 		if constexpr (Master_matrix::Option_list::has_row_access){
	// 			col.set_rows(&this->rows_);
	// 		}
	// 	}
	// }
}

template<class Master_matrix>
template<class Boundary_type>
inline std::vector<typename Master_matrix::cell_rep_type> Chain_matrix<Master_matrix>::insert_boundary(
		const Boundary_type &boundary, dimension_type dim)
{
	return insert_boundary(nextIndex_, boundary, dim);
}

template<class Master_matrix>
template<class Boundary_type>
inline std::vector<typename Master_matrix::cell_rep_type> Chain_matrix<Master_matrix>::insert_boundary(
		id_index faceID, const Boundary_type &boundary, dimension_type dim)
{
	if constexpr (!Master_matrix::Option_list::has_map_column_container){
		if (pivotToColumnIndex_.size() <= faceID){
			pivotToColumnIndex_.resize(faceID * 2 + 1, -1);
		}
	}

	//TODO: wrong position
	if constexpr (Master_matrix::Option_list::has_vine_update && Master_matrix::Option_list::has_column_pairings){
		if constexpr (Master_matrix::Option_list::has_map_column_container){
			swap_opt::CP::pivotToPosition_.try_emplace(faceID, _nextPosition());
		} else {
			if (swap_opt::CP::pivotToPosition_.size() <= faceID)
				swap_opt::CP::pivotToPosition_.resize(pivotToColumnIndex_.size(), -1);
			swap_opt::CP::pivotToPosition_[faceID] = _nextPosition();
		}
	}

	if constexpr (Master_matrix::Option_list::has_matrix_maximal_dimension_access){
		dim_opt::update_up(dim == -1 ? (boundary.size() == 0 ? 0 : boundary.size() - 1) : dim);
	}

	return _reduce_boundary(faceID, boundary, dim);
}

template<class Master_matrix>
inline typename Chain_matrix<Master_matrix>::Column_type &
Chain_matrix<Master_matrix>::get_column(index columnIndex)
{
	if constexpr (Master_matrix::Option_list::has_map_column_container){
		return matrix_.at(columnIndex);
	} else {
		return matrix_[columnIndex];
	}
}

template<class Master_matrix>
inline const typename Chain_matrix<Master_matrix>::Column_type &
Chain_matrix<Master_matrix>::get_column(index columnIndex) const
{
	if constexpr (Master_matrix::Option_list::has_map_column_container){
		return matrix_.at(columnIndex);
	} else {
		return matrix_[columnIndex];
	}
}

template<class Master_matrix>
inline void Chain_matrix<Master_matrix>::remove_maximal_face(id_index faceID)
{
	static_assert(Master_matrix::Option_list::has_removable_columns,
			"'remove_maximal_face' is not implemented for the chosen options.");
	static_assert(Master_matrix::Option_list::has_map_column_container && 
					Master_matrix::Option_list::has_vine_update &&
					Master_matrix::Option_list::has_column_pairings,
			"'remove_maximal_face' is not implemented for the chosen options.");

	//TODO: find simple test to verify that col at columnIndex is maximal even without row access.

	const auto& pivotToPosition = swap_opt::CP::pivotToPosition_;
	auto it = pivotToPosition.find(faceID);
	if (it == pivotToPosition.end()) return;	//face does not exists. put an assert instead?
	pos_index startPos = it->second;
	index startIndex = pivotToColumnIndex_.at(faceID);

	if (startPos != _nextPosition() - 1){
		std::vector<index> colToSwap;
		colToSwap.reserve(matrix_.size());

		for (auto& p : pivotToPosition){
			if (p.second > startPos) colToSwap.push_back(pivotToColumnIndex_.at(p.first));
		}
		std::sort(colToSwap.begin(), colToSwap.end(), 
					[&](index c1, index c2){ 
						return pivotToPosition.at(get_pivot(c1)) < pivotToPosition.at(get_pivot(c2));
					});

		for (index i : colToSwap){
			startIndex = swap_opt::vine_swap(startIndex, i);
		}
	}

	_remove_last(startIndex);
}

template<class Master_matrix>
inline void Chain_matrix<Master_matrix>::remove_maximal_face(id_index faceID, const std::vector<index>& columnsToSwap)
{
	static_assert(Master_matrix::Option_list::has_removable_columns,
			"'remove_maximal_face' is not implemented for the chosen options.");
	static_assert(Master_matrix::Option_list::has_map_column_container && 
					Master_matrix::Option_list::has_vine_update,
			"'remove_maximal_face' is not implemented for the chosen options.");

	//TODO: find simple test to verify that col at columnIndex is maximal even without row access.

	index startIndex = pivotToColumnIndex_.at(faceID);

	for (index i : columnsToSwap){
		startIndex = swap_opt::vine_swap(startIndex, i);
	}

	_remove_last(startIndex);
}

template<class Master_matrix>
inline void Chain_matrix<Master_matrix>::remove_last()
{
	static_assert(Master_matrix::Option_list::has_removable_columns,
			"'remove_last' is not implemented for the chosen options.");
	static_assert(Master_matrix::Option_list::has_map_column_container || 
					!Master_matrix::Option_list::has_vine_update,
			"'remove_last' is not implemented for the chosen options.");

	if constexpr (Master_matrix::Option_list::has_vine_update){
		//carefull: linear because of the search of the last index. It is better to keep track of the IDIdx index
		//of the last column while performing swaps (or the MatIdx with the return values of `vine_swap` + get_pivot) 
		//and then call `remove_maximal_face` with it and an empty `columnsToSwap`.

		id_index pivot = 0;
		index colIndex;
		for (auto& p : pivotToColumnIndex_){
			if (p.first > pivot){	//pivots have to be strictly increasing in order of filtration
				pivot = p.first;
				colIndex = p.second;
			}
		}
		_remove_last(colIndex);
	} else {
		_remove_last(nextIndex_ - 1);
	}
}

template<class Master_matrix>
inline typename Chain_matrix<Master_matrix>::index Chain_matrix<Master_matrix>::get_number_of_columns() const
{
	if constexpr (Master_matrix::Option_list::has_map_column_container){
//		return nextInsertIndex_;	//if erased columns are viewed as zero columns, otherwise use matrix size.
		return matrix_.size();
	} else {
		return nextIndex_;	//matrix could have been resized much bigger while insert
	}
}

template<class Master_matrix>
inline typename Chain_matrix<Master_matrix>::dimension_type Chain_matrix<Master_matrix>::get_column_dimension(index columnIndex) const
{
	if constexpr (Master_matrix::Option_list::has_map_column_container){
		return matrix_.at(columnIndex).get_dimension();
	} else {
		return matrix_[columnIndex].get_dimension();
	}
}

template<class Master_matrix>
inline void Chain_matrix<Master_matrix>::add_to(index sourceColumnIndex, index targetColumnIndex)
{
	if constexpr (Master_matrix::Option_list::has_map_column_container){
		auto& col = matrix_.at(targetColumnIndex);
		_add_to(col, [&](){col += matrix_.at(sourceColumnIndex);});
	} else {
		auto& col = matrix_[targetColumnIndex];
		_add_to(col, [&](){col += matrix_[sourceColumnIndex];});
	}
}

template<class Master_matrix>
inline void Chain_matrix<Master_matrix>::multiply_target_and_add_to(index sourceColumnIndex, const Field_element_type& coefficient, index targetColumnIndex)
{
	if constexpr (Master_matrix::Option_list::has_map_column_container){
		auto& col = matrix_.at(targetColumnIndex);
		_add_to(col, [&](){col.multiply_and_add(coefficient, matrix_.at(sourceColumnIndex));});
	} else {
		auto& col = matrix_[targetColumnIndex];
		_add_to(col, [&](){col.multiply_and_add(coefficient, matrix_[sourceColumnIndex]);});
	}
}

template<class Master_matrix>
inline void Chain_matrix<Master_matrix>::multiply_source_and_add_to(const Field_element_type& coefficient, index sourceColumnIndex, index targetColumnIndex)
{
	if constexpr (Master_matrix::Option_list::has_map_column_container){
		auto& col = matrix_.at(targetColumnIndex);
		_add_to(col, [&](){col.multiply_and_add(matrix_.at(sourceColumnIndex), coefficient);});
	} else {
		auto& col = matrix_[targetColumnIndex];
		_add_to(col, [&](){col.multiply_and_add(matrix_[sourceColumnIndex], coefficient);});
	}
}

template<class Master_matrix>
inline bool Chain_matrix<Master_matrix>::is_zero_cell(index columnIndex, id_index rowIndex) const
{
	if constexpr (Master_matrix::Option_list::has_map_column_container){
		return !matrix_.at(columnIndex).is_non_zero(rowIndex);
	} else {
		return !matrix_[columnIndex].is_non_zero(rowIndex);
	}
}

template<class Master_matrix>
inline bool Chain_matrix<Master_matrix>::is_zero_column(index columnIndex)
{
	if constexpr (Master_matrix::Option_list::has_map_column_container){
		return matrix_.at(columnIndex).is_empty();
	} else {
		return matrix_[columnIndex].is_empty();
	}
}

template<class Master_matrix>
inline typename Chain_matrix<Master_matrix>::index Chain_matrix<Master_matrix>::get_column_with_pivot(id_index faceID) const
{
	if constexpr (Master_matrix::Option_list::has_map_column_container){
		return pivotToColumnIndex_.at(faceID);
	} else {
		return pivotToColumnIndex_[faceID];
	}
}

template<class Master_matrix>
inline typename Chain_matrix<Master_matrix>::id_index Chain_matrix<Master_matrix>::get_pivot(index columnIndex)
{
	if constexpr (Master_matrix::Option_list::has_map_column_container){
		return matrix_.at(columnIndex).get_pivot();
	} else {
		return matrix_[columnIndex].get_pivot();
	}
}

template<class Master_matrix>
inline Chain_matrix<Master_matrix> &Chain_matrix<Master_matrix>::operator=(
		const Chain_matrix& other)
{
	dim_opt::operator=(other);
	swap_opt::operator=(other);
	pair_opt::operator=(other);
	rep_opt::operator=(other);
	matrix_.clear();
	pivotToColumnIndex_ = other.pivotToColumnIndex_;
	nextIndex_ = other.nextIndex_;
	operators_ = other.operators_;
	cellPool_ = other.cellPool_;

	matrix_.reserve(other.matrix_.size());
	if constexpr (Master_matrix::Option_list::has_row_access){
		ra_opt::operator=(other);
		if constexpr (Master_matrix::Option_list::has_map_column_container){
			for (const auto& p : other.matrix_){
				const Column_type& col = p.second;
				matrix_.try_emplace(p.first, Column_type(col, col.get_column_index(), ra_opt::rows_, operators_, cellPool_));
			}
		} else {
			for (const auto& col : other.matrix_){
				matrix_.emplace_back(col, col.get_column_index(), ra_opt::rows_, operators_, cellPool_);
			}
		}
	} else {
		if constexpr (Master_matrix::Option_list::has_map_column_container){
			for (const auto& p : other.matrix_){
				const Column_type& col = p.second;
				matrix_.try_emplace(p.first, Column_type(col, operators_, cellPool_));
			}
		} else {
			for (const auto& col : other.matrix_){
				matrix_.emplace_back(col, operators_, cellPool_);
			}
		}
	}

	return *this;
}

template<class Master_matrix>
inline void Chain_matrix<Master_matrix>::print() const
{
	std::cout << "Column Matrix:\n";
	if constexpr (!Master_matrix::Option_list::has_map_column_container){
		for (id_index i = 0; i < pivotToColumnIndex_.size() && pivotToColumnIndex_[i] != -1; ++i){
			index pos = pivotToColumnIndex_[i];
			const Column_type& col = matrix_[pos];
			for (const auto &cell : col){
				std::cout << cell.get_row_index() << " ";
			}
			std::cout << "(" << i << ", " << pos << ")\n";
		}
		if constexpr (Master_matrix::Option_list::has_row_access){
			std::cout << "\n";
			std::cout << "Row Matrix:\n";
			for (id_index i = 0; i < pivotToColumnIndex_.size() && pivotToColumnIndex_[i] != -1; ++i){
				index pos = pivotToColumnIndex_[i];
				const Row_type& row = ra_opt::get_row(pos);
				for (const auto &cell : row){
					std::cout << cell.get_column_index() << " ";
				}
				std::cout << "(" << i << ", " << pos << ")\n";
			}
		}
	} else {
		for (const auto& p : pivotToColumnIndex_){
			const Column_type& col = matrix_.at(p.second);
			for (const auto &cell : col){
				std::cout << cell.get_row_index() << " ";
			}
			std::cout << "(" << p.first << ", " << p.second << ")\n";
		}
		if constexpr (Master_matrix::Option_list::has_row_access){
			std::cout << "\n";
			std::cout << "Row Matrix:\n";
			for (const auto& p : pivotToColumnIndex_){
				const Row_type& row = ra_opt::get_row(p.first);
				for (const auto &cell : row){
					std::cout << cell.get_column_index() << " ";
				}
				std::cout << "(" << p.first << ", " << p.second << ")\n";
			}
		}
	}
	std::cout << "\n";
}

template<class Master_matrix>
template<class Boundary_type>
inline std::vector<typename Master_matrix::cell_rep_type> Chain_matrix<Master_matrix>::_reduce_boundary(
		id_index faceID, const Boundary_type& boundary, dimension_type dim)
{
	tmp_column_type column(boundary.begin(), boundary.end());
	if (dim == -1) dim = boundary.begin() == boundary.end() ? 0 : boundary.size() - 1;
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
			column.insert(faceID);
		else
			column.emplace(faceID, 1);
		_insert_chain(column, dim);
		return chainsInF;
	}

	index currentIndex = get_column_with_pivot(get_last());

	while (get_column(currentIndex).is_paired())
	{
		_reduce_by_G(column, chainsInH, currentIndex);

		if (column.empty()) {
			//produce the sum of all col_h in chains_in_H
			_build_from_H(faceID, column, chainsInH);
			//create a new cycle (in F) sigma - \sum col_h
			_insert_chain(column, dim);
			return chainsInF;
		}

		currentIndex = get_column_with_pivot(get_last());
	}

	while (!column.empty())
	{
		currentIndex = get_column_with_pivot(get_last());

		if (!get_column(currentIndex).is_paired()) {
			//only fills currentEssentialCycleIndices if Z2 coefficients, so chainsInF remains empty
			_reduce_by_F(column, chainsInF, currentIndex);
		} else {
			_reduce_by_G(column, chainsInH, currentIndex);
		}
	}

	_update_largest_death_in_F(chainsInF);

	//Compute the new column zzsh + \sum col_h, for col_h in chains_in_H
	_build_from_H(faceID, column, chainsInH);

	//Create and insert (\sum col_h) + sigma (in H, paired with chain_fp) in matrix_
	if constexpr (Master_matrix::Option_list::is_z2)
		_insert_chain(column, dim, chainsInF[0]);
	else
		_insert_chain(column, dim, chainsInF[0].first);

	return chainsInF;
}

template<class Master_matrix>
inline void Chain_matrix<Master_matrix>::_reduce_by_G(
		tmp_column_type &column,
		std::vector<cell_rep_type>& chainsInH,
		index currentIndex)
{
	Column_type& col = get_column(currentIndex);
	if constexpr (Master_matrix::Option_list::is_z2){
		_add_to(col, column, 1u);	//Reduce with the column col_g
		chainsInH.push_back(col.get_paired_chain_index());//keep the col_h with which col_g is paired
	} else {
		Field_element_type coef = col.get_pivot_value();
		coef = operators_->get_inverse(coef);
		coef = operators_->multiply(coef, operators_->get_characteristic() - column.rbegin()->second);

		_add_to(col, column, coef);	//Reduce with the column col_g
		chainsInH.emplace_back(col.get_paired_chain_index(), coef);//keep the col_h with which col_g is paired
	}
}

template<class Master_matrix>
inline void Chain_matrix<Master_matrix>::_reduce_by_F(
		tmp_column_type& column,
		std::vector<cell_rep_type>& chainsInF,
		index currentIndex)
{
	Column_type& col = get_column(currentIndex);
	if constexpr (Master_matrix::Option_list::is_z2){
		_add_to(col, column, 1u);	//Reduce with the column col_g
		chainsInF.push_back(currentIndex);
	} else {
		Field_element_type coef = col.get_pivot_value();
		coef = operators_->get_inverse(coef);
		coef = operators_->multiply(coef, operators_->get_characteristic() - column.rbegin()->second);

		_add_to(col, column, coef);	//Reduce with the column col_g
		chainsInF.emplace_back(currentIndex, operators_->get_characteristic() - coef);
	}
}

template<class Master_matrix>
inline void Chain_matrix<Master_matrix>::_build_from_H(
		id_index faceID, 
		tmp_column_type& column,
		std::vector<cell_rep_type>& chainsInH)
{
	if constexpr (Master_matrix::Option_list::is_z2){
		column.insert(faceID);
		for (index idx_h : chainsInH) {
			_add_to(get_column(idx_h), column, 1u);
		}
	} else {
		column.emplace(faceID, 1);
		for (std::pair<index,Field_element_type>& idx_h : chainsInH) {
			_add_to(get_column(idx_h.first), column, idx_h.second);
		}
	}
}

template<class Master_matrix>
inline void Chain_matrix<Master_matrix>::_update_largest_death_in_F(
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
		get_column(toUpdate) *= chainsInF[0].second;
		for (auto other_col_it = chainsInF.begin() + 1;
			other_col_it != chainsInF.end();
			 ++other_col_it)
		{
			multiply_source_and_add_to(other_col_it->second, other_col_it->first, toUpdate);
		}
	}
}

template<class Master_matrix>
inline void Chain_matrix<Master_matrix>::_insert_chain(
		const tmp_column_type &column, dimension_type dimension)
{
	id_index pivot;
	if constexpr (Master_matrix::Option_list::is_z2){
		pivot = *(column.rbegin());
	} else {
		pivot = column.rbegin()->first;
	}
	if constexpr (Master_matrix::Option_list::has_map_column_container){
		pivotToColumnIndex_.try_emplace(pivot, nextIndex_);

		if constexpr (Master_matrix::Option_list::has_row_access){
			matrix_.try_emplace(nextIndex_, Column_type(nextIndex_, column, dimension, ra_opt::rows_, operators_, cellPool_));
		} else {
			matrix_.try_emplace(nextIndex_, Column_type(column, dimension, operators_, cellPool_));
		}
	} else {
		if constexpr (Master_matrix::Option_list::has_row_access){
			matrix_.emplace_back(nextIndex_, column, dimension, ra_opt::rows_, operators_, cellPool_);
		} else {
			matrix_.emplace_back(column, dimension, operators_, cellPool_);
		}
		
		pivotToColumnIndex_[pivot] = nextIndex_;
	}

	_add_bar(dimension);
	
	++nextIndex_;
}

template<class Master_matrix>
inline void Chain_matrix<Master_matrix>::_insert_chain(
		const tmp_column_type &column, dimension_type dimension, index pair)
{
	//true when no vine updates and if nextIndex_ is updated in remove_last for special case of no vines
	//because then PosIdx == MatIdx
	pos_index pairPos = pair;

	id_index pivot;
	if constexpr (Master_matrix::Option_list::is_z2){
		pivot = *(column.rbegin());
	} else {
		pivot = column.rbegin()->first;
	}

	if constexpr (Master_matrix::Option_list::has_map_column_container){
		pivotToColumnIndex_.try_emplace(pivot, nextIndex_);

		if constexpr (Master_matrix::Option_list::has_row_access){
			matrix_.try_emplace(nextIndex_, Column_type(nextIndex_, column, dimension, ra_opt::rows_, operators_, cellPool_));
		} else {
			matrix_.try_emplace(nextIndex_, Column_type(column, dimension, operators_, cellPool_));
		}

		matrix_.at(nextIndex_).assign_paired_chain(pair);
		auto& p = matrix_.at(pair);
		p.assign_paired_chain(nextIndex_);

		if constexpr (Master_matrix::Option_list::has_column_pairings && Master_matrix::Option_list::has_vine_update){
			pairPos = swap_opt::CP::pivotToPosition_[p.get_pivot()];
		}
	} else {
		if constexpr (Master_matrix::Option_list::has_row_access){
			matrix_.emplace_back(nextIndex_, column, dimension, ra_opt::rows_, operators_, cellPool_);
		} else {
			matrix_.emplace_back(column, dimension, operators_, cellPool_);
		}

		matrix_[nextIndex_].assign_paired_chain(pair);
		matrix_[pair].assign_paired_chain(nextIndex_);
		pivotToColumnIndex_[pivot] = nextIndex_;

		if constexpr (Master_matrix::Option_list::has_column_pairings && Master_matrix::Option_list::has_vine_update){
			pairPos = swap_opt::CP::pivotToPosition_[matrix_[pair].get_pivot()];
		}
	}

	_update_barcode(pairPos);
	
	++nextIndex_;
}

template<class Master_matrix>
inline void Chain_matrix<Master_matrix>::_add_to(
		const Column_type &column,
		tmp_column_type &set,
		[[maybe_unused]] unsigned int coef)
{
	if constexpr (Master_matrix::Option_list::is_z2){
		std::pair<typename std::set<index>::iterator,bool> res_insert;
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
				p.second = operators_->multiply_and_add(p.second, coef, res_it->second);
				set.erase(res_it);
				if (p.second != Field_operators::get_additive_identity()){
					set.insert(p);
				}
			} else {
				p.second = operators_->multiply(p.second, coef);
				set.insert(p);
			}
		}
	}
}

template<class Master_matrix>
template<typename F>
inline void Chain_matrix<Master_matrix>::_add_to(Column_type& target, F&& addition)
{
	auto pivot = target.get_pivot();
	addition();

	if (pivot != target.get_pivot()){
		if constexpr (Master_matrix::Option_list::has_map_column_container){
			std::swap(pivotToColumnIndex_.at(pivot), pivotToColumnIndex_.at(target.get_pivot()));
		} else {
			std::swap(pivotToColumnIndex_[pivot], pivotToColumnIndex_[target.get_pivot()]);
		}
	}
}

template<class Master_matrix>
inline void Chain_matrix<Master_matrix>::_remove_last(index lastIndex)
{
	static_assert(Master_matrix::Option_list::has_removable_columns,
			"'_remove_last' is not implemented for the chosen options.");
	static_assert(Master_matrix::Option_list::has_map_column_container || !Master_matrix::Option_list::has_vine_update,
			"'_remove_last' is not implemented for the chosen options.");

	id_index pivot;

	if constexpr (Master_matrix::Option_list::has_map_column_container){
		auto itToErase = matrix_.find(lastIndex);
		Column_type& colToErase = itToErase->second;
		pivot = colToErase.get_pivot();

		if constexpr (Master_matrix::Option_list::has_matrix_maximal_dimension_access){
			dim_opt::update_down(colToErase.get_dimension());
		}

		if (colToErase.is_paired()) matrix_.at(colToErase.get_paired_chain_index()).unassign_paired_chain();
		pivotToColumnIndex_.erase(pivot);
		matrix_.erase(itToErase);
	} else {
		assert(lastIndex == nextIndex_ - 1 && nextIndex_ == matrix_.size() && "Indexation problem.");

		Column_type& colToErase = matrix_[lastIndex];
		pivot = colToErase.get_pivot();

		if constexpr (Master_matrix::Option_list::has_matrix_maximal_dimension_access){
			dim_opt::update_down(colToErase.get_dimension());
		}

		if (colToErase.is_paired()) matrix_.at(colToErase.get_paired_chain_index()).unassign_paired_chain();
		pivotToColumnIndex_[pivot] = -1;
		matrix_.pop_back();
		//TODO: resize matrix_ when a lot is removed? Could be not the best strategy if user inserts a lot back afterwards.
	}

	if constexpr (!Master_matrix::Option_list::has_vine_update){
		--nextIndex_;	//should not be updated when there are vine updates, as possibly lastIndex != nextIndex - 1 
	}

	if constexpr (Master_matrix::Option_list::has_column_pairings){
		auto it = _indexToBar().find(--_nextPosition());
		typename barcode_type::iterator bar = it->second;

		if (bar->death == -1) _barcode().erase(bar);
		else bar->death = -1;

		_indexToBar().erase(it);
		if constexpr (Master_matrix::Option_list::has_vine_update) swap_opt::CP::pivotToPosition_.erase(pivot);
	}

	if constexpr (Master_matrix::Option_list::has_row_access){
		assert(ra_opt::get_row(pivot).size() == 0 && "Column asked to be removed do not corresponds to a maximal simplex.");
		if constexpr (Master_matrix::Option_list::has_removable_rows){
			ra_opt::erase_row(pivot);
		}
	}
}

template<class Master_matrix>
inline void Chain_matrix<Master_matrix>::_update_barcode(pos_index birth)
{
	if constexpr (Master_matrix::Option_list::has_column_pairings){
		if constexpr (Master_matrix::Option_list::has_removable_columns){
			auto& barIt = _indexToBar().at(birth);
			barIt->death = _nextPosition();
			_indexToBar().try_emplace(_nextPosition(), barIt);	//list so iterators are stable
		} else {
			_barcode()[_indexToBar()[birth]].death = _nextPosition();
			_indexToBar().push_back(_indexToBar()[birth]);
		}
		++_nextPosition();
	}
}

template<class Master_matrix>
inline void Chain_matrix<Master_matrix>::_add_bar(dimension_type dim)
{
	if constexpr (Master_matrix::Option_list::has_column_pairings){
		_barcode().emplace_back(dim, _nextPosition(), -1);
		if constexpr (Master_matrix::Option_list::has_removable_columns){
			_indexToBar().try_emplace(_nextPosition(), --_barcode().end());
		} else {
			_indexToBar().push_back(_barcode().size() - 1);
		}
		++_nextPosition();
	}
}

template<class Master_matrix>
inline constexpr typename Chain_matrix<Master_matrix>::barcode_type &
Chain_matrix<Master_matrix>::_barcode()
{
	if constexpr (Master_matrix::Option_list::has_vine_update)
		return swap_opt::template Chain_pairing<Master_matrix>::barcode_;
	else
		return pair_opt::barcode_;
}

template<class Master_matrix>
inline constexpr typename Chain_matrix<Master_matrix>::bar_dictionnary_type &
Chain_matrix<Master_matrix>::_indexToBar()
{
	if constexpr (Master_matrix::Option_list::has_vine_update)
		return swap_opt::template Chain_pairing<Master_matrix>::indexToBar_;
	else
		return pair_opt::indexToBar_;
}

template<class Master_matrix>
inline constexpr typename Chain_matrix<Master_matrix>::pos_index &
Chain_matrix<Master_matrix>::_nextPosition()
{
	if constexpr (Master_matrix::Option_list::has_vine_update)
		return swap_opt::template Chain_pairing<Master_matrix>::nextPosition_;
	else
		return pair_opt::nextPosition_;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // PM_CHAIN_MATRIX_H

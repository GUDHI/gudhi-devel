/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022-23 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_BASE_MATRIX_H
#define PM_BASE_MATRIX_H

#include <iostream>	//print() only
// #include <stdexcept>
#include <vector>
#include <utility>	//std::swap, std::move & std::exchange

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class Base_matrix
		: public Master_matrix::template Base_swap_option<Base_matrix<Master_matrix> >, 
		  protected Master_matrix::Matrix_row_access_option
{
public:
	using index = typename Master_matrix::index;
	using dimension_type = typename Master_matrix::dimension_type;
	using Field_operators = typename Master_matrix::Field_operators;
	using Field_element_type = typename Master_matrix::element_type;
	using Column_type = typename Master_matrix::Column_type;
	using container_type = typename Master_matrix::boundary_type;
	using Row_type = typename Master_matrix::Row_type;
	using Cell_constructor = typename Master_matrix::Cell_constructor;

	Base_matrix(Field_operators* operators, Cell_constructor* cellConstructor);
	template<class Container_type = container_type>
	Base_matrix(const std::vector<Container_type>& columns, Field_operators* operators, Cell_constructor* cellConstructor);
	Base_matrix(unsigned int numberOfColumns, Field_operators* operators, Cell_constructor* cellConstructor);
	Base_matrix(const Base_matrix& matrixToCopy, Field_operators* operators = nullptr, Cell_constructor* cellConstructor = nullptr);
	Base_matrix(Base_matrix&& other) noexcept;

	template<class Container_type = container_type>
	void insert_column(const Container_type& column);
	template<class Container_type = container_type>
	void insert_column(const Container_type& column, index columnIndex);
	template<class Boundary_type>
	void insert_boundary(const Boundary_type& boundary, dimension_type dim = -1);	//same as insert_column
	Column_type& get_column(index columnIndex);
	//get_row(rowIndex) --> simplex ID (=/= columnIndex)
	Row_type& get_row(index rowIndex);
	void remove_column(index columnIndex);
	void remove_last();
	void erase_row(index rowIndex);		//assumes the row is empty, just thought as index a cleanup

	index get_number_of_columns() const;

	template<class Cell_range_or_column_index>
	void add_to(const Cell_range_or_column_index& sourceColumn, index targetColumnIndex);
	template<class Cell_range_or_column_index>
	void multiply_target_and_add_to(const Cell_range_or_column_index& sourceColumn, const Field_element_type& coefficient, index targetColumnIndex);
	template<class Cell_range_or_column_index>
	void multiply_source_and_add_to(const Field_element_type& coefficient, const Cell_range_or_column_index& sourceColumn, index targetColumnIndex);

	void zero_cell(index columnIndex, index rowIndex);
	void zero_column(index columnIndex);
	bool is_zero_cell(index columnIndex, index rowIndex) const;
	bool is_zero_column(index columnIndex);

	void reset(Field_operators* operators, Cell_constructor* cellConstructor){
		matrix_.clear();
		nextInsertIndex_ = 0;
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

	Base_matrix& operator=(const Base_matrix& other);
	friend void swap(Base_matrix& matrix1, Base_matrix& matrix2){
		swap(static_cast<typename Master_matrix::template Base_swap_option<Base_matrix<Master_matrix> >&>(matrix1),
			 static_cast<typename Master_matrix::template Base_swap_option<Base_matrix<Master_matrix> >&>(matrix2));
		matrix1.matrix_.swap(matrix2.matrix_);
		std::swap(matrix1.nextInsertIndex_, matrix2.nextInsertIndex_);
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

	void print();  //for debug

private:
	using swap_opt = typename Master_matrix::template Base_swap_option<Base_matrix<Master_matrix> >;
	using ra_opt = typename Master_matrix::Matrix_row_access_option;
	using matrix_type = typename Master_matrix::column_container_type;
	using cell_rep_type = typename std::conditional<
								Master_matrix::Option_list::is_z2,
								index,
								std::pair<index,Field_element_type>
							>::type;

	friend swap_opt;	//direct access to matrix_ to avoid row reorder.

	matrix_type matrix_;
	index nextInsertIndex_;
	Field_operators* operators_;
	Cell_constructor* cellPool_;

	template<class Container_type = container_type>
	void _insert(const Container_type& column, index columnIndex, dimension_type dim);
	void _orderRowsIfNecessary();
};

template<class Master_matrix>
inline Base_matrix<Master_matrix>::Base_matrix(Field_operators* operators, Cell_constructor* cellConstructor)
	: swap_opt(),
	  ra_opt(),
	  nextInsertIndex_(0),
	  operators_(operators),
	  cellPool_(cellConstructor)
{}

template<class Master_matrix>
template<class Container_type>
inline Base_matrix<Master_matrix>::Base_matrix(const std::vector<Container_type> &columns, Field_operators* operators, Cell_constructor* cellConstructor)
	: swap_opt(columns.size()),
	  ra_opt(columns.size()),	//not ideal if max row index is much smaller than max column index, does that happen often?
	  matrix_(!Master_matrix::Option_list::has_map_column_container && Master_matrix::Option_list::has_row_access ? 0 : columns.size()),
	  nextInsertIndex_(columns.size()),
	  operators_(operators),
	  cellPool_(cellConstructor)
{
	if constexpr (!Master_matrix::Option_list::has_map_column_container && Master_matrix::Option_list::has_row_access)
		matrix_.reserve(columns.size());

	for (index i = 0; i < columns.size(); i++){
		if constexpr (Master_matrix::Option_list::has_map_column_container){
			if constexpr (Master_matrix::Option_list::has_row_access){
				matrix_.try_emplace(i, Column_type(i, columns[i], ra_opt::rows_, operators_, cellPool_));
			} else {
				matrix_.try_emplace(i, Column_type(columns[i], operators_, cellPool_));
			}
		} else {
			if constexpr (Master_matrix::Option_list::has_row_access){
				matrix_.emplace_back(i, columns[i], ra_opt::rows_, operators_, cellPool_);
			} else {
				matrix_[i] = Column_type(columns[i], operators_, cellPool_);
			}
		}
	}
}

template<class Master_matrix>
inline Base_matrix<Master_matrix>::Base_matrix(unsigned int numberOfColumns, Field_operators* operators, Cell_constructor* cellConstructor)
	: swap_opt(numberOfColumns),
	  ra_opt(numberOfColumns),
	  matrix_(!Master_matrix::Option_list::has_map_column_container && Master_matrix::Option_list::has_row_access ? 0 : numberOfColumns),
	  nextInsertIndex_(0),
	  operators_(operators),
	  cellPool_(cellConstructor)
{
	if constexpr (!Master_matrix::Option_list::has_map_column_container && Master_matrix::Option_list::has_row_access)
		matrix_.reserve(numberOfColumns);
}

template<class Master_matrix>
inline Base_matrix<Master_matrix>::Base_matrix(const Base_matrix &matrixToCopy, Field_operators* operators, Cell_constructor* cellConstructor)
	: swap_opt(static_cast<const swap_opt&>(matrixToCopy)),
	  ra_opt(static_cast<const ra_opt&>(matrixToCopy)),
	  nextInsertIndex_(matrixToCopy.nextInsertIndex_),
	  operators_(operators == nullptr ? matrixToCopy.operators_ : operators),
	  cellPool_(cellConstructor == nullptr ? matrixToCopy.cellPool_ : cellConstructor)
{
	matrix_.reserve(matrixToCopy.matrix_.size());
	if constexpr (Master_matrix::Option_list::has_map_column_container){
		for (const auto& p : matrixToCopy.matrix_){
			const Column_type& col = p.second;
			if constexpr (Master_matrix::Option_list::has_row_access){
				matrix_.try_emplace(p.first, Column_type(col, col.get_column_index(), ra_opt::rows_, operators_, cellPool_));
			} else {
				matrix_.try_emplace(p.first, Column_type(col, operators_, cellPool_));
			}
		}
	} else {
		for (const auto& col : matrixToCopy.matrix_){
			if constexpr (Master_matrix::Option_list::has_row_access){
				matrix_.emplace_back(col, col.get_column_index(), ra_opt::rows_, operators_, cellPool_);
			} else {
				matrix_.emplace_back(col, operators_, cellPool_);
			}
		}
	}
}

template<class Master_matrix>
inline Base_matrix<Master_matrix>::Base_matrix(Base_matrix &&other) noexcept
	: swap_opt(std::move(static_cast<swap_opt&>(other))),
	  ra_opt(std::move(static_cast<ra_opt&>(other))),
	  matrix_(std::move(other.matrix_)),
	  nextInsertIndex_(std::exchange(other.nextInsertIndex_, 0)),
	  operators_(std::exchange(other.operators_, nullptr)),
	  cellPool_(std::exchange(other.cellPool_, nullptr))
{
	//TODO: not sur this is necessary, as the address of rows_ should not change from the move, no?
	// if constexpr (Master_matrix::Option_list::has_row_access){
	// 	if constexpr (Master_matrix::Option_list::has_map_column_container){
	// 		for (auto& p : matrix_){
	// 			p.second.set_rows(&this->rows_);
	// 		}
	// 	} else {
	// 		for (auto& col : matrix_){
	// 			col.set_rows(&this->rows_);
	// 		}
	// 	}
	// }
}

template<class Master_matrix>
template<class Container_type>
inline void Base_matrix<Master_matrix>::insert_column(const Container_type &column)
{
	_insert(column, nextInsertIndex_, column.size() == 0 ? 0 : column.size() - 1);
	++nextInsertIndex_;
}

template<class Master_matrix>
template<class Container_type>
inline void Base_matrix<Master_matrix>::insert_column(const Container_type &column, index columnIndex)
{
	if (columnIndex >= nextInsertIndex_) nextInsertIndex_ = columnIndex + 1;
	_insert(column, columnIndex, column.size() == 0 ? 0 : column.size() - 1);
}

template<class Master_matrix>
template<class Boundary_type>
inline void Base_matrix<Master_matrix>::insert_boundary(const Boundary_type &boundary, dimension_type dim)
{
	if (dim == -1) dim = boundary.size() == 0 ? 0 : boundary.size() - 1;
	_insert(boundary, nextInsertIndex_++, dim);
}

template<class Master_matrix>
inline typename Base_matrix<Master_matrix>::Column_type& Base_matrix<Master_matrix>::get_column(index columnIndex)
{
	_orderRowsIfNecessary();
	if constexpr (Master_matrix::Option_list::has_map_column_container){
		return matrix_.at(columnIndex);
	} else {
		return matrix_[columnIndex];
	}
}

template<class Master_matrix>
inline typename Base_matrix<Master_matrix>::Row_type &Base_matrix<Master_matrix>::get_row(index rowIndex)
{
	static_assert(Master_matrix::Option_list::has_row_access, "Row access has to be enabled for this method.");

	_orderRowsIfNecessary();
	return ra_opt::get_row(rowIndex);
}

template<class Master_matrix>
inline void Base_matrix<Master_matrix>::remove_column(index columnIndex)
{
	static_assert(Master_matrix::Option_list::has_map_column_container,
			"'remove_column' is not implemented for the chosen options.");
	
	if (columnIndex == nextInsertIndex_ - 1) --nextInsertIndex_;	//assumes somehow that the user does not do random things with the insertions at given column indices

	matrix_.erase(columnIndex);
}

template<class Master_matrix>
inline void Base_matrix<Master_matrix>::remove_last()
{
	--nextInsertIndex_;	//assumes that the columns are consecutive
	if constexpr (Master_matrix::Option_list::has_map_column_container){
		matrix_.erase(nextInsertIndex_);
	} else {
		if constexpr (Master_matrix::Option_list::has_row_access){
			assert(nextInsertIndex_ == matrix_.size() - 1 && "Indexation problem.");
			matrix_.pop_back();
		} else {
			matrix_[nextInsertIndex_].clear();
		}
	}
}

template<class Master_matrix>
inline void Base_matrix<Master_matrix>::erase_row(index rowIndex)
{
	if constexpr (Master_matrix::Option_list::has_column_and_row_swaps || Master_matrix::Option_list::has_vine_update){
		if constexpr (Master_matrix::Option_list::has_row_access && Master_matrix::Option_list::has_removable_rows){
			ra_opt::erase_row(swap_opt::indexToRow_[rowIndex]);
		}
		if constexpr (Master_matrix::Option_list::has_map_column_container){
			auto it = swap_opt::indexToRow_.find(rowIndex);
			swap_opt::rowToIndex_.erase(it->second);
			swap_opt::indexToRow_.erase(it);
		}
	} else {
		if constexpr (Master_matrix::Option_list::has_row_access && Master_matrix::Option_list::has_removable_rows){
			ra_opt::erase_row(rowIndex);
		}
	}
}

template<class Master_matrix>
inline typename Base_matrix<Master_matrix>::index Base_matrix<Master_matrix>::get_number_of_columns() const
{
	if constexpr (Master_matrix::Option_list::has_map_column_container){
//		return nextInsertIndex_;	//if erased columns are viewed as zero columns, otherwise use matrix size.
		return matrix_.size();
	} else {
		return nextInsertIndex_;	//matrix could have been resized much bigger while insert
	}
}

template<class Master_matrix>
template<class Cell_range_or_column_index>
inline void Base_matrix<Master_matrix>::add_to(const Cell_range_or_column_index& sourceColumn, index targetColumnIndex)
{
	if constexpr (std::is_integral_v<Cell_range_or_column_index>){
		if constexpr (Master_matrix::Option_list::has_map_column_container){
			matrix_.at(targetColumnIndex) += matrix_.at(sourceColumn);
		} else {
			matrix_[targetColumnIndex] += matrix_[sourceColumn];
		}
	} else {
		if constexpr (Master_matrix::Option_list::has_map_column_container){
			matrix_.at(targetColumnIndex) += sourceColumn;
		} else {
			matrix_[targetColumnIndex] += sourceColumn;
		}
	}
}

template<class Master_matrix>
template<class Cell_range_or_column_index>
inline void Base_matrix<Master_matrix>::multiply_target_and_add_to(const Cell_range_or_column_index& sourceColumn, const Field_element_type& coefficient, index targetColumnIndex)
{
	if constexpr (std::is_integral_v<Cell_range_or_column_index>){
		if constexpr (Master_matrix::Option_list::has_map_column_container){
			matrix_.at(targetColumnIndex).multiply_and_add(coefficient, matrix_.at(sourceColumn));
		} else {
			matrix_[targetColumnIndex].multiply_and_add(coefficient, matrix_[sourceColumn]);
		}
	} else {
		if constexpr (Master_matrix::Option_list::has_map_column_container){
			matrix_.at(targetColumnIndex).multiply_and_add(coefficient, sourceColumn);
		} else {
			matrix_[targetColumnIndex].multiply_and_add(coefficient, sourceColumn);
		}
	}
}

template<class Master_matrix>
template<class Cell_range_or_column_index>
inline void Base_matrix<Master_matrix>::multiply_source_and_add_to(const Field_element_type& coefficient, const Cell_range_or_column_index& sourceColumn, index targetColumnIndex)
{
	if constexpr (std::is_integral_v<Cell_range_or_column_index>){
		if constexpr (Master_matrix::Option_list::has_map_column_container){
			matrix_.at(targetColumnIndex).multiply_and_add(matrix_.at(sourceColumn), coefficient);
		} else {
			matrix_[targetColumnIndex].multiply_and_add(matrix_[sourceColumn], coefficient);
		}
	} else {
		if constexpr (Master_matrix::Option_list::has_map_column_container){
			matrix_.at(targetColumnIndex).multiply_and_add(sourceColumn, coefficient);
		} else {
			matrix_[targetColumnIndex].multiply_and_add(sourceColumn, coefficient);
		}
	}
}

template<class Master_matrix>
inline void Base_matrix<Master_matrix>::zero_cell(index columnIndex, index rowIndex)
{
	if constexpr (Master_matrix::Option_list::has_map_column_container){
		if constexpr (Master_matrix::Option_list::has_column_and_row_swaps || Master_matrix::Option_list::has_vine_update){
			matrix_.at(columnIndex).clear(swap_opt::indexToRow_.at(rowIndex));	//TODO: reorganize columns
		} else {
			matrix_.at(columnIndex).clear(rowIndex);	//TODO: reorganize columns
		}
	} else {
		if constexpr (Master_matrix::Option_list::has_column_and_row_swaps || Master_matrix::Option_list::has_vine_update){
			matrix_[columnIndex].clear(swap_opt::indexToRow_[rowIndex]);	//TODO: reorganize columns
		} else {
			matrix_[columnIndex].clear(rowIndex);	//TODO: reorganize columns
		}
	}
}

template<class Master_matrix>
inline void Base_matrix<Master_matrix>::zero_column(index columnIndex)
{
	if constexpr (Master_matrix::Option_list::has_map_column_container){
		matrix_.at(columnIndex).clear();	//TODO: reorganize columns
	} else {
		matrix_[columnIndex].clear();	//TODO: reorganize columns
	}
}

template<class Master_matrix>
inline bool Base_matrix<Master_matrix>::is_zero_cell(index columnIndex, index rowIndex) const
{
	if constexpr (Master_matrix::Option_list::has_map_column_container){
		if constexpr (Master_matrix::Option_list::has_column_and_row_swaps || Master_matrix::Option_list::has_vine_update){
			return !(matrix_.at(columnIndex).is_non_zero(swap_opt::indexToRow_.at(rowIndex)));
		} else {
			return !(matrix_.at(columnIndex).is_non_zero(rowIndex));
		}
	} else {	//operator[] non const for maps
		if constexpr (Master_matrix::Option_list::has_column_and_row_swaps || Master_matrix::Option_list::has_vine_update){
			return !(matrix_[columnIndex].is_non_zero(swap_opt::indexToRow_[rowIndex]));
		} else {
			return !(matrix_[columnIndex].is_non_zero(rowIndex));
		}
	}
}

template<class Master_matrix>
inline bool Base_matrix<Master_matrix>::is_zero_column(index columnIndex)
{
	if constexpr (Master_matrix::Option_list::has_map_column_container){
		return matrix_.at(columnIndex).is_empty();
	} else {
		return matrix_[columnIndex].is_empty();
	}
}

template<class Master_matrix>
inline Base_matrix<Master_matrix> &Base_matrix<Master_matrix>::operator=(const Base_matrix& other)
{
	swap_opt::operator=(other);
	ra_opt::operator=(other);
	matrix_.clear();
	nextInsertIndex_ = other.nextInsertIndex_;
	operators_ = other.operators_;
	cellPool_ = other.cellPool_;

	matrix_.reserve(other.matrix_.size());
	if constexpr (Master_matrix::Option_list::has_map_column_container){
		for (const auto& p : other.matrix_){
			const Column_type& col = p.second;
			if constexpr (Master_matrix::Option_list::has_row_access){
				matrix_.try_emplace(p.first, Column_type(col, col.get_column_index(), ra_opt::rows_, operators_, cellPool_));
			} else {
				matrix_.try_emplace(p.first, Column_type(col, operators_, cellPool_));
			}
		}
	} else {
		for (const auto& col : other.matrix_){
			if constexpr (Master_matrix::Option_list::has_row_access){
				matrix_.emplace_back(col, col.get_column_index(), ra_opt::rows_, operators_, cellPool_);
			} else {
				matrix_.emplace_back(col, operators_, cellPool_);
			}
		}
	}

	return *this;
}

template<class Master_matrix>
inline void Base_matrix<Master_matrix>::print()
{
	_orderRowsIfNecessary();
	std::cout << "Base_matrix:\n";
	for (index i = 0; i < nextInsertIndex_; ++i){
		const Column_type& col = matrix_[i];
		for (const auto& e : col.get_content(nextInsertIndex_)){
			if (e == 0u) std::cout << "- ";
			else std::cout << e << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";
	if constexpr (Master_matrix::Option_list::has_row_access){
		std::cout << "Row Matrix:\n";
		for (index i = 0; i < nextInsertIndex_; ++i){
			const auto& row = ra_opt::rows_[i];
			for (const auto &cell : row){
				std::cout << cell.get_column_index() << " ";
			}
			std::cout << "(" << i << ")\n";
		}
		std::cout << "\n";
	}
}

template<class Master_matrix>
template<class Container_type>
inline void Base_matrix<Master_matrix>::_insert(const Container_type& column, index columnIndex, dimension_type dim)
{
	_orderRowsIfNecessary();

	index pivot = 0;
	if (column.begin() != column.end()){
		if constexpr (Master_matrix::Option_list::is_z2){
			pivot = *std::prev(column.end());
		} else {
			pivot = std::prev(column.end())->first;
		}
		if constexpr (Master_matrix::Option_list::has_row_access && !Master_matrix::Option_list::has_removable_rows)
			if (ra_opt::rows_->size() <= pivot) ra_opt::rows_->resize(pivot + 1);
	}
	
	if constexpr (Master_matrix::Option_list::has_map_column_container){
		if constexpr (Master_matrix::Option_list::has_column_and_row_swaps || Master_matrix::Option_list::has_vine_update){
			for (auto id : column){
				index idx;
				if constexpr (Master_matrix::Option_list::is_z2){
					idx = id;
				} else {
					idx = id.first;
				}
				swap_opt::indexToRow_[idx] = idx;
				swap_opt::rowToIndex_[idx] = idx;
			}
		}

		if constexpr (Master_matrix::Option_list::has_row_access){
			matrix_.try_emplace(columnIndex, Column_type(columnIndex, column, dim, ra_opt::rows_, operators_, cellPool_));
		} else {
			matrix_.try_emplace(columnIndex, column, dim, operators_, cellPool_);
		}
	} else {
		if constexpr (Master_matrix::Option_list::has_column_and_row_swaps || Master_matrix::Option_list::has_vine_update){
			index size = swap_opt::indexToRow_.size();
			if (size <= pivot){
				for (index i = size; i <= pivot; i++){
					swap_opt::indexToRow_.push_back(i);
					swap_opt::rowToIndex_.push_back(i);
				}
			}
		}
		if constexpr (Master_matrix::Option_list::has_row_access){
			matrix_.emplace_back(columnIndex, column, dim, ra_opt::rows_, operators_, cellPool_);
		} else {
			if (matrix_.size() <= columnIndex) {
				matrix_.resize(columnIndex + 1);
			}
			matrix_[columnIndex] = Column_type(column, dim, operators_, cellPool_);
		}
	}
}

template<class Master_matrix>
inline void Base_matrix<Master_matrix>::_orderRowsIfNecessary(){
	if constexpr (Master_matrix::Option_list::has_column_and_row_swaps || Master_matrix::Option_list::has_vine_update){
		if (swap_opt::rowSwapped_) swap_opt::_orderRows();
	}
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // PM_BASE_MATRIX_H

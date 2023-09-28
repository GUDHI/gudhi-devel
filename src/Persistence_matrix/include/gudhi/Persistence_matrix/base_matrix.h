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
	using Field_element_type = typename Master_matrix::Field_type;
	using Column_type = typename Master_matrix::Column_type;
	using container_type = typename Master_matrix::boundary_type;
	using Row_type = typename Master_matrix::Row_type;

	Base_matrix();
	template<class Container_type = container_type>
	Base_matrix(const std::vector<Container_type>& columns);
	Base_matrix(unsigned int numberOfColumns);
	Base_matrix(const Base_matrix& matrixToCopy);
	Base_matrix(Base_matrix&& other) noexcept;

	template<class Container_type = container_type>
	void insert_column(const Container_type& column);
	template<class Container_type = container_type>
	void insert_column(const Container_type& column, int columnIndex);
	template<class Boundary_type>
	void insert_boundary(const Boundary_type& boundary);	//same as insert_column
	Column_type& get_column(index columnIndex);
	const Column_type& get_column(index columnIndex) const;
	//get_row(rowIndex) --> simplex ID (=/= columnIndex)
	Row_type& get_row(index rowIndex);
	const Row_type& get_row(index rowIndex) const;
	void erase_column(index columnIndex);
	void erase_row(index rowIndex);		//assumes the row is empty, just thought as index a cleanup

	unsigned int get_number_of_columns() const;

	// void add_to(index sourceColumnIndex, index targetColumnIndex);
	template<typename Index_type>
	std::enable_if_t<std::is_integral_v<Index_type> > add_to(Index_type sourceColumnIndex, Index_type targetColumnIndex);
	template<class Cell_range>
	std::enable_if_t<!std::is_integral_v<Cell_range> > add_to(const Cell_range& sourceColumn, index targetColumnIndex);
	template<class Cell_range>
	void add_to(const Cell_range& sourceColumn, const Field_element_type& coefficient, index targetColumnIndex);
	template<class Cell_range>
	void add_to(const Field_element_type& coefficient, const Cell_range& sourceColumn, index targetColumnIndex);
	//necessary because of vector columns ? TODO: base vector with naive cell clear ?
	void add_to(Column_type& sourceColumn, index targetColumnIndex);
	void add_to(Column_type& sourceColumn, const Field_element_type& coefficient, index targetColumnIndex);
	void add_to(const Field_element_type& coefficient, Column_type& sourceColumn, index targetColumnIndex);

	void zero_cell(index columnIndex, index rowIndex);
	void zero_column(index columnIndex);
	bool is_zero_cell(index columnIndex, index rowIndex) const;
	bool is_zero_column(index columnIndex);

	Base_matrix& operator=(const Base_matrix& other);
	friend void swap(Base_matrix& matrix1, Base_matrix& matrix2){
		swap(static_cast<typename Master_matrix::template Base_swap_option<Base_matrix<Master_matrix> >&>(matrix1),
			 static_cast<typename Master_matrix::template Base_swap_option<Base_matrix<Master_matrix> >&>(matrix2));
		matrix1.matrix_.swap(matrix2.matrix_);
		std::swap(matrix1.nextInsertIndex_, matrix2.nextInsertIndex_);

		if constexpr (Master_matrix::Option_list::has_row_access){
			swap(static_cast<typename Master_matrix::Matrix_row_access_option&>(matrix1),
				 static_cast<typename Master_matrix::Matrix_row_access_option&>(matrix2));
			if constexpr (Master_matrix::Option_list::has_removable_columns){
				for (auto& p : matrix1.matrix_){
					p.second.set_rows(&matrix1.rows_);
				}
				for (auto& p : matrix2.matrix_){
					p.second.set_rows(&matrix2.rows_);
				}
			} else {
				for (auto& col : matrix1.matrix_){
					col.set_rows(&matrix1.rows_);
				}
				for (auto& col : matrix2.matrix_){
					col.set_rows(&matrix2.rows_);
				}
			}
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

	void _orderRowsIfNecessary();
};

template<class Master_matrix>
inline Base_matrix<Master_matrix>::Base_matrix()
	: swap_opt(),
	  ra_opt(),
	  nextInsertIndex_(0)
{}

template<class Master_matrix>
template<class Container_type>
inline Base_matrix<Master_matrix>::Base_matrix(const std::vector<Container_type> &columns)
	: swap_opt(columns.size()),
	  ra_opt(columns.size()),	//not ideal if max row index is much smaller than max column index, does that happen often?
	  matrix_(!Master_matrix::Option_list::has_removable_columns && Master_matrix::Option_list::has_row_access ? 0 : columns.size()),
	  nextInsertIndex_(columns.size())
{
	if constexpr (!Master_matrix::Option_list::has_removable_columns && Master_matrix::Option_list::has_row_access)
		matrix_.reserve(columns.size());

	for (unsigned int i = 0; i < columns.size(); i++){
		if constexpr (Master_matrix::Option_list::has_removable_columns){
			if constexpr (Master_matrix::Option_list::has_row_access){
				matrix_.try_emplace(i, Column_type(i, columns[i], ra_opt::rows_));
			} else {
				matrix_.try_emplace(i, Column_type(columns[i]));
			}
		} else {
			if constexpr (Master_matrix::Option_list::has_row_access){
				matrix_.emplace_back(i, columns[i], ra_opt::rows_);
			} else {
				matrix_[i] = Column_type(columns[i]);
			}
		}
	}
}

template<class Master_matrix>
inline Base_matrix<Master_matrix>::Base_matrix(unsigned int numberOfColumns)
	: swap_opt(numberOfColumns),
	  ra_opt(numberOfColumns),
	  matrix_(!Master_matrix::Option_list::has_removable_columns && Master_matrix::Option_list::has_row_access ? 0 : numberOfColumns),
	  nextInsertIndex_(0)
{
	if constexpr (!Master_matrix::Option_list::has_removable_columns && Master_matrix::Option_list::has_row_access)
		matrix_.reserve(numberOfColumns);
}

template<class Master_matrix>
inline Base_matrix<Master_matrix>::Base_matrix(const Base_matrix &matrixToCopy)
	: swap_opt(static_cast<const swap_opt&>(matrixToCopy)),
	  ra_opt(static_cast<const ra_opt&>(matrixToCopy)),
	  nextInsertIndex_(matrixToCopy.nextInsertIndex_)
{
	if constexpr (Master_matrix::Option_list::has_row_access){
		matrix_.reserve(matrixToCopy.matrix_.size());
		if constexpr (Master_matrix::Option_list::has_removable_columns){
			for (const auto& p : matrixToCopy.matrix_){
				const Column_type& col = p.second;
				matrix_.try_emplace(p.first, Column_type(col, col.get_column_index(), ra_opt::rows_));
			}
		} else {
			for (const auto& col : matrixToCopy.matrix_){
				matrix_.emplace_back(col, col.get_column_index(), ra_opt::rows_);
			}
		}
	} else {
		matrix_ = matrix_type(matrixToCopy.matrix_);
	}
}

template<class Master_matrix>
inline Base_matrix<Master_matrix>::Base_matrix(Base_matrix &&other) noexcept
	: swap_opt(std::move(static_cast<swap_opt&>(other))),
	  ra_opt(std::move(static_cast<ra_opt&>(other))),
	  matrix_(std::move(other.matrix_)),
	  nextInsertIndex_(std::exchange(other.nextInsertIndex_, 0))
{
	if constexpr (Master_matrix::Option_list::has_row_access){
		if constexpr (Master_matrix::Option_list::has_removable_columns){
			for (auto& p : matrix_){
				p.second.set_rows(&this->rows_);
			}
		} else {
			for (auto& col : matrix_){
				col.set_rows(&this->rows_);
			}
		}
	}
}

template<class Master_matrix>
template<class Container_type>
inline void Base_matrix<Master_matrix>::insert_column(const Container_type &column)
{
	_orderRowsIfNecessary();

	if constexpr (Master_matrix::Option_list::has_row_access && !Master_matrix::Option_list::has_removable_rows){
		unsigned int pivot;
		if constexpr (Master_matrix::Option_list::is_z2){
			pivot = *std::prev(column.end());
		} else {
			pivot = std::prev(column.end())->first;
		}
		if (ra_opt::rows_.size() <= pivot) ra_opt::rows_.resize(pivot + 1);
	}
	
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		if constexpr (Master_matrix::Option_list::has_column_and_row_swaps){
			swap_opt::indexToRow_[nextInsertIndex_] = nextInsertIndex_;
			swap_opt::rowToIndex_[nextInsertIndex_] = nextInsertIndex_;
		}

		if constexpr (Master_matrix::Option_list::has_row_access){
			matrix_.try_emplace(nextInsertIndex_, Column_type(nextInsertIndex_, column, ra_opt::rows_));
			++nextInsertIndex_;
		} else {
			matrix_.try_emplace(nextInsertIndex_++, column);
		}
	} else {
		if constexpr (Master_matrix::Option_list::has_row_access){
			matrix_.emplace_back(nextInsertIndex_++, column, ra_opt::rows_);
		} else {
			unsigned int size = matrix_.size();
			if (size <= nextInsertIndex_) {
				if constexpr (Master_matrix::Option_list::has_column_and_row_swaps){
					for (unsigned int i = size; i <= size * 2; i++){
						swap_opt::indexToRow_.push_back(i);
						swap_opt::rowToIndex_.push_back(i);
					}
				}
				matrix_.resize(size * 2);
			}
			matrix_[nextInsertIndex_++] = Column_type(column);
		}
	}
}

template<class Master_matrix>
template<class Container_type>
inline void Base_matrix<Master_matrix>::insert_column(const Container_type &column, int columnIndex)
{
	static_assert(Master_matrix::Option_list::has_removable_columns, "Specification of the column index only possible for removable columns.");

	index id = columnIndex < 0 ? nextInsertIndex_ : columnIndex;
	assert(matrix_.find(id) == matrix_.end() && "Column already existing at given index.");
	if (columnIndex > static_cast<int>(nextInsertIndex_)) nextInsertIndex_ = columnIndex + 1;
	else if (columnIndex < 0) nextInsertIndex_++;

	if constexpr (Master_matrix::Option_list::has_column_and_row_swaps){
		if (swap_opt::rowSwapped_) swap_opt::_orderRows();
		swap_opt::indexToRow_[id] = id;
		swap_opt::rowToIndex_[id] = id;
	}

	if constexpr (Master_matrix::Option_list::has_row_access){
		if constexpr (!Master_matrix::Option_list::has_removable_rows){
			unsigned int pivot;
			if constexpr (Master_matrix::Option_list::is_z2){
				pivot = *std::prev(column.end());
			} else {
				pivot = std::prev(column.end())->first;
			}
			if (ra_opt::rows_.size() <= pivot) ra_opt::rows_.resize(pivot + 1);
		}
		matrix_.try_emplace(id, Column_type(id, column, ra_opt::rows_));
	} else {
		matrix_.try_emplace(id, column);
	}
}

template<class Master_matrix>
template<class Boundary_type>
inline void Base_matrix<Master_matrix>::insert_boundary(const Boundary_type &boundary)
{
	insert_column(boundary);
}

template<class Master_matrix>
inline typename Base_matrix<Master_matrix>::Column_type& Base_matrix<Master_matrix>::get_column(index columnIndex)
{
	_orderRowsIfNecessary();
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		return matrix_.at(columnIndex);
	} else {
		return matrix_[columnIndex];
	}
}

template<class Master_matrix>
inline const typename Base_matrix<Master_matrix>::Column_type& Base_matrix<Master_matrix>::get_column(index columnIndex) const
{
	_orderRowsIfNecessary();
	if constexpr (Master_matrix::Option_list::has_removable_columns){
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
inline const typename Base_matrix<Master_matrix>::Row_type &Base_matrix<Master_matrix>::get_row(index rowIndex) const
{
	static_assert(Master_matrix::Option_list::has_row_access, "Row access has to be enabled for this method.");

	_orderRowsIfNecessary();
	return ra_opt::get_row(rowIndex);
}

template<class Master_matrix>
inline void Base_matrix<Master_matrix>::erase_column(index columnIndex)
{
	static_assert(Master_matrix::Option_list::has_removable_columns,
			"'erase_column' is not imatrix_.erase(columnIndex);mplemented for the chosen options.");
	
	if (columnIndex == nextInsertIndex_ - 1) --nextInsertIndex_;

	matrix_.erase(columnIndex);
}

template<class Master_matrix>
inline void Base_matrix<Master_matrix>::erase_row(index rowIndex)
{
	if constexpr (Master_matrix::Option_list::has_column_and_row_swaps){
		if constexpr (Master_matrix::Option_list::has_row_access && Master_matrix::Option_list::has_removable_rows){
			ra_opt::erase_row(swap_opt::indexToRow_[rowIndex]);
		}
		
		auto it = swap_opt::indexToRow_.find(rowIndex);
		swap_opt::rowToIndex_.erase(it->second);
		swap_opt::indexToRow_.erase(it);
	} else {
		if constexpr (Master_matrix::Option_list::has_row_access && Master_matrix::Option_list::has_removable_rows){
			ra_opt::erase_row(rowIndex);
		}
	}
}

template<class Master_matrix>
inline unsigned int Base_matrix<Master_matrix>::get_number_of_columns() const
{
	if constexpr (Master_matrix::Option_list::has_removable_columns){
//		return nextInsertIndex_;	//if erased columns are viewed as zero columns, otherwise use matrix size.
		return matrix_.size();
	} else {
		return nextInsertIndex_;	//matrix could have been resized much bigger while insert
	}
}

// template<class Master_matrix>
// inline void Base_matrix<Master_matrix>::add_to(index sourceColumnIndex, index targetColumnIndex)
// {
// 	if constexpr (Master_matrix::Option_list::has_removable_columns){
// 		matrix_.at(targetColumnIndex) += matrix_.at(sourceColumnIndex);
// 	} else {
// 		matrix_[targetColumnIndex] += matrix_[sourceColumnIndex];
// 	}
// }

template<class Master_matrix>
template<typename Index_type>
inline std::enable_if_t<std::is_integral_v<Index_type> > Base_matrix<Master_matrix>::add_to(Index_type sourceColumnIndex, Index_type targetColumnIndex)
{
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		matrix_.at(targetColumnIndex) += matrix_.at(sourceColumnIndex);
	} else {
		matrix_[targetColumnIndex] += matrix_[sourceColumnIndex];
	}
}

template<class Master_matrix>
template<class Cell_range>
inline std::enable_if_t<!std::is_integral_v<Cell_range> > Base_matrix<Master_matrix>::add_to(const Cell_range& sourceColumn, index targetColumnIndex)
{
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		matrix_.at(targetColumnIndex) += sourceColumn;
	} else {
		matrix_[targetColumnIndex] += sourceColumn;
	}
}

template<class Master_matrix>
template<class Cell_range>
inline void Base_matrix<Master_matrix>::add_to(const Cell_range& sourceColumn, const Field_element_type& coefficient, index targetColumnIndex)
{
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		matrix_.at(targetColumnIndex).multiply_and_add(coefficient, sourceColumn);
	} else {
		matrix_[targetColumnIndex].multiply_and_add(coefficient, sourceColumn);
	}
}

template<class Master_matrix>
template<class Cell_range>
inline void Base_matrix<Master_matrix>::add_to(const Field_element_type& coefficient, const Cell_range& sourceColumn, index targetColumnIndex)
{
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		matrix_.at(targetColumnIndex).multiply_and_add(sourceColumn, coefficient);
	} else {
		matrix_[targetColumnIndex].multiply_and_add(sourceColumn, coefficient);
	}
}

template<class Master_matrix>
inline void Base_matrix<Master_matrix>::add_to(Column_type& sourceColumn, index targetColumnIndex)
{
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		matrix_.at(targetColumnIndex) += sourceColumn;
	} else {
		matrix_[targetColumnIndex] += sourceColumn;
	}
}

template<class Master_matrix>
inline void Base_matrix<Master_matrix>::add_to(Column_type& sourceColumn, const Field_element_type& coefficient, index targetColumnIndex)
{
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		matrix_.at(targetColumnIndex).multiply_and_add(coefficient, sourceColumn);
	} else {
		matrix_[targetColumnIndex].multiply_and_add(coefficient, sourceColumn);
	}
}

template<class Master_matrix>
inline void Base_matrix<Master_matrix>::add_to(const Field_element_type& coefficient, Column_type& sourceColumn, index targetColumnIndex)
{
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		matrix_.at(targetColumnIndex).multiply_and_add(sourceColumn, coefficient);
	} else {
		matrix_[targetColumnIndex].multiply_and_add(sourceColumn, coefficient);
	}
}

template<class Master_matrix>
inline void Base_matrix<Master_matrix>::zero_cell(index columnIndex, index rowIndex)
{
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		if constexpr (Master_matrix::Option_list::has_column_and_row_swaps){
			matrix_.at(columnIndex).clear(swap_opt::indexToRow_.at(rowIndex));	//TODO: reorganize columns
		} else {
			matrix_.at(columnIndex).clear(rowIndex);	//TODO: reorganize columns
		}
	} else {
		if constexpr (Master_matrix::Option_list::has_column_and_row_swaps){
			matrix_[columnIndex].clear(swap_opt::indexToRow_[rowIndex]);	//TODO: reorganize columns
		} else {
			matrix_[columnIndex].clear(rowIndex);	//TODO: reorganize columns
		}
	}
}

template<class Master_matrix>
inline void Base_matrix<Master_matrix>::zero_column(index columnIndex)
{
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		matrix_.at(columnIndex).clear();	//TODO: reorganize columns
	} else {
		matrix_[columnIndex].clear();	//TODO: reorganize columns
	}
}

template<class Master_matrix>
inline bool Base_matrix<Master_matrix>::is_zero_cell(index columnIndex, index rowIndex) const
{
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		if constexpr (Master_matrix::Option_list::has_column_and_row_swaps){
			return !(matrix_.at(columnIndex).is_non_zero(swap_opt::indexToRow_.at(rowIndex)));
		} else {
			return !(matrix_.at(columnIndex).is_non_zero(rowIndex));
		}
	} else {	//operator[] non const for maps
		if constexpr (Master_matrix::Option_list::has_column_and_row_swaps){
			return !(matrix_[columnIndex].is_non_zero(swap_opt::indexToRow_[rowIndex]));
		} else {
			return !(matrix_[columnIndex].is_non_zero(rowIndex));
		}
	}
}

template<class Master_matrix>
inline bool Base_matrix<Master_matrix>::is_zero_column(index columnIndex)
{
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		return matrix_.at(columnIndex).is_empty();
	} else {
		return matrix_[columnIndex].is_empty();
	}
}

template<class Master_matrix>
inline Base_matrix<Master_matrix> &Base_matrix<Master_matrix>::operator=(const Base_matrix& other)
{
	swap_opt::operator=(other);
	nextInsertIndex_ = other.nextInsertIndex_;

	if constexpr (Master_matrix::Option_list::has_row_access){
		ra_opt::operator=(other);
		matrix_.reserve(other.matrix_.size());
		if constexpr (Master_matrix::Option_list::has_removable_columns){
			for (const auto& p : other.matrix_){
				const Column_type& col = p.second;
				matrix_.try_emplace(p.first, Column_type(col, col.get_column_index(), ra_opt::rows_));
			}
		} else {
			for (const auto& col : other.matrix_){
				matrix_.emplace_back(col, col.get_column_index(), ra_opt::rows_);
			}
		}
	} else {
		matrix_ = other.matrix_;
	}

	return *this;
}

template<class Master_matrix>
inline void Base_matrix<Master_matrix>::print()
{
	std::cout << "Base_matrix:\n";
	for (unsigned int i = 0; i < nextInsertIndex_; ++i){
		const Column_type& col = matrix_[i];
		for (const auto e : col.get_content(nextInsertIndex_)){
			if (e == 0u) std::cout << "- ";
			else std::cout << e << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";
	if constexpr (Master_matrix::Option_list::has_row_access){
		std::cout << "Row Matrix:\n";
		for (unsigned int i = 0; i < nextInsertIndex_; ++i){
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
inline void Base_matrix<Master_matrix>::_orderRowsIfNecessary(){
	if constexpr (Master_matrix::Option_list::has_column_and_row_swaps){
		if (swap_opt::rowSwapped_) swap_opt::_orderRows();
	}
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // PM_BASE_MATRIX_H

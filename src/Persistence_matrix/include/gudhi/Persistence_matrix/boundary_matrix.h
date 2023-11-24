/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022-23 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_BOUNDARY_MATRIX_H
#define PM_BOUNDARY_MATRIX_H

#include <iostream>	//print() only
#include <vector>
#include <utility>	//std::swap, std::move & std::exchange

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class Boundary_matrix		//TODO: factorize/inheritate/compose with base matrix?
		: public Master_matrix::Matrix_dimension_option,
		  public Master_matrix::template Base_swap_option<Boundary_matrix<Master_matrix> >,
		  public Master_matrix::Base_pairing_option, 
		  public Master_matrix::Matrix_row_access_option
{
public:
	using index = typename Master_matrix::index;
	using dimension_type = typename Master_matrix::dimension_type;
	using Field_element_type = typename Master_matrix::Field_type;
	using Column_type = typename Master_matrix::Column_type;
	using boundary_type = typename Master_matrix::boundary_type;
	using Row_type = typename Master_matrix::Row_type;

	Boundary_matrix();
	template<class Boundary_type = boundary_type>
	Boundary_matrix(const std::vector<Boundary_type>& orderedBoundaries);
	Boundary_matrix(unsigned int numberOfColumns);
	Boundary_matrix(const Boundary_matrix& matrixToCopy);
	Boundary_matrix(Boundary_matrix&& other) noexcept;

	template<class Boundary_type = boundary_type>
	void insert_boundary(const Boundary_type& boundary, dimension_type dim = -1);	//does not update barcode as it needs reduction
	Column_type& get_column(index columnIndex);
	Row_type& get_row(index rowIndex);
	void remove_maximal_simplex(index columnIndex);		//update barcode if already computed
														//does not verify if it really was maximal

	unsigned int get_number_of_columns() const;

	dimension_type get_column_dimension(index columnIndex) const;

	//avoid calling with pairing option or make it such that it makes sense for persistence
	//=================================================================
	void add_to(index sourceColumnIndex, index targetColumnIndex);
	void add_to(index sourceColumnIndex, const Field_element_type& coefficient, index targetColumnIndex);
	void add_to(const Field_element_type& coefficient, index sourceColumnIndex, index targetColumnIndex);
	//TODO: are those other versions below really necessary for a boundary matrix?
	// void add_to(const Column_type& sourceColumn, index targetColumnIndex);
	// void add_to(const Column_type& sourceColumn, const Field_element_type& coefficient, index targetColumnIndex);
	// void add_to(const Field_element_type& coefficient, const Column_type& sourceColumn, index targetColumnIndex);

	void zero_cell(index columnIndex, index rowIndex);
	void zero_column(index columnIndex);
	//=================================================================
	bool is_zero_cell(index columnIndex, index rowIndex) const;
	bool is_zero_column(index columnIndex);

	int get_pivot(index columnIndex);

	Boundary_matrix& operator=(const Boundary_matrix& other);
	friend void swap(Boundary_matrix& matrix1, Boundary_matrix& matrix2){
		swap(static_cast<typename Master_matrix::Matrix_dimension_option&>(matrix1),
			 static_cast<typename Master_matrix::Matrix_dimension_option&>(matrix2));
		swap(static_cast<typename Master_matrix::template Base_swap_option<Boundary_matrix<Master_matrix> >&>(matrix1),
			 static_cast<typename Master_matrix::template Base_swap_option<Boundary_matrix<Master_matrix> >&>(matrix2));
		swap(static_cast<typename Master_matrix::Base_pairing_option&>(matrix1),
			 static_cast<typename Master_matrix::Base_pairing_option&>(matrix2));
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
	using dim_opt = typename Master_matrix::Matrix_dimension_option;
	using swap_opt = typename Master_matrix::template Base_swap_option<Boundary_matrix<Master_matrix> >;
	using pair_opt = typename Master_matrix::Base_pairing_option;
	using ra_opt = typename Master_matrix::Matrix_row_access_option;
	using matrix_type = typename Master_matrix::column_container_type;

	friend swap_opt;

	matrix_type matrix_;
	index nextInsertIndex_;

	static const bool activeDimOption = Master_matrix::Option_list::has_matrix_maximal_dimension_access || Master_matrix::dimensionIsNeeded;
	static const bool activeSwapOption = Master_matrix::Option_list::has_column_and_row_swaps || Master_matrix::Option_list::has_vine_update;
	static const bool activePairingOption = Master_matrix::Option_list::has_column_pairings && !Master_matrix::Option_list::has_vine_update && !Master_matrix::Option_list::can_retrieve_representative_cycles;
};

template<class Master_matrix>
inline Boundary_matrix<Master_matrix>::Boundary_matrix()
	: dim_opt(-1),
	  swap_opt(),
	  pair_opt(),
	  ra_opt(),
	  nextInsertIndex_(0)
{}

template<class Master_matrix>
template<class Boundary_type>
inline Boundary_matrix<Master_matrix>::Boundary_matrix(const std::vector<Boundary_type> &orderedBoundaries)
	: dim_opt(-1),
	  swap_opt(orderedBoundaries.size()),
	  pair_opt(),
	  ra_opt(orderedBoundaries.size()),
	  nextInsertIndex_(orderedBoundaries.size())
{
	matrix_.reserve(orderedBoundaries.size());

	for (unsigned int i = 0; i < orderedBoundaries.size(); i++){
		if constexpr (Master_matrix::Option_list::has_removable_columns){
			if constexpr (Master_matrix::Option_list::has_row_access){
				matrix_.try_emplace(i, Column_type(i, orderedBoundaries[i], ra_opt::rows_));
			} else {
				matrix_.try_emplace(i, Column_type(orderedBoundaries[i]));
			}
			if constexpr (activeDimOption){
				dim_opt::update_up(matrix_.at(i).get_dimension());
			}
		} else {
			if constexpr (Master_matrix::Option_list::has_row_access){
				matrix_.emplace_back(i, orderedBoundaries[i], ra_opt::rows_);
			} else {
				matrix_.emplace_back(orderedBoundaries[i]);
			}
			if constexpr (activeDimOption){
				dim_opt::update_up(matrix_[i].get_dimension());
			}
		}
	}
}

template<class Master_matrix>
inline Boundary_matrix<Master_matrix>::Boundary_matrix(unsigned int numberOfColumns)
	: dim_opt(-1),
	  swap_opt(numberOfColumns),
	  pair_opt(),
	  ra_opt(numberOfColumns),
	  matrix_(!Master_matrix::Option_list::has_removable_columns && Master_matrix::Option_list::has_row_access ? 0 : numberOfColumns),
	  nextInsertIndex_(0)
{
	if constexpr (!Master_matrix::Option_list::has_removable_columns && Master_matrix::Option_list::has_row_access)
		matrix_.reserve(numberOfColumns);
}

template<class Master_matrix>
inline Boundary_matrix<Master_matrix>::Boundary_matrix(const Boundary_matrix &matrixToCopy)
	: dim_opt(static_cast<const dim_opt&>(matrixToCopy)),
	  swap_opt(static_cast<const swap_opt&>(matrixToCopy)),
	  pair_opt(static_cast<const pair_opt&>(matrixToCopy)),
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
inline Boundary_matrix<Master_matrix>::Boundary_matrix(Boundary_matrix &&other) noexcept
	: dim_opt(std::move(static_cast<dim_opt&>(other))),
	  swap_opt(std::move(static_cast<swap_opt&>(other))),
	  pair_opt(std::move(static_cast<pair_opt&>(other))),
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
template<class Boundary_type>
inline void Boundary_matrix<Master_matrix>::insert_boundary(const Boundary_type &boundary, dimension_type dim)
{
	if (dim == -1) dim = boundary.size() == 0 ? 0 : boundary.size() - 1;

	if constexpr (activeSwapOption){
		if (swap_opt::rowSwapped_) swap_opt::_orderRows();
	}

	if constexpr (Master_matrix::Option_list::has_row_access && !Master_matrix::Option_list::has_removable_rows){
		unsigned int pivot;
		if constexpr (Master_matrix::Option_list::is_z2){
			pivot = *std::prev(boundary.end());
		} else {
			pivot = std::prev(boundary.end())->first;
		}
		if (ra_opt::rows_.size() <= pivot) ra_opt::rows_.resize(pivot + 1);
	}
	
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		if constexpr (activeSwapOption){
			swap_opt::indexToRow_[nextInsertIndex_] = nextInsertIndex_;
			swap_opt::rowToIndex_[nextInsertIndex_] = nextInsertIndex_;
		}

		if constexpr (Master_matrix::Option_list::has_row_access){
			matrix_.try_emplace(nextInsertIndex_, Column_type(nextInsertIndex_, boundary, dim, ra_opt::rows_));
			++nextInsertIndex_;
		} else {
			matrix_.try_emplace(nextInsertIndex_++, boundary, dim);
		}
	} else {
		if constexpr (Master_matrix::Option_list::has_row_access){
			matrix_.emplace_back(nextInsertIndex_++, boundary, dim, ra_opt::rows_);
		} else {
			unsigned int size = matrix_.size();
			if (size <= nextInsertIndex_) {
				if constexpr (activeSwapOption){
					for (unsigned int i = size; i <= size * 2; i++){
						swap_opt::indexToRow_.push_back(i);
						swap_opt::rowToIndex_.push_back(i);
					}
				}
				matrix_.resize(size * 2);
			}
			matrix_[nextInsertIndex_++] = Column_type(boundary, dim);
		}
	}

	if constexpr (activeDimOption){
		dim_opt::update_up(boundary.size() == 0 ? 0 : boundary.size() - 1);
	}
}

template<class Master_matrix>
inline typename Boundary_matrix<Master_matrix>::Column_type &Boundary_matrix<Master_matrix>::get_column(index columnIndex)
{
	if constexpr (activeSwapOption){
		if (swap_opt::rowSwapped_) swap_opt::_orderRows();
	}

	if constexpr (Master_matrix::Option_list::has_removable_columns){
		return matrix_.at(columnIndex);
	} else {
		return matrix_[columnIndex];
	}
}

template<class Master_matrix>
inline typename Boundary_matrix<Master_matrix>::Row_type& Boundary_matrix<Master_matrix>::get_row(index rowIndex)
{
	static_assert(Master_matrix::Option_list::has_row_access,
			"'get_row' is not implemented for the chosen options.");

	if constexpr (activeSwapOption){
		if (swap_opt::rowSwapped_) swap_opt::_orderRows();
	}

	return ra_opt::get_row(rowIndex);
}

template<class Master_matrix>
inline void Boundary_matrix<Master_matrix>::remove_maximal_simplex(index columnIndex)
{
	static_assert(Master_matrix::Option_list::has_removable_columns,
			"'remove_maximal_simplex' is not implemented for the chosen options.");

	if (columnIndex == nextInsertIndex_ - 1) --nextInsertIndex_;

	if constexpr (activeDimOption){
		dim_opt::update_down(matrix_.at(columnIndex).get_dimension());
	}

	matrix_.erase(columnIndex);

	if constexpr (Master_matrix::Option_list::has_row_access && Master_matrix::Option_list::has_removable_rows){
		ra_opt::erase_row(columnIndex);		//maximal, so empty
	}

	if constexpr (activeSwapOption){
		auto it = swap_opt::indexToRow_.find(columnIndex);
		swap_opt::rowToIndex_.erase(it->second);
		swap_opt::indexToRow_.erase(it);
	}
	if constexpr (activePairingOption){
		pair_opt::_remove_maximal(columnIndex);
	}
}

template<class Master_matrix>
inline unsigned int Boundary_matrix<Master_matrix>::get_number_of_columns() const
{
	if constexpr (Master_matrix::Option_list::has_removable_columns){
//		return nextInsertIndex_;	//if erased columns are viewed as zero columns, otherwise use matrix size.
		return matrix_.size();
	} else {
		return nextInsertIndex_;	//matrix could have been resized much bigger while insert
	}
}

template<class Master_matrix>
inline typename Boundary_matrix<Master_matrix>::dimension_type Boundary_matrix<Master_matrix>::get_column_dimension(index columnIndex) const
{
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		return matrix_.at(columnIndex).get_dimension();
	} else {
		return matrix_[columnIndex].get_dimension();
	}
}

template<class Master_matrix>
inline void Boundary_matrix<Master_matrix>::add_to(index sourceColumnIndex, index targetColumnIndex)
{
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		matrix_.at(targetColumnIndex) += matrix_.at(sourceColumnIndex);
	} else {
		matrix_[targetColumnIndex] += matrix_[sourceColumnIndex];
	}
}

template<class Master_matrix>
inline void Boundary_matrix<Master_matrix>::add_to(index sourceColumnIndex, const Field_element_type& coefficient, index targetColumnIndex)
{
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		matrix_.at(targetColumnIndex).multiply_and_add(coefficient, matrix_.at(sourceColumnIndex));
	} else {
		matrix_[targetColumnIndex].multiply_and_add(coefficient, matrix_[sourceColumnIndex]);
	}
}

template<class Master_matrix>
inline void Boundary_matrix<Master_matrix>::add_to(const Field_element_type& coefficient, index sourceColumnIndex, index targetColumnIndex)
{
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		matrix_.at(targetColumnIndex).multiply_and_add(matrix_.at(sourceColumnIndex), coefficient);
	} else {
		matrix_[targetColumnIndex].multiply_and_add(matrix_[sourceColumnIndex], coefficient);
	}
}

// template<class Master_matrix>
// inline void Boundary_matrix<Master_matrix>::add_to(Column_type& sourceColumn, index targetColumnIndex)
// {
// 	if constexpr (Master_matrix::Option_list::has_removable_columns){
// 		matrix_.at(targetColumnIndex) += sourceColumn;
// 	} else {
// 		matrix_[targetColumnIndex] += sourceColumn;
// 	}
// }

// template<class Master_matrix>
// inline void Boundary_matrix<Master_matrix>::add_to(Column_type& sourceColumn, const Field_element_type& coefficient, index targetColumnIndex)
// {
// 	if constexpr (Master_matrix::Option_list::has_removable_columns){
// 		matrix_.at(targetColumnIndex).multiply_and_add(coefficient, sourceColumn);
// 	} else {
// 		matrix_[targetColumnIndex].multiply_and_add(coefficient, sourceColumn);
// 	}
// }

// template<class Master_matrix>
// inline void Boundary_matrix<Master_matrix>::add_to(const Field_element_type& coefficient, Column_type& sourceColumn, index targetColumnIndex)
// {
// 	if constexpr (Master_matrix::Option_list::has_removable_columns){
// 		matrix_.at(targetColumnIndex).multiply_and_add(sourceColumn, coefficient);
// 	} else {
// 		matrix_[targetColumnIndex].multiply_and_add(sourceColumn, coefficient);
// 	}
// }

// template<class Master_matrix>
// inline void Boundary_matrix<Master_matrix>::add_to(const Column_type& sourceColumn, index targetColumnIndex)
// {
// 	if constexpr (Master_matrix::Option_list::has_removable_columns){
// 		matrix_.at(targetColumnIndex) += sourceColumn;
// 	} else {
// 		matrix_[targetColumnIndex] += sourceColumn;
// 	}
// }

// template<class Master_matrix>
// inline void Boundary_matrix<Master_matrix>::add_to(const Column_type& sourceColumn, const Field_element_type& coefficient, index targetColumnIndex)
// {
// 	if constexpr (Master_matrix::Option_list::has_removable_columns){
// 		matrix_.at(targetColumnIndex).multiply_and_add(coefficient, sourceColumn);
// 	} else {
// 		matrix_[targetColumnIndex].multiply_and_add(coefficient, sourceColumn);
// 	}
// }

// template<class Master_matrix>
// inline void Boundary_matrix<Master_matrix>::add_to(const Field_element_type& coefficient, const Column_type& sourceColumn, index targetColumnIndex)
// {
// 	if constexpr (Master_matrix::Option_list::has_removable_columns){
// 		matrix_.at(targetColumnIndex).multiply_and_add(sourceColumn, coefficient);
// 	} else {
// 		matrix_[targetColumnIndex].multiply_and_add(sourceColumn, coefficient);
// 	}
// }

template<class Master_matrix>
inline void Boundary_matrix<Master_matrix>::zero_cell(index columnIndex, index rowIndex)
{
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		if constexpr (activeSwapOption){
			matrix_.at(columnIndex).clear(swap_opt::indexToRow_[rowIndex]);
		} else {
			matrix_.at(columnIndex).clear(rowIndex);
		}
	} else {
		if constexpr (activeSwapOption){
			matrix_[columnIndex].clear(swap_opt::indexToRow_[rowIndex]);
		} else {
			matrix_[columnIndex].clear(rowIndex);
		}
	}
}

template<class Master_matrix>
inline void Boundary_matrix<Master_matrix>::zero_column(index columnIndex)
{
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		matrix_.at(columnIndex).clear();
	} else {
		matrix_[columnIndex].clear();
	}
}

template<class Master_matrix>
inline bool Boundary_matrix<Master_matrix>::is_zero_cell(index columnIndex, index rowIndex) const
{
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		if constexpr (activeSwapOption){
			return !(matrix_.at(columnIndex).is_non_zero(swap_opt::indexToRow_.at(rowIndex)));
		} else {
			return !(matrix_.at(columnIndex).is_non_zero(rowIndex));
		}
	} else {	//operator[] non const for maps
		if constexpr (activeSwapOption){
			return !(matrix_[columnIndex].is_non_zero(swap_opt::indexToRow_[rowIndex]));
		} else {
			return !(matrix_[columnIndex].is_non_zero(rowIndex));
		}
	}
}

template<class Master_matrix>
inline bool Boundary_matrix<Master_matrix>::is_zero_column(index columnIndex)
{
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		return matrix_.at(columnIndex).is_empty();
	} else {
		return matrix_[columnIndex].is_empty();
	}
}

template<class Master_matrix>
inline int Boundary_matrix<Master_matrix>::get_pivot(index columnIndex)
{
	if constexpr (activeSwapOption){
		if (swap_opt::rowSwapped_) swap_opt::_orderRows();
	}
	
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		return matrix_.at(columnIndex).get_pivot();
	} else {
		return matrix_[columnIndex].get_pivot();
	}
}

template<class Master_matrix>
inline Boundary_matrix<Master_matrix> &Boundary_matrix<Master_matrix>::operator=(const Boundary_matrix& other)
{
	dim_opt::operator=(other);
	swap_opt::operator=(other);
	pair_opt::operator=(other);

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
inline void Boundary_matrix<Master_matrix>::print()
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

} //namespace persistence_matrix
} //namespace Gudhi

#endif // PM_BOUNDARY_MATRIX_H

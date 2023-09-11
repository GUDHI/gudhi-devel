/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022-23 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_RU_MATRIX_H
#define PM_RU_MATRIX_H

#include <vector>
#include <utility>	//std::exchange
#include <iostream>	//print() only

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class RU_matrix;

template<class Master_matrix>
class RU_matrix
		: public Master_matrix::RU_pairing_option,
		  public Master_matrix::RU_vine_swap_option,
		  public Master_matrix::RU_representative_cycles_option
{
public:
	using Field_element_type = typename Master_matrix::Field_type;
	using Column_type = typename Master_matrix::Column_type;
	using Row_type = typename Master_matrix::Row_type;
	using boundary_type = typename Master_matrix::boundary_type;
	using index = typename Master_matrix::index;
	using dimension_type = typename Master_matrix::dimension_type;

	RU_matrix();
	template<class Boundary_type = boundary_type>
	RU_matrix(const std::vector<Boundary_type>& orderedBoundaries);
	RU_matrix(unsigned int numberOfColumns);
	RU_matrix(const RU_matrix& matrixToCopy);
	RU_matrix(RU_matrix&& other) noexcept;

	template<class Boundary_type = boundary_type>
	void insert_boundary(const Boundary_type& boundary);
	Column_type& get_column(index columnIndex, bool inR = true);
	const Column_type& get_column(index columnIndex, bool inR = true) const;
	//get_row(rowIndex) --> simplex ID (=/= columnIndex)
	Row_type& get_row(index rowIndex, bool inR = true);
	const Row_type& get_row(index rowIndex, bool inR = true) const;
	void remove_maximal_simplex(index columnIndex);

	dimension_type get_max_dimension() const;
	unsigned int get_number_of_columns() const;

	dimension_type get_column_dimension(index columnIndex) const;

	//avoid calling with specialized options or make it such that it makes sense for persistence
	//=================================================================
	void add_to(index sourceColumnIndex, index targetColumnIndex);
	void add_to(Column_type& sourceColumn, index targetColumnIndex);
	void add_to(Column_type& sourceColumn, const Field_element_type& coefficient, index targetColumnIndex);
	void add_to(const Field_element_type& coefficient, Column_type& sourceColumn, index targetColumnIndex);
	void add_to(const Column_type& sourceColumn, index targetColumnIndex);
	void add_to(const Column_type& sourceColumn, const Field_element_type& coefficient, index targetColumnIndex);
	void add_to(const Field_element_type& coefficient, const Column_type& sourceColumn, index targetColumnIndex);

	void zero_cell(index columnIndex, index rowIndex, bool inR = true);
	void zero_column(index columnIndex, bool inR = true);
	//=================================================================
	bool is_zero_cell(index columnIndex, index rowIndex, bool inR = true) const;
	bool is_zero_column(index columnIndex, bool inR = true);

	index get_column_with_pivot(index simplexIndex) const;	//assumes that pivot exists
	index get_pivot(index columnIndex);

	RU_matrix& operator=(RU_matrix other);
	friend void swap(RU_matrix& matrix1, RU_matrix& matrix2){
		swap(static_cast<typename Master_matrix::RU_pairing_option&>(matrix1),
			 static_cast<typename Master_matrix::RU_pairing_option&>(matrix2));
		swap(static_cast<typename Master_matrix::RU_vine_swap_option&>(matrix1),
			 static_cast<typename Master_matrix::RU_vine_swap_option&>(matrix2));
		swap(static_cast<typename Master_matrix::RU_representative_cycles_option&>(matrix1),
			 static_cast<typename Master_matrix::RU_representative_cycles_option&>(matrix2));
		swap(matrix1.reducedMatrixR_, matrix2.reducedMatrixR_);
		swap(matrix1.mirrorMatrixU_, matrix2.mirrorMatrixU_);
		matrix1.pivotToColumnIndex_.swap(matrix2.pivotToColumnIndex_);
		std::swap(matrix1.nextInsertIndex_, matrix2.nextInsertIndex_);
	}

	void print();  //for debug

private:
	using swap_opt = typename Master_matrix::RU_vine_swap_option;
	using pair_opt = typename Master_matrix::RU_pairing_option;
	using rep_opt = typename Master_matrix::RU_representative_cycles_option;
	using dictionnary_type = typename Master_matrix::template dictionnary_type<int>;
	using barcode_type = typename Master_matrix::barcode_type;
	using bar_dictionnary_type = typename Master_matrix::bar_dictionnary_type;
	using matrix_type = typename Master_matrix::Boundary_matrix_type;

	matrix_type reducedMatrixR_;
	matrix_type mirrorMatrixU_;	//make U not accessible by default and add option to enable access? Inaccessible, it needs less options and we could avoid some ifs.
	dictionnary_type pivotToColumnIndex_;
	index nextInsertIndex_;

	void _initialize_U();
	void _reduce();
	void _reduce_last_column();

	constexpr barcode_type& _barcode();
	constexpr bar_dictionnary_type& _indexToBar();
};

template<class Master_matrix>
inline RU_matrix<Master_matrix>::RU_matrix()
	: Master_matrix::RU_pairing_option(),
	  Master_matrix::RU_vine_swap_option(reducedMatrixR_, mirrorMatrixU_, pivotToColumnIndex_),
	  Master_matrix::RU_representative_cycles_option(reducedMatrixR_, mirrorMatrixU_),
	  nextInsertIndex_(0)
{}

template<class Master_matrix>
template<class Boundary_type>
inline RU_matrix<Master_matrix>::RU_matrix(const std::vector<Boundary_type> &orderedBoundaries)
	: Master_matrix::RU_pairing_option(),
	  Master_matrix::RU_vine_swap_option(reducedMatrixR_, mirrorMatrixU_, pivotToColumnIndex_),
	  Master_matrix::RU_representative_cycles_option(reducedMatrixR_, mirrorMatrixU_),
	  reducedMatrixR_(orderedBoundaries),
	  mirrorMatrixU_(orderedBoundaries.size()),
	  nextInsertIndex_(orderedBoundaries.size())
{
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		pivotToColumnIndex_.reserve(orderedBoundaries.size());
	} else {
		pivotToColumnIndex_.resize(orderedBoundaries.size(), -1);
	}

	_initialize_U();
	_reduce();
}

template<class Master_matrix>
inline RU_matrix<Master_matrix>::RU_matrix(unsigned int numberOfColumns)
	: Master_matrix::RU_pairing_option(),
	  Master_matrix::RU_vine_swap_option(reducedMatrixR_, mirrorMatrixU_, pivotToColumnIndex_),
	  Master_matrix::RU_representative_cycles_option(reducedMatrixR_, mirrorMatrixU_),
	  reducedMatrixR_(numberOfColumns),
	  mirrorMatrixU_(numberOfColumns),
	  nextInsertIndex_(0)
{
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		pivotToColumnIndex_.reserve(numberOfColumns);
	} else {
		pivotToColumnIndex_.resize(numberOfColumns, -1);
	}
}

template<class Master_matrix>
inline RU_matrix<Master_matrix>::RU_matrix(const RU_matrix &matrixToCopy)
	: Master_matrix::RU_pairing_option(matrixToCopy),
	  Master_matrix::RU_vine_swap_option(matrixToCopy),
	  Master_matrix::RU_representative_cycles_option(matrixToCopy),
	  reducedMatrixR_(matrixToCopy.reducedMatrixR_),
	  mirrorMatrixU_(matrixToCopy.mirrorMatrixU_),
	  pivotToColumnIndex_(matrixToCopy.pivotToColumnIndex_),
	  nextInsertIndex_(matrixToCopy.nextInsertIndex_)
{
	if constexpr (Master_matrix::Option_list::can_retrieve_representative_cycles){
		rep_opt::reducedMatrixR_ = &reducedMatrixR_;
		rep_opt::mirrorMatrixU_ = &mirrorMatrixU_;
	}
	if constexpr (Master_matrix::Option_list::has_vine_update){
		swap_opt::reducedMatrixR_ = &reducedMatrixR_;
		swap_opt::mirrorMatrixU_ = &mirrorMatrixU_;
		swap_opt::pivotToColumnIndex_ = &pivotToColumnIndex_;
	}
}

template<class Master_matrix>
inline RU_matrix<Master_matrix>::RU_matrix(RU_matrix &&other) noexcept
	: Master_matrix::RU_pairing_option(std::move(other)),
	  Master_matrix::RU_vine_swap_option(std::move(other)),
	  Master_matrix::RU_representative_cycles_option(std::move(other)),
	  reducedMatrixR_(std::move(other.reducedMatrixR_)),
	  mirrorMatrixU_(std::move(other.mirrorMatrixU_)),
	  pivotToColumnIndex_(std::move(other.pivotToColumnIndex_)),
	  nextInsertIndex_(std::exchange(other.nextInsertIndex_, 0))
{
	if constexpr (Master_matrix::Option_list::can_retrieve_representative_cycles){
		rep_opt::reducedMatrixR_ = &reducedMatrixR_;
		rep_opt::mirrorMatrixU_ = &mirrorMatrixU_;
	}
	if constexpr (Master_matrix::Option_list::has_vine_update){
		swap_opt::reducedMatrixR_ = &reducedMatrixR_;
		swap_opt::mirrorMatrixU_ = &mirrorMatrixU_;
		swap_opt::pivotToColumnIndex_ = &pivotToColumnIndex_;
	}
}

template<class Master_matrix>
template<class Boundary_type>
inline void RU_matrix<Master_matrix>::insert_boundary(const Boundary_type &boundary)
{
	reducedMatrixR_.insert_boundary(boundary);

	if constexpr (Master_matrix::Option_list::is_z2) {
		mirrorMatrixU_.insert_boundary({nextInsertIndex_});
	} else {
		mirrorMatrixU_.insert_boundary({{nextInsertIndex_, 1}});
	}

	if constexpr (!Master_matrix::Option_list::has_removable_columns){
		while (pivotToColumnIndex_.size() <= nextInsertIndex_)
			pivotToColumnIndex_.resize(pivotToColumnIndex_.size()*2, -1);
	}
	
	_reduce_last_column();
	++nextInsertIndex_;
}

template<class Master_matrix>
inline typename RU_matrix<Master_matrix>::Column_type &
RU_matrix<Master_matrix>::get_column(index columnIndex, bool inR)
{
	if (inR){
		return reducedMatrixR_.get_column(columnIndex);
	}
	return mirrorMatrixU_.get_column(columnIndex);
}

template<class Master_matrix>
inline const typename RU_matrix<Master_matrix>::Column_type &
RU_matrix<Master_matrix>::get_column(index columnIndex, bool inR) const
{
	if (inR){
		return reducedMatrixR_.get_column(columnIndex);
	}
	return mirrorMatrixU_.get_column(columnIndex);
}

template<class Master_matrix>
inline typename RU_matrix<Master_matrix>::Row_type&
RU_matrix<Master_matrix>::get_row(index rowIndex, bool inR)
{
	static_assert(Master_matrix::Option_list::has_row_access,
			"'get_row' is not implemented for the chosen options.");

	if (inR){
		return reducedMatrixR_.get_row(rowIndex);
	}
	return mirrorMatrixU_.get_row(rowIndex);
}

template<class Master_matrix>
inline const typename RU_matrix<Master_matrix>::Row_type&
RU_matrix<Master_matrix>::get_row(index rowIndex, bool inR) const
{
	static_assert(Master_matrix::Option_list::has_row_access,
			"'get_row' is not implemented for the chosen options.");

	if (inR){
		return reducedMatrixR_.get_row(rowIndex);
	}
	return mirrorMatrixU_.get_row(rowIndex);
}

template<class Master_matrix>
inline void RU_matrix<Master_matrix>::remove_maximal_simplex([[maybe_unused]] index columnIndex)
{
	static_assert(Master_matrix::Option_list::has_removable_columns,
			"'remove_maximal_simplex' is not implemented for the chosen options.");

	if (columnIndex == nextInsertIndex_ - 1) --nextInsertIndex_;

	if constexpr (Master_matrix::Option_list::has_column_pairings){
		typename barcode_type::iterator bar = _indexToBar().at(columnIndex);

		if (bar->death == -1) _barcode().erase(bar);
		else bar->death = -1;

		_indexToBar().erase(columnIndex);
	}

	pivotToColumnIndex_.erase(reducedMatrixR_.get_column(columnIndex).get_pivot());

	reducedMatrixR_.remove_maximal_simplex(columnIndex);
	mirrorMatrixU_.remove_maximal_simplex(columnIndex);
}

template<class Master_matrix>
inline typename RU_matrix<Master_matrix>::dimension_type RU_matrix<Master_matrix>::get_max_dimension() const
{
	return reducedMatrixR_.get_max_dimension();
}

template<class Master_matrix>
inline unsigned int RU_matrix<Master_matrix>::get_number_of_columns() const
{
	return reducedMatrixR_.get_number_of_columns();
}

template<class Master_matrix>
inline typename RU_matrix<Master_matrix>::dimension_type RU_matrix<Master_matrix>::get_column_dimension(index columnIndex) const
{
	return reducedMatrixR_.get_column_dimension(columnIndex);
}

template<class Master_matrix>
inline void RU_matrix<Master_matrix>::add_to(index sourceColumnIndex, index targetColumnIndex)
{
	reducedMatrixR_.add_to(sourceColumnIndex, targetColumnIndex);
	mirrorMatrixU_.add_to(sourceColumnIndex, targetColumnIndex);
}

template<class Master_matrix>
inline void RU_matrix<Master_matrix>::add_to(Column_type& sourceColumn, index targetColumnIndex)
{
	reducedMatrixR_.add_to(sourceColumn, targetColumnIndex);
	mirrorMatrixU_.add_to(sourceColumn, targetColumnIndex);
}

template<class Master_matrix>
inline void RU_matrix<Master_matrix>::add_to(Column_type& sourceColumn, const Field_element_type& coefficient, index targetColumnIndex)
{
	reducedMatrixR_.add_to(sourceColumn, coefficient, targetColumnIndex);
	mirrorMatrixU_.add_to(sourceColumn, coefficient, targetColumnIndex);
}

template<class Master_matrix>
inline void RU_matrix<Master_matrix>::add_to(const Field_element_type& coefficient, Column_type& sourceColumn, index targetColumnIndex)
{
	reducedMatrixR_.add_to(coefficient, sourceColumn, targetColumnIndex);
	mirrorMatrixU_.add_to(coefficient, sourceColumn, targetColumnIndex);
}

template<class Master_matrix>
inline void RU_matrix<Master_matrix>::add_to(const Column_type& sourceColumn, index targetColumnIndex)
{
	reducedMatrixR_.add_to(sourceColumn, targetColumnIndex);
	mirrorMatrixU_.add_to(sourceColumn, targetColumnIndex);
}

template<class Master_matrix>
inline void RU_matrix<Master_matrix>::add_to(const Column_type& sourceColumn, const Field_element_type& coefficient, index targetColumnIndex)
{
	reducedMatrixR_.add_to(sourceColumn, coefficient, targetColumnIndex);
	mirrorMatrixU_.add_to(sourceColumn, coefficient, targetColumnIndex);
}

template<class Master_matrix>
inline void RU_matrix<Master_matrix>::add_to(const Field_element_type& coefficient, const Column_type& sourceColumn, index targetColumnIndex)
{
	reducedMatrixR_.add_to(coefficient, sourceColumn, targetColumnIndex);
	mirrorMatrixU_.add_to(coefficient, sourceColumn, targetColumnIndex);
}

template<class Master_matrix>
inline void RU_matrix<Master_matrix>::zero_cell(index columnIndex, index rowIndex, bool inR)
{
	if (inR){
		return reducedMatrixR_.zero_cell(columnIndex, rowIndex);
	}
	return mirrorMatrixU_.zero_cell(columnIndex, rowIndex);
}

template<class Master_matrix>
inline void RU_matrix<Master_matrix>::zero_column(index columnIndex, bool inR)
{
	if (inR){
		return reducedMatrixR_.zero_column(columnIndex);
	}
	return mirrorMatrixU_.zero_column(columnIndex);
}

template<class Master_matrix>
inline bool RU_matrix<Master_matrix>::is_zero_cell(index columnIndex, index rowIndex, bool inR) const
{
	if (inR){
		return reducedMatrixR_.is_zero_cell(columnIndex, rowIndex);
	}
	return mirrorMatrixU_.is_zero_cell(columnIndex, rowIndex);
}

template<class Master_matrix>
inline bool RU_matrix<Master_matrix>::is_zero_column(index columnIndex, bool inR)
{
	if (inR){
		return reducedMatrixR_.is_zero_column(columnIndex);
	}
	return mirrorMatrixU_.is_zero_column(columnIndex);
}

template<class Master_matrix>
inline typename RU_matrix<Master_matrix>::index RU_matrix<Master_matrix>::get_column_with_pivot(index simplexIndex) const
{
	if constexpr (Master_matrix::Option_list::has_removable_columns){
		return pivotToColumnIndex_.at(simplexIndex);
	} else {
		return pivotToColumnIndex_[simplexIndex];
	}
}

template<class Master_matrix>
inline typename RU_matrix<Master_matrix>::index RU_matrix<Master_matrix>::get_pivot(index columnIndex)
{
	return reducedMatrixR_.get_column(columnIndex).get_pivot();
}

template<class Master_matrix>
inline RU_matrix<Master_matrix> &RU_matrix<Master_matrix>::operator=(RU_matrix other)
{
	swap_opt::operator=(other);
	pair_opt::operator=(other);
	rep_opt::operator=(other);
	swap(reducedMatrixR_, other.reducedMatrixR_);
	swap(mirrorMatrixU_, other.mirrorMatrixU_);
	pivotToColumnIndex_.swap(other.pivotToColumnIndex_);
	std::swap(nextInsertIndex_, other.nextInsertIndex_);
	return *this;
}

template<class Master_matrix>
inline void RU_matrix<Master_matrix>::print()
{
	std::cout << "R_matrix:\n";
	reducedMatrixR_.print();
	std::cout << "U_matrix:\n";
	mirrorMatrixU_.print();
}

template<class Master_matrix>
inline void RU_matrix<Master_matrix>::_initialize_U()
{
	typename std::conditional<Master_matrix::Option_list::is_z2, index, std::pair<index,Field_element_type> >::type id;
	if constexpr (!Master_matrix::Option_list::is_z2) id.second = 1;

	for (unsigned int i = 0; i < reducedMatrixR_.get_number_of_columns(); i++){
		if constexpr (Master_matrix::Option_list::is_z2) id = i;
		else id.first = i;
		mirrorMatrixU_.insert_boundary({id});
	}
}

template<class Master_matrix>
inline void RU_matrix<Master_matrix>::_reduce()
{
	auto get_column_with_pivot_ = [&](int pivot)->int{
		if (pivot == -1) return -1;
		if constexpr (Master_matrix::Option_list::has_removable_columns){
			auto it = pivotToColumnIndex_.find(pivot);
			if (it == pivotToColumnIndex_.end()) return -1;
			else return it->second;
		} else {
			return pivotToColumnIndex_[pivot];
		}
	};

	if constexpr (Master_matrix::Option_list::has_column_pairings){
		if constexpr (Master_matrix::Option_list::has_removable_columns){
			_indexToBar().reserve(reducedMatrixR_.get_number_of_columns());
		} else {
			_indexToBar().resize(reducedMatrixR_.get_number_of_columns(), -1);
		}
	}

	for (unsigned int i = 0; i < reducedMatrixR_.get_number_of_columns(); i++){
		if (!(reducedMatrixR_.is_zero_column(i)))
		{
			Column_type &curr = reducedMatrixR_.get_column(i);
			int pivot = curr.get_pivot();
			index currIndex = get_column_with_pivot_(pivot);

			while (pivot != -1 && currIndex != -1){
				if constexpr (Master_matrix::Option_list::is_z2){
					curr += reducedMatrixR_.get_column(currIndex);
					mirrorMatrixU_.get_column(i) += mirrorMatrixU_.get_column(currIndex);
				} else {
					Column_type &toadd = reducedMatrixR_.get_column(currIndex);
					typename Master_matrix::Field_type coef = curr.get_pivot_value();
					coef = coef.get_inverse();
					coef *= (Master_matrix::Field_type::get_characteristic() - static_cast<unsigned int>(toadd.get_pivot_value()));

					curr.multiply_and_add(coef, toadd);
					mirrorMatrixU_.get_column(i).multiply_and_add(coef, mirrorMatrixU_.get_column(currIndex));
				}

				pivot = curr.get_pivot();
				currIndex = get_column_with_pivot_(pivot);
			}

			if (pivot != -1){
				if constexpr (Master_matrix::Option_list::has_removable_columns){
					pivotToColumnIndex_.try_emplace(pivot, i);
					if constexpr (Master_matrix::Option_list::has_column_pairings){
						_indexToBar().at(pivot)->death = i;
						_indexToBar().try_emplace(i, _indexToBar().at(pivot));
					}
				} else {
					pivotToColumnIndex_[pivot] = i;
					if constexpr (Master_matrix::Option_list::has_column_pairings){
						_barcode()[_indexToBar()[pivot]].death = i;
						_indexToBar()[i] = _indexToBar()[pivot];
					}
				}
			} else if constexpr (Master_matrix::Option_list::has_column_pairings){
				if constexpr (Master_matrix::Option_list::has_removable_columns){
					_barcode().emplace_back(get_column_dimension(i), i, -1);
					_indexToBar().try_emplace(i, --_barcode().end());
				} else {
					_barcode().emplace_back(get_column_dimension(i), i, -1);
					_indexToBar()[i] = _barcode().size() - 1;
				}
			}
		} else if constexpr (Master_matrix::Option_list::has_column_pairings){
			if constexpr (Master_matrix::Option_list::has_removable_columns){
				_barcode().emplace_back(0, i, -1);
				_indexToBar().try_emplace(i, --_barcode().end());
			} else {
				_barcode().emplace_back(0, i, -1);
				_indexToBar()[i] = _barcode().size() - 1;
			}
		}
	}
}

template<class Master_matrix>
inline void RU_matrix<Master_matrix>::_reduce_last_column()
{
	auto get_column_with_pivot_ = [&](int pivot)->int{
		if (pivot == -1) return -1;
		if constexpr (Master_matrix::Option_list::has_removable_columns){
			auto it = pivotToColumnIndex_.find(pivot);
			if (it == pivotToColumnIndex_.end()) return -1;
			else return it->second;
		} else {
			return pivotToColumnIndex_[pivot];
		}
	};

	Column_type &curr = reducedMatrixR_.get_column(nextInsertIndex_);
	if (curr.is_empty()) {
		if constexpr (Master_matrix::Option_list::has_removable_columns){
			if constexpr (Master_matrix::Option_list::has_column_pairings){
				_barcode().emplace_back(0, nextInsertIndex_, -1);
				_indexToBar().try_emplace(nextInsertIndex_, --_barcode().end());
			}
		} else {
			if constexpr (Master_matrix::Option_list::has_column_pairings){
				_barcode().emplace_back(0, nextInsertIndex_, -1);
				_indexToBar().push_back(_barcode().size() - 1);
			}
		}
		return;
	}

	int pivot = curr.get_pivot();
	index currIndex = get_column_with_pivot_(pivot);

	while (pivot != -1 && currIndex != -1){
		if constexpr (Master_matrix::Option_list::is_z2){
			curr += reducedMatrixR_.get_column(currIndex);
			mirrorMatrixU_.get_column(nextInsertIndex_) += mirrorMatrixU_.get_column(currIndex);
		} else {
			Column_type &toadd = reducedMatrixR_.get_column(currIndex);
			typename Master_matrix::Field_type coef = curr.get_pivot_value();
			coef = coef.get_inverse();
			coef *= (Master_matrix::Field_type::get_characteristic() - static_cast<unsigned int>(toadd.get_pivot_value()));
			curr *= coef;
			curr += toadd;
			mirrorMatrixU_.get_column(nextInsertIndex_) *= coef;
			mirrorMatrixU_.get_column(nextInsertIndex_) += mirrorMatrixU_.get_column(currIndex);
		}

		pivot = curr.get_pivot();
		currIndex = get_column_with_pivot_(pivot);
	}

	if (pivot != -1){
		if constexpr (Master_matrix::Option_list::has_removable_columns){
			pivotToColumnIndex_.try_emplace(pivot, nextInsertIndex_);
			if constexpr (Master_matrix::Option_list::has_column_pairings){
				_indexToBar().at(pivot)->death = nextInsertIndex_;
				_indexToBar().try_emplace(nextInsertIndex_, _indexToBar().at(pivot));
			}
		} else {
			pivotToColumnIndex_[pivot] = nextInsertIndex_;
			if constexpr (Master_matrix::Option_list::has_column_pairings){
				_barcode()[_indexToBar()[pivot]].death = nextInsertIndex_;
				_indexToBar().push_back(_indexToBar()[pivot]);
			}
		}
	} else if constexpr (Master_matrix::Option_list::has_column_pairings){
		if constexpr (Master_matrix::Option_list::has_removable_columns){
			_barcode().emplace_back(get_column_dimension(nextInsertIndex_), nextInsertIndex_, -1);
			_indexToBar().try_emplace(nextInsertIndex_, --_barcode().end());
		} else {
			_barcode().emplace_back(get_column_dimension(nextInsertIndex_), nextInsertIndex_, -1);
			_indexToBar().push_back(_barcode().size() - 1);
		}
	}
}

template<class Master_matrix>
inline constexpr typename RU_matrix<Master_matrix>::barcode_type &RU_matrix<Master_matrix>::_barcode()
{
	if constexpr (Master_matrix::Option_list::has_vine_update)
		return swap_opt::template RU_pairing<Master_matrix>::barcode_;
	else
		return pair_opt::barcode_;
}

template<class Master_matrix>
inline constexpr typename RU_matrix<Master_matrix>::bar_dictionnary_type &RU_matrix<Master_matrix>::_indexToBar()
{
	if constexpr (Master_matrix::Option_list::has_vine_update)
		return swap_opt::template RU_pairing<Master_matrix>::indexToBar_;
	else
		return pair_opt::indexToBar_;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // PM_RU_MATRIX_H

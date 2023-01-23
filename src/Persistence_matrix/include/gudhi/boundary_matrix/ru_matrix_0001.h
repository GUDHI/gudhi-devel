/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef RU_MATRIX_0001_H
#define RU_MATRIX_0001_H

#include "../utilities/utilities.h"
#include "boundary_matrix_0001.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class RU_matrix_with_removals
		: public Master_matrix::RU_pairing_option,
		  public Master_matrix::RU_vine_swap_option,
		  public Master_matrix::RU_representative_cycles_option
{
public:
	using Column_type = typename Master_matrix::Column_type;
	using Row_type = void;
	using boundary_type = typename Master_matrix::boundary_type;

	RU_matrix_with_removals();
	template<class Boundary_type = boundary_type>
	RU_matrix_with_removals(const std::vector<Boundary_type>& orderedBoundaries);
	RU_matrix_with_removals(unsigned int numberOfColumns);
	RU_matrix_with_removals(const RU_matrix_with_removals& matrixToCopy);
	RU_matrix_with_removals(RU_matrix_with_removals&& other) noexcept;

	template<class Boundary_type = boundary_type>
	void insert_boundary(const Boundary_type& boundary);
	Column_type& get_column(index columnIndex, bool inR = true);
	const Column_type& get_column(index columnIndex, bool inR = true) const;
	Row_type get_row(index rowIndex) const;
	void erase_last();

	dimension_type get_max_dimension() const;
	unsigned int get_number_of_columns() const;

	dimension_type get_column_dimension(index columnIndex) const;

	void add_to(index sourceColumnIndex, index targetColumnIndex);

	void zero_cell(index columnIndex, index rowIndex);
	void zero_column(index columnIndex);
	bool is_zero_cell(index columnIndex, index rowIndex, bool inR = true) const;
	bool is_zero_column(index columnIndex, bool inR = true);

	index get_column_with_pivot(index simplexIndex) const;
	index get_pivot(index columnIndex);

	RU_matrix_with_removals& operator=(RU_matrix_with_removals other);
	friend void swap(RU_matrix_with_removals& matrix1, RU_matrix_with_removals& matrix2){
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

	Boundary_matrix_with_removals<Master_matrix> reducedMatrixR_;
	Boundary_matrix_with_removals<Master_matrix> mirrorMatrixU_;	//make U not accessible by default and add option to enable access? Inaccessible, it needs less options and we could avoid some ifs.
	dictionnary_type pivotToColumnIndex_;
	index nextInsertIndex_;

	void _initialize_U();
	void _reduce();
	void _reduce_last_column();

	static constexpr bool _barcode_option_is_active();
	constexpr barcode_type& _barcode();
	constexpr bar_dictionnary_type& _indexToBar();
};

template<class Master_matrix>
inline RU_matrix_with_removals<Master_matrix>::RU_matrix_with_removals()
	: Master_matrix::RU_pairing_option(),
	  Master_matrix::RU_vine_swap_option(reducedMatrixR_, mirrorMatrixU_, pivotToColumnIndex_),
	  Master_matrix::RU_representative_cycles_option(reducedMatrixR_, mirrorMatrixU_),
	  nextInsertIndex_(0)
{}

template<class Master_matrix>
template<class Boundary_type>
inline RU_matrix_with_removals<Master_matrix>::RU_matrix_with_removals(const std::vector<Boundary_type> &orderedBoundaries)
	: Master_matrix::RU_pairing_option(),
	  Master_matrix::RU_vine_swap_option(reducedMatrixR_, mirrorMatrixU_, pivotToColumnIndex_),
	  Master_matrix::RU_representative_cycles_option(reducedMatrixR_, mirrorMatrixU_),
	  reducedMatrixR_(orderedBoundaries),
	  mirrorMatrixU_(orderedBoundaries.size()),
	  pivotToColumnIndex_(orderedBoundaries.size()),
	  nextInsertIndex_(orderedBoundaries.size())
{
	_initialize_U();
	_reduce();
}

template<class Master_matrix>
inline RU_matrix_with_removals<Master_matrix>::RU_matrix_with_removals(unsigned int numberOfColumns)
	: Master_matrix::RU_pairing_option(),
	  Master_matrix::RU_vine_swap_option(reducedMatrixR_, mirrorMatrixU_, pivotToColumnIndex_),
	  Master_matrix::RU_representative_cycles_option(reducedMatrixR_, mirrorMatrixU_),
	  reducedMatrixR_(numberOfColumns),
	  mirrorMatrixU_(numberOfColumns),
	  pivotToColumnIndex_(numberOfColumns),
	  nextInsertIndex_(0)
{}

template<class Master_matrix>
inline RU_matrix_with_removals<Master_matrix>::RU_matrix_with_removals(const RU_matrix_with_removals &matrixToCopy)
	: Master_matrix::RU_pairing_option(matrixToCopy),
	  Master_matrix::RU_vine_swap_option(matrixToCopy),
	  Master_matrix::RU_representative_cycles_option(matrixToCopy),
	  reducedMatrixR_(matrixToCopy.reducedMatrixR_),
	  mirrorMatrixU_(matrixToCopy.mirrorMatrixU_),
	  pivotToColumnIndex_(matrixToCopy.pivotToColumnIndex_),
	  nextInsertIndex_(matrixToCopy.nextInsertIndex_)
{
	if constexpr (rep_opt::isActive_){
		rep_opt::reducedMatrixR_ = &reducedMatrixR_;
		rep_opt::mirrorMatrixU_ = &mirrorMatrixU_;
	}
	if constexpr (swap_opt::isActive_){
		swap_opt::reducedMatrixR_ = &reducedMatrixR_;
		swap_opt::mirrorMatrixU_ = &mirrorMatrixU_;
		swap_opt::pivotToColumnIndex_ = &pivotToColumnIndex_;
	}
}

template<class Master_matrix>
inline RU_matrix_with_removals<Master_matrix>::RU_matrix_with_removals(RU_matrix_with_removals &&other) noexcept
	: Master_matrix::RU_pairing_option(std::move(other)),
	  Master_matrix::RU_vine_swap_option(std::move(other)),
	  Master_matrix::RU_representative_cycles_option(std::move(other)),
	  reducedMatrixR_(std::move(other.reducedMatrixR_)),
	  mirrorMatrixU_(std::move(other.mirrorMatrixU_)),
	  pivotToColumnIndex_(std::move(other.pivotToColumnIndex_)),
	  nextInsertIndex_(std::exchange(other.nextInsertIndex_, 0))
{
	if constexpr (rep_opt::isActive_){
		rep_opt::reducedMatrixR_ = &reducedMatrixR_;
		rep_opt::mirrorMatrixU_ = &mirrorMatrixU_;
	}
	if constexpr (swap_opt::isActive_){
		swap_opt::reducedMatrixR_ = &reducedMatrixR_;
		swap_opt::mirrorMatrixU_ = &mirrorMatrixU_;
		swap_opt::pivotToColumnIndex_ = &pivotToColumnIndex_;
	}
}

template<class Master_matrix>
template<class Boundary_type>
inline void RU_matrix_with_removals<Master_matrix>::insert_boundary(const Boundary_type &boundary)
{
	reducedMatrixR_.insert_boundary(boundary);

	boundary_type id(1);
	if constexpr (Master_matrix::Field_type::get_characteristic() == 2) {
		id[0] = nextInsertIndex_;
	} else {
		id[0].first = nextInsertIndex_;
		id[0].second = 1;
	}
	mirrorMatrixU_.insert_boundary(id);

	_reduce_last_column();
	++nextInsertIndex_;
}

template<class Master_matrix>
inline typename RU_matrix_with_removals<Master_matrix>::Column_type &
RU_matrix_with_removals<Master_matrix>::get_column(index columnIndex, bool inR)
{
	if (inR){
		return reducedMatrixR_.get_column(columnIndex);
	}
	return mirrorMatrixU_.get_column(columnIndex);
}

template<class Master_matrix>
inline const typename RU_matrix_with_removals<Master_matrix>::Column_type &
RU_matrix_with_removals<Master_matrix>::get_column(index columnIndex, bool inR) const
{
	if (inR){
		return reducedMatrixR_.get_column(columnIndex);
	}
	return mirrorMatrixU_.get_column(columnIndex);
}

template<class Master_matrix>
inline typename RU_matrix_with_removals<Master_matrix>::Row_type
RU_matrix_with_removals<Master_matrix>::get_row(index rowIndex) const
{
	static_assert(static_cast<int>(Master_matrix::Field_type::get_characteristic()) == -1,
			"'get_row' is not implemented for the chosen options.");
}

template<class Master_matrix>
inline void RU_matrix_with_removals<Master_matrix>::erase_last()
{
	--nextInsertIndex_;

	if constexpr (_barcode_option_is_active()){
		typename barcode_type::iterator bar = _indexToBar().at(nextInsertIndex_);

		if (bar->death == -1) _barcode().erase(bar);
		else bar->death = -1;

		_indexToBar().erase(nextInsertIndex_);
	}

	pivotToColumnIndex_.erase(reducedMatrixR_.get_column(nextInsertIndex_).get_pivot());

	reducedMatrixR_.erase_last();
	mirrorMatrixU_.erase_last();
}

template<class Master_matrix>
inline dimension_type RU_matrix_with_removals<Master_matrix>::get_max_dimension() const
{
	return reducedMatrixR_.get_max_dimension();
}

template<class Master_matrix>
inline unsigned int RU_matrix_with_removals<Master_matrix>::get_number_of_columns() const
{
	return reducedMatrixR_.get_number_of_columns();
}

template<class Master_matrix>
inline dimension_type RU_matrix_with_removals<Master_matrix>::get_column_dimension(index columnIndex) const
{
	return reducedMatrixR_.get_column_dimension(columnIndex);
}

template<class Master_matrix>
inline void RU_matrix_with_removals<Master_matrix>::add_to(index sourceColumnIndex, index targetColumnIndex)
{
	reducedMatrixR_.add_to(sourceColumnIndex, targetColumnIndex);
	mirrorMatrixU_.add_to(sourceColumnIndex, targetColumnIndex);
}

template<class Master_matrix>
inline void RU_matrix_with_removals<Master_matrix>::zero_cell(index columnIndex, index rowIndex)
{
	static_assert(static_cast<int>(Master_matrix::Field_type::get_characteristic()) == -1,
			"'zero_cell' is not implemented for the chosen options.");
}

template<class Master_matrix>
inline void RU_matrix_with_removals<Master_matrix>::zero_column(index columnIndex)
{
	static_assert(static_cast<int>(Master_matrix::Field_type::get_characteristic()) == -1,
			"'zero_column' is not implemented for the chosen options.");
}

template<class Master_matrix>
inline bool RU_matrix_with_removals<Master_matrix>::is_zero_cell(index columnIndex, index rowIndex, bool inR) const
{
	if (inR){
		return reducedMatrixR_.is_zero_cell(columnIndex, rowIndex);
	}
	return mirrorMatrixU_.is_zero_cell(columnIndex, rowIndex);
}

template<class Master_matrix>
inline bool RU_matrix_with_removals<Master_matrix>::is_zero_column(index columnIndex, bool inR)
{
	if (inR){
		return reducedMatrixR_.is_zero_column(columnIndex);
	}
	return mirrorMatrixU_.is_zero_column(columnIndex);
}

template<class Master_matrix>
inline index RU_matrix_with_removals<Master_matrix>::get_column_with_pivot(index simplexIndex) const
{
	return pivotToColumnIndex_.at(simplexIndex);
}

template<class Master_matrix>
inline index RU_matrix_with_removals<Master_matrix>::get_pivot(index columnIndex)
{
	return reducedMatrixR_.get_column(columnIndex).get_pivot();
}

template<class Master_matrix>
inline RU_matrix_with_removals<Master_matrix> &RU_matrix_with_removals<Master_matrix>::operator=(RU_matrix_with_removals other)
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
inline void RU_matrix_with_removals<Master_matrix>::print()
{
	std::cout << "R_matrix:\n";
	reducedMatrixR_.print();
	std::cout << "U_matrix:\n";
	mirrorMatrixU_.print();
}

template<class Master_matrix>
inline void RU_matrix_with_removals<Master_matrix>::_initialize_U()
{
	boundary_type id(1);
	if constexpr (Master_matrix::Field_type::get_characteristic() != 2) id[0].second = 1;

	for (unsigned int i = 0; i < reducedMatrixR_.get_number_of_columns(); i++){
		if constexpr (Master_matrix::Field_type::get_characteristic() == 2) id[0] = i;
		else id[0].first = i;
		mirrorMatrixU_.insert_boundary(id);
	}
}

template<class Master_matrix>
inline void RU_matrix_with_removals<Master_matrix>::_reduce()
{
	if constexpr (_barcode_option_is_active()){
		_indexToBar().reserve(reducedMatrixR_.get_number_of_columns());
	}

	for (unsigned int i = 0; i < reducedMatrixR_.get_number_of_columns(); i++){
		if (!(reducedMatrixR_.is_zero_column(i)))
		{
			Column_type &curr = reducedMatrixR_.get_column(i);
			int pivot = curr.get_pivot();

			while (pivot != -1 && pivotToColumnIndex_.find(pivot) != pivotToColumnIndex_.end()){
				if constexpr (Master_matrix::Field_type::get_characteristic() == 2){
					curr += reducedMatrixR_.get_column(pivotToColumnIndex_.at(pivot));
					mirrorMatrixU_.get_column(i) += mirrorMatrixU_.get_column(pivotToColumnIndex_.at(pivot));
				} else {
					Column_type &toadd = reducedMatrixR_.get_column(pivotToColumnIndex_.at(pivot));
					typename Master_matrix::Field_type coef = curr.get_pivot_value();
					coef = coef.get_inverse();
					coef *= (Master_matrix::Field_type::get_characteristic() - static_cast<unsigned int>(toadd.get_pivot_value()));
					curr *= coef;
					curr += toadd;
					mirrorMatrixU_.get_column(i) *= coef;
					mirrorMatrixU_.get_column(i) += mirrorMatrixU_.get_column(pivotToColumnIndex_.at(pivot));
				}

				pivot = curr.get_pivot();
			}

			if (pivot != -1){
				pivotToColumnIndex_.emplace(pivot, i);
				if constexpr (_barcode_option_is_active()){
					_indexToBar().at(pivot)->death = i;
					_indexToBar().emplace(i, _indexToBar().at(pivot));
				}
			} else if constexpr (_barcode_option_is_active()){
				_barcode().emplace_back(get_column_dimension(i), i, -1);
				_indexToBar().emplace(i, --_barcode().end());
			}
		} else if constexpr (_barcode_option_is_active()){
			_barcode().emplace_back(0, i, -1);
			_indexToBar().emplace(i, --_barcode().end());
		}
	}
}

template<class Master_matrix>
inline void RU_matrix_with_removals<Master_matrix>::_reduce_last_column()
{
	Column_type &curr = reducedMatrixR_.get_column(nextInsertIndex_);
	if (curr.is_empty()) {
		if constexpr (_barcode_option_is_active()){
			_barcode().emplace_back(0, nextInsertIndex_, -1);
			_indexToBar().emplace(nextInsertIndex_, --_barcode().end());
		}
		return;
	}

	int pivot = curr.get_pivot();

	while (pivot != -1 && pivotToColumnIndex_.find(pivot) != pivotToColumnIndex_.end()){
		if constexpr (Master_matrix::Field_type::get_characteristic() == 2){
			curr += reducedMatrixR_.get_column(pivotToColumnIndex_.at(pivot));
			mirrorMatrixU_.get_column(nextInsertIndex_) += mirrorMatrixU_.get_column(pivotToColumnIndex_.at(pivot));
		} else {
			Column_type &toadd = reducedMatrixR_.get_column(pivotToColumnIndex_.at(pivot));
			typename Master_matrix::Field_type coef = curr.get_pivot_value();
			coef = coef.get_inverse();
			coef *= (Master_matrix::Field_type::get_characteristic() - static_cast<unsigned int>(toadd.get_pivot_value()));
			curr *= coef;
			curr += toadd;
			mirrorMatrixU_.get_column(nextInsertIndex_) *= coef;
			mirrorMatrixU_.get_column(nextInsertIndex_) += mirrorMatrixU_.get_column(pivotToColumnIndex_.at(pivot));
		}

		pivot = curr.get_pivot();
	}

	if (pivot != -1){
		pivotToColumnIndex_.emplace(pivot, nextInsertIndex_);
		if constexpr (_barcode_option_is_active()){
			_indexToBar().at(pivot)->death = nextInsertIndex_;
			_indexToBar().emplace(nextInsertIndex_, _indexToBar().at(pivot));
		}
	} else if constexpr (_barcode_option_is_active()){
		_barcode().emplace_back(get_column_dimension(nextInsertIndex_), nextInsertIndex_, -1);
		_indexToBar().emplace(nextInsertIndex_, --_barcode().end());
	}
}

template<class Master_matrix>
inline constexpr bool RU_matrix_with_removals<Master_matrix>::_barcode_option_is_active()
{
	return swap_opt::isActive_ || pair_opt::isActive_;
}

template<class Master_matrix>
inline constexpr typename RU_matrix_with_removals<Master_matrix>::barcode_type &RU_matrix_with_removals<Master_matrix>::_barcode()
{
	if constexpr (swap_opt::isActive_)
		return swap_opt::template RU_pairing<Master_matrix>::barcode_;
	else
		return pair_opt::barcode_;
}

template<class Master_matrix>
inline constexpr typename RU_matrix_with_removals<Master_matrix>::bar_dictionnary_type &RU_matrix_with_removals<Master_matrix>::_indexToBar()
{
	if constexpr (swap_opt::isActive_)
		return swap_opt::template RU_pairing<Master_matrix>::indexToBar_;
	else
		return pair_opt::indexToBar_;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // RU_MATRIX_0001_H

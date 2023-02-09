/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef BASE_MATRIX_1000_H
#define BASE_MATRIX_1000_H

#include <iostream>
#include <vector>
#include <unordered_map>

//#include <boost/pending/disjoint_sets.hpp>

#include "../utilities/utilities.h"
#include "../utilities/union_find.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class Base_matrix_with_column_compression				//swap possible, but needs to be readapted: TODO later
{
public:
	using Field_element_type = typename Master_matrix::Field_type;
	using Column_type = typename Master_matrix::Column_type;
	using boundary_type = typename Master_matrix::boundary_type;
	using Row_type = void;

	Base_matrix_with_column_compression();
	template<class Container_type>
	Base_matrix_with_column_compression(const std::vector<Container_type>& columns);
	Base_matrix_with_column_compression(unsigned int numberOfColumns);
	Base_matrix_with_column_compression(const Base_matrix_with_column_compression& matrixToCopy);
	Base_matrix_with_column_compression(Base_matrix_with_column_compression&& other) noexcept;

	template<class Container_type>
	void insert_column(const Container_type& column);
	template<class Boundary_type>
	void insert_boundary(const Boundary_type& boundary);
	const Column_type& get_column(index columnIndex);
	Row_type get_row(index rowIndex) const;
	void erase_column(index columnIndex);
	void erase_row(index rowIndex);

	unsigned int get_number_of_columns() const;

	void add_to(index sourceColumnIndex, index targetColumnIndex);	//adds classes and not individual columns
	void add_to(const Column_type& sourceColumn, index targetColumnIndex);
	void add_to(const Column_type& sourceColumn, const Field_element_type& coefficient, index targetColumnIndex);
	void add_to(const Field_element_type& coefficient, const Column_type& sourceColumn, index targetColumnIndex);

	void zero_cell(index columnIndex, index rowIndex);
	void zero_column(index columnIndex);
	bool is_zero_cell(index columnIndex, index rowIndex);
	bool is_zero_column(index columnIndex);

	Base_matrix_with_column_compression& operator=(Base_matrix_with_column_compression other);
	friend void swap(Base_matrix_with_column_compression& matrix1, Base_matrix_with_column_compression& matrix2){
		matrix1.matrix_.swap(matrix2.matrix_);
		swap(matrix1.columnClasses_, matrix2.columnClasses_);
		matrix1.repToColumn_.swap(matrix2.repToColumn_);
		std::swap(matrix1.nextColumnIndex_, matrix2.nextColumnIndex_);
	}

	void print();  //for debug

private:
	std::unordered_map<Column_type, index> matrix_;
	Union_find columnClasses_;
	std::vector<const Column_type*> repToColumn_;
	unsigned int nextColumnIndex_;

	void _insert_column(Column_type& column, index columnIndex);
};

template<class Master_matrix>
inline Base_matrix_with_column_compression<Master_matrix>::Base_matrix_with_column_compression()
	: nextColumnIndex_(0)
{}

template<class Master_matrix>
template<class Container_type>
inline Base_matrix_with_column_compression<Master_matrix>::Base_matrix_with_column_compression(const std::vector<Container_type> &columns)
	: matrix_(columns.size()), columnClasses_(columns.size()), repToColumn_(columns.size(), nullptr), nextColumnIndex_(0)
{
	for (const Container_type& c : columns){
		insert_column(c);
	}
}

template<class Master_matrix>
inline Base_matrix_with_column_compression<Master_matrix>::Base_matrix_with_column_compression(unsigned int numberOfColumns)
	: matrix_(numberOfColumns), columnClasses_(numberOfColumns), repToColumn_(numberOfColumns, nullptr), nextColumnIndex_(0)
{}

template<class Master_matrix>
inline Base_matrix_with_column_compression<Master_matrix>::Base_matrix_with_column_compression(const Base_matrix_with_column_compression &matrixToCopy)
	: matrix_(matrixToCopy.matrix_),
	  columnClasses_(matrixToCopy.columnClasses_),
	  repToColumn_(matrixToCopy.repToColumn_),
	  nextColumnIndex_(matrixToCopy.nextColumnIndex_)
{
	for (unsigned int i = 0; i < repToColumn_.size(); ++i){
		if (repToColumn_[i] != nullptr){
			repToColumn_[i] = &(matrix_.find(*repToColumn_[i])->first);
		}
	}
}

template<class Master_matrix>
inline Base_matrix_with_column_compression<Master_matrix>::Base_matrix_with_column_compression(Base_matrix_with_column_compression &&other) noexcept
	: matrix_(std::move(other.matrix_)),
	  columnClasses_(std::move(other.columnClasses_)),
	  repToColumn_(std::move(other.repToColumn_)),
	  nextColumnIndex_(std::exchange(other.nextColumnIndex_, 0))
{}

template<class Master_matrix>
template<class Container_type>
inline void Base_matrix_with_column_compression<Master_matrix>::insert_column(const Container_type &column)
{
	columnClasses_.initialize(nextColumnIndex_);
	if (repToColumn_.size() <= nextColumnIndex_) repToColumn_.resize(nextColumnIndex_ * 2, nullptr);
	Column_type col(column);
	_insert_column(col, nextColumnIndex_);
	nextColumnIndex_++;
}

template<class Master_matrix>
template<class Boundary_type>
inline void Base_matrix_with_column_compression<Master_matrix>::insert_boundary(const Boundary_type &boundary)
{
	insert_column(boundary);
}

template<class Master_matrix>
inline const typename Base_matrix_with_column_compression<Master_matrix>::Column_type&
Base_matrix_with_column_compression<Master_matrix>::get_column(index columnIndex)
{
	return *repToColumn_[columnClasses_.find(columnIndex)];
}

template<class Master_matrix>
inline typename Base_matrix_with_column_compression<Master_matrix>::Row_type
Base_matrix_with_column_compression<Master_matrix>::get_row(index rowIndex) const
{
	static_assert(Master_matrix::Option_list::has_row_access,
			"'get_row' is not implemented for the chosen options.");
}

template<class Master_matrix>
inline void Base_matrix_with_column_compression<Master_matrix>::erase_column(index columnIndex)
{
	static_assert(Master_matrix::Option_list::has_row_access,
			"'erase_column' is not implemented for the chosen options.");
}

template<class Master_matrix>
inline void Base_matrix_with_column_compression<Master_matrix>::erase_row(index rowIndex)
{
	static_assert(Master_matrix::Option_list::has_row_access,
			"'erase_row' is not implemented for the chosen options.");
}

template<class Master_matrix>
inline unsigned int Base_matrix_with_column_compression<Master_matrix>::get_number_of_columns() const
{
	return nextColumnIndex_;
}

template<class Master_matrix>
inline void Base_matrix_with_column_compression<Master_matrix>::add_to(index sourceColumnIndex, index targetColumnIndex)
{
	//handle case where targetRep == sourceRep?
	Column_type target(get_column(targetColumnIndex));
	auto res = matrix_.find(target);
	index targetRep = res->second;
	matrix_.erase(res);
	target += get_column(sourceColumnIndex);
	_insert_column(target, targetRep);
}

template<class Master_matrix>
inline void Base_matrix_with_column_compression<Master_matrix>::add_to(const Column_type& sourceColumn, index targetColumnIndex)
{
	Column_type target(get_column(targetColumnIndex));
	auto res = matrix_.find(target);
	index targetRep = res->second;
	matrix_.erase(res);
	target += sourceColumn;
	_insert_column(target, targetRep);
}

template<class Master_matrix>
inline void Base_matrix_with_column_compression<Master_matrix>::add_to(const Column_type& sourceColumn, const Field_element_type& coefficient, index targetColumnIndex)
{
	Column_type target(get_column(targetColumnIndex));
	auto res = matrix_.find(target);
	index targetRep = res->second;
	matrix_.erase(res);
	target.multiply_and_add(coefficient, sourceColumn);
	_insert_column(target, targetRep);
}

template<class Master_matrix>
inline void Base_matrix_with_column_compression<Master_matrix>::add_to(const Field_element_type& coefficient, const Column_type& sourceColumn, index targetColumnIndex)
{
	Column_type target(get_column(targetColumnIndex));
	auto res = matrix_.find(target);
	index targetRep = res->second;
	matrix_.erase(res);
	target.multiply_and_add(sourceColumn, coefficient);
	_insert_column(target, targetRep);
}

template<class Master_matrix>
inline void Base_matrix_with_column_compression<Master_matrix>::zero_cell(index columnIndex, index rowIndex)
{
	static_assert(Master_matrix::Option_list::has_row_access,
			"'zero_cell' is not implemented for the chosen options.");
}

template<class Master_matrix>
inline void Base_matrix_with_column_compression<Master_matrix>::zero_column(index columnIndex)
{
	static_assert(Master_matrix::Option_list::has_row_access,
			"'zero_column' is not implemented for the chosen options.");
}

template<class Master_matrix>
inline bool Base_matrix_with_column_compression<Master_matrix>::is_zero_cell(index columnIndex, index rowIndex)
{
	return !(repToColumn_[columnClasses_.find(columnIndex)]->is_non_zero(rowIndex));
}

template<class Master_matrix>
inline bool Base_matrix_with_column_compression<Master_matrix>::is_zero_column(index columnIndex)
{
	return repToColumn_[columnClasses_.find(columnIndex)]->is_empty();
}

template<class Master_matrix>
inline Base_matrix_with_column_compression<Master_matrix> &Base_matrix_with_column_compression<Master_matrix>::operator=(Base_matrix_with_column_compression other)
{
	matrix_.swap(other.matrix_);
	swap(columnClasses_, other.columnClasses_);
	repToColumn_.swap(other.repToColumn_);
	std::swap(nextColumnIndex_, other.nextColumnIndex_);
	return *this;
}

template<class Master_matrix>
inline void Base_matrix_with_column_compression<Master_matrix>::print()
{
	std::cout << "Compressed_matrix:\n";
	for (auto& p : matrix_){
		Column_type& col = p.first;
		for (auto e : col.get_content(nextColumnIndex_)){
			if (e == 0u) std::cout << "- ";
			else std::cout << e << " ";
		}
		std::cout << "(";
		for (unsigned int i = 0; i < nextColumnIndex_; ++i){
			if (columnClasses_.find(i) == p.second)
				std::cout << i << " ";
		}
		std::cout << ")\n";
	}
	std::cout << "\n";
}

template<class Master_matrix>
inline void Base_matrix_with_column_compression<Master_matrix>::_insert_column(Column_type &column, index columnIndex)
{
	auto res = matrix_.emplace(column, columnIndex);
	repToColumn_[columnIndex] = &(res.first->first);
	if (!res.second){
		columnClasses_.merge(columnIndex, res.first->second);
		repToColumn_[res.first->second] = &(res.first->first);
		res.first->second = columnClasses_.find(columnIndex);
	}
}

//	struct Union_find{
//	public:
//		Union_find() : sets_(ranks_.data(), parents_.data()){};
//		Union_find(const Union_find& toCopy)
//			: ranks_(toCopy.ranks_),
//			  parents_(toCopy.parents_),
//			  sets_(ranks_.data(), parents_.data()){};
//		Union_find(Union_find&& other) noexcept
//			: ranks_(std::move(other.ranks_)),
//			  parents_(std::move(other.parents_)),
//			  sets_(ranks_.data(), parents_.data()){};	//sets_ stores pointers, it seems that there is no guarantee they are preserved?

//		void initialize(index id){
//			if (id >= parents_.size()){
//				parents_.resize(id + 1);
//				ranks_.resize(id + 1);
//			}
//			sets_.make_set(id);
//		}

//		void merge(index id1, index id2){
//			sets_.link(id1, id2);
//		}

////		void find_and_merge(index id1, index id2){
////			sets_.union_set(id1, id2);
////		}

//		index find(index id){
//			return sets_.find_set(id);
//		}

//		friend void swap(Union_find& uf1, Union_find& uf2){
//			uf1.ranks_.swap(uf2.ranks_);
//			uf1.parents_.swap(uf2.parents_);
////			uf1.sets_ = boost::disjoint_sets<int*, index*>(uf1.ranks_.data(), uf1.parents_.data());
////			uf2.sets_ = boost::disjoint_sets<int*, index*>(uf2.ranks_.data(), uf2.parents_.data());
//			std::swap(uf1.sets_, uf2.sets_);
//		}

//	private:
//		std::vector<int> ranks_;
//		std::vector<index> parents_;
//		boost::disjoint_sets<int*, index*> sets_;
//	};

} //namespace persistence_matrix
} //namespace Gudhi

#endif // BASE_MATRIX_1000_H

/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef BASE_MATRIX_1010_H
#define BASE_MATRIX_1010_H

#include <iostream>
#include <vector>
#include <unordered_set>
#include <set>

//#include <boost/pending/disjoint_sets.hpp>

#include "../utilities/utilities.h"
#include "../utilities/union_find.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class Base_matrix_with_column_compression_with_row_access				//swap possible, but needs to be readapted: TODO later
{
public:
	using Column_type = typename Master_matrix::Column_type;
	using Row_type = typename Master_matrix::Row_type;
	using Field_element_type = typename Master_matrix::Field_type;

	Base_matrix_with_column_compression_with_row_access();
	template<class Container_type>
	Base_matrix_with_column_compression_with_row_access(const std::vector<Container_type>& columns);
	Base_matrix_with_column_compression_with_row_access(unsigned int numberOfColumns);
	Base_matrix_with_column_compression_with_row_access(const Base_matrix_with_column_compression_with_row_access& matrixToCopy);
	Base_matrix_with_column_compression_with_row_access(Base_matrix_with_column_compression_with_row_access&& other) noexcept;

	template<class Container_type>
	void insert_column(const Container_type& column);
	template<class Boundary_type>
	void insert_boundary(const Boundary_type& boundary);
	Column_type& get_column(index columnIndex);
	Row_type& get_row(index rowIndex);
	const Row_type& get_row(index rowIndex) const;
	void erase_column(index columnIndex);
	void erase_row(index rowIndex);

	unsigned int get_number_of_columns() const;

	void add_to(index sourceColumnIndex, index targetColumnIndex);
	void add_to(const Column_type& sourceColumn, index targetColumnIndex);
	void add_to(const Column_type& sourceColumn, const Field_element_type& coefficient, index targetColumnIndex);
	void add_to(const Field_element_type& coefficient, const Column_type& sourceColumn, index targetColumnIndex);

	void zero_cell(index columnIndex, index rowIndex);
	void zero_column(index columnIndex);
	bool is_zero_cell(index columnIndex, index rowIndex);
	bool is_zero_column(index columnIndex);

	Base_matrix_with_column_compression_with_row_access& operator=(const Base_matrix_with_column_compression_with_row_access& other);
	friend void swap(Base_matrix_with_column_compression_with_row_access& matrix1, Base_matrix_with_column_compression_with_row_access& matrix2){
		matrix1.rows_.swap(matrix2.rows_);
		matrix1.columnToRep_.swap(matrix2.columnToRep_);
		swap(matrix1.columnClasses_, matrix2.columnClasses_);
		matrix1.repToColumn_.swap(matrix2.repToColumn_);
		std::swap(matrix1.nextColumnIndex_, matrix2.nextColumnIndex_);
		for (auto& p : matrix1.repToColumn_){
			Column_type& col = p.second;
			col.set_rows(&matrix1.rows_);
		}
		for (auto& p : matrix2.repToColumn_){
			Column_type& col = p.second;
			col.set_rows(&matrix2.rows_);
		}
	}

	void print();  //for debug

private:
	using rows_type = typename Master_matrix::row_container_type;
	using cell_rep_type = typename std::conditional<
								Master_matrix::Option_list::is_z2,
								index,
								std::pair<index,typename Master_matrix::Field_type>
							>::type;
	using tmp_column_type = typename std::conditional<
								Master_matrix::Option_list::is_z2,
								std::set<index>,
								std::set<std::pair<index,Field_element_type>,CellPairComparator<Field_element_type> >
							>::type;

	struct ColumnPointerHash {
		size_t operator()(const Column_type* c) const
		{
			return std::hash<Column_type>()(*c);
		}
	};

	struct ColumnPointerEqual : std::binary_function<const Column_type*, const Column_type*, bool> {
		bool operator()(const Column_type* c1, const Column_type* c2) const
		{
			return *c1 == *c2;
		}
	};

	rows_type rows_;	//has to be destroyed after repToColumn_
	std::unordered_map<Column_type*, index, ColumnPointerHash, ColumnPointerEqual> columnToRep_;
	Union_find columnClasses_;
	std::unordered_map<index,Column_type> repToColumn_;
	unsigned int nextColumnIndex_;

	void _insert_column(typename std::unordered_map<index,Column_type>::iterator& it, index columnIndex);
//	void _add_to(const Column_type& column, tmp_column_type& set);
};

template<class Master_matrix>
inline Base_matrix_with_column_compression_with_row_access<Master_matrix>::Base_matrix_with_column_compression_with_row_access()
	: nextColumnIndex_(0)
{}

template<class Master_matrix>
template<class Container_type>
inline Base_matrix_with_column_compression_with_row_access<Master_matrix>::Base_matrix_with_column_compression_with_row_access(const std::vector<Container_type> &columns)
	: rows_(columns.size()),
	  columnToRep_(columns.size()),
	  columnClasses_(columns.size()),
	  repToColumn_(columns.size()),
	  nextColumnIndex_(0)
{
	for (const Container_type& c : columns){
		insert_column(c);
	}
}

template<class Master_matrix>
inline Base_matrix_with_column_compression_with_row_access<Master_matrix>::Base_matrix_with_column_compression_with_row_access(unsigned int numberOfColumns)
	: rows_(numberOfColumns),
	  columnToRep_(numberOfColumns),
	  columnClasses_(numberOfColumns),
	  repToColumn_(numberOfColumns),
	  nextColumnIndex_(0)
{}

template<class Master_matrix>
inline Base_matrix_with_column_compression_with_row_access<Master_matrix>::Base_matrix_with_column_compression_with_row_access(const Base_matrix_with_column_compression_with_row_access &matrixToCopy)
	: rows_(matrixToCopy.rows_.size()),
	  columnToRep_(matrixToCopy.columnToRep_.size()),
	  columnClasses_(matrixToCopy.columnClasses_),
	  repToColumn_(matrixToCopy.repToColumn_.size()),
	  nextColumnIndex_(matrixToCopy.nextColumnIndex_)
{
	for (const auto& p : matrixToCopy.repToColumn_){
		const Column_type& col = p.second;
		std::vector<cell_rep_type> tmp(col.begin(), col.end());
		auto res = repToColumn_.emplace(p.first, Column_type(col.get_column_index(), tmp, col.get_dimension(), rows_));
		columnToRep_.emplace(&(res.first->second), p.first);
	}
}

template<class Master_matrix>
inline Base_matrix_with_column_compression_with_row_access<Master_matrix>::Base_matrix_with_column_compression_with_row_access(Base_matrix_with_column_compression_with_row_access &&other) noexcept
	: rows_(std::move(other.rows_)),
	  columnToRep_(std::move(other.columnToRep_)),
	  columnClasses_(std::move(other.columnClasses_)),
	  repToColumn_(std::move(other.repToColumn_)),
	  nextColumnIndex_(std::exchange(other.nextColumnIndex_, 0))
{
	for (auto& p : repToColumn_){
		Column_type& col = p.second;
		col.set_rows(&rows_);
	}
}

template<class Master_matrix>
template<class Container_type>
inline void Base_matrix_with_column_compression_with_row_access<Master_matrix>::insert_column(const Container_type &column)
{
	columnClasses_.initialize(nextColumnIndex_);
	if (rows_.size() <= nextColumnIndex_) rows_.resize(nextColumnIndex_ + 1);
	auto p = repToColumn_.emplace(nextColumnIndex_, Column_type(nextColumnIndex_, column, rows_));
	_insert_column(p.first, nextColumnIndex_);
	nextColumnIndex_++;
}

template<class Master_matrix>
template<class Boundary_type>
inline void Base_matrix_with_column_compression_with_row_access<Master_matrix>::insert_boundary(const Boundary_type &boundary)
{
	insert_column(boundary);
}

template<class Master_matrix>
inline typename Base_matrix_with_column_compression_with_row_access<Master_matrix>::Column_type &Base_matrix_with_column_compression_with_row_access<Master_matrix>::get_column(index columnIndex)
{
	return repToColumn_.at(columnClasses_.find(columnIndex));
}

template<class Master_matrix>
inline typename Base_matrix_with_column_compression_with_row_access<Master_matrix>::Row_type& Base_matrix_with_column_compression_with_row_access<Master_matrix>::get_row(index rowIndex)
{
	return rows_[rowIndex];
}

template<class Master_matrix>
inline const typename Base_matrix_with_column_compression_with_row_access<Master_matrix>::Row_type& Base_matrix_with_column_compression_with_row_access<Master_matrix>::get_row(index rowIndex) const
{
	return rows_[rowIndex];
}

template<class Master_matrix>
inline void Base_matrix_with_column_compression_with_row_access<Master_matrix>::erase_column(index columnIndex)
{
	static_assert(Master_matrix::Option_list::has_removable_columns,
			"'erase_column' is not implemented for the chosen options.");
}

template<class Master_matrix>
inline void Base_matrix_with_column_compression_with_row_access<Master_matrix>::erase_row(index rowIndex)
{
	static_assert(Master_matrix::Option_list::has_removable_columns,
			"'erase_row' is not implemented for the chosen options.");
}

template<class Master_matrix>
inline unsigned int Base_matrix_with_column_compression_with_row_access<Master_matrix>::get_number_of_columns() const
{
	return nextColumnIndex_;
}

template<class Master_matrix>
inline void Base_matrix_with_column_compression_with_row_access<Master_matrix>::add_to(index sourceColumnIndex, index targetColumnIndex)
{
	//handle case where targetRep == sourceRep?
	index targetRep = columnClasses_.find(targetColumnIndex);
	auto itTarget = repToColumn_.find(targetRep);
	Column_type& target = itTarget->second;
	columnToRep_.erase(&target);
	target += get_column(sourceColumnIndex);
	_insert_column(itTarget, targetRep);
}

template<class Master_matrix>
inline void Base_matrix_with_column_compression_with_row_access<Master_matrix>::add_to(const Column_type& sourceColumn, index targetColumnIndex)
{
	index targetRep = columnClasses_.find(targetColumnIndex);
	auto itTarget = repToColumn_.find(targetRep);
	Column_type& target = itTarget->second;
	columnToRep_.erase(&target);
	target += sourceColumn;
	_insert_column(itTarget, targetRep);
}

template<class Master_matrix>
inline void Base_matrix_with_column_compression_with_row_access<Master_matrix>::add_to(const Column_type& sourceColumn, const Field_element_type& coefficient, index targetColumnIndex)
{
	index targetRep = columnClasses_.find(targetColumnIndex);
	auto itTarget = repToColumn_.find(targetRep);
	Column_type& target = itTarget->second;
	columnToRep_.erase(&target);
	target.multiply_and_add(coefficient, sourceColumn);
	_insert_column(itTarget, targetRep);
}

template<class Master_matrix>
inline void Base_matrix_with_column_compression_with_row_access<Master_matrix>::add_to(const Field_element_type& coefficient, const Column_type& sourceColumn, index targetColumnIndex)
{
	index targetRep = columnClasses_.find(targetColumnIndex);
	auto itTarget = repToColumn_.find(targetRep);
	Column_type& target = itTarget->second;
	columnToRep_.erase(&target);
	target.multiply_and_add(sourceColumn, coefficient);
	_insert_column(itTarget, targetRep);
}

template<class Master_matrix>
inline void Base_matrix_with_column_compression_with_row_access<Master_matrix>::zero_cell(index columnIndex, index rowIndex)
{
	static_assert(Master_matrix::Option_list::has_removable_columns,
			"'zero_cell' is not implemented for the chosen options.");
}

template<class Master_matrix>
inline void Base_matrix_with_column_compression_with_row_access<Master_matrix>::zero_column(index columnIndex)
{
	static_assert(Master_matrix::Option_list::has_removable_columns,
			"'zero_column' is not implemented for the chosen options.");
}

template<class Master_matrix>
inline bool Base_matrix_with_column_compression_with_row_access<Master_matrix>::is_zero_cell(index columnIndex, index rowIndex)
{
	return !(repToColumn_.at(columnClasses_.find(columnIndex)).is_non_zero(rowIndex));
}

template<class Master_matrix>
inline bool Base_matrix_with_column_compression_with_row_access<Master_matrix>::is_zero_column(index columnIndex)
{
	return repToColumn_.at(columnClasses_.find(columnIndex)).is_empty();
}

template<class Master_matrix>
inline Base_matrix_with_column_compression_with_row_access<Master_matrix> &
Base_matrix_with_column_compression_with_row_access<Master_matrix>::operator=(const Base_matrix_with_column_compression_with_row_access& other)
{
	rows_.resize(other.rows_.size());
	columnToRep_.reserve(other.columnToRep_.size());
	repToColumn_.reserve(other.repToColumn_.size());
	nextColumnIndex_ = other.nextColumnIndex_;
	for (const auto& p : other.repToColumn_){
		const Column_type& col = p.second;
		std::vector<cell_rep_type> tmp(col.begin(), col.end());
		auto res = repToColumn_.emplace(p.first, Column_type(col.get_column_index(), tmp, col.get_dimension(), rows_));
		columnToRep_.emplace(&(*res.first), p.first);
	}
	return *this;
}

template<class Master_matrix>
inline void Base_matrix_with_column_compression_with_row_access<Master_matrix>::print()
{
	std::cout << "Compressed_matrix:\n";
	for (auto& p : columnToRep_){
		Column_type* col = p.first;
		for (auto e : col->get_content(nextColumnIndex_)){
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
	std::cout << "Row Matrix:\n";
	for (unsigned int i = 0; i < rows_.size(); ++i){
		const Row_type& row = rows_[i];
		for (const auto &cell : row){
			std::cout << cell.get_column_index() << " ";
		}
		std::cout << "(" << i << ")\n";
	}
	std::cout << "\n";
}

template<class Master_matrix>
inline void Base_matrix_with_column_compression_with_row_access<Master_matrix>::_insert_column(typename std::unordered_map<index,Column_type>::iterator& it, index columnIndex)
{
	Column_type* col = &(it->second);
	auto res = columnToRep_.emplace(col, columnIndex);
	if (!res.second){
		index doubleRep = res.first->second;
		index newRep = columnClasses_.merge(columnIndex, doubleRep);
		if (newRep == columnIndex){
			auto pos = columnToRep_.erase(res.first);
			repToColumn_.erase(doubleRep);
			columnToRep_.emplace_hint(pos, col, columnIndex);
		} else {
			repToColumn_.erase(it);
		}
	}
}

//template<class Master_matrix>
//inline void Base_matrix_with_column_compression_with_row_access<Master_matrix>::_add_to(
//		const Column_type &column,
//		tmp_column_type &set)
//{
//	if constexpr (Master_matrix::Option_list::is_z2){
//		std::pair<std::set<index>::iterator,bool> res_insert;
//		for (const auto &cell : column) {
//			res_insert = set.insert(cell.get_row_index());
//			if (!res_insert.second) {
//				set.erase(res_insert.first);
//			}
//		}
//	} else {
//		for (const auto &cell : column) {
//			std::pair<index,Field_element_type> p(cell.get_row_index(), cell.get_element());
//			auto res_it = set.find(p);

//			if (res_it != set.end()){
//				p.second += res_it->second;
//				set.erase(res_it);
//				if (p.second != 0){
//					set.insert(p);
//				}
//			} else {
//				set.insert(p);
//			}
//		}
//	}
//}

//	struct Union_find{
//	public:
//		Union_find() : sets_(ranks_.data(), parents_.data()){};
//		Union_find(unsigned int numberOfElements)
//			: ranks_(numberOfElements),
//			  parents_(numberOfElements),
//			  sets_(ranks_.data(), parents_.data()){};
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
//				sets_ = boost::disjoint_sets<int*, index*>(ranks_.data(), parents_.data());
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

#endif // BASE_MATRIX_1010_H

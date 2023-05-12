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
#include <unordered_map>
//#include <map>

#include <gudhi/Simple_object_pool.h>
#include <boost/intrusive/set.hpp>

#include "../utilities/utilities.h"
#include "../utilities/union_find.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class Base_matrix_with_column_compression_with_row_access				//swap possible, but needs to be readapted: TODO later
{
public:
//	using Column_type = typename Master_matrix::Column_type;
	using Row_type = typename Master_matrix::Row_type;
	using Field_element_type = typename Master_matrix::Field_type;

	class Column_type : public Master_matrix::Column_type,
			public boost::intrusive::set_base_hook<boost::intrusive::link_mode<boost::intrusive::normal_link> >
	{
	public:
		using Base = typename Master_matrix::Column_type;

		Column_type()
			: Base()
		{}
		template<class Container_type>
		Column_type(const Container_type& nonZeroRowIndices)
			: Base(nonZeroRowIndices)
		{}
		template<class Container_type>
		Column_type(const Container_type& nonZeroRowIndices, dimension_type dimension)
			: Base(nonZeroRowIndices, dimension)
		{}
		template<class Row_container_type>
		Column_type(index columnIndex, Row_container_type &rowContainer)
			: Base(columnIndex, rowContainer)
		{}
		template<class Container_type, class Row_container_type>
		Column_type(index columnIndex, const Container_type& nonZeroRowIndices, Row_container_type &rowContainer)
			: Base(columnIndex, nonZeroRowIndices, rowContainer)
		{}
		template<class Container_type, class Row_container_type>
		Column_type(index columnIndex, const Container_type& nonZeroRowIndices, dimension_type dimension, Row_container_type &rowContainer)
			: Base(columnIndex, nonZeroRowIndices, dimension, rowContainer)
		{}
		Column_type(const Column_type& column)
			: Base(static_cast<const Base&>(column))
		{}
		Column_type(const Column_type& column, index columnIndex)
			: Base(static_cast<const Base&>(column), columnIndex)
		{}
		template<class Row_container_type>
		Column_type(const Column_type& column, index columnIndex, Row_container_type &rowContainer)
			: Base(static_cast<const Base&>(column), columnIndex, rowContainer)
		{}
		Column_type(Column_type&& column) noexcept
			: Base(std::move(static_cast<Base&&>(column)))
		{}

		index get_rep() const{
			return rep_;
		}
		void set_rep(const index &rep){
			rep_ = rep;
		}

		struct Hasher {
			size_t operator()(const Column_type& c) const
			{
				return std::hash<Base>()(c);
			}
		};

	private:
		index rep_;
	};

	Base_matrix_with_column_compression_with_row_access();
	template<class Container_type>
	Base_matrix_with_column_compression_with_row_access(const std::vector<Container_type>& columns);
	Base_matrix_with_column_compression_with_row_access(unsigned int numberOfColumns);
	Base_matrix_with_column_compression_with_row_access(const Base_matrix_with_column_compression_with_row_access& matrixToCopy);
	Base_matrix_with_column_compression_with_row_access(Base_matrix_with_column_compression_with_row_access&& other) noexcept;
	~Base_matrix_with_column_compression_with_row_access();

	template<class Container_type>
	void insert_column(const Container_type& column);
//	void insert_column(index rowIndex);
//	void insert_column(index rowIndex, const Field_element_type& coefficient);
	template<class Boundary_type>
	void insert_boundary(const Boundary_type& boundary);
	const Column_type& get_column(index columnIndex);
//	Row_type& get_row(index rowIndex);
	const Row_type& get_row(index rowIndex) const;
	void erase_column(index columnIndex);
	void erase_row(index rowIndex);

	unsigned int get_number_of_columns() const;

	void add_to(index sourceColumnIndex, index targetColumnIndex);
//	void add_to(const Column_type& sourceColumn, index targetColumnIndex);
//	void add_to(const Column_type& sourceColumn, const Field_element_type& coefficient, index targetColumnIndex);
//	void add_to(const Field_element_type& coefficient, const Column_type& sourceColumn, index targetColumnIndex);
	template<class Cell_range>
	void add_to(const Cell_range& sourceColumn, index targetColumnIndex);
	template<class Cell_range>
	void add_to(const Cell_range& sourceColumn, const Field_element_type& coefficient, index targetColumnIndex);
	template<class Cell_range>
	void add_to(const Field_element_type& coefficient, const Cell_range& sourceColumn, index targetColumnIndex);
//	void add_to(const Field_element_type& coefficient, const typename Master_matrix::Column_type& sourceColumn, index targetColumnIndex);

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
//		for (auto& p : matrix1.columnToRep_){
//			p.first->set_rows(&matrix1.rows_);
//		}
//		for (auto& p : matrix2.columnToRep_){
//			p.first->set_rows(&matrix2.rows_);
//		}
		for (auto& col : matrix1.columnToRep_){
			col.set_rows(&matrix1.rows_);
		}
		for (auto& col : matrix2.columnToRep_){
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

//	struct ColumnPointerHash {
//		size_t operator()(const Column_type* c) const
//		{
//			return std::hash<Column_type>()(*c);
//		}
//	};

//	struct ColumnPointerEqual {
//		bool operator()(const Column_type* c1, const Column_type* c2) const
//		{
//			return *c1 == *c2;
//		}
//	};

//	struct ColumnPointerStrictlySmallerThan {
//		bool operator()(const Column_type* c1, const Column_type* c2) const
//		{
//			return *c1 < *c2;
//		}
//	};

	//The disposer object function for boost intrusive container
	struct delete_disposer
	{
		void operator()(Column_type *delete_this){
			columnPool_.destroy(delete_this);
		}
	};

//	using col_dict_type = std::unordered_map<Column_type*, index, ColumnPointerHash, ColumnPointerEqual>;
	using col_dict_type = boost::intrusive::set<Column_type, boost::intrusive::constant_time_size<false> >;

	rows_type rows_;	//has to be destroyed after repToColumn_
	col_dict_type columnToRep_;
//	boost::unordered_flat_map<Column_type*, index, ColumnPointerHash, ColumnPointerEqual> columnToRep_;
//	absl::flat_hash_map<Column_type*, index, ColumnPointerHash, ColumnPointerEqual> columnToRep_;
//	boost::container::flat_map<Column_type*, index, ColumnPointerStrictlySmallerThan> columnToRep_;
//	std::map<Column_type*, index, ColumnPointerStrictlySmallerThan> columnToRep_;
	Union_find columnClasses_;
	std::vector<Column_type*> repToColumn_;
	unsigned int nextColumnIndex_;
	inline static Simple_object_pool<Column_type> columnPool_;
	inline static const Column_type empty_column_;

	void _insert_column(index columnIndex);
	void _insert_double_column(index columnIndex, typename col_dict_type::iterator& doubleIt);
};

template<class Master_matrix>
inline Base_matrix_with_column_compression_with_row_access<Master_matrix>::Base_matrix_with_column_compression_with_row_access()
	: nextColumnIndex_(0)
{}

template<class Master_matrix>
template<class Container_type>
inline Base_matrix_with_column_compression_with_row_access<Master_matrix>::Base_matrix_with_column_compression_with_row_access(const std::vector<Container_type> &columns)
	: /*columnToRep_(columns.size()),*/
	  columnClasses_(columns.size()),
	  repToColumn_(columns.size(), nullptr),
	  nextColumnIndex_(0)
{
	if constexpr (!Master_matrix::Option_list::has_removable_rows){
		rows_.resize(columns.size());
	}

	for (const Container_type& c : columns){
		insert_column(c);
	}
}

template<class Master_matrix>
inline Base_matrix_with_column_compression_with_row_access<Master_matrix>::Base_matrix_with_column_compression_with_row_access(unsigned int numberOfColumns)
	: /*columnToRep_(numberOfColumns),*/
	  columnClasses_(numberOfColumns),
	  repToColumn_(numberOfColumns, nullptr),
	  nextColumnIndex_(0)
{
	if constexpr (!Master_matrix::Option_list::has_removable_rows){
		rows_.resize(numberOfColumns);
	}
}

template<class Master_matrix>
inline Base_matrix_with_column_compression_with_row_access<Master_matrix>::Base_matrix_with_column_compression_with_row_access(const Base_matrix_with_column_compression_with_row_access &matrixToCopy)
	: /*columnToRep_(matrixToCopy.columnToRep_.size()),*/
	  columnClasses_(matrixToCopy.columnClasses_),
	  repToColumn_(matrixToCopy.repToColumn_.size(), nullptr),
	  nextColumnIndex_(0)
{
	if constexpr (!Master_matrix::Option_list::has_removable_rows){
		rows_.resize(matrixToCopy.rows_.size());
	}

	for (const Column_type* col : matrixToCopy.repToColumn_){
		if (col != nullptr){
//			repToColumn_[nextColumnIndex_] = new Column_type(*col, col->get_column_index(), rows_);
			repToColumn_[nextColumnIndex_] = columnPool_.construct(*col, col->get_column_index(), rows_);
//			columnToRep_.emplace_hint(columnToRep_.end(), repToColumn_[nextColumnIndex_], nextColumnIndex_);
			columnToRep_.insert(columnToRep_.end(), *repToColumn_[nextColumnIndex_]);
			repToColumn_[nextColumnIndex_]->set_rep(nextColumnIndex_);
		}
		nextColumnIndex_++;
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
//	for (auto& p : columnToRep_){
//		p.first->set_rows(&rows_);
//	}
	for (Column_type& col : columnToRep_){
		col.set_rows(&rows_);
	}
}

//int count = 0;
template<class Master_matrix>
inline Base_matrix_with_column_compression_with_row_access<Master_matrix>::~Base_matrix_with_column_compression_with_row_access()
{
//	for (auto col : repToColumn_){
//		if (col != nullptr)
////			delete col;
//			columnPool_.destroy(col);
//	}
	columnToRep_.clear_and_dispose(delete_disposer());
}

template<class Master_matrix>
template<class Container_type>
inline void Base_matrix_with_column_compression_with_row_access<Master_matrix>::insert_column(const Container_type &column)
{
//	columnClasses_.initialize(nextColumnIndex_);
//	if (rows_.size() <= nextColumnIndex_) rows_.resize(nextColumnIndex_ + 1);
//	repToColumn_.emplace_back(nextColumnIndex_, column, rows_);
//	_insert_column(nextColumnIndex_);
////	auto p = repToColumn_.emplace(nextColumnIndex_, Column_type(nextColumnIndex_, column, rows_));
////	_insert_column(p.first, nextColumnIndex_);
//	nextColumnIndex_++;

	columnClasses_.initialize(nextColumnIndex_);
	if constexpr (!Master_matrix::Option_list::has_removable_rows) {
		if (rows_.size() <= nextColumnIndex_) rows_.resize(nextColumnIndex_ + 1);
	}
//	repToColumn_[nextColumnIndex_] = new Column_type(nextColumnIndex_, column, rows_);
	if (repToColumn_.size() == nextColumnIndex_){
		repToColumn_.push_back(columnPool_.construct(nextColumnIndex_, column, rows_));
	} else {
		repToColumn_[nextColumnIndex_] = columnPool_.construct(nextColumnIndex_, column, rows_);
	}
	_insert_column(nextColumnIndex_);
//	Column_type& col = *repToColumn_[nextColumnIndex_];

////	auto res = columnToRep_.emplace_hint(columnToRep_.end(), &col, nextColumnIndex_);
////	if (res->second != nextColumnIndex_){
////		_insert_double_column(nextColumnIndex_, res);
////	}
//	columnToRep_.insert(columnToRep_.end(), col);
//	repToColumn_[nextColumnIndex_]->set_rep(nextColumnIndex_);
	nextColumnIndex_++;
}

//void insert_column(index rowIndex);
//void insert_column(index rowIndex, const Field_element_type& coefficient);

template<class Master_matrix>
template<class Boundary_type>
inline void Base_matrix_with_column_compression_with_row_access<Master_matrix>::insert_boundary(const Boundary_type &boundary)
{
	insert_column(boundary);
}

template<class Master_matrix>
inline const typename Base_matrix_with_column_compression_with_row_access<Master_matrix>::Column_type &Base_matrix_with_column_compression_with_row_access<Master_matrix>::get_column(index columnIndex)
{
	auto col = repToColumn_[columnClasses_.find(columnIndex)];
	if (col == nullptr) return empty_column_;
	return *col;
//	return *repToColumn_[columnClasses_.find(columnIndex)];
}

//template<class Master_matrix>
//inline typename Base_matrix_with_column_compression_with_row_access<Master_matrix>::Row_type& Base_matrix_with_column_compression_with_row_access<Master_matrix>::get_row(index rowIndex)
//{
//	if constexpr (Master_matrix::Option_list::has_removable_rows) {
//		return rows_[rowIndex];
//	} else {
//		return rows_.at(rowIndex);
//	}
//}

template<class Master_matrix>
inline const typename Base_matrix_with_column_compression_with_row_access<Master_matrix>::Row_type& Base_matrix_with_column_compression_with_row_access<Master_matrix>::get_row(index rowIndex) const
{
	if constexpr (!Master_matrix::Option_list::has_removable_rows) {
		return rows_[rowIndex];
	} else {
		return rows_.at(rowIndex);
	}
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
	rows_.erase(rowIndex);
}

template<class Master_matrix>
inline unsigned int Base_matrix_with_column_compression_with_row_access<Master_matrix>::get_number_of_columns() const
{
	return nextColumnIndex_;
//	return columnToRep_.size();
}

template<class Master_matrix>
inline void Base_matrix_with_column_compression_with_row_access<Master_matrix>::add_to(index sourceColumnIndex, index targetColumnIndex)
{
	//handle case where targetRep == sourceRep?
	index targetRep = columnClasses_.find(targetColumnIndex);
	Column_type& target = *repToColumn_[targetRep];
//	columnToRep_.erase(&target);
	columnToRep_.erase(target);
	target += get_column(sourceColumnIndex);
	_insert_column(targetRep);
}

template<class Master_matrix>
template<class Cell_range>
//inline void Base_matrix_with_column_compression_with_row_access<Master_matrix>::add_to(const Column_type& sourceColumn, index targetColumnIndex)
inline void Base_matrix_with_column_compression_with_row_access<Master_matrix>::add_to(const Cell_range& sourceColumn, index targetColumnIndex)
{
	index targetRep = columnClasses_.find(targetColumnIndex);
	Column_type& target = *repToColumn_[targetRep];
//	columnToRep_.erase(&target);
	columnToRep_.erase(target);
	target += sourceColumn;
	_insert_column(targetRep);
}

template<class Master_matrix>
template<class Cell_range>
//inline void Base_matrix_with_column_compression_with_row_access<Master_matrix>::add_to(const Column_type& sourceColumn, const Field_element_type& coefficient, index targetColumnIndex)
inline void Base_matrix_with_column_compression_with_row_access<Master_matrix>::add_to(const Cell_range& sourceColumn, const Field_element_type& coefficient, index targetColumnIndex)
{
	index targetRep = columnClasses_.find(targetColumnIndex);
	Column_type& target = *repToColumn_[targetRep];
//	columnToRep_.erase(&target);
	columnToRep_.erase(target);
	target.multiply_and_add(coefficient, sourceColumn);
	_insert_column(targetRep);
}

template<class Master_matrix>
template<class Cell_range>
//inline void Base_matrix_with_column_compression_with_row_access<Master_matrix>::add_to(const Field_element_type& coefficient, const Column_type& sourceColumn, index targetColumnIndex)
//inline void Base_matrix_with_column_compression_with_row_access<Master_matrix>::add_to(const Field_element_type& coefficient, const typename Master_matrix::Column_type& sourceColumn, index targetColumnIndex)
inline void Base_matrix_with_column_compression_with_row_access<Master_matrix>::add_to(const Field_element_type& coefficient, const Cell_range& sourceColumn, index targetColumnIndex)
{
//	index targetRep = columnClasses_.find(targetColumnIndex);
////	auto itTarget = repToColumn_.find(targetRep);
////	Column_type& target = itTarget->second;
////	columnToRep_.erase(&target);
////	target.multiply_and_add(sourceColumn, coefficient);
////	_insert_column(itTarget, targetRep);
//	Column_type& target = *repToColumn_[targetRep];
//	columnToRep_.erase(&target);
//	target.multiply_and_add(sourceColumn, coefficient);

	index targetRep = columnClasses_.find(targetColumnIndex);
	Column_type& target = *repToColumn_[targetRep];
//	columnToRep_.erase(&target);
	columnToRep_.erase(target);
	target.multiply_and_add(sourceColumn, coefficient);
	_insert_column(targetRep);

//	Column_type& col = *repToColumn_[targetRep];
//	if (col.is_empty()){
//		auto res = _insert_try_empty(targetRep, col);
//		if (res->second != targetRep){
//			_insert_double(res, targetRep, col);
//		}
//		return;
//	}
//	auto res = _insert_try(targetRep, col);
//	if (res->second != targetRep){
//		_insert_double(res, targetRep, col);
//	}
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
	auto col = repToColumn_[columnClasses_.find(columnIndex)];
	if (col == nullptr) return true;
	return !col->is_non_zero(rowIndex);
//	return !(repToColumn_[columnClasses_.find(columnIndex)]->is_non_zero(rowIndex));
}

template<class Master_matrix>
inline bool Base_matrix_with_column_compression_with_row_access<Master_matrix>::is_zero_column(index columnIndex)
{
	auto col = repToColumn_[columnClasses_.find(columnIndex)];
	if (col == nullptr) return true;
	return col->is_empty();
//	return repToColumn_[columnClasses_.find(columnIndex)]->is_empty();
}

template<class Master_matrix>
inline Base_matrix_with_column_compression_with_row_access<Master_matrix> &
Base_matrix_with_column_compression_with_row_access<Master_matrix>::operator=(const Base_matrix_with_column_compression_with_row_access& other)
{
	for (auto col : repToColumn_){
		if (col != nullptr){
//			delete col;
			columnPool_.destroy(col);
			col = nullptr;
		}
	}
	columnClasses_ = other.columnClasses_;
	if constexpr (!Master_matrix::Option_list::has_removable_rows)
			rows_.resize(other.rows_.size());
	columnToRep_.reserve(other.columnToRep_.size());
	repToColumn_.resize(other.repToColumn_.size(), nullptr);
	nextColumnIndex_ = 0;
	for (const Column_type* col : other.repToColumn_){
//		repToColumn_[nextColumnIndex_] = new Column_type(*col, col->get_column_index(), rows_);
		repToColumn_[nextColumnIndex_] = columnPool_.construct(*col, col->get_column_index(), rows_);
//		columnToRep_.emplace_hint(columnToRep_.end(), repToColumn_[nextColumnIndex_], nextColumnIndex_);
		columnToRep_.insert(columnToRep_.end(), *repToColumn_[nextColumnIndex_]);
		repToColumn_[nextColumnIndex_]->set_rep(nextColumnIndex_);
		nextColumnIndex_++;
	}
	return *this;
}

template<class Master_matrix>
inline void Base_matrix_with_column_compression_with_row_access<Master_matrix>::print()
{
	std::cout << "Compressed_matrix:\n";
//	for (auto& p : columnToRep_){
	for (Column_type& col : columnToRep_){
//		Column_type* col = p.first;
		for (auto e : col->get_content(nextColumnIndex_)){
			if (e == 0u) std::cout << "- ";
			else std::cout << e << " ";
		}
		std::cout << "(";
		for (unsigned int i = 0; i < nextColumnIndex_; ++i){
//			if (columnClasses_.find(i) == p.second)
			if (columnClasses_.find(i) == col.get_rep())
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
inline void Base_matrix_with_column_compression_with_row_access<Master_matrix>::_insert_column(index columnIndex)
{
	Column_type& col = *repToColumn_[columnIndex];

	if (col.is_empty()){
//		auto res = columnToRep_.emplace_hint(columnToRep_.begin(), &col, columnIndex);
//		if (res->second != columnIndex){
//			_insert_double_column(columnIndex, res);
//		}
		columnPool_.destroy(&col);  // delete curr_col;
		repToColumn_[columnIndex] = nullptr;
		return;
	}

//	auto res = columnToRep_.emplace(&col, columnIndex);
	repToColumn_[columnIndex]->set_rep(columnIndex);
	auto res = columnToRep_.insert(col);
//	if (res.first->second != columnIndex){
//		_insert_double_column(columnIndex, res.first);
//	}
	if (res.first->get_rep() != columnIndex){
		_insert_double_column(columnIndex, res.first);
	}

//		auto res = columnToRep_.emplace_hint(columnToRep_.end(), &col, columnIndex);
//		if (res->second != columnIndex){
//			_insert_double_column(columnIndex, res);
//		}
}

template<class Master_matrix>
inline void Base_matrix_with_column_compression_with_row_access<Master_matrix>::_insert_double_column(
		index columnIndex, typename col_dict_type::iterator& doubleIt)
{
//	index doubleRep = doubleIt->second;
	index doubleRep = doubleIt->get_rep();
	index newRep = columnClasses_.merge(columnIndex, doubleRep);

//	delete repToColumn_[columnIndex];
	columnPool_.destroy(repToColumn_[columnIndex]);
	repToColumn_[columnIndex] = nullptr;

	if (newRep == columnIndex){
		std::swap(repToColumn_[doubleRep], repToColumn_[columnIndex]);
//		doubleIt->second = columnIndex;
		doubleIt->set_rep(columnIndex);
	}
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // BASE_MATRIX_1010_H

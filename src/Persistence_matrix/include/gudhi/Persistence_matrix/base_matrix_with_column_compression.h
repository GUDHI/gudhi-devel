/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022-23 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_BASE_MATRIX_COMPRESSION_H
#define PM_BASE_MATRIX_COMPRESSION_H

#include <iostream>	//print() only
#include <vector>
#include <set>
#include <utility>	//std::swap, std::move & std::exchange

#include <boost/intrusive/set.hpp>
#include <boost/pending/disjoint_sets.hpp>

#include <gudhi/Simple_object_pool.h>
// #include <gudhi/Persistence_matrix/union_find.h>

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class Base_matrix_with_column_compression : protected Master_matrix::Matrix_row_access_option
{
public:
	using index = typename Master_matrix::index;
	using dimension_type = typename Master_matrix::dimension_type;
	using Field_element_type = typename Master_matrix::Field_type;
	using Row_type = typename Master_matrix::Row_type;
	// using Column_type = typename Master_matrix::Column_type;

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
		template<class Container_type, class Row_container_type>
		Column_type(index columnIndex, const Container_type& nonZeroRowIndices, Row_container_type &rowContainer)
			: Base(columnIndex, nonZeroRowIndices, rowContainer)
		{}
		template<class Container_type>
		Column_type(const Container_type& nonZeroRowIndices, dimension_type dimension)
			: Base(nonZeroRowIndices, dimension)
		{}
		template<class Container_type, class Row_container_type>
		Column_type(index columnIndex, const Container_type& nonZeroRowIndices, dimension_type dimension, Row_container_type &rowContainer)
			: Base(columnIndex, nonZeroRowIndices, dimension, rowContainer)
		{}
		Column_type(const Column_type& column)
			: Base(static_cast<const Base&>(column))
		{}
		template<class Row_container_type>
		Column_type(const Column_type& column, index columnIndex, Row_container_type &rowContainer)
			: Base(static_cast<const Base&>(column), columnIndex, rowContainer)
		{}
		Column_type(Column_type&& column) noexcept
			: Base(std::move(static_cast<Base&>(column)))
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
				return std::hash<Base>()(static_cast<Base>(c));
			}
		};

	private:
		index rep_;
	};

	Base_matrix_with_column_compression();
	template<class Container_type>
	Base_matrix_with_column_compression(const std::vector<Container_type>& columns);
	Base_matrix_with_column_compression(unsigned int numberOfColumns);
	Base_matrix_with_column_compression(const Base_matrix_with_column_compression& matrixToCopy);
	Base_matrix_with_column_compression(Base_matrix_with_column_compression&& other) noexcept;
	~Base_matrix_with_column_compression();

	template<class Container_type>
	void insert_column(const Container_type& column);
	template<class Boundary_type>
	void insert_boundary(const Boundary_type& boundary, dimension_type dim = -1);	//same as insert_column
	const Column_type& get_column(index columnIndex);	//non const because of path compression in union-find
	//get_row(rowIndex) --> simplex ID (=/= columnIndex)
	const Row_type& get_row(index rowIndex) const;
	void erase_row(index rowIndex);		//assumes the row is empty, just thought as index a cleanup

	unsigned int get_number_of_columns() const;

	template<class Cell_range_or_column_index>
	void add_to(const Cell_range_or_column_index& sourceColumn, index targetColumnIndex);
	template<class Cell_range_or_column_index>
	void add_to(const Cell_range_or_column_index& sourceColumn, const Field_element_type& coefficient, index targetColumnIndex);
	template<class Cell_range_or_column_index>
	void add_to(const Field_element_type& coefficient, const Cell_range_or_column_index& sourceColumn, index targetColumnIndex);

	bool is_zero_cell(index columnIndex, index rowIndex);
	bool is_zero_column(index columnIndex);

	Base_matrix_with_column_compression& operator=(const Base_matrix_with_column_compression& other);
	friend void swap(Base_matrix_with_column_compression& matrix1, Base_matrix_with_column_compression& matrix2){
		matrix1.columnToRep_.swap(matrix2.columnToRep_);
		swap(matrix1.columnClasses_, matrix2.columnClasses_);
		matrix1.repToColumn_.swap(matrix2.repToColumn_);
		std::swap(matrix1.nextColumnIndex_, matrix2.nextColumnIndex_);

		if constexpr (Master_matrix::Option_list::has_row_access){
			swap(static_cast<typename Master_matrix::Matrix_row_access_option&>(matrix1),
				 static_cast<typename Master_matrix::Matrix_row_access_option&>(matrix2));
			for (auto& col : matrix1.columnToRep_){
				col.set_rows(&matrix1.rows_);
			}
			for (auto& col : matrix2.columnToRep_){
				col.set_rows(&matrix2.rows_);
			}
		}
	}

	void print();  //for debug

private:
	//The disposer object function for boost intrusive container
	struct delete_disposer
	{
		void operator()(Column_type *delete_this){
			columnPool_.destroy(delete_this);
		}
	};

	using ra_opt = typename Master_matrix::Matrix_row_access_option;
	using cell_rep_type = typename std::conditional<
								Master_matrix::Option_list::is_z2,
								index,
								std::pair<index,Field_element_type>
							>::type;
	using tmp_column_type = typename std::conditional<
								Master_matrix::Option_list::is_z2,
								std::set<index>,
								std::set<std::pair<index,Field_element_type>,typename Master_matrix::CellPairComparator>
							>::type;
	using col_dict_type = boost::intrusive::set<Column_type, boost::intrusive::constant_time_size<false> >;

	col_dict_type columnToRep_;
	boost::disjoint_sets_with_storage<> columnClasses_;
	// Union_find columnClasses_;
	std::vector<Column_type*> repToColumn_;
	unsigned int nextColumnIndex_;
	inline static Simple_object_pool<Column_type> columnPool_;
	inline static const Column_type empty_column_;

	void _insert_column(index columnIndex);
	void _insert_double_column(index columnIndex, typename col_dict_type::iterator& doubleIt);
};

template<class Master_matrix>
inline Base_matrix_with_column_compression<Master_matrix>::Base_matrix_with_column_compression()
	: ra_opt(),
	  nextColumnIndex_(0)
{}

template<class Master_matrix>
template<class Container_type>
inline Base_matrix_with_column_compression<Master_matrix>::Base_matrix_with_column_compression(const std::vector<Container_type> &columns)
	: ra_opt(columns.size()),
	  columnClasses_(columns.size()),
	  repToColumn_(columns.size(), nullptr),
	  nextColumnIndex_(0)
{
	for (const Container_type& c : columns){
		insert_column(c);
	}
}

template<class Master_matrix>
inline Base_matrix_with_column_compression<Master_matrix>::Base_matrix_with_column_compression(unsigned int numberOfColumns)
	: ra_opt(numberOfColumns),
	  columnClasses_(numberOfColumns),
	  repToColumn_(numberOfColumns, nullptr),
	  nextColumnIndex_(0)
{}

template<class Master_matrix>
inline Base_matrix_with_column_compression<Master_matrix>::Base_matrix_with_column_compression(const Base_matrix_with_column_compression &matrixToCopy)
	: ra_opt(static_cast<const ra_opt&>(matrixToCopy)),
	  columnClasses_(matrixToCopy.columnClasses_),
	  repToColumn_(matrixToCopy.repToColumn_.size(), nullptr),
	  nextColumnIndex_(0)
{
	for (const Column_type* col : matrixToCopy.repToColumn_){
		if (col != nullptr){
			if constexpr (Master_matrix::Option_list::has_row_access){
				repToColumn_[nextColumnIndex_] = columnPool_.construct(*col, col->get_column_index(), ra_opt::rows_);
			} else {
				repToColumn_[nextColumnIndex_] = columnPool_.construct(*col);
			}
			columnToRep_.insert(columnToRep_.end(), *repToColumn_[nextColumnIndex_]);
			repToColumn_[nextColumnIndex_]->set_rep(nextColumnIndex_);
		}
		nextColumnIndex_++;
	}
}

template<class Master_matrix>
inline Base_matrix_with_column_compression<Master_matrix>::Base_matrix_with_column_compression(Base_matrix_with_column_compression &&other) noexcept
	: ra_opt(std::move(static_cast<ra_opt&>(other))),
	  columnToRep_(std::move(other.columnToRep_)),
	  columnClasses_(std::move(other.columnClasses_)),
	  repToColumn_(std::move(other.repToColumn_)),
	  nextColumnIndex_(std::exchange(other.nextColumnIndex_, 0))
{
	if constexpr (Master_matrix::Option_list::has_row_access){
		for (Column_type& col : columnToRep_){
			col.set_rows(&this->rows_);
		}
	}
}

template<class Master_matrix>
inline Base_matrix_with_column_compression<Master_matrix>::~Base_matrix_with_column_compression()
{
	columnToRep_.clear_and_dispose(delete_disposer());
}

template<class Master_matrix>
template<class Container_type>
inline void Base_matrix_with_column_compression<Master_matrix>::insert_column(const Container_type &column)
{
	insert_boundary(column);
}

template<class Master_matrix>
template<class Boundary_type>
inline void Base_matrix_with_column_compression<Master_matrix>::insert_boundary(const Boundary_type &boundary, dimension_type dim)
{
	if (dim == -1) dim = boundary.size() == 0 ? 0 : boundary.size() - 1;

	if constexpr (Master_matrix::Option_list::has_row_access && !Master_matrix::Option_list::has_removable_rows){
		if (boundary.begin() != boundary.end()){
			unsigned int pivot;
			if constexpr (Master_matrix::Option_list::is_z2){
				pivot = *std::prev(boundary.end());
			} else {
				pivot = std::prev(boundary.end())->first;
			}
			if (ra_opt::rows_.size() <= pivot) ra_opt::rows_.resize(pivot + 1);
		}
	}

	// columnClasses_.initialize(nextColumnIndex_);
	if (repToColumn_.size() == nextColumnIndex_){
		columnClasses_.link(nextColumnIndex_, nextColumnIndex_);	//could perhaps be avoided, if find_set returns something special when it does not find
		if constexpr (Master_matrix::Option_list::has_row_access){
			repToColumn_.push_back(columnPool_.construct(nextColumnIndex_, boundary, dim, ra_opt::rows_));
		} else {
			repToColumn_.push_back(columnPool_.construct(boundary, dim));
		}
	} else {
		if constexpr (Master_matrix::Option_list::has_row_access){
			repToColumn_[nextColumnIndex_] = columnPool_.construct(nextColumnIndex_, boundary, dim, ra_opt::rows_);
		} else {
			repToColumn_[nextColumnIndex_] = columnPool_.construct(boundary, dim);
		}
	}
	_insert_column(nextColumnIndex_);

	nextColumnIndex_++;
}

template<class Master_matrix>
inline const typename Base_matrix_with_column_compression<Master_matrix>::Column_type& 
Base_matrix_with_column_compression<Master_matrix>::get_column(index columnIndex)
{
	// auto col = repToColumn_[columnClasses_.find(columnIndex)];
	auto col = repToColumn_[columnClasses_.find_set(columnIndex)];
	if (col == nullptr) return empty_column_;
	return *col;
}

template<class Master_matrix>
inline const typename Base_matrix_with_column_compression<Master_matrix>::Row_type&
Base_matrix_with_column_compression<Master_matrix>::get_row(index columnIndex) const
{
	static_assert(Master_matrix::Option_list::has_row_access, "Row access has to be enabled for this method.");

	return ra_opt::get_row(columnIndex);
}

template<class Master_matrix>
inline void Base_matrix_with_column_compression<Master_matrix>::erase_row(index rowIndex)
{	
	if constexpr (Master_matrix::Option_list::has_row_access && Master_matrix::Option_list::has_removable_rows){
		ra_opt::erase_row(rowIndex);
	}
}

template<class Master_matrix>
inline unsigned int Base_matrix_with_column_compression<Master_matrix>::get_number_of_columns() const
{
	return nextColumnIndex_;
}

template<class Master_matrix>
template<class Cell_range_or_column_index>
inline void Base_matrix_with_column_compression<Master_matrix>::add_to(const Cell_range_or_column_index& sourceColumn, index targetColumnIndex)
{
	//handle case where targetRep == sourceRep?
	// index targetRep = columnClasses_.find(targetColumnIndex);
	index targetRep = columnClasses_.find_set(targetColumnIndex);
	Column_type& target = *repToColumn_[targetRep];
	columnToRep_.erase(target);
	if constexpr (std::is_integral_v<Cell_range_or_column_index>){
		target += get_column(sourceColumn);
	} else {
		target += sourceColumn;
	}
	_insert_column(targetRep);
}

template<class Master_matrix>
template<class Cell_range_or_column_index>
inline void Base_matrix_with_column_compression<Master_matrix>::add_to(const Cell_range_or_column_index& sourceColumn, const Field_element_type& coefficient, index targetColumnIndex)
{
	//handle case where targetRep == sourceRep?
	// index targetRep = columnClasses_.find(targetColumnIndex);
	index targetRep = columnClasses_.find_set(targetColumnIndex);
	Column_type& target = *repToColumn_[targetRep];
	columnToRep_.erase(target);
	if constexpr (std::is_integral_v<Cell_range_or_column_index>){
		target.multiply_and_add(coefficient, get_column(sourceColumn));
	} else {
		target.multiply_and_add(coefficient, sourceColumn);
	}
	_insert_column(targetRep);
}

template<class Master_matrix>
template<class Cell_range_or_column_index>
inline void Base_matrix_with_column_compression<Master_matrix>::add_to(const Field_element_type& coefficient, const Cell_range_or_column_index& sourceColumn, index targetColumnIndex)
{
	//handle case where targetRep == sourceRep?
	// index targetRep = columnClasses_.find(targetColumnIndex);
	index targetRep = columnClasses_.find_set(targetColumnIndex);
	Column_type& target = *repToColumn_[targetRep];
	columnToRep_.erase(target);
	if constexpr (std::is_integral_v<Cell_range_or_column_index>){
		target.multiply_and_add(get_column(sourceColumn), coefficient);
	} else {
		target.multiply_and_add(sourceColumn, coefficient);
	}
	_insert_column(targetRep);
}

template<class Master_matrix>
inline bool Base_matrix_with_column_compression<Master_matrix>::is_zero_cell(index columnIndex, index rowIndex)
{
	// auto col = repToColumn_[columnClasses_.find(columnIndex)];
	auto col = repToColumn_[columnClasses_.find_set(columnIndex)];
	if (col == nullptr) return true;
	return !col->is_non_zero(rowIndex);
}

template<class Master_matrix>
inline bool Base_matrix_with_column_compression<Master_matrix>::is_zero_column(index columnIndex)
{
	// auto col = repToColumn_[columnClasses_.find(columnIndex)];
	auto col = repToColumn_[columnClasses_.find_set(columnIndex)];
	if (col == nullptr) return true;
	return col->is_empty();
}

template<class Master_matrix>
inline Base_matrix_with_column_compression<Master_matrix>&
Base_matrix_with_column_compression<Master_matrix>::operator=(const Base_matrix_with_column_compression& other)
{
	for (auto col : repToColumn_){
		if (col != nullptr){
			columnPool_.destroy(col);
			col = nullptr;
		}
	}
	ra_opt::operator=(other);
	columnClasses_ = other.columnClasses_;
	columnToRep_.reserve(other.columnToRep_.size());
	repToColumn_.resize(other.repToColumn_.size(), nullptr);
	nextColumnIndex_ = 0;
	for (const Column_type* col : other.repToColumn_){
		if constexpr (Master_matrix::Option_list::has_row_access){
			repToColumn_[nextColumnIndex_] = columnPool_.construct(*col, col->get_column_index(), ra_opt::rows_);
		} else {
			repToColumn_[nextColumnIndex_] = columnPool_.construct(*col);
		}
		columnToRep_.insert(columnToRep_.end(), *repToColumn_[nextColumnIndex_]);
		repToColumn_[nextColumnIndex_]->set_rep(nextColumnIndex_);
		nextColumnIndex_++;
	}
	return *this;
}

template<class Master_matrix>
inline void Base_matrix_with_column_compression<Master_matrix>::print()
{
	std::cout << "Compressed_matrix:\n";
	for (Column_type& col : columnToRep_){
		for (auto e : col->get_content(nextColumnIndex_)){
			if (e == 0u) std::cout << "- ";
			else std::cout << e << " ";
		}
		std::cout << "(";
		for (unsigned int i = 0; i < nextColumnIndex_; ++i){
			// if (columnClasses_.find(i) == col.get_rep())
			if (columnClasses_.find_set(i) == col.get_rep())
				std::cout << i << " ";
		}
		std::cout << ")\n";
	}
	std::cout << "\n";
	std::cout << "Row Matrix:\n";
	for (unsigned int i = 0; i < ra_opt::rows_.size(); ++i){
		const Row_type& row = ra_opt::rows_[i];
		for (const auto &cell : row){
			std::cout << cell.get_column_index() << " ";
		}
		std::cout << "(" << i << ")\n";
	}
	std::cout << "\n";
}

template<class Master_matrix>
inline void Base_matrix_with_column_compression<Master_matrix>::_insert_column(index columnIndex)
{
	Column_type& col = *repToColumn_[columnIndex];

	if (col.is_empty()){
		columnPool_.destroy(&col);  // delete curr_col;
		repToColumn_[columnIndex] = nullptr;
		return;
	}

	col.set_rep(columnIndex);
	auto res = columnToRep_.insert(col);
	if (res.first->get_rep() != columnIndex){
		_insert_double_column(columnIndex, res.first);
	}
}

template<class Master_matrix>
inline void Base_matrix_with_column_compression<Master_matrix>::_insert_double_column(
		index columnIndex, typename col_dict_type::iterator& doubleIt)
{
	index doubleRep = doubleIt->get_rep();
	// index newRep = columnClasses_.merge(columnIndex, doubleRep);
	columnClasses_.link(columnIndex, doubleRep);	//both should be representatives
	index newRep = columnClasses_.find_set(columnIndex);

	columnPool_.destroy(repToColumn_[columnIndex]);
	repToColumn_[columnIndex] = nullptr;

	if (newRep == columnIndex){
		std::swap(repToColumn_[doubleRep], repToColumn_[columnIndex]);
		doubleIt->set_rep(columnIndex);
	}
}


} //namespace persistence_matrix
} //namespace Gudhi

#endif // PM_BASE_MATRIX_COMPRESSION_H

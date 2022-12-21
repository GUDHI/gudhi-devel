/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef Z2_INTRUSIVE_LIST_COLUMN_H
#define Z2_INTRUSIVE_LIST_COLUMN_H

#include <iostream>
#include <vector>

#include <boost/intrusive/list.hpp>

#include "../utilities/utilities.h"
#include "cell.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Cell_type, class Column_pairing_option, class Row_access_option>
class Z2_intrusive_list_column : public Column_pairing_option, public Row_access_option
{
public:
//	using Cell = Z2_intrusive_list_cell;
	using Cell = Cell_type;
	using Column_type = boost::intrusive::list <
							Cell
						  , boost::intrusive::constant_time_size<false>
						  , boost::intrusive::base_hook< base_hook_matrix_list_column >  >;
	using iterator = typename Column_type::iterator;
	using const_iterator = typename Column_type::const_iterator;

	Z2_intrusive_list_column();
	template<class Container_type>
	Z2_intrusive_list_column(const Container_type& nonZeroRowIndices);
	template<class Container_type>
	Z2_intrusive_list_column(const Container_type& nonZeroRowIndices, dimension_type dimension);
	template<class Row_container_type>
	Z2_intrusive_list_column(index columnIndex, Row_container_type &rowContainer);
	template<class Container_type, class Row_container_type>
	Z2_intrusive_list_column(index columnIndex, const Container_type& nonZeroRowIndices, Row_container_type &rowContainer);
	template<class Container_type, class Row_container_type>
	Z2_intrusive_list_column(index columnIndex, const Container_type& nonZeroRowIndices, dimension_type dimension, Row_container_type &rowContainer);
	Z2_intrusive_list_column(const Z2_intrusive_list_column& column);
	Z2_intrusive_list_column(const Z2_intrusive_list_column& column, index columnIndex);
	Z2_intrusive_list_column(Z2_intrusive_list_column&& column) noexcept;
	~Z2_intrusive_list_column();

	std::vector<bool> get_content(unsigned int columnLength) const;
	bool is_non_zero(index rowIndex) const;
	bool is_empty() const;
	dimension_type get_dimension() const;

	iterator begin() noexcept;
	const_iterator begin() const noexcept;
	iterator end() noexcept;
	const_iterator end() const noexcept;

	Z2_intrusive_list_column& operator+=(Z2_intrusive_list_column const &column);
	friend Z2_intrusive_list_column operator+(Z2_intrusive_list_column column1, Z2_intrusive_list_column const& column2){
		column1 += column2;
		return column1;
	}

	Z2_intrusive_list_column& operator*=(unsigned int v);
	friend Z2_intrusive_list_column operator*(Z2_intrusive_list_column column, unsigned int const& v){
		column *= v;
		return column;
	}
	friend Z2_intrusive_list_column operator*(unsigned int const& v, Z2_intrusive_list_column column){
		column *= v;
		return column;
	}

	Z2_intrusive_list_column& operator=(const Z2_intrusive_list_column& other);

	friend void swap(Z2_intrusive_list_column& col1, Z2_intrusive_list_column& col2){
		swap(static_cast<Column_pairing_option&>(col1),
			 static_cast<Column_pairing_option&>(col2));
		swap(static_cast<Row_access_option&>(col1),
			 static_cast<Row_access_option&>(col2));
		std::swap(col1.dim_, col2.dim_);
		col1.column_.swap(col2.column_);
	}

protected:
	dimension_type dim_;
	Column_type column_;

	void _delete_cell(iterator& it);
	void _insert_cell(index rowIndex, const iterator& position);

private:
	//Cloner object function
	struct new_cloner
	{
	   Cell *operator()(const Cell &clone_this)
	   {  return new Cell(clone_this);  }
	};

	//The disposer object function
	struct delete_disposer
	{
	   void operator()(Cell *delete_this)
	   {  delete delete_this;  }
	};
};

template<class Cell_type, class Column_pairing_option, class Row_access_option>
inline Z2_intrusive_list_column<Cell_type,Column_pairing_option,Row_access_option>::Z2_intrusive_list_column() : dim_(0)
{
	static_assert(!Row_access_option::isActive_, "When row access option enabled, a row container has to be provided.");
}

template<class Cell_type, class Column_pairing_option, class Row_access_option>
template<class Container_type>
inline Z2_intrusive_list_column<Cell_type,Column_pairing_option,Row_access_option>::Z2_intrusive_list_column(const Container_type &nonZeroRowIndices)
	: dim_(nonZeroRowIndices.size() == 0 ? 0 : nonZeroRowIndices.size() - 1)
{
	static_assert(!Row_access_option::isActive_, "When row access option enabled, a row container has to be provided.");

	for (index id : nonZeroRowIndices){
		_insert_cell(id, column_.end());
	}
}

template<class Cell_type, class Column_pairing_option, class Row_access_option>
template<class Container_type>
inline Z2_intrusive_list_column<Cell_type,Column_pairing_option,Row_access_option>::Z2_intrusive_list_column(const Container_type &nonZeroRowIndices, dimension_type dimension)
	: dim_(dimension)
{
	static_assert(!Row_access_option::isActive_, "When row access option enabled, a row container has to be provided.");

	for (index id : nonZeroRowIndices){
		_insert_cell(id, column_.end());
	}
}

template<class Cell_type, class Column_pairing_option, class Row_access_option>
template<class Row_container_type>
inline Z2_intrusive_list_column<Cell_type,Column_pairing_option,Row_access_option>::Z2_intrusive_list_column(
		index columnIndex, Row_container_type &rowContainer)
	: Row_access_option(columnIndex, rowContainer), dim_(0)
{}

template<class Cell_type, class Column_pairing_option, class Row_access_option>
template<class Container_type, class Row_container_type>
inline Z2_intrusive_list_column<Cell_type,Column_pairing_option,Row_access_option>::Z2_intrusive_list_column(
		index columnIndex, const Container_type &nonZeroRowIndices, Row_container_type &rowContainer)
	: Row_access_option(columnIndex, rowContainer), dim_(nonZeroRowIndices.size() == 0 ? 0 : nonZeroRowIndices.size() - 1)
{
	for (index id : nonZeroRowIndices){
		_insert_cell(id, column_.end());
	}
}

template<class Cell_type, class Column_pairing_option, class Row_access_option>
template<class Container_type, class Row_container_type>
inline Z2_intrusive_list_column<Cell_type,Column_pairing_option,Row_access_option>::Z2_intrusive_list_column(
		index columnIndex, const Container_type &nonZeroRowIndices, dimension_type dimension, Row_container_type &rowContainer)
	: Row_access_option(columnIndex, rowContainer), dim_(dimension)
{
	for (index id : nonZeroRowIndices){
		_insert_cell(id, column_.end());
	}
}

template<class Cell_type, class Column_pairing_option, class Row_access_option>
inline Z2_intrusive_list_column<Cell_type,Column_pairing_option,Row_access_option>::Z2_intrusive_list_column(const Z2_intrusive_list_column &column)
	: Column_pairing_option(column),
	  dim_(column.dim_)
{
	static_assert(!Row_access_option::isActive_,
			"Copy constructor not available when row access option enabled.");

	column_.clone_from(column.column_, new_cloner(), delete_disposer());
}

template<class Cell_type, class Column_pairing_option, class Row_access_option>
inline Z2_intrusive_list_column<Cell_type,Column_pairing_option,Row_access_option>::Z2_intrusive_list_column(
		const Z2_intrusive_list_column &column, index columnIndex)
	: Column_pairing_option(column),
	  Row_access_option(columnIndex, *column.rows_),
	  dim_(column.dim_)
{
	for (const Cell& cell : column.column_){
		_insert_cell(cell.get_row_index(), column_.end());
	}
}

template<class Cell_type, class Column_pairing_option, class Row_access_option>
inline Z2_intrusive_list_column<Cell_type,Column_pairing_option,Row_access_option>::Z2_intrusive_list_column(Z2_intrusive_list_column &&column) noexcept
	: Column_pairing_option(std::move(column)),
	  Row_access_option(std::move(column)),
	  dim_(std::exchange(column.dim_, 0)),
	  column_(std::move(column.column_))
{}

template<class Cell_type, class Column_pairing_option, class Row_access_option>
inline Z2_intrusive_list_column<Cell_type,Column_pairing_option,Row_access_option>::~Z2_intrusive_list_column()
{
	for (iterator c_it = column_.begin(); c_it != column_.end(); ){
		_delete_cell(c_it);
	}
}

template<class Cell_type, class Column_pairing_option, class Row_access_option>
inline std::vector<bool> Z2_intrusive_list_column<Cell_type,Column_pairing_option,Row_access_option>::get_content(unsigned int columnLength) const
{
	std::vector<bool> container(columnLength);
	for (auto it = column_.begin(); it != column_.end() && it->get_row_index() < columnLength; ++it){
		container[it->get_row_index()] = 1;
	}
	return container;
}

template<class Cell_type, class Column_pairing_option, class Row_access_option>
inline bool Z2_intrusive_list_column<Cell_type,Column_pairing_option,Row_access_option>::is_non_zero(index rowIndex) const
{
	for (const Cell& cell : column_)
		if (cell.get_row_index() == rowIndex) return true;

	return false;
}

template<class Cell_type, class Column_pairing_option, class Row_access_option>
inline bool Z2_intrusive_list_column<Cell_type,Column_pairing_option,Row_access_option>::is_empty() const
{
	return column_.empty();
}

template<class Cell_type, class Column_pairing_option, class Row_access_option>
inline dimension_type Z2_intrusive_list_column<Cell_type,Column_pairing_option,Row_access_option>::get_dimension() const
{
	return dim_;
}

template<class Cell_type, class Column_pairing_option, class Row_access_option>
inline typename Z2_intrusive_list_column<Cell_type,Column_pairing_option,Row_access_option>::iterator
Z2_intrusive_list_column<Cell_type,Column_pairing_option,Row_access_option>::begin() noexcept
{
	return column_.begin();
}

template<class Cell_type, class Column_pairing_option, class Row_access_option>
inline typename Z2_intrusive_list_column<Cell_type,Column_pairing_option,Row_access_option>::const_iterator
Z2_intrusive_list_column<Cell_type,Column_pairing_option,Row_access_option>::begin() const noexcept
{
	return column_.begin();
}

template<class Cell_type, class Column_pairing_option, class Row_access_option>
inline typename Z2_intrusive_list_column<Cell_type,Column_pairing_option,Row_access_option>::iterator
Z2_intrusive_list_column<Cell_type,Column_pairing_option,Row_access_option>::end() noexcept
{
	return column_.end();
}

template<class Cell_type, class Column_pairing_option, class Row_access_option>
inline typename Z2_intrusive_list_column<Cell_type,Column_pairing_option,Row_access_option>::const_iterator
Z2_intrusive_list_column<Cell_type,Column_pairing_option,Row_access_option>::end() const noexcept
{
	return column_.end();
}

template<class Cell_type, class Column_pairing_option, class Row_access_option>
inline Z2_intrusive_list_column<Cell_type,Column_pairing_option,Row_access_option> &Z2_intrusive_list_column<Cell_type,Column_pairing_option,Row_access_option>::operator+=(Z2_intrusive_list_column const &column)
{
	Column_type& tc = column_;
	const Column_type& sc = column.column_;

	auto it1 = tc.begin();
	auto it2 = sc.begin();
	while (it1 != tc.end() && it2 != sc.end())
	{
		if (it1->get_row_index() < it2->get_row_index()) {
			++it1;
		} else {
			if (it1->get_row_index() > it2->get_row_index()) {
				_insert_cell(it2->get_row_index(), it1);
			} else {
				_delete_cell(it1);
			}
			++it2;
		}
	}

	while (it2 != sc.end()) {
		_insert_cell(it2->get_row_index(), tc.end());
		++it2;
	}

	return *this;
}

template<class Cell_type, class Column_pairing_option, class Row_access_option>
inline Z2_intrusive_list_column<Cell_type,Column_pairing_option,Row_access_option> &Z2_intrusive_list_column<Cell_type,Column_pairing_option,Row_access_option>::operator*=(unsigned int v)
{
	if (v % 2 == 0){
		auto it = column_.begin();
		while (it != column_.end()){
			_delete_cell(it);
		}
	}

	return *this;
}

template<class Cell_type, class Column_pairing_option, class Row_access_option>
inline Z2_intrusive_list_column<Cell_type,Column_pairing_option,Row_access_option> &Z2_intrusive_list_column<Cell_type,Column_pairing_option,Row_access_option>::operator=(const Z2_intrusive_list_column& other)
{
	static_assert (!Row_access_option::isActive_, "= assignement not enabled with row access option.");

	Column_pairing_option::operator=(other);
	dim_ = other.dim_;
	column_.clone_from(other.column_, new_cloner(), delete_disposer());
	return *this;
}

template<class Cell_type, class Column_pairing_option, class Row_access_option>
inline void Z2_intrusive_list_column<Cell_type,Column_pairing_option,Row_access_option>::_delete_cell(iterator &it)
{
	iterator tmp_it = it;
	++it;
	Cell* tmp_ptr = &(*tmp_it);
	if constexpr (Row_access_option::isActive_)
		Row_access_option::unlink(tmp_ptr);
	column_.erase(tmp_it);
	delete tmp_ptr;
}

template<class Cell_type, class Column_pairing_option, class Row_access_option>
inline void Z2_intrusive_list_column<Cell_type,Column_pairing_option,Row_access_option>::_insert_cell(
		index rowIndex, const iterator &position)
{
	if constexpr (Row_access_option::isActive_){
		Cell *new_cell = new Cell(Row_access_option::columnIndex_, rowIndex);
		column_.insert(position, *new_cell);
		Row_access_option::insert_cell(rowIndex, new_cell);
	} else {
		Cell *new_cell = new Cell(rowIndex);
		column_.insert(position, *new_cell);
	}
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // Z2_INTRUSIVE_LIST_COLUMN_H

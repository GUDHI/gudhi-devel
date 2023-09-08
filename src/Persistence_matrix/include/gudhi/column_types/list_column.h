/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef LIST_COLUMN_H
#define LIST_COLUMN_H

#include <initializer_list>
#include <iostream>
#include <list>
#include <utility>
#include <vector>

#include "../utilities/utilities.h"
#include "cell.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Field_element_type, class Cell_type, class Row_access_option>
class List_column : public Row_access_option
{
public:
//	using Cell = Base_cell<Field_element_type>;
	using Cell = Cell_type;
	using Column_type = std::list<Cell>;
	using iterator = typename Column_type::iterator;
	using const_iterator = typename Column_type::const_iterator;
	using reverse_iterator = typename Column_type::reverse_iterator;
	using const_reverse_iterator = typename Column_type::const_reverse_iterator;

	List_column();
	template<class Container_type>
	List_column(const Container_type& nonZeroRowIndices);
	template<class Container_type>
	List_column(const Container_type& nonZeroRowIndices, dimension_type dimension);
	template<class Row_container_type>
	List_column(index columnIndex, Row_container_type &rowContainer);
	template<class Container_type, class Row_container_type>
	List_column(index columnIndex, const Container_type& nonZeroRowIndices, Row_container_type &rowContainer);
	template<class Container_type, class Row_container_type>
	List_column(index columnIndex, const Container_type& nonZeroRowIndices, dimension_type dimension, Row_container_type &rowContainer);
	List_column(const List_column& column);
	List_column(const List_column& column, index columnIndex);
	template<class Row_container_type>
	List_column(const List_column& column, index columnIndex, Row_container_type &rowContainer);
	List_column(List_column&& column) noexcept;
	~List_column();

	std::vector<Field_element_type> get_content(int columnLength = -1) const;
	bool is_non_zero(index rowIndex) const;
	bool is_empty() const;
	dimension_type get_dimension() const;
	template<class Map_type>
	void reorder(Map_type& valueMap);
	void clear();

	iterator begin() noexcept;
	const_iterator begin() const noexcept;
	iterator end() noexcept;
	const_iterator end() const noexcept;
	reverse_iterator rbegin() noexcept;
	const_reverse_iterator rbegin() const noexcept;
	reverse_iterator rend() noexcept;
	const_reverse_iterator rend() const noexcept;

	template<class Cell_range>
	List_column& operator+=(Cell_range const &column);
	friend List_column operator+(List_column column1, List_column const& column2){
		column1 += column2;
		return column1;
	}

	List_column& operator*=(unsigned int v);
	friend List_column operator*(List_column column, unsigned int const& v){
		column *= v;
		return column;
	}
	friend List_column operator*(unsigned int const& v, List_column column){
		column *= v;
		return column;
	}

	//this = v * this + column
	template<class Cell_range>
	List_column& multiply_and_add(const Field_element_type& v, const Cell_range& column);
	//this = this + column * v
	template<class Cell_range>
	List_column& multiply_and_add(const Cell_range& column, const Field_element_type& v);

	friend bool operator==(const List_column& c1, const List_column& c2){
		if (&c1 == &c2) return true;

		auto it1 = c1.column_.begin();
		auto it2 = c2.column_.begin();
		if (c1.column_.size() != c2.column_.size()) return false;
		while (it1 != c1.column_.end() && it2 != c2.column_.end()) {
			if (it1->get_row_index() != it2->get_row_index() || it1->get_element() != it2->get_element())
				return false;
			++it1; ++it2;
		}
		return true;
	}
	friend bool operator<(const List_column& c1, const List_column& c2){
		if (&c1 == &c2) return false;

		auto it1 = c1.column_.begin();
		auto it2 = c2.column_.begin();
		while (it1 != c1.column_.end() && it2 != c2.column_.end()) {
			if (it1->get_row_index() != it2->get_row_index())
				return it1->get_row_index() < it2->get_row_index();
			if (it1->get_element() != it2->get_element())
				return it1->get_element() < it2->get_element();
			++it1; ++it2;
		}
		return it2 != c2.column_.end();
	}

	List_column& operator=(List_column other);

	friend void swap(List_column& col1, List_column& col2){
		swap(static_cast<Row_access_option&>(col1),
			 static_cast<Row_access_option&>(col2));
		std::swap(col1.dim_, col2.dim_);
		col1.column_.swap(col2.column_);
	}

protected:
	dimension_type dim_;
	Column_type column_;

	void _delete_cell(iterator& it);
	void _insert_cell(const Field_element_type& value, index rowIndex, const iterator& position);
	void _update_cell(const Field_element_type& value, index rowIndex, const iterator& position);
};

template<class Field_element_type, class Cell_type, class Row_access_option>
inline List_column<Field_element_type,Cell_type,Row_access_option>::List_column() : dim_(0)
{
//	static_assert(!Row_access_option::isActive_, "When row access option enabled, a row container has to be provided.");
}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Container_type>
inline List_column<Field_element_type,Cell_type,Row_access_option>::List_column(const Container_type &nonZeroRowIndices)
	: dim_(nonZeroRowIndices.size() == 0 ? 0 : nonZeroRowIndices.size() - 1), column_(nonZeroRowIndices.size())
{
//	static_assert(!Row_access_option::isActive_, "When row access option enabled, a row container has to be provided.");

	auto it = column_.begin();
	for (const auto& p : nonZeroRowIndices){
		_update_cell(p.second, p.first, it++);
	}
}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Container_type>
inline List_column<Field_element_type,Cell_type,Row_access_option>::List_column(const Container_type &nonZeroRowIndices, dimension_type dimension)
	: dim_(dimension), column_(nonZeroRowIndices.size())
{
//	static_assert(!Row_access_option::isActive_, "When row access option enabled, a row container has to be provided.");

	auto it = column_.begin();
	for (const auto& p : nonZeroRowIndices){
		_update_cell(p.second, p.first, it++);
	}
}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Row_container_type>
inline List_column<Field_element_type,Cell_type,Row_access_option>::List_column(
		index columnIndex, Row_container_type &rowContainer)
	: Row_access_option(columnIndex, rowContainer), dim_(0)
{}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Container_type, class Row_container_type>
inline List_column<Field_element_type,Cell_type,Row_access_option>::List_column(
		index columnIndex, const Container_type &nonZeroRowIndices, Row_container_type &rowContainer)
	: Row_access_option(columnIndex, rowContainer), dim_(nonZeroRowIndices.size() == 0 ? 0 : nonZeroRowIndices.size() - 1), column_(nonZeroRowIndices.size())
{
	auto it = column_.begin();
	for (const auto& p : nonZeroRowIndices){
		_update_cell(p.second, p.first, it++);
	}
}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Container_type, class Row_container_type>
inline List_column<Field_element_type,Cell_type,Row_access_option>::List_column(
		index columnIndex, const Container_type &nonZeroRowIndices, dimension_type dimension, Row_container_type &rowContainer)
	: Row_access_option(columnIndex, rowContainer), dim_(dimension), column_(nonZeroRowIndices.size())
{
	auto it = column_.begin();
	for (const auto& p : nonZeroRowIndices){
		_update_cell(p.second, p.first, it++);
	}
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline List_column<Field_element_type,Cell_type,Row_access_option>::List_column(const List_column &column)
	: dim_(column.dim_),
	  column_(column.column_)
{
	static_assert(!Row_access_option::isActive_,
			"Copy constructor not available when row access option enabled.");
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline List_column<Field_element_type,Cell_type,Row_access_option>::List_column(
		const List_column &column, index columnIndex)
	: Row_access_option(columnIndex, *column.rows_),
	  dim_(column.dim_),
	  column_(column.column_.size())
{
	auto it = column_.begin();
	for (const Cell& cell : column.column_){
		_update_cell(cell.get_element(), cell.get_row_index(), it++);
	}
}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Row_container_type>
inline List_column<Field_element_type,Cell_type,Row_access_option>::List_column(
		const List_column &column, index columnIndex, Row_container_type &rowContainer)
	: Row_access_option(columnIndex, rowContainer),
	  dim_(column.dim_),
	  column_(column.column_.size())
{
	auto it = column_.begin();
	for (const Cell& cell : column.column_){
		_update_cell(cell.get_element(), cell.get_row_index(), it++);
	}
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline List_column<Field_element_type,Cell_type,Row_access_option>::List_column(List_column &&column) noexcept
	: Row_access_option(std::move(column)),
	  dim_(std::exchange(column.dim_, 0)),
	  column_(std::move(column.column_))
{}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline List_column<Field_element_type,Cell_type,Row_access_option>::~List_column()
{
	if constexpr (Row_access_option::isActive_){
		for (Cell& cell : column_)
			Row_access_option::unlink(&cell);
	}
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline std::vector<Field_element_type> List_column<Field_element_type,Cell_type,Row_access_option>::get_content(int columnLength) const
{
	if (columnLength < 0) columnLength = column_.back().get_row_index() + 1;

	std::vector<Field_element_type> container(columnLength);
	for (auto it = column_.begin(); it != column_.end() && it->get_row_index() < static_cast<index>(columnLength); ++it){
		container[it->get_row_index()] = it->get_element();
	}
	return container;
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline bool List_column<Field_element_type,Cell_type,Row_access_option>::is_non_zero(index rowIndex) const
{
	for (const Cell& v : column_){
		if (v.get_row_index() == rowIndex) return true;
	}
	return false;
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline bool List_column<Field_element_type,Cell_type,Row_access_option>::is_empty() const
{
	return column_.empty();
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline dimension_type List_column<Field_element_type,Cell_type,Row_access_option>::get_dimension() const
{
	return dim_;
}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Map_type>
inline void List_column<Field_element_type,Cell_type,Row_access_option>::reorder(Map_type &valueMap)
{
	iterator it = column_.begin();
	while (it != column_.end()) {
		Cell* cell = &(*it);
		if constexpr (Row_access_option::isActive_) Row_access_option::unlink(cell);
		cell->set_row_index(valueMap[cell->get_row_index()]);
		it++;
	}
	//all cells have to be deleted first, to avoid problem with insertion when row is a set
	if constexpr (Row_access_option::isActive_){
		for (auto it = column_.begin(); it != column_.end(); ++it) {
			Cell* cell = &(*it);
			Row_access_option::insert_cell(cell->get_row_index(), cell);
		}
	}
	column_.sort();
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline typename List_column<Field_element_type,Cell_type,Row_access_option>::iterator
List_column<Field_element_type,Cell_type,Row_access_option>::begin() noexcept
{
	return column_.begin();
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline typename List_column<Field_element_type,Cell_type,Row_access_option>::const_iterator
List_column<Field_element_type,Cell_type,Row_access_option>::begin() const noexcept
{
	return column_.begin();
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline typename List_column<Field_element_type,Cell_type,Row_access_option>::iterator
List_column<Field_element_type,Cell_type,Row_access_option>::end() noexcept
{
	return column_.end();
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline typename List_column<Field_element_type,Cell_type,Row_access_option>::const_iterator
List_column<Field_element_type,Cell_type,Row_access_option>::end() const noexcept
{
	return column_.end();
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline typename List_column<Field_element_type,Cell_type,Row_access_option>::reverse_iterator
List_column<Field_element_type,Cell_type,Row_access_option>::rbegin() noexcept
{
	return column_.rbegin();
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline typename List_column<Field_element_type,Cell_type,Row_access_option>::const_reverse_iterator
List_column<Field_element_type,Cell_type,Row_access_option>::rbegin() const noexcept
{
	return column_.rbegin();
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline typename List_column<Field_element_type,Cell_type,Row_access_option>::reverse_iterator
List_column<Field_element_type,Cell_type,Row_access_option>::rend() noexcept
{
	return column_.rend();
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline typename List_column<Field_element_type,Cell_type,Row_access_option>::const_reverse_iterator
List_column<Field_element_type,Cell_type,Row_access_option>::rend() const noexcept
{
	return column_.rend();
}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Cell_range>
inline List_column<Field_element_type,Cell_type,Row_access_option> &List_column<Field_element_type,Cell_type,Row_access_option>::operator+=(Cell_range const &column)
{
	if (column.begin() == column.end()) return *this;
	if (column_.empty()){
		if constexpr (Row_access_option::isActive_){
//			column_.resize(column.column_.size());
//			auto it = column_.begin();
			for (const Cell& cell : column)
				_insert_cell(cell.get_element(), cell.get_row_index(), column_.end());
//				_update_cell(cell.get_element(), cell.get_row_index(), it++);
		} else {
			std::copy(column.begin(), column.end(), std::back_inserter(column_));
		}
		return *this;
	}

	const_iterator itToAdd = column.begin();
	iterator itTarget = column_.begin();
	index curRowToAdd = itToAdd->get_row_index();
	index curRowTarget = itTarget->get_row_index();

	while (itToAdd != column.end() && itTarget != column_.end())
	{
		if (curRowToAdd == curRowTarget){
			itTarget->get_element() += itToAdd->get_element();
			if (itTarget->get_element() == Field_element_type::get_additive_identity()) _delete_cell(itTarget);
			else {
				if constexpr (Row_access_option::isActive_)
						Row_access_option::update_cell(*itTarget);
				itTarget++;
			}
			itToAdd++;
		} else if (curRowToAdd < curRowTarget){
			_insert_cell(itToAdd->get_element(), curRowToAdd, itTarget);
			itToAdd++;
		} else {
			itTarget++;
		}

		if (itToAdd != column.end()) curRowToAdd = itToAdd->get_row_index();
		if (itTarget != column_.end()) curRowTarget = itTarget->get_row_index();
	}

	while (itToAdd != column.end()){
		curRowToAdd = itToAdd->get_row_index();
		if (itToAdd != column.end()){
			_insert_cell(itToAdd->get_element(), curRowToAdd, column_.end());
			itToAdd++;
		}
	}

	return *this;
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline List_column<Field_element_type,Cell_type,Row_access_option> &List_column<Field_element_type,Cell_type,Row_access_option>::operator*=(unsigned int v)
{
//	v %= Field_element_type::get_characteristic();
	Field_element_type val(v);

	if (val == 0u) {
		clear();
		return *this;
	}

	if (val == 1u) return *this;

	for (Cell& cell : column_){
		cell.get_element() *= val;
		if constexpr (Row_access_option::isActive_)
				Row_access_option::update_cell(cell);
	}

	return *this;
}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Cell_range>
inline List_column<Field_element_type,Cell_type,Row_access_option> &
List_column<Field_element_type,Cell_type,Row_access_option>::multiply_and_add(const Field_element_type& val, const Cell_range& column)
{
	if (val == 0u) {
		clear();
	}

	const_iterator itToAdd = column.begin();
	iterator itTarget = column_.begin();
	index curRowToAdd = itToAdd->get_row_index();
	index curRowTarget = itTarget->get_row_index();

	while (itToAdd != column.end() && itTarget != column_.end())
	{
		if (curRowToAdd == curRowTarget){
			itTarget->get_element() *= val;
			itTarget->get_element() += itToAdd->get_element();
			if (itTarget->get_element() == Field_element_type::get_additive_identity()) _delete_cell(itTarget);
			else {
				if constexpr (Row_access_option::isActive_)
						Row_access_option::update_cell(*itTarget);
				itTarget++;
			}
			itToAdd++;
		} else if (curRowToAdd < curRowTarget){
			_insert_cell(itToAdd->get_element(), curRowToAdd, itTarget);
			itToAdd++;
		} else {
			itTarget->get_element() *= val;
			if constexpr (Row_access_option::isActive_)
					Row_access_option::update_cell(*itTarget);
			itTarget++;
		}

		if (itToAdd != column.end()) curRowToAdd = itToAdd->get_row_index();
		if (itTarget != column_.end()) curRowTarget = itTarget->get_row_index();
	}

	while (itTarget != column_.end()){
		itTarget->get_element() *= val;
		if constexpr (Row_access_option::isActive_)
				Row_access_option::update_cell(*itTarget);
		++itTarget;
	}

	while (itToAdd != column.end()){
		curRowToAdd = itToAdd->get_row_index();
		if (itToAdd != column.end()){
			_insert_cell(itToAdd->get_element(), curRowToAdd, column_.end());
			itToAdd++;
		}
	}

	return *this;
}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Cell_range>
inline List_column<Field_element_type,Cell_type,Row_access_option> &
List_column<Field_element_type,Cell_type,Row_access_option>::multiply_and_add(const Cell_range& column, const Field_element_type& val)
{
	if (val == 0u) {
		return *this;
	}

	const_iterator itToAdd = column.begin();
	iterator itTarget = column_.begin();
	index curRowToAdd = itToAdd->get_row_index();
	index curRowTarget = itTarget->get_row_index();

	while (itToAdd != column.end() && itTarget != column_.end())
	{
		if (curRowToAdd == curRowTarget){
			itTarget->get_element() += (itToAdd->get_element() * val);
			if (itTarget->get_element() == Field_element_type::get_additive_identity())
				_delete_cell(itTarget);
			else {
				if constexpr (Row_access_option::isActive_)
						Row_access_option::update_cell(*itTarget);
				itTarget++;
			}
			itToAdd++;
		} else if (curRowToAdd < curRowTarget){
			_insert_cell(itToAdd->get_element() * val, curRowToAdd, itTarget);
			itToAdd++;
		} else {
			itTarget++;
		}

		if (itToAdd != column.end()) curRowToAdd = itToAdd->get_row_index();
		if (itTarget != column_.end()) curRowTarget = itTarget->get_row_index();
	}

	while (itToAdd != column.end()){
		curRowToAdd = itToAdd->get_row_index();
		if (itToAdd != column.end()){
			_insert_cell(itToAdd->get_element() * val, curRowToAdd, column_.end());
			itToAdd++;
		}
	}

	return *this;
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline List_column<Field_element_type,Cell_type,Row_access_option> &List_column<Field_element_type,Cell_type,Row_access_option>::operator=(List_column other)
{
	static_assert (!Row_access_option::isActive_, "= assignement not enabled with row access option.");

	std::swap(dim_, other.dim_);
	column_.swap(other.column_);
	return *this;
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline void List_column<Field_element_type,Cell_type,Row_access_option>::_delete_cell(iterator &it)
{
	if constexpr (Row_access_option::isActive_)
		Row_access_option::unlink(&(*it));
	column_.erase(it++);
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline void List_column<Field_element_type,Cell_type,Row_access_option>::_insert_cell(
		const Field_element_type &value, index rowIndex, const iterator &position)
{
	if constexpr (Row_access_option::isActive_){
		auto it = column_.insert(position, Cell(value, Row_access_option::columnIndex_, rowIndex));
		Row_access_option::insert_cell(rowIndex, &(*it));
	} else {
		column_.emplace(position, value, rowIndex);
	}
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline void List_column<Field_element_type,Cell_type,Row_access_option>::_update_cell(
		const Field_element_type &value, index rowIndex, const iterator &position)
{
	if constexpr (Row_access_option::isActive_){
		*position = Cell(value, Row_access_option::columnIndex_, rowIndex);
		Row_access_option::insert_cell(rowIndex, &(*position));
	} else {
		*position = Cell(value, rowIndex);
	}
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline void List_column<Field_element_type,Cell_type,Row_access_option>::clear()
{
	if constexpr (Row_access_option::isActive_){
		for (Cell& cell : column_)
			Row_access_option::unlink(&cell);
	}
	column_.clear();
}

} //namespace persistence_matrix
} //namespace Gudhi

template<class Field_element_type, class Cell_type, class Row_access_option>
struct std::hash<Gudhi::persistence_matrix::List_column<Field_element_type,Cell_type,Row_access_option> >
{
	size_t operator()(const Gudhi::persistence_matrix::List_column<Field_element_type,Cell_type,Row_access_option>& column) const
	{
		std::size_t seed = 0;
		for (auto& cell : column){
			seed ^= std::hash<unsigned int>()(cell.get_row_index() * static_cast<unsigned int>(cell.get_element())) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		}
		return seed;
	}
};

#endif // LIST_COLUMN_H

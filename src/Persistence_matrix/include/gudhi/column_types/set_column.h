/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef SET_COLUMN_H
#define SET_COLUMN_H

#include <iostream>
#include <list>
#include <set>

#include "../utilities/utilities.h"
#include "../utilities/Zp_field.h"
#include "cell.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Field_element_type, class Cell_type, class Row_access_option>
class Set_column : public Row_access_option
{
public:
//	using Cell = Base_cell<Field_element_type>;
	using Cell = Cell_type;
	using Column_type = std::set<Cell>;
	using iterator = typename Column_type::iterator;
	using const_iterator = typename Column_type::const_iterator;
	using reverse_iterator = typename Column_type::reverse_iterator;
	using const_reverse_iterator = typename Column_type::const_reverse_iterator;

	Set_column();
	template<class Container_type>
	Set_column(const Container_type& nonZeroRowIndices);
	template<class Container_type>
	Set_column(const Container_type& nonZeroRowIndices, dimension_type dimension);
	template<class Row_container_type>
	Set_column(index columnIndex, Row_container_type &rowContainer);
	template<class Container_type, class Row_container_type>
	Set_column(index columnIndex, const Container_type& nonZeroRowIndices, Row_container_type &rowContainer);
	template<class Container_type, class Row_container_type>
	Set_column(index columnIndex, const Container_type& nonZeroRowIndices, dimension_type dimension, Row_container_type &rowContainer);
	Set_column(const Set_column& column);
	Set_column(const Set_column& column, index columnIndex);
	template<class Row_container_type>
	Set_column(const Set_column& column, index columnIndex, Row_container_type &rowContainer);
	Set_column(Set_column&& column) noexcept;
	~Set_column();

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
	Set_column& operator+=(Cell_range const &column);
	friend Set_column operator+(Set_column column1, Set_column const& column2){
		column1 += column2;
		return column1;
	}
	Set_column& operator*=(unsigned int v);
	friend Set_column operator*(Set_column column, unsigned int const& v){
		column *= v;
		return column;
	}
	friend Set_column operator*(unsigned int const& v, Set_column column){
		column *= v;
		return column;
	}

	//this = v * this + column
	template<class Cell_range>
	Set_column& multiply_and_add(const Field_element_type& v, const Cell_range& column);
	//this = this + column * v
	template<class Cell_range>
	Set_column& multiply_and_add(const Cell_range& column, const Field_element_type& v);

	friend bool operator==(const Set_column& c1, const Set_column& c2){
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
	friend bool operator<(const Set_column& c1, const Set_column& c2){
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

	Set_column& operator=(Set_column other);

	friend void swap(Set_column& col1, Set_column& col2){
		swap(static_cast<Row_access_option&>(col1),
			 static_cast<Row_access_option&>(col2));
		std::swap(col1.dim_, col2.dim_);
		col1.column_.swap(col2.column_);
	}

protected:
	int dim_;
	Column_type column_;

	void _delete_cell(iterator& it);
	void _insert_cell(const Field_element_type& value, index rowIndex, const iterator &position);
};

template<class Field_element_type, class Cell_type, class Row_access_option>
inline Set_column<Field_element_type,Cell_type,Row_access_option>::Set_column() : dim_(0)
{
//	static_assert(!Row_access_option::isActive_, "When row access option enabled, a row container has to be provided.");
}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Container_type>
inline Set_column<Field_element_type,Cell_type,Row_access_option>::Set_column(const Container_type &nonZeroRowIndices)
	: dim_(nonZeroRowIndices.size() == 0 ? 0 : nonZeroRowIndices.size() - 1)
{
//	static_assert(!Row_access_option::isActive_, "When row access option enabled, a row container has to be provided.");

	for (const auto& p : nonZeroRowIndices){
		_insert_cell(p.second, p.first, column_.end());
	}
}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Container_type>
inline Set_column<Field_element_type,Cell_type,Row_access_option>::Set_column(const Container_type &nonZeroRowIndices, dimension_type dimension)
	: dim_(dimension)
{
//	static_assert(!Row_access_option::isActive_, "When row access option enabled, a row container has to be provided.");

	for (const auto& p : nonZeroRowIndices){
		_insert_cell(p.second, p.first, column_.end());
	}
}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Row_container_type>
inline Set_column<Field_element_type,Cell_type,Row_access_option>::Set_column(
		index columnIndex, Row_container_type &rowContainer)
	: Row_access_option(columnIndex, rowContainer), dim_(0)
{}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Container_type, class Row_container_type>
inline Set_column<Field_element_type,Cell_type,Row_access_option>::Set_column(
		index columnIndex, const Container_type &nonZeroRowIndices, Row_container_type &rowContainer)
	: Row_access_option(columnIndex, rowContainer), dim_(nonZeroRowIndices.size() == 0 ? 0 : nonZeroRowIndices.size() - 1)
{
	for (const auto& p : nonZeroRowIndices){
		_insert_cell(p.second, p.first, column_.end());
	}
}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Container_type, class Row_container_type>
inline Set_column<Field_element_type,Cell_type,Row_access_option>::Set_column(
		index columnIndex, const Container_type &nonZeroRowIndices, dimension_type dimension, Row_container_type &rowContainer)
	: Row_access_option(columnIndex, rowContainer), dim_(dimension)
{
	for (const auto& p : nonZeroRowIndices){
		_insert_cell(p.second, p.first, column_.end());
	}
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline Set_column<Field_element_type,Cell_type,Row_access_option>::Set_column(const Set_column &column)
	: dim_(column.dim_),
	  column_(column.column_)
{
	static_assert(!Row_access_option::isActive_,
			"Copy constructor not available when row access option enabled.");
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline Set_column<Field_element_type,Cell_type,Row_access_option>::Set_column(
		const Set_column &column, index columnIndex)
	: Row_access_option(columnIndex, *column.rows_),
	  dim_(column.dim_)
{
	for (const Cell& cell : column.column_){
		_insert_cell(cell.get_element(), cell.get_row_index(), column_.end());
	}
}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Row_container_type>
inline Set_column<Field_element_type,Cell_type,Row_access_option>::Set_column(
		const Set_column &column, index columnIndex, Row_container_type &rowContainer)
	: Row_access_option(columnIndex, rowContainer),
	  dim_(column.dim_)
{
	for (const Cell& cell : column.column_){
		_insert_cell(cell.get_element(), cell.get_row_index(), column_.end());
	}
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline Set_column<Field_element_type,Cell_type,Row_access_option>::Set_column(Set_column &&column) noexcept
	: Row_access_option(std::move(column)),
	  dim_(std::exchange(column.dim_, 0)),
	  column_(std::move(column.column_))
{}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline Set_column<Field_element_type,Cell_type,Row_access_option>::~Set_column()
{
	if constexpr (Row_access_option::isActive_){
		for (const Cell& cell : column_)
			Row_access_option::unlink(cell);
	}
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline std::vector<Field_element_type> Set_column<Field_element_type,Cell_type,Row_access_option>::get_content(int columnLength) const
{
	if (columnLength < 0) columnLength = column_.rbegin()->get_row_index() + 1;

	std::vector<Field_element_type> container(columnLength);
	for (auto it = column_.begin(); it != column_.end() && it->get_row_index() < static_cast<index>(columnLength); ++it){
		container[it->get_row_index()] = it->get_element();
	}
	return container;
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline bool Set_column<Field_element_type,Cell_type,Row_access_option>::is_non_zero(index rowIndex) const
{
	if constexpr (Row_access_option::isActive_){
		return column_.find(Cell(0, Row_access_option::columnIndex_, rowIndex)) != column_.end();
	} else {
		return column_.find(Cell(0, rowIndex)) != column_.end();
	}
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline bool Set_column<Field_element_type,Cell_type,Row_access_option>::is_empty() const
{
	return column_.empty();
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline dimension_type Set_column<Field_element_type,Cell_type,Row_access_option>::get_dimension() const
{
	return dim_;
}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Map_type>
inline void Set_column<Field_element_type,Cell_type,Row_access_option>::reorder(Map_type &valueMap)
{
	Column_type newSet;
	for (auto it = column_.begin(); it != column_.end(); ) {
		if constexpr (Row_access_option::isActive_) {
			newSet.emplace_hint(newSet.end(), it->get_element(), Row_access_option::columnIndex_, valueMap[it->get_row_index()]);
			auto ittemp = it;
			++it;
			_delete_cell(ittemp);
		} else {
			newSet.emplace_hint(newSet.end(), it->get_element(), valueMap[it->get_row_index()]);
			++it;
		}
	}
	//all cells have to be deleted first, to avoid problem with insertion when row is a set
	if constexpr (Row_access_option::isActive_) {
		for (auto it = newSet.begin(); it != newSet.end(); ++it) {
			Row_access_option::insert_cell(it->get_row_index(), *it);
		}
	}
	column_.swap(newSet);
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline typename Set_column<Field_element_type,Cell_type,Row_access_option>::iterator
Set_column<Field_element_type,Cell_type,Row_access_option>::begin() noexcept
{
	return column_.begin();
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline typename Set_column<Field_element_type,Cell_type,Row_access_option>::const_iterator
Set_column<Field_element_type,Cell_type,Row_access_option>::begin() const noexcept
{
	return column_.begin();
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline typename Set_column<Field_element_type,Cell_type,Row_access_option>::iterator
Set_column<Field_element_type,Cell_type,Row_access_option>::end() noexcept
{
	return column_.end();
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline typename Set_column<Field_element_type,Cell_type,Row_access_option>::const_iterator
Set_column<Field_element_type,Cell_type,Row_access_option>::end() const noexcept
{
	return column_.end();
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline typename Set_column<Field_element_type,Cell_type,Row_access_option>::reverse_iterator
Set_column<Field_element_type,Cell_type,Row_access_option>::rbegin() noexcept
{
	return column_.rbegin();
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline typename Set_column<Field_element_type,Cell_type,Row_access_option>::const_reverse_iterator
Set_column<Field_element_type,Cell_type,Row_access_option>::rbegin() const noexcept
{
	return column_.rbegin();
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline typename Set_column<Field_element_type,Cell_type,Row_access_option>::reverse_iterator
Set_column<Field_element_type,Cell_type,Row_access_option>::rend() noexcept
{
	return column_.rend();
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline typename Set_column<Field_element_type,Cell_type,Row_access_option>::const_reverse_iterator
Set_column<Field_element_type,Cell_type,Row_access_option>::rend() const noexcept
{
	return column_.rend();
}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Cell_range>
inline Set_column<Field_element_type,Cell_type,Row_access_option> &Set_column<Field_element_type,Cell_type,Row_access_option>::operator+=(Cell_range const &column)
{
	for (const Cell& v : column){
		auto c = column_.find(v);
		if (c != column_.end()){
			index r = c->get_row_index();
			Field_element_type coef = c->get_element();
			coef += v.get_element();
			_delete_cell(c);
			if (coef != Field_element_type::get_additive_identity())
				_insert_cell(coef, r, column_.end());
		} else
			_insert_cell(v.get_element(), v.get_row_index(), column_.end());
	}

	return *this;
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline Set_column<Field_element_type,Cell_type,Row_access_option> &Set_column<Field_element_type,Cell_type,Row_access_option>::operator*=(unsigned int v)
{
//	v %= Field_element_type::get_characteristic();
	Field_element_type val(v);

	if (val == Field_element_type::get_multiplicative_identity()) return *this;

	if (val == 0u) {
		clear();
		return *this;
	}

	Column_type newColumn;

	for (const Cell& cell : column_){
		if constexpr (Row_access_option::isActive_){
			Row_access_option::unlink(cell);
		}

		Cell newCell(cell);
		newCell.get_element() *= val;

		if constexpr (Row_access_option::isActive_){
			auto it = newColumn.insert(newCell);
			Row_access_option::insert_cell(newCell.get_row_index(), *it.first);
		} else {
			newColumn.insert(newCell);
		}
	}

	column_.swap(newColumn);

	return *this;
}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Cell_range>
inline Set_column<Field_element_type,Cell_type,Row_access_option> &
Set_column<Field_element_type,Cell_type,Row_access_option>::multiply_and_add(const Field_element_type& val, const Cell_range& column)
{
	if (val == 0u) {
		clear();
	}

	Column_type newColumn;
	auto emplace_new_cell = [&](Field_element_type value, index row_index){
		if constexpr (Row_access_option::isActive_){
			auto it = newColumn.emplace_hint(newColumn.end(), value, Row_access_option::columnIndex_, row_index);
			Row_access_option::insert_cell(row_index, *it);
		} else {
			newColumn.emplace_hint(newColumn.end(), value, row_index);
		}
	};

	auto itTarget = column_.begin();
	auto itSource = column.begin();
	while (itTarget != column_.end() && itSource != column.end())
	{
		if (itTarget->get_row_index() < itSource->get_row_index()) {
			if constexpr (Row_access_option::isActive_){
				Row_access_option::unlink(*itTarget);
			}
			emplace_new_cell(itTarget->get_element() * val, itTarget->get_row_index());
			++itTarget;
		} else if (itTarget->get_row_index() > itSource->get_row_index()) {
			emplace_new_cell(itSource->get_element(), itSource->get_row_index());
			++itSource;
		} else {
			if constexpr (Row_access_option::isActive_){
				Row_access_option::unlink(*itTarget);
			}
			Field_element_type f = itTarget->get_element() * val + itSource->get_element();
			if (f != Field_element_type::get_additive_identity()){
				emplace_new_cell(f, itTarget->get_row_index());
			}
			++itTarget;
			++itSource;
		}
	}

	while (itTarget != column_.end()){
		if constexpr (Row_access_option::isActive_){
			Row_access_option::unlink(*itTarget);
		}
		emplace_new_cell(itTarget->get_element() * val, itTarget->get_row_index());
		++itTarget;
	}

	while (itSource != column.end()) {
		emplace_new_cell(itSource->get_element(), itSource->get_row_index());
		++itSource;
	}

	column_.swap(newColumn);

	return *this;
}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Cell_range>
inline Set_column<Field_element_type,Cell_type,Row_access_option> &
Set_column<Field_element_type,Cell_type,Row_access_option>::multiply_and_add(const Cell_range& column, const Field_element_type& val)
{
	if (val == 0u) {
		return *this;
	}

	for (const Cell& v : column){
		auto c = column_.find(v);
		if (c != column_.end()){
			index r = c->get_row_index();
			Field_element_type coef(c->get_element());
			coef += (v.get_element() * val);
			_delete_cell(c);
			if (coef != Field_element_type::get_additive_identity())
				_insert_cell(coef, r, column_.end());
		} else
			_insert_cell(v.get_element() * val, v.get_row_index(), column_.end());
	}

	return *this;
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline Set_column<Field_element_type,Cell_type,Row_access_option> &Set_column<Field_element_type,Cell_type,Row_access_option>::operator=(Set_column other)
{
	static_assert (!Row_access_option::isActive_, "= assignement not enabled with row access option.");

	std::swap(dim_, other.dim_);
	column_.swap(other.column_);
	return *this;
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline void Set_column<Field_element_type,Cell_type,Row_access_option>::_delete_cell(iterator &it)
{
	if constexpr (Row_access_option::isActive_)
		Row_access_option::unlink(*it);
	column_.erase(it);
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline void Set_column<Field_element_type,Cell_type,Row_access_option>::_insert_cell(
		const Field_element_type &value, index rowIndex, const iterator &position)
{
	if constexpr (Row_access_option::isActive_){
		auto it = column_.emplace_hint(position, value, Row_access_option::columnIndex_, rowIndex);
		Row_access_option::insert_cell(rowIndex, *it);
	} else {
		column_.emplace_hint(position, value, rowIndex);
	}
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline void Set_column<Field_element_type,Cell_type,Row_access_option>::clear()
{
	if constexpr (Row_access_option::isActive_){
		for (const Cell& cell : column_)
			Row_access_option::unlink(cell);
	}
	column_.clear();
}

} //namespace persistence_matrix
} //namespace Gudhi

template<class Field_element_type, class Cell_type, class Row_access_option>
struct std::hash<Gudhi::persistence_matrix::Set_column<Field_element_type,Cell_type,Row_access_option> >
{
	size_t operator()(const Gudhi::persistence_matrix::Set_column<Field_element_type,Cell_type,Row_access_option>& column) const
	{
		std::size_t seed = 0;
		for (auto& cell : column){
			seed ^= std::hash<unsigned int>()(cell.get_row_index() * static_cast<unsigned int>(cell.get_element())) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		}
		return seed;
	}
};

#endif // SET_COLUMN_H

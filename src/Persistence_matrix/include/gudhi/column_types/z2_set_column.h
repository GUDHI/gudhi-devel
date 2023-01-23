/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef Z2_SET_COLUMN_H
#define Z2_SET_COLUMN_H

#include <iostream>
#include <list>
#include <set>

#include "../utilities/utilities.h"
#include "cell.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Cell_type, class Row_access_option>
class Z2_set_column : public Row_access_option
{
public:
//	using Cell = Z2_base_cell;
	using Cell = Cell_type;
	using Column_type = std::set<Cell>;
	using iterator = typename Column_type::iterator;
	using const_iterator = typename Column_type::const_iterator;
	using reverse_iterator = typename Column_type::reverse_iterator;
	using const_reverse_iterator = typename Column_type::const_reverse_iterator;

	Z2_set_column();
	template<class Container_type>
	Z2_set_column(const Container_type& nonZeroRowIndices);
	template<class Container_type>
	Z2_set_column(const Container_type& nonZeroRowIndices, dimension_type dimension);
	template<class Row_container_type>
	Z2_set_column(index columnIndex, Row_container_type &rowContainer);
	template<class Container_type, class Row_container_type>
	Z2_set_column(index columnIndex, const Container_type& nonZeroRowIndices, Row_container_type &rowContainer);
	template<class Container_type, class Row_container_type>
	Z2_set_column(index columnIndex, const Container_type& nonZeroRowIndices, dimension_type dimension, Row_container_type &rowContainer);
	Z2_set_column(const Z2_set_column& column);
	Z2_set_column(const Z2_set_column& column, index columnIndex);
	Z2_set_column(Z2_set_column&& column) noexcept;
	~Z2_set_column();

	std::vector<bool> get_content(int columnLength = -1) const;
	bool is_non_zero(index rowIndex) const;
	bool is_empty() const;
	dimension_type get_dimension() const;
	template<class Map_type>
	void reorder(Map_type& valueMap);

	iterator begin() noexcept;
	const_iterator begin() const noexcept;
	iterator end() noexcept;
	const_iterator end() const noexcept;
	reverse_iterator rbegin() noexcept;
	const_reverse_iterator rbegin() const noexcept;
	reverse_iterator rend() noexcept;
	const_reverse_iterator rend() const noexcept;

	Z2_set_column& operator+=(Z2_set_column const &column);
	friend Z2_set_column operator+(Z2_set_column column1, Z2_set_column const& column2){
		column1 += column2;
		return column1;
	}

	Z2_set_column& operator*=(unsigned int v);
	friend Z2_set_column operator*(Z2_set_column column, unsigned int const& v){
		column *= v;
		return column;
	}
	friend Z2_set_column operator*(unsigned int const& v, Z2_set_column column){
		column *= v;
		return column;
	}

	friend bool operator==(const Z2_set_column& c1, const Z2_set_column& c2){
		return c1.column_ == c2.column_;
	}
	friend bool operator<(const Z2_set_column& c1, const Z2_set_column& c2){
		return c1.column_ < c2.column_;
	}

	Z2_set_column& operator=(Z2_set_column other);

	friend void swap(Z2_set_column& col1, Z2_set_column& col2){
		swap(static_cast<Row_access_option&>(col1),
			 static_cast<Row_access_option&>(col2));
		std::swap(col1.dim_, col2.dim_);
		col1.column_.swap(col2.column_);
	}

protected:
	int dim_;
	Column_type column_;

	void _delete_cell(iterator& it);
	void _insert_cell(index rowIndex);
};

template<class Cell_type, class Row_access_option>
inline Z2_set_column<Cell_type,Row_access_option>::Z2_set_column() : dim_(0)
{
	static_assert(!Row_access_option::isActive_, "When row access option enabled, a row container has to be provided.");
}

template<class Cell_type, class Row_access_option>
template<class Container_type>
inline Z2_set_column<Cell_type,Row_access_option>::Z2_set_column(const Container_type &nonZeroRowIndices)
	: dim_(nonZeroRowIndices.size() == 0 ? 0 : nonZeroRowIndices.size() - 1),
	  column_(nonZeroRowIndices.begin(), nonZeroRowIndices.end())
{
	static_assert(!Row_access_option::isActive_, "When row access option enabled, a row container has to be provided.");
}

template<class Cell_type, class Row_access_option>
template<class Container_type>
inline Z2_set_column<Cell_type,Row_access_option>::Z2_set_column(const Container_type &nonZeroRowIndices, dimension_type dimension)
	: dim_(dimension),
	  column_(nonZeroRowIndices.begin(), nonZeroRowIndices.end())
{
	static_assert(!Row_access_option::isActive_, "When row access option enabled, a row container has to be provided.");
}

template<class Cell_type, class Row_access_option>
template<class Row_container_type>
inline Z2_set_column<Cell_type,Row_access_option>::Z2_set_column(
		index columnIndex, Row_container_type &rowContainer)
	: Row_access_option(columnIndex, rowContainer), dim_(0)
{}

template<class Cell_type, class Row_access_option>
template<class Container_type, class Row_container_type>
inline Z2_set_column<Cell_type,Row_access_option>::Z2_set_column(
		index columnIndex, const Container_type &nonZeroRowIndices, Row_container_type &rowContainer)
	: Row_access_option(columnIndex, rowContainer), dim_(nonZeroRowIndices.size() == 0 ? 0 : nonZeroRowIndices.size() - 1)
{
	for (index id : nonZeroRowIndices){
		_insert_cell(id);
	}
}

template<class Cell_type, class Row_access_option>
template<class Container_type, class Row_container_type>
inline Z2_set_column<Cell_type,Row_access_option>::Z2_set_column(
		index columnIndex, const Container_type &nonZeroRowIndices, dimension_type dimension, Row_container_type &rowContainer)
	: Row_access_option(columnIndex, rowContainer), dim_(dimension)
{
	for (index id : nonZeroRowIndices){
		_insert_cell(id);
	}
}

template<class Cell_type, class Row_access_option>
inline Z2_set_column<Cell_type,Row_access_option>::Z2_set_column(const Z2_set_column &column)
	: dim_(column.dim_),
	  column_(column.column_)
{
	static_assert(!Row_access_option::isActive_,
			"Copy constructor not available when row access option enabled.");
}

template<class Cell_type, class Row_access_option>
inline Z2_set_column<Cell_type,Row_access_option>::Z2_set_column(
		const Z2_set_column &column, index columnIndex)
	: Row_access_option(columnIndex, *column.rows_),
	  dim_(column.dim_)
{
	for (const Cell& cell : column.column_){
		_insert_cell(cell.get_row_index());
	}
}

template<class Cell_type, class Row_access_option>
inline Z2_set_column<Cell_type,Row_access_option>::Z2_set_column(Z2_set_column &&column) noexcept
	: Row_access_option(std::move(column)),
	  dim_(std::exchange(column.dim_, 0)),
	  column_(std::move(column.column_))
{}

template<class Cell_type, class Row_access_option>
inline Z2_set_column<Cell_type,Row_access_option>::~Z2_set_column()
{
	if constexpr (Row_access_option::isActive_){
		for (const Cell& cell : column_)
			Row_access_option::unlink(cell);
	}
}

template<class Cell_type, class Row_access_option>
inline std::vector<bool> Z2_set_column<Cell_type,Row_access_option>::get_content(int columnLength) const
{
	if (columnLength < 0) columnLength = column_.rbegin()->get_row_index() + 1;

	std::vector<bool> container(columnLength, 0);
	for (auto it = column_.begin(); it != column_.end() && it->get_row_index() < columnLength; ++it){
		container[it->get_row_index()] = 1;
	}
	return container;
}

template<class Cell_type, class Row_access_option>
inline bool Z2_set_column<Cell_type,Row_access_option>::is_non_zero(index rowIndex) const
{
	if constexpr (Row_access_option::isActive_){
		return column_.find(Cell(Row_access_option::columnIndex_, rowIndex)) != column_.end();
	} else {
		return column_.find(rowIndex) != column_.end();
	}
}

template<class Cell_type, class Row_access_option>
inline bool Z2_set_column<Cell_type,Row_access_option>::is_empty() const
{
	return column_.empty();
}

template<class Cell_type, class Row_access_option>
inline dimension_type Z2_set_column<Cell_type,Row_access_option>::get_dimension() const
{
	return dim_;
}

template<class Cell_type, class Row_access_option>
template<class Map_type>
inline void Z2_set_column<Cell_type,Row_access_option>::reorder(Map_type &valueMap)
{
	Column_type newSet;
	for (auto it = column_.begin(); it != column_.end(); ) {
		if constexpr (Row_access_option::isActive_) {
			newSet.emplace_hint(newSet.end(), Row_access_option::columnIndex_, valueMap[it->get_row_index()]);
			auto ittemp = it;
			++it;
			_delete_cell(ittemp);
		} else {
			newSet.emplace_hint(newSet.end(), valueMap[it->get_row_index()]);
			++it;
		}
	}
	//all cells have to be deleted first, to avoid problem with insertion when row is a set
	if constexpr (Row_access_option::isActive_) {
		for (auto cell : newSet) {
			Row_access_option::insert_cell(cell.get_row_index(), cell);
		}
	}
	column_.swap(newSet);
}

template<class Cell_type, class Row_access_option>
inline typename Z2_set_column<Cell_type,Row_access_option>::iterator
Z2_set_column<Cell_type,Row_access_option>::begin() noexcept
{
	return column_.begin();
}

template<class Cell_type, class Row_access_option>
inline typename Z2_set_column<Cell_type,Row_access_option>::const_iterator
Z2_set_column<Cell_type,Row_access_option>::begin() const noexcept
{
	return column_.begin();
}

template<class Cell_type, class Row_access_option>
inline typename Z2_set_column<Cell_type,Row_access_option>::iterator
Z2_set_column<Cell_type,Row_access_option>::end() noexcept
{
	return column_.end();
}

template<class Cell_type, class Row_access_option>
inline typename Z2_set_column<Cell_type,Row_access_option>::const_iterator
Z2_set_column<Cell_type,Row_access_option>::end() const noexcept
{
	return column_.end();
}

template<class Cell_type, class Row_access_option>
inline typename Z2_set_column<Cell_type,Row_access_option>::reverse_iterator
Z2_set_column<Cell_type,Row_access_option>::rbegin() noexcept
{
	return column_.rbegin();
}

template<class Cell_type, class Row_access_option>
inline typename Z2_set_column<Cell_type,Row_access_option>::const_reverse_iterator
Z2_set_column<Cell_type,Row_access_option>::rbegin() const noexcept
{
	return column_.rbegin();
}

template<class Cell_type, class Row_access_option>
inline typename Z2_set_column<Cell_type,Row_access_option>::reverse_iterator
Z2_set_column<Cell_type,Row_access_option>::rend() noexcept
{
	return column_.rend();
}

template<class Cell_type, class Row_access_option>
inline typename Z2_set_column<Cell_type,Row_access_option>::const_reverse_iterator
Z2_set_column<Cell_type,Row_access_option>::rend() const noexcept
{
	return column_.rend();
}

template<class Cell_type, class Row_access_option>
inline Z2_set_column<Cell_type,Row_access_option> &Z2_set_column<Cell_type,Row_access_option>::operator+=(Z2_set_column const &column)
{
	for (const Cell& v : column.column_){
		auto it = column_.find(v);
		if (it != column_.end()){
			_delete_cell(it);
		} else {
			_insert_cell(v.get_row_index());
		}
	}

	return *this;
}

template<class Cell_type, class Row_access_option>
inline Z2_set_column<Cell_type,Row_access_option> &Z2_set_column<Cell_type,Row_access_option>::operator*=(unsigned int v)
{
	if (v % 2 == 0){
		if constexpr (Row_access_option::isActive_){
			for (const Cell& cell : column_)
				Row_access_option::unlink(cell);
		}
		column_.clear();
	}

	return *this;
}

template<class Cell_type, class Row_access_option>
inline Z2_set_column<Cell_type,Row_access_option> &Z2_set_column<Cell_type,Row_access_option>::operator=(Z2_set_column other)
{
	static_assert (!Row_access_option::isActive_, "= assignement not enabled with row access option.");

	std::swap(dim_, other.dim_);
	column_.swap(other.column_);
	return *this;
}

template<class Cell_type, class Row_access_option>
inline void Z2_set_column<Cell_type,Row_access_option>::_delete_cell(iterator &it)
{
	if constexpr (Row_access_option::isActive_)
		Row_access_option::unlink(*it);
	column_.erase(it);
}

template<class Cell_type, class Row_access_option>
inline void Z2_set_column<Cell_type,Row_access_option>::_insert_cell(index rowIndex)
{
	if constexpr (Row_access_option::isActive_){
		auto it = column_.emplace_hint(column_.end(), Row_access_option::columnIndex_, rowIndex);
		Row_access_option::insert_cell(rowIndex, *it);
	} else {
		column_.emplace_hint(column_.end(), rowIndex);
	}
}

} //namespace persistence_matrix
} //namespace Gudhi

template<class Cell_type, class Row_access_option>
struct std::hash<Gudhi::persistence_matrix::Z2_set_column<Cell_type,Row_access_option> >
{
	size_t operator()(const Gudhi::persistence_matrix::Z2_set_column<Cell_type,Row_access_option>& column) const
	{
		std::size_t seed = 0;
		for (auto& cell : column){
			seed ^= std::hash<unsigned int>()(cell.get_row_index()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		}
		return seed;
	}
};

#endif // Z2_SET_COLUMN_H

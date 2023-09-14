/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef INTRUSIVE_LIST_COLUMN_H
#define INTRUSIVE_LIST_COLUMN_H

#include <iostream>
#include <vector>

#include <boost/intrusive/list.hpp>
#include <gudhi/Simple_object_pool.h>

#include "../utilities/utilities.h"
// #include "cell.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Field_element_type, class Cell_type, class Row_access_option>
class Intrusive_list_column : public Row_access_option
{
public:
	using base_hook_matrix_list_column = boost::intrusive::list_base_hook <
				boost::intrusive::tag < matrix_column_tag >
			, boost::intrusive::link_mode < boost::intrusive::safe_link >
			>;
//	using Cell = Intrusive_list_cell<Field_element_type>;
	using Cell = Cell_type;
	using Column_type = boost::intrusive::list <
							Cell
						  , boost::intrusive::constant_time_size<false>
						  , boost::intrusive::base_hook< base_hook_matrix_list_column >  >;
	using iterator = typename Column_type::iterator;
	using const_iterator = typename Column_type::const_iterator;
	using reverse_iterator = typename Column_type::reverse_iterator;
	using const_reverse_iterator = typename Column_type::const_reverse_iterator;

	Intrusive_list_column();
	template<class Container_type>
	Intrusive_list_column(const Container_type& nonZeroRowIndices);
	template<class Container_type>
	Intrusive_list_column(const Container_type& nonZeroRowIndices, dimension_type dimension);
	template<class Row_container_type>
	Intrusive_list_column(index columnIndex, Row_container_type &rowContainer);
	template<class Container_type, class Row_container_type>
	Intrusive_list_column(index columnIndex, const Container_type& nonZeroRowIndices, Row_container_type &rowContainer);
	template<class Container_type, class Row_container_type>
	Intrusive_list_column(index columnIndex, const Container_type& nonZeroRowIndices, dimension_type dimension, Row_container_type &rowContainer);
	Intrusive_list_column(const Intrusive_list_column& column);
	Intrusive_list_column(const Intrusive_list_column& column, index columnIndex);
	template<class Row_container_type>
	Intrusive_list_column(const Intrusive_list_column& column, index columnIndex, Row_container_type &rowContainer);
	Intrusive_list_column(Intrusive_list_column&& column) noexcept;
	~Intrusive_list_column();

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
	Intrusive_list_column& operator+=(Cell_range const &column);
//	Intrusive_list_column& operator+=(Intrusive_list_column const &column);
	friend Intrusive_list_column operator+(Intrusive_list_column column1, Intrusive_list_column const& column2){
		column1 += column2;
		return column1;
	}

	Intrusive_list_column& operator*=(unsigned int v);
	friend Intrusive_list_column operator*(Intrusive_list_column column, unsigned int const& v){
		column *= v;
		return column;
	}
	friend Intrusive_list_column operator*(unsigned int const& v, Intrusive_list_column column){
		column *= v;
		return column;
	}

	//this = v * this + column
	template<class Cell_range>
	Intrusive_list_column& multiply_and_add(const Field_element_type& v, const Cell_range& column);
//	Intrusive_list_column& multiply_and_add(const Field_element_type& v, const Intrusive_list_column& column);
	//this = this + column * v
	template<class Cell_range>
	Intrusive_list_column& multiply_and_add(const Cell_range& column, const Field_element_type& v);
//	Intrusive_list_column& multiply_and_add(const Intrusive_list_column& column, const Field_element_type& v);

	friend bool operator==(const Intrusive_list_column& c1, const Intrusive_list_column& c2){
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
	friend bool operator<(const Intrusive_list_column& c1, const Intrusive_list_column& c2){
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

	//Disabled with row access.
	Intrusive_list_column& operator=(const Intrusive_list_column& other);

	friend void swap(Intrusive_list_column& col1, Intrusive_list_column& col2){
		swap(static_cast<Row_access_option&>(col1),
			 static_cast<Row_access_option&>(col2));
		std::swap(col1.dim_, col2.dim_);
		col1.column_.swap(col2.column_);
	}

protected:
	dimension_type dim_;
	Column_type column_;
//	inline static Simple_object_pool<Column_type> columnPool_;
	inline static Simple_object_pool<Cell> cellPool_;

	void _delete_cell(iterator& it);
	void _insert_cell(const Field_element_type& value, index rowIndex, const iterator& position);

private:
	//Cloner object function for boost intrusive container
	struct new_cloner
	{
		Cell *operator()(const Cell &clone_this)
		{
//			return new Cell(clone_this);
			return cellPool_.construct(clone_this);
		}
	};

	//The disposer object function for boost intrusive container
	struct delete_disposer
	{
		delete_disposer(){};
		delete_disposer(Intrusive_list_column* col) : col_(col)
		{};

		void operator()(Cell *delete_this)
		{
			if constexpr (Row_access_option::isActive_)
				col_->unlink(delete_this);
//			delete delete_this;
			cellPool_.destroy(delete_this);
		}

		Intrusive_list_column* col_;
	};
};

template<class Field_element_type, class Cell_type, class Row_access_option>
inline Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::Intrusive_list_column()
	: dim_(0)
{
//	static_assert(!Row_access_option::isActive_, "When row access option enabled, a row container has to be provided.");
}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Container_type>
inline Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::Intrusive_list_column(
		const Container_type &nonZeroRowIndices)
	: dim_(nonZeroRowIndices.size() == 0 ? 0 : nonZeroRowIndices.size() - 1)
{
//	static_assert(!Row_access_option::isActive_, "When row access option enabled, a row container has to be provided.");

	for (const auto& p : nonZeroRowIndices){
		_insert_cell(p.second, p.first, column_.end());
	}
}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Container_type>
inline Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::Intrusive_list_column(
		const Container_type &nonZeroRowIndices, dimension_type dimension)
	: dim_(dimension)
{
//	static_assert(!Row_access_option::isActive_, "When row access option enabled, a row container has to be provided.");

	for (const auto& p : nonZeroRowIndices){
		_insert_cell(p.second, p.first, column_.end());
	}
}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Row_container_type>
inline Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::Intrusive_list_column(
		index columnIndex, Row_container_type &rowContainer)
	: Row_access_option(columnIndex, rowContainer), dim_(0)
{}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Container_type, class Row_container_type>
inline Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::Intrusive_list_column(
		index columnIndex, const Container_type &nonZeroRowIndices, Row_container_type &rowContainer)
	: Row_access_option(columnIndex, rowContainer), dim_(nonZeroRowIndices.size() == 0 ? 0 : nonZeroRowIndices.size() - 1)
{
	for (const auto& p : nonZeroRowIndices){
		_insert_cell(p.second, p.first, column_.end());
	}
}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Container_type, class Row_container_type>
inline Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::Intrusive_list_column(
		index columnIndex, const Container_type &nonZeroRowIndices, dimension_type dimension, Row_container_type &rowContainer)
	: Row_access_option(columnIndex, rowContainer), dim_(dimension)
{
	for (const auto& p : nonZeroRowIndices){
		_insert_cell(p.second, p.first, column_.end());
	}
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::Intrusive_list_column(
		const Intrusive_list_column &column)
	: dim_(column.dim_)
{
	static_assert(!Row_access_option::isActive_,
			"Copy constructor not available when row access option enabled.");

	column_.clone_from(column.column_, new_cloner(), delete_disposer());
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::Intrusive_list_column(
		const Intrusive_list_column &column, index columnIndex)
	: Row_access_option(columnIndex, *column.rows_),
	  dim_(column.dim_)
{
	for (const Cell& cell : column.column_){
		_insert_cell(cell.get_element(), cell.get_row_index(), column_.end());
	}
}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Row_container_type>
inline Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::Intrusive_list_column(
		const Intrusive_list_column &column, index columnIndex, Row_container_type &rowContainer)
	: Row_access_option(columnIndex, rowContainer),
	  dim_(column.dim_)
{
	for (const Cell& cell : column.column_){
		_insert_cell(cell.get_element(), cell.get_row_index(), column_.end());
	}
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::Intrusive_list_column(
		Intrusive_list_column &&column) noexcept
	: Row_access_option(std::move(column)),
	  dim_(std::exchange(column.dim_, 0)),
	  column_(std::move(column.column_))
{}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::~Intrusive_list_column()
{
//	for (iterator c_it = column_.begin(); c_it != column_.end(); ){
//		_delete_cell(c_it);
//	}
	clear();
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline std::vector<Field_element_type> Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::get_content(int columnLength) const
{
	if (columnLength < 0) columnLength = column_.back().get_row_index() + 1;

	std::vector<Field_element_type> container(columnLength);
	for (auto it = column_.begin(); it != column_.end() && it->get_row_index() < static_cast<index>(columnLength); ++it){
		container[it->get_row_index()] = it->get_element();
	}
	return container;
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline bool Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::is_non_zero(index rowIndex) const
{
	for (const Cell& cell : column_)
		if (cell.get_row_index() == rowIndex) return true;

	return false;
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline bool Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::is_empty() const
{
	return column_.empty();
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline dimension_type Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::get_dimension() const
{
	return dim_;
}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Map_type>
inline void Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::reorder(Map_type &valueMap)
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
inline typename Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::iterator
Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::begin() noexcept
{
	return column_.begin();
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline typename Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::const_iterator
Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::begin() const noexcept
{
	return column_.begin();
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline typename Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::iterator
Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::end() noexcept
{
	return column_.end();
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline typename Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::const_iterator
Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::end() const noexcept
{
	return column_.end();
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline typename Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::reverse_iterator
Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::rbegin() noexcept
{
	return column_.rbegin();
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline typename Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::const_reverse_iterator
Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::rbegin() const noexcept
{
	return column_.rbegin();
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline typename Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::reverse_iterator
Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::rend() noexcept
{
	return column_.rend();
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline typename Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::const_reverse_iterator
Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::rend() const noexcept
{
	return column_.rend();
}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Cell_range>
inline Intrusive_list_column<Field_element_type,Cell_type,Row_access_option> &
Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::operator+=(Cell_range const &column)
//Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::operator+=(Intrusive_list_column const &column)
{
//	const Column_type& sc = column.column_;

	auto itTarget = column_.begin();
	auto itSource = column.begin();
	while (itTarget != column_.end() && itSource != column.end())
	{
		if (itTarget->get_row_index() < itSource->get_row_index()) {
			++itTarget;
		} else if (itTarget->get_row_index() > itSource->get_row_index()) {
			_insert_cell(itSource->get_element(), itSource->get_row_index(), itTarget);
			++itSource;
		} else {
			itTarget->get_element() += itSource->get_element();
			if (itTarget->get_element() == Field_element_type::get_additive_identity()){
				_delete_cell(itTarget);
			} else {
				if constexpr (Row_access_option::isActive_)
						Row_access_option::update_cell(*itTarget);
				++itTarget;
			}
			++itSource;
		}
	}

	while (itSource != column.end()) {
		_insert_cell(itSource->get_element(), itSource->get_row_index(), column_.end());
		++itSource;
	}

	return *this;
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline Intrusive_list_column<Field_element_type,Cell_type,Row_access_option> &
Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::operator*=(unsigned int v)
{
//	v %= Field_element_type::get_characteristic();		//don't work because of multifields...
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
inline Intrusive_list_column<Field_element_type,Cell_type,Row_access_option> &
Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::multiply_and_add(const Field_element_type& val, const Cell_range& column)
{
	if (val == 0u) {
		clear();
	}

	auto itTarget = column_.begin();
	auto itSource = column.begin();
	while (itTarget != column_.end() && itSource != column.end())
	{
		if (itTarget->get_row_index() < itSource->get_row_index()) {
			itTarget->get_element() *= val;
			if constexpr (Row_access_option::isActive_)
					Row_access_option::update_cell(*itTarget);
			++itTarget;
		} else if (itTarget->get_row_index() > itSource->get_row_index()) {
			_insert_cell(itSource->get_element(), itSource->get_row_index(), itTarget);
			++itSource;
		} else {
			itTarget->get_element() *= val;
			itTarget->get_element() += itSource->get_element();
			if (itTarget->get_element() == Field_element_type::get_additive_identity()){
				_delete_cell(itTarget);
			} else {
				if constexpr (Row_access_option::isActive_)
						Row_access_option::update_cell(*itTarget);
				++itTarget;
			}
			++itSource;
		}
	}

	while (itTarget != column_.end()){
		itTarget->get_element() *= val;
		if constexpr (Row_access_option::isActive_)
				Row_access_option::update_cell(*itTarget);
		itTarget++;
	}

	while (itSource != column.end()) {
		_insert_cell(itSource->get_element(), itSource->get_row_index(), column_.end());
		++itSource;
	}

	return *this;
}

template<class Field_element_type, class Cell_type, class Row_access_option>
template<class Cell_range>
inline Intrusive_list_column<Field_element_type,Cell_type,Row_access_option> &
Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::multiply_and_add(const Cell_range& column, const Field_element_type& val)
//Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::multiply_and_add(const Intrusive_list_column& column, const Field_element_type& val)
{
	if (val == 0u) {
		return *this;
	}

//	const Column_type& sc = column.column_;
	auto itTarget = column_.begin();
	auto itSource = column.begin();
	while (itTarget != column_.end() && itSource != column.end())
	{
		if (itTarget->get_row_index() < itSource->get_row_index()) {
			++itTarget;
		} else if (itTarget->get_row_index() > itSource->get_row_index()) {
			_insert_cell(itSource->get_element() * val, itSource->get_row_index(), itTarget);
			++itSource;
		} else {
			itTarget->get_element() += (itSource->get_element() * val);
			if (itTarget->get_element() == Field_element_type::get_additive_identity()){
				_delete_cell(itTarget);
			} else {
				if constexpr (Row_access_option::isActive_)
						Row_access_option::update_cell(*itTarget);
				++itTarget;
			}
			++itSource;
		}
	}

	while (itSource != column.end()) {
		_insert_cell(itSource->get_element() * val, itSource->get_row_index(), column_.end());
		++itSource;
	}

	return *this;
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline Intrusive_list_column<Field_element_type,Cell_type,Row_access_option> &
Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::operator=(const Intrusive_list_column& other)
{
	static_assert (!Row_access_option::isActive_, "= assignement not enabled with row access option.");

	dim_ = other.dim_;
	column_.clone_from(other.column_, new_cloner(), delete_disposer());
	return *this;
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline void Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::_delete_cell(iterator &it)
{
//	iterator tmp_it = it;
//	++it;
//	Cell* tmp_ptr = &(*tmp_it);
//	if constexpr (Row_access_option::isActive_)
//		Row_access_option::unlink(tmp_ptr);
//	column_.erase(tmp_it);
//	delete tmp_ptr;
	it = column_.erase_and_dispose(it, delete_disposer(this));
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline void Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::_insert_cell(
		const Field_element_type &value, index rowIndex, const iterator &position)
{
	if constexpr (Row_access_option::isActive_){
//		Cell *new_cell = new Cell(value, Row_access_option::columnIndex_, rowIndex);
		Cell *new_cell = cellPool_.construct(value, Row_access_option::columnIndex_, rowIndex);
		column_.insert(position, *new_cell);
		Row_access_option::insert_cell(rowIndex, new_cell);
	} else {
//		Cell *new_cell = new Cell(value, rowIndex);
		Cell *new_cell = cellPool_.construct(value, rowIndex);
		column_.insert(position, *new_cell);
	}
}

template<class Field_element_type, class Cell_type, class Row_access_option>
inline void Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>::clear()
{
	column_.clear_and_dispose(delete_disposer(this));
}

} //namespace persistence_matrix
} //namespace Gudhi

template<class Field_element_type, class Cell_type, class Row_access_option>
struct std::hash<Gudhi::persistence_matrix::Intrusive_list_column<Field_element_type,Cell_type,Row_access_option> >
{
	size_t operator()(const Gudhi::persistence_matrix::Intrusive_list_column<Field_element_type,Cell_type,Row_access_option>& column) const
	{
		std::size_t seed = 0;
		for (auto& cell : column){
			seed ^= std::hash<unsigned int>()(cell.get_row_index() * static_cast<unsigned int>(cell.get_element())) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		}
		return seed;
	}
};

#endif // INTRUSIVE_LIST_COLUMN_H

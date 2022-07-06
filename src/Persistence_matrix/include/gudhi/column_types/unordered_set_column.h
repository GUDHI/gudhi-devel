/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef UNORDEREDSETCOLUMN_H
#define UNORDEREDSETCOLUMN_H

#include <iostream>
#include <list>
#include <unordered_set>
#include <algorithm>

#include "../utilities.h"
#include "../Zp_field.h"
#include "cell.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Field_element_type = Zp_field_element<11> >
class Unordered_set_column
{
public:
	using Cell = Base_cell<Field_element_type>;

	Unordered_set_column();
	template<class Boundary_type>
	Unordered_set_column(Boundary_type& boundary);
	Unordered_set_column(Unordered_set_column& column);
	Unordered_set_column(Unordered_set_column&& column) noexcept;

	std::vector<Field_element_type> get_content(unsigned int columnLength);
	bool is_non_zero(index rowIndex) const;
	bool is_empty();
	dimension_type get_dimension() const;
	int get_pivot();
	Field_element_type get_pivot_value();
	void clear();
	void clear(index rowIndex);
	void reorder(std::vector<index>& valueMap);

	Unordered_set_column& operator+=(Unordered_set_column const &column);
	template<class Friend_field_element_type>
	friend Unordered_set_column<Friend_field_element_type> operator+(Unordered_set_column<Friend_field_element_type> column1, Unordered_set_column<Friend_field_element_type> const& column2);
	Unordered_set_column& operator*=(unsigned int const &v);
	template<class Friend_field_element_type>
	friend Unordered_set_column<Friend_field_element_type> operator*(Unordered_set_column<Friend_field_element_type> column, unsigned int const& v);
	template<class Friend_field_element_type>
	friend Unordered_set_column<Friend_field_element_type> operator*(unsigned int const& v, Unordered_set_column<Friend_field_element_type> const column);

	Unordered_set_column& operator=(Unordered_set_column other);

	template<class Friend_field_element_type>
	friend void swap(Unordered_set_column<Friend_field_element_type>& col1,
					 Unordered_set_column<Friend_field_element_type>& col2);

private:
	int dim_;
	std::unordered_set<Cell> column_;
	bool pivotChanged_;
	Cell pivot_;
};

template<class Field_element_type>
inline Unordered_set_column<Field_element_type>::Unordered_set_column()
	: dim_(0), pivotChanged_(false)
{}

template<class Field_element_type>
template<class Boundary_type>
inline Unordered_set_column<Field_element_type>::Unordered_set_column(Boundary_type &boundary)
	: dim_(boundary.size() == 0 ? 0 : boundary.size() - 1),
	  column_(boundary.size()),
	  pivotChanged_(false),
	  pivot_(boundary.empty() ? 0 : Cell(boundary.rbegin()->second, boundary.rbegin()->first))
{
	for (std::pair<index,Field_element_type>& p : boundary){
		column_.insert(Cell(p.second, p.first));
	}
}

template<class Field_element_type>
inline Unordered_set_column<Field_element_type>::Unordered_set_column(Unordered_set_column &column)
	: dim_(column.dim_),
	  column_(column.column_),
	  pivotChanged_(column.pivotChanged_),
	  pivot_(column.pivot_)
{}

template<class Field_element_type>
inline Unordered_set_column<Field_element_type>::Unordered_set_column(Unordered_set_column &&column) noexcept
	: dim_(std::exchange(column.dim_, 0)),
	  column_(std::move(column.column_)),
	  pivotChanged_(std::exchange(column.pivotChanged_, 0)),
	  pivot_(std::move(column.pivot_))
{}

template<class Field_element_type>
inline std::vector<Field_element_type> Unordered_set_column<Field_element_type>::get_content(unsigned int columnLength)
{
	std::vector<Field_element_type> container(columnLength);
	for (auto it = column_.begin(); it != column_.end() && it->get_row_index() < columnLength; ++it){
		container[it->get_row_index()] = it->element();
	}
	return container;
}

template<class Field_element_type>
inline bool Unordered_set_column<Field_element_type>::is_non_zero(index rowIndex) const
{
	return column_.find(Cell(0, rowIndex)) != column_.end();
}

template<class Field_element_type>
inline bool Unordered_set_column<Field_element_type>::is_empty()
{
	return column_.empty();
}

template<class Field_element_type>
inline dimension_type Unordered_set_column<Field_element_type>::get_dimension() const
{
	return dim_;
}

template<class Field_element_type>
inline int Unordered_set_column<Field_element_type>::get_pivot()
{
	if (pivotChanged_){
		pivot_ = column_.size() == 0 ?
					0
				  : *std::max_element(column_.begin(), column_.end());
		pivotChanged_ = false;
	}

	if (pivot_.get_element() == 0) return -1;
	return pivot_.get_row_index();
}

template<class Field_element_type>
inline Field_element_type Unordered_set_column<Field_element_type>::get_pivot_value()
{
	if (pivotChanged_){
		pivot_ = column_.size() == 0 ?
					-1
				  : *std::max_element(column_.begin(), column_.end());
		pivotChanged_ = false;
	}

	return pivot_.get_element();
}

template<class Field_element_type>
inline void Unordered_set_column<Field_element_type>::clear()
{
	column_.clear();
	pivot_ = 0;
	pivotChanged_ = false;
}

template<class Field_element_type>
inline void Unordered_set_column<Field_element_type>::clear(index rowIndex)
{
	column_.erase(Cell(0, rowIndex));
	if (rowIndex == pivot_.get_row_index()) pivotChanged_ = true;
}

template<class Field_element_type>
inline void Unordered_set_column<Field_element_type>::reorder(std::vector<index> &valueMap)
{
	std::unordered_set<Cell> newSet;
	for (const Cell& v : column_) newSet.insert(Cell(v.get_element(), valueMap.at(v.get_row_index())));
	column_.swap(newSet);
	pivotChanged_ = true;
}

template<class Field_element_type>
inline Unordered_set_column<Field_element_type> &Unordered_set_column<Field_element_type>::operator+=(Unordered_set_column const &column)
{
	for (const Cell& v : column.column_){
		auto c = column_.find(v);
		if (c != column_.end()){
			c->get_element() += v.get_element();
			if (c->get_element() == 0) column_.erase(c);
			if (v == pivot_) pivotChanged_ = true;
		} else {
			column_.insert(v);
			if (pivot_ < v){
				pivot_ = v;
				pivotChanged_ = false;
			}
		}
	}

	return *this;
}

template<class Field_element_type>
inline Unordered_set_column<Field_element_type> &Unordered_set_column<Field_element_type>::operator*=(unsigned int const &v)
{
	v %= Field_element_type::get_characteristic();

	if (v == 0) {
		column_.clear();
		pivot_ = 0;
		pivotChanged_ = false;
		return *this;
	}

	for (Cell& cell : column_){
		cell.get_element() *= v;
	}

	pivot_.get_element() *= v;

	return *this;
}

template<class Field_element_type>
inline Unordered_set_column<Field_element_type> &Unordered_set_column<Field_element_type>::operator=(Unordered_set_column other)
{
	std::swap(dim_, other.dim_);
	std::swap(column_, other.column_);
	std::swap(pivotChanged_, other.pivotChanged_);
	std::swap(pivot_, other.pivot_);
	return *this;
}

template<class Friend_field_element_type>
Unordered_set_column<Friend_field_element_type> operator+(
		Unordered_set_column<Friend_field_element_type> column1,
		Unordered_set_column<Friend_field_element_type> const& column2)
{
	column1 += column2;
	return column1;
}

template<class Friend_field_element_type>
Unordered_set_column<Friend_field_element_type> operator*(
		Unordered_set_column<Friend_field_element_type> column, unsigned int const& v)
{
	column *= v;
	return column;
}

template<class Friend_field_element_type>
Unordered_set_column<Friend_field_element_type> operator*(
		unsigned int const& v, Unordered_set_column<Friend_field_element_type> column)
{
	column *= v;
	return column;
}

template<class Friend_field_element_type>
inline void swap(Unordered_set_column<Friend_field_element_type>& col1,
				 Unordered_set_column<Friend_field_element_type>& col2)
{
	std::swap(col1.dim_, col2.dim_);
	col1.column_.swap(col2.column_);
	std::swap(col1.pivotChanged_, col2.pivotChanged_);
	std::swap(col1.pivot_, col2.pivot_);
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // UNORDEREDSETCOLUMN_H

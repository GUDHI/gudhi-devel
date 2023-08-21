/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CELL_H
#define CELL_H

#include <boost/intrusive/list.hpp>
#include <boost/intrusive/set.hpp>

#include "../utilities/utilities.h"

namespace Gudhi {
namespace persistence_matrix {

struct matrix_row_tag;
struct matrix_column_tag;

using base_hook_matrix_row = boost::intrusive::list_base_hook<
			boost::intrusive::tag < matrix_row_tag >
		  , boost::intrusive::link_mode < boost::intrusive::auto_unlink >
	>;
using base_hook_matrix_list_column = boost::intrusive::list_base_hook <
			boost::intrusive::tag < matrix_column_tag >
		  , boost::intrusive::link_mode < boost::intrusive::safe_link >
		>;
using base_hook_matrix_set_column = boost::intrusive::set_base_hook <
			boost::intrusive::tag < matrix_column_tag >
		  , boost::intrusive::link_mode < boost::intrusive::safe_link >
		>;

class Z2_base_cell
{
public:
	Z2_base_cell() : rowIndex_(0){};
	Z2_base_cell(index rowIndex) : rowIndex_(rowIndex){};
	Z2_base_cell(const Z2_base_cell& cell) : rowIndex_(cell.rowIndex_){};
	Z2_base_cell(Z2_base_cell&& cell) noexcept : rowIndex_(std::exchange(cell.rowIndex_, 0)){};

	index get_row_index() const{
		return rowIndex_;
	};

	void set_row_index(const index &rowIndex){
		rowIndex_ = rowIndex;
	};

	Z2_base_cell& operator=(Z2_base_cell other){
		std::swap(rowIndex_, other.rowIndex_);
		return *this;
	};

	friend bool operator<(const Z2_base_cell& c1, const Z2_base_cell& c2){
		return c1.get_row_index() < c2.get_row_index();
	}
	friend bool operator==(const Z2_base_cell& c1, const Z2_base_cell& c2){
		return c1.get_row_index() == c2.get_row_index();
	}

	operator unsigned int() const{
		return rowIndex_;
	}

private:
	index rowIndex_;
};

class Z2_row_cell : public Z2_base_cell
{
public:
	Z2_row_cell()
		: Z2_base_cell(), columnIndex_(0){};
	Z2_row_cell(index columnIndex, index rowIndex)
		: Z2_base_cell(rowIndex), columnIndex_(columnIndex){};
	Z2_row_cell(const Z2_row_cell& cell)
		: Z2_base_cell(cell), columnIndex_(cell.columnIndex_){};
	Z2_row_cell(Z2_row_cell&& cell) noexcept
		: Z2_base_cell(std::move(cell)), columnIndex_(std::exchange(cell.columnIndex_, 0)){};

	index get_column_index() const{
		return columnIndex_;
	};
//	void set_column_index(index newIndex){
//		columnIndex_ = newIndex;
//	};

	Z2_row_cell& operator=(Z2_row_cell other){
		std::swap(columnIndex_, other.columnIndex_);
		Z2_base_cell::operator=(other);
		return *this;
	};

//	friend bool operator<(const Z2_row_cell& c1, const Z2_row_cell& c2){
//		if (c1.get_row_index() == c2.get_row_index())
//			return c1.get_column_index() < c2.get_column_index();
//		return c1.get_row_index() < c2.get_row_index();
//	}
//	friend bool operator==(const Z2_row_cell& c1, const Z2_row_cell& c2){
//		return c1.get_row_index() == c2.get_row_index() &&
//				c1.get_column_index() == c2.get_column_index();
//	}

private:
	index columnIndex_;
};

struct Z2_intrusive_row_cell : public Z2_row_cell, public base_hook_matrix_row
{
	Z2_intrusive_row_cell()
		: Z2_row_cell(){};
	Z2_intrusive_row_cell(index columnIndex, index rowIndex)
		: Z2_row_cell(columnIndex, rowIndex){};
	Z2_intrusive_row_cell(const Z2_intrusive_row_cell& cell)
		: Z2_row_cell(cell){};
	Z2_intrusive_row_cell(Z2_intrusive_row_cell&& cell) noexcept
		: Z2_row_cell(std::move(cell)){};

	Z2_intrusive_row_cell& operator=(Z2_intrusive_row_cell other){
		Z2_row_cell::operator=(other);
		return *this;
	};

	// using base_hook_matrix_row = base_hook_matrix_row;		//why ?????
};

template<class Field_element_type>
class Base_cell : public Z2_base_cell
{
public:
	Base_cell()
		: Z2_base_cell(), element_(0){};
	Base_cell(const Field_element_type& element, index rowIndex)
		: Z2_base_cell(rowIndex), element_(element){};
	Base_cell(const Base_cell& cell)
		: Z2_base_cell(cell), element_(cell.element_){};
	Base_cell(Base_cell&& cell) noexcept
		: Z2_base_cell(std::move(cell)), element_(std::move(cell.element_)){};

	Field_element_type& get_element(){
		return element_;
	};

	const Field_element_type& get_element() const{
		return element_;
	};

	Base_cell& operator=(Base_cell other){
		swap(element_, other.element_);
		Z2_base_cell::operator=(other);
		return *this;
	};

	void set_element(const Field_element_type &element){
		element_ = element;
	}

	operator std::pair<unsigned int,Field_element_type>() const{
		return {Z2_base_cell::get_row_index(), element_};
	}

private:
	Field_element_type element_;
};

template<class Field_element_type>
class Row_cell : public Z2_row_cell
{
public:
	using Field_type = Field_element_type;

	Row_cell()
		: Z2_row_cell(), element_(0){};
	Row_cell(const Field_element_type& element, index columnIndex, index rowIndex)
		: Z2_row_cell(columnIndex, rowIndex), element_(element){};
	Row_cell(const Row_cell& cell)
		: Z2_row_cell(cell), element_(cell.element_){};
	Row_cell(Row_cell&& cell) noexcept
		: Z2_row_cell(std::move(cell)), element_(std::move(cell.element_)){};

	Field_element_type& get_element(){
		return element_;
	};

	const Field_element_type& get_element() const{
		return element_;
	};

	void set_element(const Field_element_type &element){
		element_ = element;
	}

	Row_cell& operator=(Row_cell other){
		swap(element_, other.element_);
		Z2_row_cell::operator=(other);
		return *this;
	};

	operator std::pair<unsigned int,Field_element_type>() const{
		return {Z2_base_cell::get_row_index(), element_};
	}

private:
	Field_element_type element_;
};

template<class Field_element_type>
struct Intrusive_row_cell : public Row_cell<Field_element_type>, public base_hook_matrix_row
{
	Intrusive_row_cell()
		: Row_cell<Field_element_type>(){};
	Intrusive_row_cell(const Field_element_type& element, index columnIndex, index rowIndex)
		: Row_cell<Field_element_type>(element, columnIndex, rowIndex){};
	Intrusive_row_cell(const Intrusive_row_cell& cell)
		: Row_cell<Field_element_type>(cell){};
	Intrusive_row_cell(Intrusive_row_cell&& cell) noexcept
		: Row_cell<Field_element_type>(std::move(cell)){};

	Intrusive_row_cell& operator=(Intrusive_row_cell other){
		Row_cell<Field_element_type>::operator=(other);
		return *this;
	};

	// using base_hook_matrix_row = base_hook_matrix_row;		//why ?????
};

struct Z2_intrusive_list_cell : public Z2_base_cell, public base_hook_matrix_list_column
{
	Z2_intrusive_list_cell()
		: Z2_base_cell(){};
	Z2_intrusive_list_cell(index rowIndex)
		: Z2_base_cell(rowIndex){};
};

template<class Field_element_type>
struct Intrusive_list_cell : public Base_cell<Field_element_type>, public base_hook_matrix_list_column
{
	Intrusive_list_cell()
		: Base_cell<Field_element_type>(){};
	Intrusive_list_cell(const Field_element_type& element, index rowIndex)
		: Base_cell<Field_element_type>(element, rowIndex){};
};

template<class Base_cell_type>
struct Z2_intrusive_list_row_cell : public Base_cell_type, public base_hook_matrix_list_column
{
	Z2_intrusive_list_row_cell()
		: Base_cell_type(){};
	Z2_intrusive_list_row_cell(index columnIndex, index rowIndex)
		: Base_cell_type(columnIndex, rowIndex){};
};

template<class Base_cell_type>
struct Intrusive_list_row_cell : public Base_cell_type, public base_hook_matrix_list_column
{
	Intrusive_list_row_cell()
		: Base_cell_type(){};
	Intrusive_list_row_cell(const typename Base_cell_type::Field_type& element, index columnIndex, index rowIndex)
		: Base_cell_type(element, columnIndex, rowIndex){};
};

struct Z2_intrusive_set_cell : public Z2_base_cell, public base_hook_matrix_set_column
{
	Z2_intrusive_set_cell()
		: Z2_base_cell(){};
	Z2_intrusive_set_cell(index rowIndex)
		: Z2_base_cell(rowIndex){};
};

template<class Field_element_type>
struct Intrusive_set_cell : public Base_cell<Field_element_type>, public base_hook_matrix_set_column
{
	Intrusive_set_cell()
		: Base_cell<Field_element_type>(){};
	Intrusive_set_cell(const Field_element_type& element, index rowIndex)
		: Base_cell<Field_element_type>(element, rowIndex){};
};

template<class Base_cell_type>
struct Z2_intrusive_set_row_cell : public Base_cell_type, public base_hook_matrix_set_column
{
	Z2_intrusive_set_row_cell()
		: Base_cell_type(){};
	Z2_intrusive_set_row_cell(index columnIndex, index rowIndex)
		: Base_cell_type(columnIndex, rowIndex){};
};

template<class Base_cell_type>
struct Intrusive_set_row_cell : public Base_cell_type, public base_hook_matrix_set_column
{
	Intrusive_set_row_cell()
		: Base_cell_type(){};
	Intrusive_set_row_cell(const typename Base_cell_type::Field_type& element, index columnIndex, index rowIndex)
		: Base_cell_type(element, columnIndex, rowIndex){};
};

} //namespace persistence_matrix
} //namespace Gudhi

template<>
struct std::hash<Gudhi::persistence_matrix::Z2_base_cell>
{
	size_t operator()(const Gudhi::persistence_matrix::Z2_base_cell& cell) const
	{
		return std::hash<unsigned int>()(cell.get_row_index());
	}
};

template<>
struct std::hash<Gudhi::persistence_matrix::Z2_row_cell>
{
	size_t operator()(const Gudhi::persistence_matrix::Z2_row_cell& cell) const
	{
//		std::size_t seed = 0;
//		seed ^= std::hash<unsigned int>()(cell.get_row_index()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
//		seed ^= std::hash<unsigned int>()(cell.get_column_index()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
//		return seed;
		return std::hash<unsigned int>()(cell.get_row_index());
	}
};

template<class Field_element_type>
struct std::hash<Gudhi::persistence_matrix::Base_cell<Field_element_type> >
{
	size_t operator()(const Gudhi::persistence_matrix::Base_cell<Field_element_type>& cell) const
	{
		return std::hash<unsigned int>()(cell.get_row_index());
	}
};

template<class Field_element_type>
struct std::hash<Gudhi::persistence_matrix::Row_cell<Field_element_type> >
{
	size_t operator()(const Gudhi::persistence_matrix::Row_cell<Field_element_type>& cell) const
	{
//		std::size_t seed = 0;
//		seed ^= std::hash<unsigned int>()(cell.get_row_index()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
//		seed ^= std::hash<unsigned int>()(cell.get_column_index()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
//		return seed;
		return std::hash<unsigned int>()(cell.get_row_index());
	}
};

#endif // CELL_H

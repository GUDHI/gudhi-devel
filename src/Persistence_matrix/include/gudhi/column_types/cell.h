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

#include "../options.h"
#include "../utilities/utilities.h"
#include "../utilities/Zp_field.h"

namespace Gudhi {
namespace persistence_matrix {

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

	void setRowIndex(const index &rowIndex){
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

	Z2_row_cell& operator=(Z2_row_cell other){
		std::swap(columnIndex_, other.columnIndex_);
		Z2_base_cell::operator=(other);
		return *this;
	};

	friend bool operator==(const Z2_row_cell& c1, const Z2_row_cell& c2)
	{
		return c1.get_row_index() == c2.get_row_index() &&
				c1.get_column_index() == c2.get_column_index();
	}

private:
	index columnIndex_;
};

template<class Field_element_type = Zp_field_element<11> >
class Base_cell : public Z2_base_cell
{
public:
	Base_cell()
		: Z2_base_cell(), element_(0){};
	Base_cell(unsigned int element, index rowIndex)
		: Z2_base_cell(rowIndex), element_(element){};
	Base_cell(const Base_cell& cell)
		: Z2_base_cell(cell), element_(cell.element_){};
	Base_cell(Base_cell&& cell) noexcept
		: Z2_base_cell(std::move(cell)), element_(std::move(cell.element_)){};

	Field_element_type& get_element(){
		return element_;
	};

	unsigned int get_element_value() const{
		return element_;
	};

	Base_cell& operator=(Base_cell other){
		std::swap(element_, other.element_);
		Z2_base_cell::operator=(other);
		return *this;
	};

private:
	Field_element_type element_;
};

template<class Field_element_type = Zp_field_element<11> >
class Row_cell : public Z2_row_cell
{
public:
	Row_cell()
		: Z2_row_cell(), element_(0){};
	Row_cell(unsigned int element, index columnIndex, index rowIndex)
		: Z2_row_cell(columnIndex, rowIndex), element_(element){};
	Row_cell(const Row_cell& cell)
		: Z2_row_cell(cell), element_(cell.element_){};
	Row_cell(Row_cell&& cell) noexcept
		: Z2_row_cell(std::move(cell)), element_(std::move(cell.element_)){};

	Field_element_type& get_element(){
		return element_;
	};

	Field_element_type get_element_value() const{
		return element_;
	};

	Row_cell& operator=(Row_cell other){
		std::swap(element_, other.element_);
		Z2_row_cell::operator=(other);
		return *this;
	};

private:
	Field_element_type element_;
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
		return std::hash<unsigned int>()(cell.get_row_index());
	}
};

#endif // CELL_H

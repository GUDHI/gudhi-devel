/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022-23 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_MATRIX_CELL_H
#define PM_MATRIX_CELL_H

#include <utility>		//std::swap, std::exchange & std::move
#include <functional>	//std::hash

namespace Gudhi {
namespace persistence_matrix {

struct Dummy_cell_column_index_mixin{
	Dummy_cell_column_index_mixin(){}
	template<typename index>
	Dummy_cell_column_index_mixin([[maybe_unused]] index columnIndex){}
	Dummy_cell_column_index_mixin([[maybe_unused]] const Dummy_cell_column_index_mixin& cell){};
	Dummy_cell_column_index_mixin([[maybe_unused]] Dummy_cell_column_index_mixin&& cell){};

	Dummy_cell_column_index_mixin& operator=([[maybe_unused]] const Dummy_cell_column_index_mixin& other){
		return *this;
	};
};

struct Dummy_cell_field_element_mixin{
	Dummy_cell_field_element_mixin(){}
	template<class Field_element_type>
	Dummy_cell_field_element_mixin([[maybe_unused]] Field_element_type t){}
	Dummy_cell_field_element_mixin([[maybe_unused]] const Dummy_cell_field_element_mixin& cell){};
	Dummy_cell_field_element_mixin([[maybe_unused]] Dummy_cell_field_element_mixin&& cell){};

	Dummy_cell_field_element_mixin& operator=([[maybe_unused]] const Dummy_cell_field_element_mixin& other){
		return *this;
	};
};

template<typename index>
class Cell_column_index
{
public:
	Cell_column_index() : columnIndex_(-1) {};
	Cell_column_index(index columnIndex) : columnIndex_(columnIndex)
	{};
	Cell_column_index(const Cell_column_index& cell) : columnIndex_(cell.columnIndex_)
	{};
	Cell_column_index(Cell_column_index&& cell) noexcept 
		: columnIndex_(std::exchange(cell.columnIndex_, 0))
	{};

	index get_column_index() const{
		return columnIndex_;
	};

	Cell_column_index& operator=(Cell_column_index other){
		std::swap(columnIndex_, other.columnIndex_);
		return *this;
	};

private:
	index columnIndex_;
};

template<class Field_element_type>
class Cell_field_element
{
public:
	Cell_field_element() : element_(0) {};
	Cell_field_element(Field_element_type element) : element_(element)
	{};
	Cell_field_element(const Cell_field_element& cell) : element_(cell.element_)
	{};
	Cell_field_element(Cell_field_element&& cell) noexcept 
		: element_(std::move(cell.element_))
	{};

	Field_element_type& get_element(){
		return element_;
	};
	const Field_element_type& get_element() const{
		return element_;
	};
	void set_element(const Field_element_type &element){
		element_ = element;
	}

	Cell_field_element& operator=(Cell_field_element other){
		swap(element_, other.element_);
		return *this;
	};

private:
	Field_element_type element_;
};

template<class Master_matrix>
class Cell : public Master_matrix::Cell_column_index_option, 
			 public Master_matrix::Cell_field_element_option,
			 public Master_matrix::row_hook_type, 
			 public Master_matrix::column_hook_type
{
private:
	using col_opt = typename Master_matrix::Cell_column_index_option;
	using field_opt = typename Master_matrix::Cell_field_element_option;

public:
	using index = typename Master_matrix::index;
	using Field_element_type = typename Master_matrix::Field_type;

	// using base_hook_matrix_row = typename Master_matrix::row_hook_type;	//temporary during factorization process

	Cell(){};
	Cell(index rowIndex) 
		: col_opt(),
		  field_opt(),
		  rowIndex_(rowIndex)
	{};
	Cell(index columnIndex, index rowIndex) 
		: col_opt(columnIndex),
		  field_opt(),
		  rowIndex_(rowIndex)
	{};
	Cell(Field_element_type element, index rowIndex) 
		: col_opt(),
		  field_opt(element),
		  rowIndex_(rowIndex)
	{};
	Cell(Field_element_type element, index columnIndex, index rowIndex) 
		: col_opt(columnIndex),
		  field_opt(element),
		  rowIndex_(rowIndex)
	{};
	Cell(const Cell& cell) 
		: col_opt(static_cast<const col_opt&>(cell)),
		  field_opt(static_cast<const field_opt&>(cell)),
		  rowIndex_(cell.rowIndex_)
	{};
	Cell(Cell&& cell) noexcept 
		: col_opt(std::move(static_cast<col_opt&>(cell))),
		  field_opt(std::move(static_cast<field_opt&>(cell))),
		  rowIndex_(std::exchange(cell.rowIndex_, 0))
	{};

	index get_row_index() const{
		return rowIndex_;
	};
	void set_row_index(const index &rowIndex){
		rowIndex_ = rowIndex;
	};

	Cell& operator=(Cell other){
		col_opt::operator=(other);
		field_opt::operator=(other);
		std::swap(rowIndex_, other.rowIndex_);
		return *this;
	};

	friend bool operator<(const Cell& c1, const Cell& c2){
		return c1.get_row_index() < c2.get_row_index();
	}
	friend bool operator==(const Cell& c1, const Cell& c2){
		return c1.get_row_index() == c2.get_row_index();
	}

	operator unsigned int() const{
		return rowIndex_;
	}
	operator std::pair<unsigned int,Field_element_type>() const{
		if constexpr (Master_matrix::Option_list::is_z2){
			return {rowIndex_, 1};
		} else {
			return {rowIndex_, field_opt::element_};
		}
	}

private:
	index rowIndex_;
};

} //namespace persistence_matrix
} //namespace Gudhi

template<class Master_matrix>
struct std::hash<Gudhi::persistence_matrix::Cell<Master_matrix> >
{
	size_t operator()(const Gudhi::persistence_matrix::Cell<Master_matrix>& cell) const
	{
		return std::hash<unsigned int>()(cell.get_row_index());
	}
};

#endif // PM_MATRIX_CELL_H

/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022-23 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_INTRUSIVE_LIST_COLUMN_H
#define PM_INTRUSIVE_LIST_COLUMN_H

#include <vector>
#include <stdexcept>
#include <type_traits>
#include <utility>	//std::swap, std::move & std::exchange

#include <boost/intrusive/list.hpp>

#include <gudhi/Simple_object_pool.h>

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class Intrusive_list_column : public Master_matrix::Row_access_option, 
							  public Master_matrix::Column_dimension_option,
							  public Master_matrix::Chain_column_option
{
public:
	using Master = Master_matrix;
	using Field_element_type = typename std::conditional<
								  Master_matrix::Option_list::is_z2,
								  bool,
								  typename Master_matrix::Field_type
							   >::type;
	using index = typename Master_matrix::index;
	using dimension_type = typename Master_matrix::dimension_type;

	using Cell = typename Master_matrix::Cell_type;
	using Column_type = boost::intrusive::list <
							Cell, 
							boost::intrusive::constant_time_size<false>, 
							boost::intrusive::base_hook< typename Master_matrix::base_hook_matrix_list_column > 
						>;
	using iterator = typename Column_type::iterator;
	using const_iterator = typename Column_type::const_iterator;
	using reverse_iterator = typename Column_type::reverse_iterator;
	using const_reverse_iterator = typename Column_type::const_reverse_iterator;

	Intrusive_list_column();
	template<class Container_type = typename Master_matrix::boundary_type>
	Intrusive_list_column(const Container_type& nonZeroRowIndices);	//has to be a boundary for boundary, has no sense for chain if dimension is needed
	template<class Container_type = typename Master_matrix::boundary_type, class Row_container_type>
	Intrusive_list_column(index columnIndex, const Container_type& nonZeroRowIndices, Row_container_type &rowContainer);	//has to be a boundary for boundary, has no sense for chain if dimension is needed
	template<class Container_type = typename Master_matrix::boundary_type>
	Intrusive_list_column(const Container_type& nonZeroChainRowIndices, dimension_type dimension);	//dimension gets ignored for base
	template<class Container_type = typename Master_matrix::boundary_type, class Row_container_type>
	Intrusive_list_column(index columnIndex, const Container_type& nonZeroChainRowIndices, dimension_type dimension, Row_container_type &rowContainer);	//dimension gets ignored for base
	Intrusive_list_column(const Intrusive_list_column& column);
	template<class Row_container_type>
	Intrusive_list_column(const Intrusive_list_column& column, index columnIndex, Row_container_type &rowContainer);
	Intrusive_list_column(Intrusive_list_column&& column) noexcept;
	~Intrusive_list_column();

	std::vector<Field_element_type> get_content(int columnLength = -1) const;
	bool is_non_zero(index rowIndex) const;
	bool is_empty() const;
	std::size_t size() const;

	//****************
	//only for base and boundary
	template<class Map_type>
	void reorder(const Map_type& valueMap);	//used for lazy row swaps
	void clear();
	void clear(index rowIndex);
	//****************

	//****************
	//only for chain and boundary
	int get_pivot() const;
	Field_element_type get_pivot_value() const;
	//****************

	iterator begin() noexcept;
	const_iterator begin() const noexcept;
	iterator end() noexcept;
	const_iterator end() const noexcept;
	reverse_iterator rbegin() noexcept;
	const_reverse_iterator rbegin() const noexcept;
	reverse_iterator rend() noexcept;
	const_reverse_iterator rend() const noexcept;

	template<class Cell_range>
	Intrusive_list_column& operator+=(const Cell_range& column);	//for base & boundary except vector
	Intrusive_list_column& operator+=(Intrusive_list_column &column);	//for chain and vector
	friend Intrusive_list_column operator+(Intrusive_list_column column1, Intrusive_list_column& column2){
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
	Intrusive_list_column& multiply_and_add(const Field_element_type& val, const Cell_range& column);	//for base & boundary except vector
	Intrusive_list_column& multiply_and_add(const Field_element_type& val, Intrusive_list_column& column);	//for chain and vector
	//this = this + column * v
	template<class Cell_range>
	Intrusive_list_column& multiply_and_add(const Cell_range& column, const Field_element_type& val);	//for base & boundary except vector
	Intrusive_list_column& multiply_and_add(Intrusive_list_column& column, const Field_element_type& val);	//for chain and vector

	friend bool operator==(const Intrusive_list_column& c1, const Intrusive_list_column& c2){
		if (&c1 == &c2) return true;

		if constexpr (Master_matrix::Option_list::is_z2){
			return c1.column_ == c2.column_;
		} else {
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
	}
	friend bool operator<(const Intrusive_list_column& c1, const Intrusive_list_column& c2){
		if (&c1 == &c2) return false;

		if constexpr (Master_matrix::Option_list::is_z2){
			return c1.column_ < c2.column_;
		} else {
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
	}

	//Disabled with row access.
	Intrusive_list_column& operator=(const Intrusive_list_column& other);

	friend void swap(Intrusive_list_column& col1, Intrusive_list_column& col2){
		swap(static_cast<typename Master_matrix::Row_access_option&>(col1),
			 static_cast<typename Master_matrix::Row_access_option&>(col2));
		swap(static_cast<typename Master_matrix::Column_dimension_option&>(col1),
			 static_cast<typename Master_matrix::Column_dimension_option&>(col2));
		swap(static_cast<typename Master_matrix::Chain_column_option&>(col1),
			 static_cast<typename Master_matrix::Chain_column_option&>(col2));
		col1.column_.swap(col2.column_);
	}

protected:
	Column_type column_;
	inline static Simple_object_pool<Cell> cellPool_;

	void _delete_cell(iterator& it);
	void _insert_cell(const Field_element_type& value, index rowIndex, const iterator& position);
	void _insert_cell(index rowIndex, const iterator& position);
	template<class Cell_range>
	bool _add(const Cell_range& column);
	template<class Cell_range>
	bool _multiply_and_add(const Field_element_type& val, const Cell_range& column);
	template<class Cell_range>
	bool _multiply_and_add(const Cell_range& column, const Field_element_type& val);

private:
	using ra_opt = typename Master_matrix::Row_access_option;
	using dim_opt = typename Master_matrix::Column_dimension_option;
	using chain_opt = typename Master_matrix::Chain_column_option;

	//Cloner object function for boost intrusive container
	struct new_cloner
	{
		Cell *operator()(const Cell& clone_this)
		{
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
			if constexpr (Master_matrix::Option_list::has_row_access)
				col_->unlink(delete_this);
			cellPool_.destroy(delete_this);
		}

		Intrusive_list_column* col_;
	};
};

template<class Master_matrix>
inline Intrusive_list_column<Master_matrix>::Intrusive_list_column() : ra_opt(), dim_opt(), chain_opt()
{}

template<class Master_matrix>
template<class Container_type>
inline Intrusive_list_column<Master_matrix>::Intrusive_list_column(const Container_type &nonZeroRowIndices)
	: ra_opt(), 
	  dim_opt(nonZeroRowIndices.size() == 0 ? 0 : nonZeroRowIndices.size() - 1), 
	  chain_opt()
{
	static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type, 
						"Constructor not available for chain columns, please specify the dimension of the chain.");

	if constexpr (Master_matrix::Option_list::is_z2){
		for (index id : nonZeroRowIndices){
			_insert_cell(id, column_.end());
		}
	} else {
		for (const auto& p : nonZeroRowIndices){
			_insert_cell(p.second, p.first, column_.end());
		}
	}
}

template<class Master_matrix>
template<class Container_type, class Row_container_type>
inline Intrusive_list_column<Master_matrix>::Intrusive_list_column(
	index columnIndex, const Container_type &nonZeroRowIndices, Row_container_type &rowContainer) 
	: ra_opt(columnIndex, rowContainer), 
	  dim_opt(nonZeroRowIndices.size() == 0 ? 0 : nonZeroRowIndices.size() - 1),
	  chain_opt([&]{
			if constexpr (Master_matrix::Option_list::is_z2){
				return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : *std::prev(nonZeroRowIndices.end());
			} else {
				return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : std::prev(nonZeroRowIndices.end())->first;
			}
		}())
{
	static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type, 
						"Constructor not available for chain columns, please specify the dimension of the chain.");

	if constexpr (Master_matrix::Option_list::is_z2){
		for (index id : nonZeroRowIndices){
			_insert_cell(id, column_.end());
		}
	} else {
		for (const auto& p : nonZeroRowIndices){
			_insert_cell(p.second, p.first, column_.end());
		}
	}
}

template<class Master_matrix>
template<class Container_type>
inline Intrusive_list_column<Master_matrix>::Intrusive_list_column(
	const Container_type &nonZeroRowIndices, dimension_type dimension) 
	: ra_opt(), 
	  dim_opt(dimension),
	  chain_opt([&]{
			if constexpr (Master_matrix::Option_list::is_z2){
				return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : *std::prev(nonZeroRowIndices.end());
			} else {
				return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : std::prev(nonZeroRowIndices.end())->first;
			}
		}())
{
	if constexpr (Master_matrix::Option_list::is_z2){
		for (index id : nonZeroRowIndices){
			_insert_cell(id, column_.end());
		}
	} else {
		for (const auto& p : nonZeroRowIndices){
			_insert_cell(p.second, p.first, column_.end());
		}
	}
}

template<class Master_matrix>
template<class Container_type, class Row_container_type>
inline Intrusive_list_column<Master_matrix>::Intrusive_list_column(
	index columnIndex, const Container_type &nonZeroRowIndices, dimension_type dimension, Row_container_type &rowContainer) 
	: ra_opt(columnIndex, rowContainer), 
	  dim_opt(dimension),
	  chain_opt([&]{
			if constexpr (Master_matrix::Option_list::is_z2){
				return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : *std::prev(nonZeroRowIndices.end());
			} else {
				return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : std::prev(nonZeroRowIndices.end())->first;
			}
		}())
{
	if constexpr (Master_matrix::Option_list::is_z2){
		for (index id : nonZeroRowIndices){
			_insert_cell(id, column_.end());
		}
	} else {
		for (const auto& p : nonZeroRowIndices){
			_insert_cell(p.second, p.first, column_.end());
		}
	}
}

template<class Master_matrix>
inline Intrusive_list_column<Master_matrix>::Intrusive_list_column(const Intrusive_list_column &column) 
	: ra_opt(), 
	  dim_opt(static_cast<const dim_opt&>(column)), 
	  chain_opt(static_cast<const chain_opt&>(column))
{
	static_assert(!Master_matrix::Option_list::has_row_access,
			"Simple copy constructor not available when row access option enabled. Please specify the new column index and the row container.");

	column_.clone_from(column.column_, new_cloner(), delete_disposer(this));
}

template<class Master_matrix>
template<class Row_container_type>
inline Intrusive_list_column<Master_matrix>::Intrusive_list_column(
	const Intrusive_list_column &column, index columnIndex, Row_container_type &rowContainer) 
	: ra_opt(columnIndex, rowContainer), 
	  dim_opt(static_cast<const dim_opt&>(column)), 
	  chain_opt(static_cast<const chain_opt&>(column))
{
	for (const Cell& cell : column.column_){
		if constexpr (Master_matrix::Option_list::is_z2){
			_insert_cell(cell.get_row_index(), column_.end());
		} else {
			_insert_cell(cell.get_element(), cell.get_row_index(), column_.end());
		}
	}
}

template<class Master_matrix>
inline Intrusive_list_column<Master_matrix>::Intrusive_list_column(Intrusive_list_column &&column) noexcept 
	: ra_opt(std::move(static_cast<ra_opt&>(column))), 
	  dim_opt(std::move(static_cast<dim_opt&>(column))), 
	  chain_opt(std::move(static_cast<chain_opt&>(column))), 
	  column_(std::move(column.column_))
{}

template<class Master_matrix>
inline Intrusive_list_column<Master_matrix>::~Intrusive_list_column()
{
	column_.clear_and_dispose(delete_disposer(this));
}

template<class Master_matrix>
inline std::vector<typename Intrusive_list_column<Master_matrix>::Field_element_type> 
Intrusive_list_column<Master_matrix>::get_content(int columnLength) const
{
	if (columnLength < 0 && column_.size() > 0) columnLength = column_.back().get_row_index() + 1;
	else if (columnLength < 0) return std::vector<Field_element_type>();

	std::vector<Field_element_type> container(columnLength);
	for (auto it = column_.begin(); it != column_.end() && it->get_row_index() < static_cast<index>(columnLength); ++it){
		if constexpr (Master_matrix::Option_list::is_z2){
			container[it->get_row_index()] = 1;
		} else {
			container[it->get_row_index()] = it->get_element();
		}
	}
	return container;
}

template<class Master_matrix>
inline bool Intrusive_list_column<Master_matrix>::is_non_zero(index rowIndex) const
{
	//could be changed to dichotomic search as column is ordered by row index, 
	//but I am not sure if it is really worth it as there is no random access
	//and the columns should not be that long anyway.
	for (const Cell& cell : column_)
		if (cell.get_row_index() == rowIndex) return true;

	return false;
}

template<class Master_matrix>
inline bool Intrusive_list_column<Master_matrix>::is_empty() const
{
	return column_.empty();
}

template<class Master_matrix>
inline std::size_t Intrusive_list_column<Master_matrix>::size() const{
	return column_.size();
}

template<class Master_matrix>
template<class Map_type>
inline void Intrusive_list_column<Master_matrix>::reorder(const Map_type &valueMap)
{
	static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type, 
						"Method not available for chain columns.");

	for (auto it = column_.begin(); it != column_.end(); ++it) {
		Cell* cell = &(*it);
		if constexpr (Master_matrix::Option_list::has_row_access) ra_opt::unlink(cell);
		cell->set_row_index(valueMap.at(cell->get_row_index()));
		if constexpr (Master_matrix::Option_list::has_intrusive_rows && Master_matrix::Option_list::has_row_access) 
			ra_opt::insert_cell(cell->get_row_index(), cell);
	}

	//all cells have to be deleted first, to avoid problem with insertion when row is a set
	if constexpr (!Master_matrix::Option_list::has_intrusive_rows && Master_matrix::Option_list::has_row_access){
		for (auto it = column_.begin(); it != column_.end(); ++it) {
			Cell* cell = &(*it);
			ra_opt::insert_cell(cell->get_row_index(), cell);
		}
	}

	column_.sort();
}

template<class Master_matrix>
inline void Intrusive_list_column<Master_matrix>::clear()
{
	static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type, 
						"Method not available for chain columns as a base element should not be empty.");

	column_.clear_and_dispose(delete_disposer(this));
}

template<class Master_matrix>
inline void Intrusive_list_column<Master_matrix>::clear(index rowIndex)
{
	static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type, 
						"Method not available for chain columns.");

	auto it = column_.begin();
	while (it != column_.end() && it->get_row_index() != rowIndex) it++;
	if (it != column_.end()) _delete_cell(it);
}

template<class Master_matrix>
inline int Intrusive_list_column<Master_matrix>::get_pivot() const
{
	static_assert(Master_matrix::isNonBasic, "Method not available for base columns.");	//could technically be, but is the notion usefull then?

	if constexpr (Master_matrix::Option_list::is_of_boundary_type){
		if (column_.empty()) return -1;
		return column_.back().get_row_index();
	} else {
		return chain_opt::get_pivot();
	}
}

template<class Master_matrix>
inline typename Intrusive_list_column<Master_matrix>::Field_element_type Intrusive_list_column<Master_matrix>::get_pivot_value() const
{
	static_assert(Master_matrix::isNonBasic, "Method not available for base columns.");	//could technically be, but is the notion usefull then?

	if constexpr (Master_matrix::Option_list::is_z2){
		return 1;
	} else {
		if constexpr (Master_matrix::Option_list::is_of_boundary_type){
			if (column_.empty()) return 0;
			return column_.back().get_element();
		} else {
			if (chain_opt::get_pivot() == -1) return Field_element_type();
			for (const Cell& cell : column_){
				if (cell.get_row_index() == static_cast<unsigned int>(chain_opt::get_pivot())) return cell.get_element();
			}
			return Field_element_type();	//should never happen if chain column is used properly
		}
	}
}

template<class Master_matrix>
inline typename Intrusive_list_column<Master_matrix>::iterator
Intrusive_list_column<Master_matrix>::begin() noexcept
{
	return column_.begin();
}

template<class Master_matrix>
inline typename Intrusive_list_column<Master_matrix>::const_iterator
Intrusive_list_column<Master_matrix>::begin() const noexcept
{
	return column_.begin();
}

template<class Master_matrix>
inline typename Intrusive_list_column<Master_matrix>::iterator
Intrusive_list_column<Master_matrix>::end() noexcept
{
	return column_.end();
}

template<class Master_matrix>
inline typename Intrusive_list_column<Master_matrix>::const_iterator
Intrusive_list_column<Master_matrix>::end() const noexcept
{
	return column_.end();
}

template<class Master_matrix>
inline typename Intrusive_list_column<Master_matrix>::reverse_iterator
Intrusive_list_column<Master_matrix>::rbegin() noexcept
{
	return column_.rbegin();
}

template<class Master_matrix>
inline typename Intrusive_list_column<Master_matrix>::const_reverse_iterator
Intrusive_list_column<Master_matrix>::rbegin() const noexcept
{
	return column_.rbegin();
}

template<class Master_matrix>
inline typename Intrusive_list_column<Master_matrix>::reverse_iterator
Intrusive_list_column<Master_matrix>::rend() noexcept
{
	return column_.rend();
}

template<class Master_matrix>
inline typename Intrusive_list_column<Master_matrix>::const_reverse_iterator
Intrusive_list_column<Master_matrix>::rend() const noexcept
{
	return column_.rend();
}

template<class Master_matrix>
template<class Cell_range>
inline Intrusive_list_column<Master_matrix> &
Intrusive_list_column<Master_matrix>::operator+=(const Cell_range &column)
{
	static_assert((!Master_matrix::isNonBasic || std::is_same_v<Cell_range, Intrusive_list_column>), 
					"For boundary columns, the range has to be a column of same type to help ensure the validity of the base element.");	//could be removed, if we give the responsability to the user.
	static_assert((!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type), 
					"For chain columns, the given column cannot be constant.");

	_add(column);

	return *this;
}

template<class Master_matrix>
inline Intrusive_list_column<Master_matrix> &
Intrusive_list_column<Master_matrix>::operator+=(Intrusive_list_column &column)
{
	if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type){
		//assumes that the addition never zeros out this column. 
		if (_add(column)) chain_opt::swap_pivots(column);
	} else {
		_add(column);
	}

	return *this;
}

template<class Master_matrix>
inline Intrusive_list_column<Master_matrix> &
Intrusive_list_column<Master_matrix>::operator*=(unsigned int v)
{
	if constexpr (Master_matrix::Option_list::is_z2){
		if (v % 2 == 0){
			if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type){
				throw std::invalid_argument("A chain column should not be multiplied by 0.");
			} else {
				clear();
			}
		}
	} else {
	//	v %= Field_element_type::get_characteristic();		//don't work because of multifields...
		Field_element_type val(v);

		if (val == 0u) {
			if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type){
				throw std::invalid_argument("A chain column should not be multiplied by 0.");
			} else {
				clear();
			}
			return *this;
		}

		if (val == 1u) return *this;

		for (Cell& cell : column_){
			cell.get_element() *= val;
			if constexpr (Master_matrix::Option_list::has_row_access)
				ra_opt::update_cell(cell);
		}
	}

	return *this;
}

template<class Master_matrix>
template<class Cell_range>
inline Intrusive_list_column<Master_matrix> &
Intrusive_list_column<Master_matrix>::multiply_and_add(const Field_element_type& val, const Cell_range& column)
{
	static_assert((!Master_matrix::isNonBasic || std::is_same_v<Cell_range, Intrusive_list_column>), 
					"For boundary columns, the range has to be a column of same type to help ensure the validity of the base element.");	//could be removed, if we give the responsability to the user.
	static_assert((!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type), 
					"For chain columns, the given column cannot be constant.");

	if constexpr (Master_matrix::Option_list::is_z2){
		if (val){
			_add(column);
		} else {
			clear();
			_add(column);
		}
	} else {
		_multiply_and_add(val, column);
	}

	return *this;
}

template<class Master_matrix>
inline Intrusive_list_column<Master_matrix> &
Intrusive_list_column<Master_matrix>::multiply_and_add(const Field_element_type& val, Intrusive_list_column& column)
{
	if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type){
		//assumes that the addition never zeros out this column. 
		if constexpr (Master_matrix::Option_list::is_z2){
			if (val){
				if (_add(column)) chain_opt::swap_pivots(column);
			} else {
				throw std::invalid_argument("A chain column should not be multiplied by 0.");
			}
		} else {
			if (_multiply_and_add(val, column)) chain_opt::swap_pivots(column);
		}
	} else {
		if constexpr (Master_matrix::Option_list::is_z2){
			if (val){
				_add(column);
			} else {
				clear();
				_add(column);
			}
		} else {
			_multiply_and_add(val, column);
		}
	}

	return *this;
}

template<class Master_matrix>
template<class Cell_range>
inline Intrusive_list_column<Master_matrix> &
Intrusive_list_column<Master_matrix>::multiply_and_add(const Cell_range& column, const Field_element_type& val)
{
	static_assert((!Master_matrix::isNonBasic || std::is_same_v<Cell_range, Intrusive_list_column>), 
					"For boundary columns, the range has to be a column of same type to help ensure the validity of the base element.");	//could be removed, if we give the responsability to the user.
	static_assert((!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type), 
					"For chain columns, the given column cannot be constant.");

	if constexpr (Master_matrix::Option_list::is_z2){
		if (val){
			_add(column);
		}
	} else {
		_multiply_and_add(column, val);
	}

	return *this;
}

template<class Master_matrix>
inline Intrusive_list_column<Master_matrix> &
Intrusive_list_column<Master_matrix>::multiply_and_add(Intrusive_list_column& column, const Field_element_type& val)
{
	if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type){
		//assumes that the addition never zeros out this column. 
		if constexpr (Master_matrix::Option_list::is_z2){
			if (val){
				if (_add(column)) chain_opt::swap_pivots(column);
			}
		} else {
			if (_multiply_and_add(column, val)) chain_opt::swap_pivots(column);
		}
	} else {
		if constexpr (Master_matrix::Option_list::is_z2){
			if (val){
				_add(column);
			}
		} else {
			_multiply_and_add(column, val);
		}
	}

	return *this;
}

template<class Master_matrix>
inline Intrusive_list_column<Master_matrix> &
Intrusive_list_column<Master_matrix>::operator=(const Intrusive_list_column& other)
{
	static_assert(!Master_matrix::Option_list::has_row_access, "= assignement not enabled with row access option.");

	dim_opt::operator=(other);
	chain_opt::operator=(other);

	column_.clear_and_dispose(delete_disposer(this));
	column_.clone_from(other.column_, new_cloner(), delete_disposer(this));
	return *this;
}

template<class Master_matrix>
inline void Intrusive_list_column<Master_matrix>::_delete_cell(iterator &it)
{
	it = column_.erase_and_dispose(it, delete_disposer(this));
}

template<class Master_matrix>
inline void Intrusive_list_column<Master_matrix>::_insert_cell(
		const Field_element_type &value, index rowIndex, const iterator &position)
{
	if constexpr (Master_matrix::Option_list::has_row_access){
		Cell *new_cell = cellPool_.construct(value, ra_opt::columnIndex_, rowIndex);
		column_.insert(position, *new_cell);
		ra_opt::insert_cell(rowIndex, new_cell);
	} else {
		Cell *new_cell = cellPool_.construct(value, rowIndex);
		column_.insert(position, *new_cell);
	}
}

template<class Master_matrix>
inline void Intrusive_list_column<Master_matrix>::_insert_cell(
		index rowIndex, const iterator &position)
{
	if constexpr (Master_matrix::Option_list::has_row_access){
		Cell *new_cell = cellPool_.construct(ra_opt::columnIndex_, rowIndex);
		column_.insert(position, *new_cell);
		ra_opt::insert_cell(rowIndex, new_cell);
	} else {
		Cell *new_cell = cellPool_.construct(rowIndex);
		column_.insert(position, *new_cell);
	}
}

template<class Master_matrix>
template<class Cell_range>
inline bool Intrusive_list_column<Master_matrix>::_add(const Cell_range &column)
{
	auto itTarget = column_.begin();
	auto itSource = column.begin();
	bool pivotIsZeroed = false;

	while (itTarget != column_.end() && itSource != column.end())
	{
		if (itTarget->get_row_index() < itSource->get_row_index()) {
			++itTarget;
		} else if (itTarget->get_row_index() > itSource->get_row_index()) {
			if constexpr (Master_matrix::Option_list::is_z2){
				_insert_cell(itSource->get_row_index(), itTarget);
			} else {
				_insert_cell(itSource->get_element(), itSource->get_row_index(), itTarget);
			}
			++itSource;
		} else {
			if constexpr (Master_matrix::Option_list::is_z2){
				if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type){
					if (static_cast<int>(itTarget->get_row_index()) == chain_opt::get_pivot()) pivotIsZeroed = true;
				}
				_delete_cell(itTarget);
			} else {
				itTarget->get_element() += itSource->get_element();
				if (itTarget->get_element() == Field_element_type::get_additive_identity()){
					if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type){
						if (static_cast<int>(itTarget->get_row_index()) == chain_opt::get_pivot()) pivotIsZeroed = true;
					}
					_delete_cell(itTarget);
				} else {
					if constexpr (Master_matrix::Option_list::has_row_access)
						ra_opt::update_cell(*itTarget);
					++itTarget;
				}
			}
			++itSource;
		}
	}

	while (itSource != column.end()) {
		if constexpr (Master_matrix::Option_list::is_z2){
			_insert_cell(itSource->get_row_index(), column_.end());
		} else {
			_insert_cell(itSource->get_element(), itSource->get_row_index(), column_.end());
		}
		++itSource;
	}

	return pivotIsZeroed;
}

template<class Master_matrix>
template<class Cell_range>
inline bool Intrusive_list_column<Master_matrix>::_multiply_and_add(const Field_element_type& val, const Cell_range& column)
{
	bool pivotIsZeroed = false;

	if (val == 0u) {
		if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type){
			throw std::invalid_argument("A chain column should not be multiplied by 0.");
			//this would not only mess up the base, but also the pivots stored.
		} else {
			clear();
		}
	}

	auto itTarget = column_.begin();
	auto itSource = column.begin();
	while (itTarget != column_.end() && itSource != column.end())
	{
		if (itTarget->get_row_index() < itSource->get_row_index()) {
			itTarget->get_element() *= val;
			if constexpr (Master_matrix::Option_list::has_row_access)
				ra_opt::update_cell(*itTarget);
			++itTarget;
		} else if (itTarget->get_row_index() > itSource->get_row_index()) {
			_insert_cell(itSource->get_element(), itSource->get_row_index(), itTarget);
			++itSource;
		} else {
			itTarget->get_element() *= val;
			itTarget->get_element() += itSource->get_element();
			if (itTarget->get_element() == Field_element_type::get_additive_identity()){
				if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type){
					if (static_cast<int>(itTarget->get_row_index()) == chain_opt::get_pivot()) pivotIsZeroed = true;
				}
				_delete_cell(itTarget);
			} else {
				if constexpr (Master_matrix::Option_list::has_row_access)
					ra_opt::update_cell(*itTarget);
				++itTarget;
			}
			++itSource;
		}
	}

	while (itTarget != column_.end()){
		itTarget->get_element() *= val;
		if constexpr (Master_matrix::Option_list::has_row_access)
			ra_opt::update_cell(*itTarget);
		itTarget++;
	}

	while (itSource != column.end()) {
		_insert_cell(itSource->get_element(), itSource->get_row_index(), column_.end());
		++itSource;
	}

	return pivotIsZeroed;
}

template<class Master_matrix>
template<class Cell_range>
inline bool Intrusive_list_column<Master_matrix>::_multiply_and_add(const Cell_range& column, const Field_element_type& val)
{
	if (val == 0u) {
		return false;
	}

	bool pivotIsZeroed = false;

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
				if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type){
					if (static_cast<int>(itTarget->get_row_index()) == chain_opt::get_pivot()) pivotIsZeroed = true;
				}
				_delete_cell(itTarget);
			} else {
				if constexpr (Master_matrix::Option_list::has_row_access)
						ra_opt::update_cell(*itTarget);
				++itTarget;
			}
			++itSource;
		}
	}

	while (itSource != column.end()) {
		_insert_cell(itSource->get_element() * val, itSource->get_row_index(), column_.end());
		++itSource;
	}

	return pivotIsZeroed;
}

} //namespace persistence_matrix
} //namespace Gudhi

template<class Master_matrix>
struct std::hash<Gudhi::persistence_matrix::Intrusive_list_column<Master_matrix> >
{
	size_t operator()(const Gudhi::persistence_matrix::Intrusive_list_column<Master_matrix>& column) const
	{
		std::size_t seed = 0;
		for (auto& cell : column){
			seed ^= std::hash<unsigned int>()(cell.get_row_index() * static_cast<unsigned int>(cell.get_element())) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		}
		return seed;
	}
};

#endif // PM_INTRUSIVE_LIST_COLUMN_H

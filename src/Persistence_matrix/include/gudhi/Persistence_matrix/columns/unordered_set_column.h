/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022-23 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_UNORDERED_SET_COLUMN_H
#define PM_UNORDERED_SET_COLUMN_H

#include <vector>
#include <stdexcept>
#include <type_traits>
#include <unordered_set>
#include <set>
#include <utility>	//std::swap, std::move & std::exchange

#include <boost/iterator/indirect_iterator.hpp>

#include <gudhi/Simple_object_pool.h>

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class Unordered_set_column : public Master_matrix::Row_access_option, 
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

	struct CellPointerHash {
		size_t operator()(const Cell* c) const
		{
			return std::hash<Cell>()(*c);
		}
	};
	struct CellPointerEq {
		bool operator()(const Cell* c1, const Cell* c2) const
		{
			return *c1 == *c2;
		}
	};
	struct CellPointerComp {
		bool operator()(const Cell* c1, const Cell* c2) const
		{
			return *c1 < *c2;
		}
	};

	using Column_type = std::unordered_set<Cell*,CellPointerHash,CellPointerEq>;
	using iterator = boost::indirect_iterator<typename Column_type::iterator>;
	using const_iterator = boost::indirect_iterator<typename Column_type::const_iterator>;

	Unordered_set_column();
	template<class Container_type = typename Master_matrix::boundary_type>
	Unordered_set_column(const Container_type& nonZeroRowIndices);	//has to be a boundary for boundary, has no sense for chain if dimension is needed
	template<class Container_type = typename Master_matrix::boundary_type, class Row_container_type>
	Unordered_set_column(index columnIndex, const Container_type& nonZeroRowIndices, Row_container_type &rowContainer);	//has to be a boundary for boundary, has no sense for chain if dimension is needed
	template<class Container_type = typename Master_matrix::boundary_type>
	Unordered_set_column(const Container_type& nonZeroChainRowIndices, dimension_type dimension);	//dimension gets ignored for base
	template<class Container_type = typename Master_matrix::boundary_type, class Row_container_type>
	Unordered_set_column(index columnIndex, const Container_type& nonZeroChainRowIndices, dimension_type dimension, Row_container_type &rowContainer);	//dimension gets ignored for base
	Unordered_set_column(const Unordered_set_column& column);
	template<class Row_container_type>
	Unordered_set_column(const Unordered_set_column& column, index columnIndex, Row_container_type &rowContainer);
	Unordered_set_column(Unordered_set_column&& column) noexcept;
	~Unordered_set_column();

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

	template<class Cell_range>
	Unordered_set_column& operator+=(const Cell_range& column);	//for base & boundary except vector
	Unordered_set_column& operator+=(Unordered_set_column &column);	//for chain and vector
	friend Unordered_set_column operator+(Unordered_set_column column1, Unordered_set_column& column2){
		column1 += column2;
		return column1;
	}

	Unordered_set_column& operator*=(unsigned int v);
	friend Unordered_set_column operator*(Unordered_set_column column, unsigned int const& v){
		column *= v;
		return column;
	}
	friend Unordered_set_column operator*(unsigned int const& v, Unordered_set_column column){
		column *= v;
		return column;
	}

	//this = v * this + column
	template<class Cell_range>
	Unordered_set_column& multiply_and_add(const Field_element_type& val, const Cell_range& column);	//for base & boundary except vector
	Unordered_set_column& multiply_and_add(const Field_element_type& val, Unordered_set_column& column);	//for chain and vector
	//this = this + column * v
	template<class Cell_range>
	Unordered_set_column& multiply_and_add(const Cell_range& column, const Field_element_type& val);	//for base & boundary except vector
	Unordered_set_column& multiply_and_add(Unordered_set_column& column, const Field_element_type& val);	//for chain and vector

	friend bool operator==(const Unordered_set_column& c1, const Unordered_set_column& c2){
		if (&c1 == &c2) return true;
		if (c1.column_.size() != c2.column_.size()) return false;

		using rep_type = typename std::conditional<
								Master_matrix::Option_list::is_z2,
								unsigned int,
								std::pair<unsigned int,unsigned int>
							>::type;

		auto it1 = c1.column_.begin();
		auto it2 = c2.column_.begin();
		std::set<rep_type> cells1, cells2;
		while (it1 != c1.column_.end()) {
			if constexpr (Master_matrix::Option_list::is_z2){
				cells1.insert((*it1)->get_row_index());
				cells2.insert((*it2)->get_row_index());
			} else {
				cells1.emplace((*it1)->get_row_index(), (*it1)->get_element());
				cells2.emplace((*it2)->get_row_index(), (*it2)->get_element());
			}
			++it1; ++it2;
		}
		return cells1 == cells2;
	}
	friend bool operator<(const Unordered_set_column& c1, const Unordered_set_column& c2){
		if (&c1 == &c2) return false;

		using rep_type = typename std::conditional<
								Master_matrix::Option_list::is_z2,
								unsigned int,
								std::pair<unsigned int,unsigned int>
							>::type;

		auto it1 = c1.column_.begin();
		auto it2 = c2.column_.begin();
		std::set<rep_type> cells1, cells2;
		while (it1 != c1.column_.end() && it2 != c2.column_.end()) {
			if constexpr (Master_matrix::Option_list::is_z2){
				cells1.insert((*it1)->get_row_index());
				cells2.insert((*it2)->get_row_index());
			} else {
				cells1.emplace((*it1)->get_row_index(), (*it1)->get_element());
				cells2.emplace((*it2)->get_row_index(), (*it2)->get_element());
			}
			++it1; ++it2;
		}
		while (it1 != c1.column_.end()){
			if constexpr (Master_matrix::Option_list::is_z2){
				cells1.insert((*it1)->get_row_index());
			} else {
				cells1.emplace((*it1)->get_row_index(), (*it1)->get_element());
			}
			++it1;
		}
		while (it2 != c2.column_.end()){
			if constexpr (Master_matrix::Option_list::is_z2){
				cells2.insert((*it2)->get_row_index());
			} else {
				cells2.emplace((*it2)->get_row_index(), (*it2)->get_element());
			}
			++it2;
		}
		return cells1 < cells2;
	}

	//Disabled with row access.
	Unordered_set_column& operator=(const Unordered_set_column& other);

	friend void swap(Unordered_set_column& col1, Unordered_set_column& col2){
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

	void _delete_cell(typename Column_type::iterator& it);
	void _insert_cell(const Field_element_type& value, index rowIndex, const typename Column_type::iterator& position);
	void _insert_cell(index rowIndex, const typename Column_type::iterator& position);
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
};

template<class Master_matrix>
inline Unordered_set_column<Master_matrix>::Unordered_set_column() : ra_opt(), dim_opt(), chain_opt()
{}

template<class Master_matrix>
template<class Container_type>
inline Unordered_set_column<Master_matrix>::Unordered_set_column(const Container_type &nonZeroRowIndices)
	: ra_opt(), 
	  dim_opt(nonZeroRowIndices.size() == 0 ? 0 : nonZeroRowIndices.size() - 1), 
	  chain_opt(), 
	  column_(nonZeroRowIndices.size())
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
inline Unordered_set_column<Master_matrix>::Unordered_set_column(
	index columnIndex, const Container_type &nonZeroRowIndices, Row_container_type &rowContainer) 
	: ra_opt(columnIndex, rowContainer), 
	  dim_opt(nonZeroRowIndices.size() == 0 ? 0 : nonZeroRowIndices.size() - 1),
	  chain_opt([&]{
			if constexpr (Master_matrix::Option_list::is_z2){
				return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : *std::prev(nonZeroRowIndices.end());
			} else {
				return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : std::prev(nonZeroRowIndices.end())->first;
			}
		}()), 
	  column_(nonZeroRowIndices.size())
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
inline Unordered_set_column<Master_matrix>::Unordered_set_column(
	const Container_type &nonZeroRowIndices, dimension_type dimension) 
	: ra_opt(), 
	  dim_opt(dimension),
	  chain_opt([&]{
			if constexpr (Master_matrix::Option_list::is_z2){
				return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : *std::prev(nonZeroRowIndices.end());
			} else {
				return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : std::prev(nonZeroRowIndices.end())->first;
			}
		}()), 
	  column_(nonZeroRowIndices.size())
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
inline Unordered_set_column<Master_matrix>::Unordered_set_column(
	index columnIndex, const Container_type &nonZeroRowIndices, dimension_type dimension, Row_container_type &rowContainer) 
	: ra_opt(columnIndex, rowContainer), 
	  dim_opt(dimension),
	  chain_opt([&]{
			if constexpr (Master_matrix::Option_list::is_z2){
				return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : *std::prev(nonZeroRowIndices.end());
			} else {
				return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : std::prev(nonZeroRowIndices.end())->first;
			}
		}()), 
	  column_(nonZeroRowIndices.size())
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
inline Unordered_set_column<Master_matrix>::Unordered_set_column(const Unordered_set_column &column) 
	: ra_opt(), 
	  dim_opt(static_cast<const dim_opt&>(column)), 
	  chain_opt(static_cast<const chain_opt&>(column)),
	  column_(column.column_.bucket_count())
{
	static_assert(!Master_matrix::Option_list::has_row_access,
			"Simple copy constructor not available when row access option enabled. Please specify the new column index and the row container.");

	for (const Cell* cell : column.column_){
		if constexpr (Master_matrix::Option_list::is_z2){
			_insert_cell(cell->get_row_index(), column_.end());
		} else {
			_insert_cell(cell->get_element(), cell->get_row_index(), column_.end());
		}
	}
}

template<class Master_matrix>
template<class Row_container_type>
inline Unordered_set_column<Master_matrix>::Unordered_set_column(
	const Unordered_set_column &column, index columnIndex, Row_container_type &rowContainer) 
	: ra_opt(columnIndex, rowContainer), 
	  dim_opt(static_cast<const dim_opt&>(column)), 
	  chain_opt(static_cast<const chain_opt&>(column)),
	  column_(column.column_.bucket_count())
{
	for (const Cell* cell : column.column_){
		if constexpr (Master_matrix::Option_list::is_z2){
			_insert_cell(cell->get_row_index(), column_.end());
		} else {
			_insert_cell(cell->get_element(), cell->get_row_index(), column_.end());
		}
	}
}

template<class Master_matrix>
inline Unordered_set_column<Master_matrix>::Unordered_set_column(Unordered_set_column &&column) noexcept 
	: ra_opt(std::move(static_cast<ra_opt&>(column))), 
	  dim_opt(std::move(static_cast<dim_opt&>(column))), 
	  chain_opt(std::move(static_cast<chain_opt&>(column))), 
	  column_(std::move(column.column_))
{}

template<class Master_matrix>
inline Unordered_set_column<Master_matrix>::~Unordered_set_column()
{
	for (auto* cell : column_){
		if constexpr (Master_matrix::Option_list::has_row_access) ra_opt::unlink(cell);
		cellPool_.destroy(cell);
	}
}

template<class Master_matrix>
inline std::vector<typename Unordered_set_column<Master_matrix>::Field_element_type> 
Unordered_set_column<Master_matrix>::get_content(int columnLength) const
{
	if (columnLength < 0 && column_.size() > 0) columnLength = (*std::max_element(column_.begin(), column_.end(), CellPointerComp()))->get_row_index() + 1;
	else if (columnLength < 0) return std::vector<Field_element_type>();

	std::vector<Field_element_type> container(columnLength, 0);
	for (auto it = column_.begin(); it != column_.end(); ++it){
		if ((*it)->get_row_index() < static_cast<index>(columnLength)){
			if constexpr (Master_matrix::Option_list::is_z2){
				container[(*it)->get_row_index()] = 1;
			} else {
				container[(*it)->get_row_index()] = (*it)->get_element();
			}
		}
	}
	return container;
}

template<class Master_matrix>
inline bool Unordered_set_column<Master_matrix>::is_non_zero(index rowIndex) const
{
	return column_.find(cellPool_.construct(rowIndex)) != column_.end();	//cell gets destroyed with the pool at the end, but I don't know if that's a good solution
}

template<class Master_matrix>
inline bool Unordered_set_column<Master_matrix>::is_empty() const
{
	return column_.empty();
}

template<class Master_matrix>
inline std::size_t Unordered_set_column<Master_matrix>::size() const{
	return column_.size();
}

template<class Master_matrix>
template<class Map_type>
inline void Unordered_set_column<Master_matrix>::reorder(const Map_type &valueMap)
{
	static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type, 
						"Method not available for chain columns.");

	Column_type newSet;

	for (Cell* cell : column_) {
		if constexpr (Master_matrix::Option_list::has_row_access) ra_opt::unlink(cell);
		cell->set_row_index(valueMap.at(cell->get_row_index()));
		newSet.insert(cell);
		if constexpr (Master_matrix::Option_list::has_row_access && Master_matrix::Option_list::has_intrusive_rows)	//intrusive list
			ra_opt::insert_cell(cell->get_row_index(), cell);
	}

	//when row is a set, all cells have to be deleted first, to avoid colliding when inserting
	if constexpr (Master_matrix::Option_list::has_row_access && !Master_matrix::Option_list::has_intrusive_rows){	//set
		for (Cell* cell : newSet) {
			ra_opt::insert_cell(cell->get_row_index(), cell);
		}
	}

	column_.swap(newSet);
}

template<class Master_matrix>
inline void Unordered_set_column<Master_matrix>::clear()
{
	static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type, 
						"Method not available for chain columns as a base element should not be empty.");

	for (auto* cell : column_){
		if constexpr (Master_matrix::Option_list::has_row_access) ra_opt::unlink(cell);
		cellPool_.destroy(cell);
	}

	column_.clear();
}

template<class Master_matrix>
inline void Unordered_set_column<Master_matrix>::clear(index rowIndex)
{
	static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type, 
						"Method not available for chain columns.");

	auto it = column_.find(cellPool_.construct(rowIndex));	//cell gets destroyed with the pool at the end, but I don't know if that's a good solution
	if (it != column_.end()){
		_delete_cell(it);
	}
}

template<class Master_matrix>
inline int Unordered_set_column<Master_matrix>::get_pivot() const
{
	static_assert(Master_matrix::isNonBasic, "Method not available for base columns.");	//could technically be, but is the notion usefull then?

	if constexpr (Master_matrix::Option_list::is_of_boundary_type){
		if (column_.empty()) return -1;
		//linear search could be avoided with storing the pivot. But even then, some modifications of the column requires the max, so not clear how much it is worth it.
		return (*std::max_element(column_.begin(), column_.end(), CellPointerComp()))->get_row_index();
	} else {
		return chain_opt::get_pivot();
	}
}

template<class Master_matrix>
inline typename Unordered_set_column<Master_matrix>::Field_element_type Unordered_set_column<Master_matrix>::get_pivot_value() const
{
	static_assert(Master_matrix::isNonBasic, "Method not available for base columns.");	//could technically be, but is the notion usefull then?

	if constexpr (Master_matrix::Option_list::is_z2){
		return 1;
	} else {
		if constexpr (Master_matrix::Option_list::is_of_boundary_type){
			if (column_.empty()) return 0;
			return (*std::max_element(column_.begin(), column_.end(), CellPointerComp()))->get_element();
		} else {
			if (chain_opt::get_pivot() == -1) return Field_element_type();
			for (const Cell* cell : column_){
				if (cell->get_row_index() == static_cast<unsigned int>(chain_opt::get_pivot())) return cell->get_element();
			}
			return Field_element_type();	//should never happen if chain column is used properly
		}
	}
}

template<class Master_matrix>
inline typename Unordered_set_column<Master_matrix>::iterator
Unordered_set_column<Master_matrix>::begin() noexcept
{
	return column_.begin();
}

template<class Master_matrix>
inline typename Unordered_set_column<Master_matrix>::const_iterator
Unordered_set_column<Master_matrix>::begin() const noexcept
{
	return column_.begin();
}

template<class Master_matrix>
inline typename Unordered_set_column<Master_matrix>::iterator
Unordered_set_column<Master_matrix>::end() noexcept
{
	return column_.end();
}

template<class Master_matrix>
inline typename Unordered_set_column<Master_matrix>::const_iterator
Unordered_set_column<Master_matrix>::end() const noexcept
{
	return column_.end();
}

template<class Master_matrix>
template<class Cell_range>
inline Unordered_set_column<Master_matrix> &
Unordered_set_column<Master_matrix>::operator+=(const Cell_range &column)
{
	static_assert((!Master_matrix::isNonBasic || std::is_same_v<Cell_range, Unordered_set_column>), 
					"For boundary columns, the range has to be a column of same type to help ensure the validity of the base element.");	//could be removed, if we give the responsability to the user.
	static_assert((!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type), 
					"For chain columns, the given column cannot be constant.");

	_add(column);

	return *this;
}

template<class Master_matrix>
inline Unordered_set_column<Master_matrix> &
Unordered_set_column<Master_matrix>::operator+=(Unordered_set_column &column)
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
inline Unordered_set_column<Master_matrix> &
Unordered_set_column<Master_matrix>::operator*=(unsigned int v)
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

		for (Cell* cell : column_){
			cell->get_element() *= val;
			if constexpr (Master_matrix::Option_list::has_row_access)
				ra_opt::update_cell(*cell);
		}
	}

	return *this;
}

template<class Master_matrix>
template<class Cell_range>
inline Unordered_set_column<Master_matrix> &
Unordered_set_column<Master_matrix>::multiply_and_add(const Field_element_type& val, const Cell_range& column)
{
	static_assert((!Master_matrix::isNonBasic || std::is_same_v<Cell_range, Unordered_set_column>), 
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
inline Unordered_set_column<Master_matrix> &
Unordered_set_column<Master_matrix>::multiply_and_add(const Field_element_type& val, Unordered_set_column& column)
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
inline Unordered_set_column<Master_matrix> &
Unordered_set_column<Master_matrix>::multiply_and_add(const Cell_range& column, const Field_element_type& val)
{
	static_assert((!Master_matrix::isNonBasic || std::is_same_v<Cell_range, Unordered_set_column>), 
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
inline Unordered_set_column<Master_matrix> &
Unordered_set_column<Master_matrix>::multiply_and_add(Unordered_set_column& column, const Field_element_type& val)
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
inline Unordered_set_column<Master_matrix> &
Unordered_set_column<Master_matrix>::operator=(const Unordered_set_column& other)
{
	static_assert(!Master_matrix::Option_list::has_row_access, "= assignement not enabled with row access option.");

	dim_opt::operator=(other);
	chain_opt::operator=(other);
	
	for (auto* cell : column_){
		if constexpr (Master_matrix::Option_list::has_row_access) ra_opt::unlink(cell);
		cellPool_.destroy(cell);
	}
	column_.clear();
	
	for (const Cell* cell : other.column_){
		if constexpr (Master_matrix::Option_list::is_z2){
			_insert_cell(cell->get_row_index(), column_.end());
		} else {
			_insert_cell(cell->get_element(), cell->get_row_index(), column_.end());
		}
	}
	
	return *this;
}

template<class Master_matrix>
inline void Unordered_set_column<Master_matrix>::_delete_cell(typename Column_type::iterator &it)
{
	if constexpr (Master_matrix::Option_list::has_row_access) ra_opt::unlink(*it);
	cellPool_.destroy(*it);
	it = column_.erase(it);
}

template<class Master_matrix>
inline void Unordered_set_column<Master_matrix>::_insert_cell(
		const Field_element_type &value, index rowIndex, const typename Column_type::iterator &position)
{
	if constexpr (Master_matrix::Option_list::has_row_access){
		Cell *new_cell = cellPool_.construct(value, ra_opt::columnIndex_, rowIndex);
		column_.insert(position, new_cell);
		ra_opt::insert_cell(rowIndex, new_cell);
	} else {
		Cell *new_cell = cellPool_.construct(value, rowIndex);
		column_.insert(position, new_cell);
	}
}

template<class Master_matrix>
inline void Unordered_set_column<Master_matrix>::_insert_cell(
		index rowIndex, const typename Column_type::iterator &position)
{
	if constexpr (Master_matrix::Option_list::has_row_access){
		Cell *new_cell = cellPool_.construct(ra_opt::columnIndex_, rowIndex);
		column_.insert(position, new_cell);
		ra_opt::insert_cell(rowIndex, new_cell);
	} else {
		Cell *new_cell = cellPool_.construct(rowIndex);
		column_.insert(position, new_cell);
	}
}

template<class Master_matrix>
template<class Cell_range>
inline bool Unordered_set_column<Master_matrix>::_add(const Cell_range &column)
{
	bool pivotIsZeroed = false;

	for (const Cell &cell : column) {
		auto it1 = column_.find(const_cast<Cell*>(&cell));
		if (it1 != column_.end()) {
			if constexpr (Master_matrix::Option_list::is_z2){
				if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type){
					if (static_cast<int>((*it1)->get_row_index()) == chain_opt::get_pivot()) pivotIsZeroed = true;
				}
				_delete_cell(it1);
			} else {
				(*it1)->get_element() += cell.get_element();
				if ((*it1)->get_element() == Field_element_type::get_additive_identity()){
					if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type){
						if (static_cast<int>((*it1)->get_row_index()) == chain_opt::get_pivot()) pivotIsZeroed = true;
					}
					_delete_cell(it1);
				} else {
					if constexpr (Master_matrix::Option_list::has_row_access)
						ra_opt::update_cell(**it1);
				}
			}
		} else {
			if constexpr (Master_matrix::Option_list::is_z2){
				_insert_cell(cell.get_row_index(), column_.end());
			} else {
				_insert_cell(cell.get_element(), cell.get_row_index(), column_.end());
			}
		}
	}

	return pivotIsZeroed;
}

template<class Master_matrix>
template<class Cell_range>
inline bool Unordered_set_column<Master_matrix>::_multiply_and_add(const Field_element_type& val, const Cell_range& column)
{
	if (val == 0u) {
		if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type){
			throw std::invalid_argument("A chain column should not be multiplied by 0.");
			//this would not only mess up the base, but also the pivots stored.
		} else {
			clear();
			for (const Cell& v : column){
				_insert_cell(v.get_element(), v.get_row_index(), column_.end());
			}
			return true;
		}
	}

	//because the column is unordered, I don't see a way to do both operations in one go...
	operator*=(val);
	return _add(column);
}

template<class Master_matrix>
template<class Cell_range>
inline bool Unordered_set_column<Master_matrix>::_multiply_and_add(const Cell_range& column, const Field_element_type& val)
{
	if (val == 0u) {
		return false;
	}

	bool pivotIsZeroed = false;

	for (const Cell &cell : column) {
		auto it1 = column_.find(const_cast<Cell*>(&cell));
		if (it1 != column_.end()) {
			(*it1)->get_element() += (cell.get_element() * val);
			if ((*it1)->get_element() == Field_element_type::get_additive_identity()){
				if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type){
					if (static_cast<int>((*it1)->get_row_index()) == chain_opt::get_pivot()) pivotIsZeroed = true;
				}
				_delete_cell(it1);
			} else {
				if constexpr (Master_matrix::Option_list::has_row_access)
					ra_opt::update_cell(**it1);
			}
		} else {
			_insert_cell(cell.get_element() * val, cell.get_row_index(), column_.end());
		}
	}

	return pivotIsZeroed;
}

} //namespace persistence_matrix
} //namespace Gudhi

template<class Master_matrix>
struct std::hash<Gudhi::persistence_matrix::Unordered_set_column<Master_matrix> >
{
	size_t operator()(const Gudhi::persistence_matrix::Unordered_set_column<Master_matrix>& column) const
	{
		std::size_t seed = 0;
		for (const auto& cell : column){
			seed ^= std::hash<unsigned int>()(cell.get_row_index() * static_cast<unsigned int>(cell.get_element())) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		}
		return seed;
	}
};

#endif // PM_UNORDERED_SET_COLUMN_H

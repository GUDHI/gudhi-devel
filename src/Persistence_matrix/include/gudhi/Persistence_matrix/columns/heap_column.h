/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022-23 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_HEAP_COLUMN_H
#define PM_HEAP_COLUMN_H

#include <vector>
#include <stdexcept>
#include <type_traits>
#include <algorithm>	//binary_search
#include <utility>	//std::swap, std::move & std::exchange

#include <boost/iterator/indirect_iterator.hpp>

#include <gudhi/Simple_object_pool.h>

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class Heap_column : public Master_matrix::Column_dimension_option,
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
	using Column_type = std::vector<Cell*>;
	using iterator = boost::indirect_iterator<typename Column_type::iterator>;
	using const_iterator = boost::indirect_iterator<typename Column_type::const_iterator>;
	using reverse_iterator = boost::indirect_iterator<typename Column_type::reverse_iterator>;
	using const_reverse_iterator = boost::indirect_iterator<typename Column_type::const_reverse_iterator>;

	Heap_column();
	template<class Container_type = typename Master_matrix::boundary_type>
	Heap_column(const Container_type& nonZeroRowIndices);	//has to be a boundary for boundary, has no sense for chain if dimension is needed
	template<class Container_type = typename Master_matrix::boundary_type>
	Heap_column(const Container_type& nonZeroChainRowIndices, dimension_type dimension);	//dimension gets ignored for base
	Heap_column(const Heap_column& column);
	Heap_column(Heap_column&& column) noexcept;
	~Heap_column();

	//just for the sake of the interface
	//row containers and column index are ignored as row access is not implemented for heap columns
	template<class Container_type = typename Master_matrix::boundary_type, class Row_container_type>
	Heap_column(index columnIndex, const Container_type& nonZeroRowIndices, Row_container_type &rowContainer);
	template<class Container_type = typename Master_matrix::boundary_type, class Row_container_type>
	Heap_column(index columnIndex, const Container_type& nonZeroChainRowIndices, dimension_type dimension, Row_container_type &rowContainer);
	template<class Row_container_type>
	Heap_column(const Heap_column& column, index columnIndex, Row_container_type &rowContainer);

	std::vector<Field_element_type> get_content(int columnLength = -1) const;
	bool is_non_zero(index rowIndex) const;
	bool is_empty();
	std::size_t size() const;	//size of the underlying vector, not the number of non zero row cells.

	//****************
	//only for base and boundary
	template<class Map_type>
	void reorder(const Map_type& valueMap);	//used for lazy row swaps
	void clear();
	void clear(index rowIndex);
	//****************

	//****************
	//only for chain and boundary
	int get_pivot();
	Field_element_type get_pivot_value();
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
	Heap_column& operator+=(const Cell_range& column);	//for base & boundary except vector	//may not work if Cell_range = Heap_column<Other>
	Heap_column& operator+=(Heap_column &column);	//for chain and vector
	friend Heap_column operator+(Heap_column column1, Heap_column& column2){
		column1 += column2;
		return column1;
	}

	Heap_column& operator*=(unsigned int v);
	friend Heap_column operator*(Heap_column column, unsigned int const& v){
		column *= v;
		return column;
	}
	friend Heap_column operator*(unsigned int const& v, Heap_column column){
		column *= v;
		return column;
	}

	//this = v * this + column
	template<class Cell_range>
	Heap_column& multiply_and_add(const Field_element_type& val, const Cell_range& column);	//for base & boundary except vector	//may not work if Cell_range = Heap_column<Other>
	Heap_column& multiply_and_add(const Field_element_type& val, Heap_column& column);	//for chain and vector
	//this = this + column * v
	template<class Cell_range>
	Heap_column& multiply_and_add(const Cell_range& column, const Field_element_type& val);	//for base & boundary except vector	//may not work if Cell_range = Heap_column<Other>
	Heap_column& multiply_and_add(Heap_column& column, const Field_element_type& val);	//for chain and vector

	friend bool operator==(const Heap_column& c1, const Heap_column& c2){
		if (&c1 == &c2) return true;
		return c1.get_content() == c2.get_content();
	}
	friend bool operator<(const Heap_column& c1, const Heap_column& c2){
		if (&c1 == &c2) return false;

		auto cont1 = c1.get_content();
		auto cont2 = c2.get_content();

		auto it1 = cont1.begin();
		auto it2 = cont2.begin();
		while (it1 != cont1.end() && it2 != cont2.end()) {
			if (*it1 == 0u && *it2 != 0u) return false;
			if (*it1 != 0u && *it2 == 0u) return true;
			if (*it1 != *it2) return *it1 < *it2;
			++it1; ++it2;
		}
		return it2 != cont2.end();
	}

	//Disabled with row access.
	Heap_column& operator=(const Heap_column& other);

	friend void swap(Heap_column& col1, Heap_column& col2){
		swap(static_cast<typename Master_matrix::Column_dimension_option&>(col1),
			 static_cast<typename Master_matrix::Column_dimension_option&>(col2));
		swap(static_cast<typename Master_matrix::Chain_column_option&>(col1),
			 static_cast<typename Master_matrix::Chain_column_option&>(col2));
		col1.column_.swap(col2.column_);
	}

protected:
	struct{
		bool operator()(const Cell* c1, const Cell* c2) const { return *c1 < *c2; }
	} cellPointerComp_;

	Column_type column_;
	unsigned int insertsSinceLastPrune_;
	inline static Simple_object_pool<Cell> cellPool_;

	void _prune();
	Cell* _pop_pivot();
	template<class Cell_range>
	bool _add(const Cell_range& column);
	template<class Cell_range>
	bool _multiply_and_add(const Field_element_type& val, const Cell_range& column);
	template<class Cell_range>
	bool _multiply_and_add(const Cell_range& column, const Field_element_type& val);

private:
	using dim_opt = typename Master_matrix::Column_dimension_option;
	using chain_opt = typename Master_matrix::Chain_column_option;
};

template<class Master_matrix>
inline Heap_column<Master_matrix>::Heap_column() : dim_opt(), chain_opt(), insertsSinceLastPrune_(0)
{}

template<class Master_matrix>
template<class Container_type>
inline Heap_column<Master_matrix>::Heap_column(const Container_type &nonZeroRowIndices)
	: dim_opt(nonZeroRowIndices.size() == 0 ? 0 : nonZeroRowIndices.size() - 1), 
	  chain_opt(), 
	  column_(nonZeroRowIndices.size(), nullptr), 
	  insertsSinceLastPrune_(0)
{
	static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type, 
						"Constructor not available for chain columns, please specify the dimension of the chain.");

	unsigned int i = 0;
	if constexpr (Master_matrix::Option_list::is_z2){
		for (index id : nonZeroRowIndices){
			column_[i++] = cellPool_.construct(id);
		}
	} else {
		for (const auto& p : nonZeroRowIndices){
			column_[i++] = cellPool_.construct(p.second, p.first);
		}
	}
	std::make_heap(column_.begin(), column_.end(), cellPointerComp_);
}

template<class Master_matrix>
template<class Container_type>
inline Heap_column<Master_matrix>::Heap_column(
	const Container_type &nonZeroRowIndices, dimension_type dimension) 
	: dim_opt(dimension),
	  chain_opt([&]{
			if constexpr (Master_matrix::Option_list::is_z2){
				return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : *std::prev(nonZeroRowIndices.end());
			} else {
				return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : std::prev(nonZeroRowIndices.end())->first;
			}
		}()), 
	  column_(nonZeroRowIndices.size(), nullptr), 
	  insertsSinceLastPrune_(0)
{
	unsigned int i = 0;
	if constexpr (Master_matrix::Option_list::is_z2){
		for (index id : nonZeroRowIndices){
			column_[i++] = cellPool_.construct(id);
		}
	} else {
		for (const auto& p : nonZeroRowIndices){
			column_[i++] = cellPool_.construct(p.second, p.first);
		}
	}
	std::make_heap(column_.begin(), column_.end(), cellPointerComp_);
}

template<class Master_matrix>
inline Heap_column<Master_matrix>::Heap_column(const Heap_column &column) 
	: dim_opt(static_cast<const dim_opt&>(column)), 
	  chain_opt(static_cast<const chain_opt&>(column)),
	  column_(column.column_.size(), nullptr), 
	  insertsSinceLastPrune_(0)
{
	static_assert(!Master_matrix::Option_list::has_row_access,
			"Simple copy constructor not available when row access option enabled. Please specify the new column index and the row container.");

	unsigned int i = 0;
	for (const Cell* cell : column.column_){
		if constexpr (Master_matrix::Option_list::is_z2){
			column_[i++] = cellPool_.construct(cell->get_row_index());
		} else {
			column_[i++] = cellPool_.construct(cell->get_element(), cell->get_row_index());
		}
	}
	//column.column_ already ordered as a heap, so no need of make_heap.
}

template<class Master_matrix>
inline Heap_column<Master_matrix>::Heap_column(Heap_column &&column) noexcept 
	: dim_opt(std::move(static_cast<dim_opt&>(column))), 
	  chain_opt(std::move(static_cast<chain_opt&>(column))), 
	  column_(std::move(column.column_)), 
	  insertsSinceLastPrune_(std::exchange(column.insertsSinceLastPrune_, 0))
{}

template<class Master_matrix>
template<class Container_type, class Row_container_type>
inline Heap_column<Master_matrix>::Heap_column(
	index columnIndex, const Container_type &nonZeroRowIndices, Row_container_type &rowContainer) 
	: dim_opt(nonZeroRowIndices.size() == 0 ? 0 : nonZeroRowIndices.size() - 1),
	  chain_opt([&]{
			if constexpr (Master_matrix::Option_list::is_z2){
				return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : *std::prev(nonZeroRowIndices.end());
			} else {
				return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : std::prev(nonZeroRowIndices.end())->first;
			}
		}()), 
	  column_(nonZeroRowIndices.size(), nullptr), 
	  insertsSinceLastPrune_(0)
{
	static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type, 
						"Constructor not available for chain columns, please specify the dimension of the chain.");

	unsigned int i = 0;
	if constexpr (Master_matrix::Option_list::is_z2){
		for (index id : nonZeroRowIndices){
			column_[i++] = cellPool_.construct(id);
		}
	} else {
		for (const auto& p : nonZeroRowIndices){
			column_[i++] = cellPool_.construct(p.second, p.first);
		}
	}
	std::make_heap(column_.begin(), column_.end(), cellPointerComp_);
}

template<class Master_matrix>
template<class Container_type, class Row_container_type>
inline Heap_column<Master_matrix>::Heap_column(
	index columnIndex, const Container_type &nonZeroRowIndices, dimension_type dimension, Row_container_type &rowContainer) 
	: dim_opt(dimension),
	  chain_opt([&]{
			if constexpr (Master_matrix::Option_list::is_z2){
				return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : *std::prev(nonZeroRowIndices.end());
			} else {
				return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : std::prev(nonZeroRowIndices.end())->first;
			}
		}()), 
	  column_(nonZeroRowIndices.size(), nullptr), 
	  insertsSinceLastPrune_(0)
{
	unsigned int i = 0;
	if constexpr (Master_matrix::Option_list::is_z2){
		for (index id : nonZeroRowIndices){
			column_[i++] = cellPool_.construct(id);
		}
	} else {
		for (const auto& p : nonZeroRowIndices){
			column_[i++] = cellPool_.construct(p.second, p.first);
		}
	}
	std::make_heap(column_.begin(), column_.end(), cellPointerComp_);
}

template<class Master_matrix>
template<class Row_container_type>
inline Heap_column<Master_matrix>::Heap_column(
	const Heap_column &column, index columnIndex, Row_container_type &rowContainer) 
	: dim_opt(static_cast<const dim_opt&>(column)), 
	  chain_opt(static_cast<const chain_opt&>(column)),
	  column_(column.column_.size(), nullptr), 
	  insertsSinceLastPrune_(0)
{
	unsigned int i = 0;
	for (const Cell* cell : column.column_){
		if constexpr (Master_matrix::Option_list::is_z2){
			column_[i++] = cellPool_.construct(cell->get_row_index());
		} else {
			column_[i++] = cellPool_.construct(cell->get_element(), cell->get_row_index());
		}
	}
	//column.column_ already ordered as a heap, so no need of make_heap.
}

template<class Master_matrix>
inline Heap_column<Master_matrix>::~Heap_column()
{
	for (auto* cell : column_){
		cellPool_.destroy(cell);
	}
}

template<class Master_matrix>
inline std::vector<typename Heap_column<Master_matrix>::Field_element_type> 
Heap_column<Master_matrix>::get_content(int columnLength) const
{
	bool pivotLength = (columnLength < 0);
	if (columnLength < 0 && column_.size() > 0) columnLength = column_.front()->get_row_index() + 1;
	else if (columnLength < 0) return std::vector<Field_element_type>();

	std::vector<Field_element_type> container(columnLength, 0);
	for (auto it = column_.begin(); it != column_.end(); ++it){
		if ((*it)->get_row_index() < static_cast<index>(columnLength)){
			if constexpr (Master_matrix::Option_list::is_z2){
				container[(*it)->get_row_index()] = !container[(*it)->get_row_index()];
			} else {
				container[(*it)->get_row_index()] += (*it)->get_element();
			}
		}
	}

	if (pivotLength){
		while (!container.empty() && container.back() == 0u) container.pop_back();
	}

	return container;
}

template<class Master_matrix>
inline bool Heap_column<Master_matrix>::is_non_zero(index rowIndex) const
{
	if constexpr (Master_matrix::Option_list::is_z2){
		bool c = false;
		for (const Cell* cell : column_){
			if (cell->get_row_index() == rowIndex) c = !c;
		}
		return c;
	} else {
		Field_element_type c(0);
		for (const Cell* cell : column_){
			if (cell->get_row_index() == rowIndex) c += cell->get_element();
		}
		return c != Field_element_type::get_additive_identity();
	}
}

template<class Master_matrix>
inline bool Heap_column<Master_matrix>::is_empty()
{
	Cell* pivot = _pop_pivot();
	if (pivot != nullptr){
		column_.push_back(pivot);
		std::push_heap(column_.begin(), column_.end(), cellPointerComp_);
		return false;
	}
	return true;
}

template<class Master_matrix>
inline std::size_t Heap_column<Master_matrix>::size() const {
	return column_.size();
}

template<class Master_matrix>
template<class Map_type>
inline void Heap_column<Master_matrix>::reorder(const Map_type &valueMap)
{
	static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type, 
						"Method not available for chain columns.");

	Column_type tempCol;
	Cell* pivot = _pop_pivot();
	while (pivot != nullptr) {
		pivot->set_row_index(valueMap.at(pivot->get_row_index()));
		tempCol.push_back(pivot);
		pivot = _pop_pivot();
	}
	column_.swap(tempCol);
	std::make_heap(column_.begin(), column_.end(), cellPointerComp_);

	insertsSinceLastPrune_ = 0;
}

template<class Master_matrix>
inline void Heap_column<Master_matrix>::clear()
{
	static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type, 
						"Method not available for chain columns as a base element should not be empty.");

	for (auto* cell : column_){
		cellPool_.destroy(cell);
	}

	column_.clear();
	insertsSinceLastPrune_ = 0;
}

template<class Master_matrix>
inline void Heap_column<Master_matrix>::clear(index rowIndex)
{
	static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type, 
						"Method not available for chain columns.");

	column_.erase(std::remove_if(column_.begin(), 
								 column_.end(), 
								 [rowIndex](const Cell* c){ return c->get_row_index() == rowIndex; }), 
				  column_.end());
	std::make_heap(column_.begin(), column_.end(), cellPointerComp_);
}

template<class Master_matrix>
inline int Heap_column<Master_matrix>::get_pivot()
{
	static_assert(Master_matrix::isNonBasic, "Method not available for base columns.");	//could technically be, but is the notion usefull then?

	if constexpr (Master_matrix::Option_list::is_of_boundary_type){
		Cell* pivot = _pop_pivot();
		if (pivot != nullptr){
			column_.push_back(pivot);
			std::push_heap(column_.begin(), column_.end(), cellPointerComp_);
			return pivot->get_row_index();
		}
		return -1;
	} else {
		return chain_opt::get_pivot();
	}
}

template<class Master_matrix>
inline typename Heap_column<Master_matrix>::Field_element_type Heap_column<Master_matrix>::get_pivot_value()
{
	static_assert(Master_matrix::isNonBasic, "Method not available for base columns.");	//could technically be, but is the notion usefull then?

	if constexpr (Master_matrix::Option_list::is_z2){
		return 1;
	} else {
		if constexpr (Master_matrix::Option_list::is_of_boundary_type){
			Cell* pivot = _pop_pivot();
			if (pivot != nullptr){
				column_.push_back(pivot);
				std::push_heap(column_.begin(), column_.end(), cellPointerComp_);
				return pivot->get_element();
			}
			return 0;
		} else {
			Field_element_type sum(0);
			if (chain_opt::get_pivot() == -1) return sum;
			for (const Cell* cell : column_){
				if (cell->get_row_index() == static_cast<unsigned int>(chain_opt::get_pivot())) sum += cell->get_element();
			}
			return sum;	//should not be 0 if properly used.
		}
	}
}

template<class Master_matrix>
inline typename Heap_column<Master_matrix>::iterator
Heap_column<Master_matrix>::begin() noexcept
{
	return column_.begin();
}

template<class Master_matrix>
inline typename Heap_column<Master_matrix>::const_iterator
Heap_column<Master_matrix>::begin() const noexcept
{
	return column_.begin();
}

template<class Master_matrix>
inline typename Heap_column<Master_matrix>::iterator
Heap_column<Master_matrix>::end() noexcept
{
	return column_.end();
}

template<class Master_matrix>
inline typename Heap_column<Master_matrix>::const_iterator
Heap_column<Master_matrix>::end() const noexcept
{
	return column_.end();
}

template<class Master_matrix>
inline typename Heap_column<Master_matrix>::reverse_iterator
Heap_column<Master_matrix>::rbegin() noexcept
{
	return column_.rbegin();
}

template<class Master_matrix>
inline typename Heap_column<Master_matrix>::const_reverse_iterator
Heap_column<Master_matrix>::rbegin() const noexcept
{
	return column_.rbegin();
}

template<class Master_matrix>
inline typename Heap_column<Master_matrix>::reverse_iterator
Heap_column<Master_matrix>::rend() noexcept
{
	return column_.rend();
}

template<class Master_matrix>
inline typename Heap_column<Master_matrix>::const_reverse_iterator
Heap_column<Master_matrix>::rend() const noexcept
{
	return column_.rend();
}

template<class Master_matrix>
template<class Cell_range>
inline Heap_column<Master_matrix> &
Heap_column<Master_matrix>::operator+=(const Cell_range &column)
{
	static_assert((!Master_matrix::isNonBasic || std::is_same_v<Cell_range, Heap_column>), 
					"For boundary columns, the range has to be a column of same type to help ensure the validity of the base element.");	//could be removed, if we give the responsability to the user.
	static_assert((!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type), 
					"For chain columns, the given column cannot be constant.");

	_add(column);

	return *this;
}

template<class Master_matrix>
inline Heap_column<Master_matrix> &
Heap_column<Master_matrix>::operator+=(Heap_column &column)
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
inline Heap_column<Master_matrix> &
Heap_column<Master_matrix>::operator*=(unsigned int v)
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
		}
	}

	return *this;
}

template<class Master_matrix>
template<class Cell_range>
inline Heap_column<Master_matrix> &
Heap_column<Master_matrix>::multiply_and_add(const Field_element_type& val, const Cell_range& column)
{
	static_assert((!Master_matrix::isNonBasic || std::is_same_v<Cell_range, Heap_column>), 
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
inline Heap_column<Master_matrix> &
Heap_column<Master_matrix>::multiply_and_add(const Field_element_type& val, Heap_column& column)
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
inline Heap_column<Master_matrix> &
Heap_column<Master_matrix>::multiply_and_add(const Cell_range& column, const Field_element_type& val)
{
	static_assert((!Master_matrix::isNonBasic || std::is_same_v<Cell_range, Heap_column>), 
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
inline Heap_column<Master_matrix> &
Heap_column<Master_matrix>::multiply_and_add(Heap_column& column, const Field_element_type& val)
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
inline Heap_column<Master_matrix> &
Heap_column<Master_matrix>::operator=(const Heap_column& other)
{
	static_assert(!Master_matrix::Option_list::has_row_access, "= assignement not enabled with row access option.");

	dim_opt::operator=(other);
	chain_opt::operator=(other);

	while (column_.size() > other.column_.size()) {
		cellPool_.destroy(column_.back());
		column_.pop_back();
	}
	
	column_.resize(other.column_.size(), nullptr);
	unsigned int i = 0;
	for (const Cell* cell : other.column_){
		if (column_[i] != nullptr){
			cellPool_.destroy(column_[i]);
		}
		if constexpr (Master_matrix::Option_list::is_z2){
			column_[i++] = cellPool_.construct(cell->get_row_index());
		} else {
			column_[i++] = cellPool_.construct(cell->get_element(), cell->get_row_index());
		}
	}
	insertsSinceLastPrune_ = other.insertsSinceLastPrune_;
	
	return *this;
}

template<class Master_matrix>
inline void Heap_column<Master_matrix>::_prune()
{
	if (insertsSinceLastPrune_ == 0) return;

	Column_type tempCol;
	Cell* pivot = _pop_pivot();
	while (pivot != nullptr) {
		tempCol.push_back(pivot);
		pivot = _pop_pivot();
	}
	column_.swap(tempCol);
	std::make_heap(column_.begin(), column_.end(), cellPointerComp_);

	insertsSinceLastPrune_ = 0;
}

template<class Master_matrix>
inline typename Heap_column<Master_matrix>::Cell* Heap_column<Master_matrix>::_pop_pivot()
{
	if (column_.empty()) {
		return nullptr;
	}

	Cell* pivot = column_.front();
	std::pop_heap(column_.begin(), column_.end(), cellPointerComp_);
	column_.pop_back();
	if constexpr (Master_matrix::Option_list::is_z2){
		while (!column_.empty() && column_.front()->get_row_index() == pivot->get_row_index())
		{
			std::pop_heap(column_.begin(), column_.end(), cellPointerComp_);
			column_.pop_back();

			if (column_.empty()) {
				return nullptr;
			}
			pivot = column_.front();
			std::pop_heap(column_.begin(), column_.end(), cellPointerComp_);
			column_.pop_back();
		}
	} else {
		while (!column_.empty() && column_.front()->get_row_index() == pivot->get_row_index())
		{
			pivot->get_element() += column_.front()->get_element();
			std::pop_heap(column_.begin(), column_.end(), cellPointerComp_);
			column_.pop_back();
		}

		if (pivot->get_element() == Field_element_type::get_additive_identity()) return _pop_pivot();
	}

	return pivot;
}

template<class Master_matrix>
template<class Cell_range>
inline bool Heap_column<Master_matrix>::_add(const Cell_range &column)
{
	if (column.begin() == column.end()) return false;
	if (column_.empty()){	//chain should never enter here.
		column_.resize(column.size());
		unsigned int i = 0;
		for (const Cell& cell : column){
			if constexpr (Master_matrix::Option_list::is_z2){
				column_[i++] = cellPool_.construct(cell.get_row_index());
			} else {
				column_[i++] = cellPool_.construct(cell.get_element(), cell.get_row_index());
			}
		}
		insertsSinceLastPrune_ = column_.size();
		return true;
	}

	Field_element_type pivotVal(1);
	
	if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type) 
		pivotVal = get_pivot_value();

	for (const Cell& cell : column) {
		++insertsSinceLastPrune_;
		if constexpr (Master_matrix::Option_list::is_z2){
			if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type){
				if (static_cast<int>(cell.get_row_index()) == chain_opt::get_pivot()) pivotVal = !pivotVal;
			}
			column_.push_back(cellPool_.construct(cell.get_row_index()));
		} else {
			if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type){
				if (static_cast<int>(cell.get_row_index()) == chain_opt::get_pivot()) pivotVal += cell.get_element();
			}
			column_.push_back(cellPool_.construct(cell.get_element(), cell.get_row_index()));
		}
		std::push_heap(column_.begin(), column_.end(), cellPointerComp_);
	}

	if (2 * insertsSinceLastPrune_ > column_.size()) _prune();

	return pivotVal == 0u;
}

template<class Master_matrix>
template<class Cell_range>
inline bool Heap_column<Master_matrix>::_multiply_and_add(const Field_element_type& val, const Cell_range& column)
{
	if (val == 0u) {
		if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type){
			throw std::invalid_argument("A chain column should not be multiplied by 0.");
			//this would not only mess up the base, but also the pivots stored.
		} else {
			clear();
		}
	}
	if (column_.empty()){	//chain should never enter here.
		column_.resize(column.size());
		unsigned int i = 0;
		for (const Cell& cell : column){
			if constexpr (Master_matrix::Option_list::is_z2){
				column_[i++] = cellPool_.construct(cell.get_row_index());
			} else {
				column_[i++] = cellPool_.construct(cell.get_element(), cell.get_row_index());
			}
		}
		insertsSinceLastPrune_ = column_.size();
		return true;
	}

	Field_element_type pivotVal(0);

	for (Cell* cell : column_){
		cell->get_element() *= val;
		if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type){
			if (static_cast<int>(cell->get_row_index()) == chain_opt::get_pivot()) pivotVal += cell->get_element();
		}
	}

	for (const Cell& cell : column) {
		++insertsSinceLastPrune_;
		if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type){
			if (static_cast<int>(cell.get_row_index()) == chain_opt::get_pivot()) pivotVal += cell.get_element();
		}
		column_.push_back(cellPool_.construct(cell.get_element(), cell.get_row_index()));
		std::push_heap(column_.begin(), column_.end(), cellPointerComp_);
	}

	if (2 * insertsSinceLastPrune_ > column_.size()) _prune();

	return pivotVal == 0u;
}

template<class Master_matrix>
template<class Cell_range>
inline bool Heap_column<Master_matrix>::_multiply_and_add(const Cell_range& column, const Field_element_type& val)
{
	if (val == 0u || column.begin() == column.end()) {
		return false;
	}
	if (column_.empty()){	//chain should never enter here.
		column_.resize(column.size());
		unsigned int i = 0;
		for (const Cell& cell : column){
			if constexpr (Master_matrix::Option_list::is_z2){
				column_[i++] = cellPool_.construct(cell.get_row_index());
			} else {
				column_[i++] = cellPool_.construct(cell.get_element(), cell.get_row_index());
			}
		}
		insertsSinceLastPrune_ = column_.size();
		return true;
	}

	Field_element_type pivotVal(1);
	
	if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type) 
		pivotVal = get_pivot_value();

	for (const Cell& cell : column) {
		++insertsSinceLastPrune_;
		if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type){
			if (static_cast<int>(cell.get_row_index()) == chain_opt::get_pivot()) pivotVal += (cell.get_element() * val);
		}
		column_.push_back(cellPool_.construct(cell.get_element() * val, cell.get_row_index()));
		std::push_heap(column_.begin(), column_.end(), cellPointerComp_);
	}

	if (2 * insertsSinceLastPrune_ > column_.size()) _prune();

	return pivotVal == 0u;
}

} //namespace persistence_matrix
} //namespace Gudhi

template<class Master_matrix>
struct std::hash<Gudhi::persistence_matrix::Heap_column<Master_matrix> >
{
	size_t operator()(const Gudhi::persistence_matrix::Heap_column<Master_matrix>& column) const
	{
		std::size_t seed = 0;
		unsigned int i = 0;
		for (bool val : column.get_content()){
			seed ^= std::hash<unsigned int>()(i++ * val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		}
		return seed;
	}
};

#endif // PM_HEAP_COLUMN_H

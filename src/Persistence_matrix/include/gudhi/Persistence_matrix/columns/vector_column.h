/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022-23 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_VECTOR_COLUMN_H
#define PM_VECTOR_COLUMN_H

#include <vector>
#include <stdexcept>
#include <type_traits>
#include <algorithm>	//binary_search
#include <unordered_set>
#include <utility>	//std::swap, std::move & std::exchange

#include <boost/iterator/indirect_iterator.hpp>

#include <gudhi/Persistence_matrix/columns/cell_constructors.h>

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix, class Cell_constructor = New_cell_constructor<typename Master_matrix::Cell_type> >
class Vector_column : public Master_matrix::Row_access_option, 
					  public Master_matrix::Column_dimension_option,
					  public Master_matrix::Chain_column_option
{
public:
	using Master = Master_matrix;
	using Field_operators = typename Master_matrix::Field_operators;
	using Field_element_type = typename Master_matrix::element_type;
	using index = typename Master_matrix::index;
	using id_index = typename Master_matrix::id_index;
	using dimension_type = typename Master_matrix::dimension_type;

	using Cell = typename Master_matrix::Cell_type;
	using Column_type = std::vector<Cell*>;
	using iterator = boost::indirect_iterator<typename Column_type::iterator>;
	using const_iterator = boost::indirect_iterator<typename Column_type::const_iterator>;
	using reverse_iterator = boost::indirect_iterator<typename Column_type::reverse_iterator>;
	using const_reverse_iterator = boost::indirect_iterator<typename Column_type::const_reverse_iterator>;

	Vector_column(Field_operators* operators = nullptr, Cell_constructor* cellConstructor = nullptr);
	template<class Container_type = typename Master_matrix::boundary_type>
	Vector_column(const Container_type& nonZeroRowIndices, Field_operators* operators, Cell_constructor* cellConstructor);	//has to be a boundary for boundary, has no sense for chain if dimension is needed
	template<class Container_type = typename Master_matrix::boundary_type, class Row_container_type>
	Vector_column(index columnIndex, const Container_type& nonZeroRowIndices, Row_container_type *rowContainer, Field_operators* operators, Cell_constructor* cellConstructor);	//has to be a boundary for boundary, has no sense for chain if dimension is needed
	template<class Container_type = typename Master_matrix::boundary_type>
	Vector_column(const Container_type& nonZeroChainRowIndices, dimension_type dimension, Field_operators* operators, Cell_constructor* cellConstructor);	//dimension gets ignored for base
	template<class Container_type = typename Master_matrix::boundary_type, class Row_container_type>
	Vector_column(index columnIndex, const Container_type& nonZeroChainRowIndices, dimension_type dimension, Row_container_type *rowContainer, Field_operators* operators, Cell_constructor* cellConstructor);	//dimension gets ignored for base
	Vector_column(const Vector_column& column, Field_operators* operators = nullptr, Cell_constructor* cellConstructor = nullptr);
	template<class Row_container_type>
	Vector_column(const Vector_column& column, index columnIndex, Row_container_type *rowContainer, Field_operators* operators = nullptr, Cell_constructor* cellConstructor = nullptr);
	Vector_column(Vector_column&& column) noexcept;
	~Vector_column();

	std::vector<Field_element_type> get_content(int columnLength = -1) const;
	bool is_non_zero(id_index rowIndex) const;
	bool is_empty() const;
	std::size_t size() const;

	//****************
	//only for base and boundary
	template<class Map_type>
	void reorder(const Map_type& valueMap);	//used for lazy row swaps
	void clear();
	void clear(id_index rowIndex);				//do not clear a cell to 0 if the cell was already 0, otherwise size/is_empty will be wrong.
	//****************

	//****************
	//only for chain and boundary
	id_index get_pivot();
	Field_element_type get_pivot_value();
	//****************

	iterator begin() noexcept;						//potentally does not ignore erased values
	const_iterator begin() const noexcept;			//potentally does not ignore erased values
	iterator end() noexcept;						//potentally does not ignore erased values
	const_iterator end() const noexcept;			//potentally does not ignore erased values
	reverse_iterator rbegin() noexcept;				//potentally does not ignore erased values
	const_reverse_iterator rbegin() const noexcept;	//potentally does not ignore erased values
	reverse_iterator rend() noexcept;				//potentally does not ignore erased values
	const_reverse_iterator rend() const noexcept;	//potentally does not ignore erased values

	template<class Cell_range>
	Vector_column& operator+=(const Cell_range& column);	//for base & boundary except vector	//may not work if Cell_range = Vector_column<Other>
	Vector_column& operator+=(Vector_column &column);	//for chain and vector

	Vector_column& operator*=(unsigned int v);

	//this = v * this + column
	template<class Cell_range>
	Vector_column& multiply_and_add(const Field_element_type& val, const Cell_range& column);	//for base & boundary except vector	//may not work if Cell_range = Vector_column<Other>
	Vector_column& multiply_and_add(const Field_element_type& val, Vector_column& column);	//for chain and vector
	//this = this + column * v
	template<class Cell_range>
	Vector_column& multiply_and_add(const Cell_range& column, const Field_element_type& val);	//for base & boundary except vector	//may not work if Cell_range = Vector_column<Other>
	Vector_column& multiply_and_add(Vector_column& column, const Field_element_type& val);	//for chain and vector

	friend bool operator==(const Vector_column& c1, const Vector_column& c2){
		if (&c1 == &c2) return true;
		if (c1.erasedValues_.empty() && c2.erasedValues_.empty() && c1.column_.size() != c2.column_.size()) 
			return false;

		auto it1 = c1.column_.begin();
		auto it2 = c2.column_.begin();
		while (it1 != c1.column_.end() && it2 != c2.column_.end()) {
			if constexpr (!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type){
				while (it1 != c1.column_.end() && c1.erasedValues_.find((*it1)->get_row_index()) != c1.erasedValues_.end()) ++it1;
				while (it2 != c2.column_.end() && c2.erasedValues_.find((*it2)->get_row_index()) != c2.erasedValues_.end()) ++it2;
				if (it1 == c1.column_.end() || it2 == c2.column_.end()) break;
			}
			if constexpr (Master_matrix::Option_list::is_z2){
				if ((*it1)->get_row_index() != (*it2)->get_row_index())
					return false;
			} else {
				if ((*it1)->get_row_index() != (*it2)->get_row_index() || (*it1)->get_element() != (*it2)->get_element())
					return false;
			}
			++it1; ++it2;
		}

		if constexpr (!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type){
			while (it1 != c1.column_.end() && c1.erasedValues_.find((*it1)->get_row_index()) != c1.erasedValues_.end()) ++it1;
			while (it2 != c2.column_.end() && c2.erasedValues_.find((*it2)->get_row_index()) != c2.erasedValues_.end()) ++it2;
			return it2 == c2.column_.end() && it1 == c1.column_.end();
		} else {
			return true;
		}
	}
	friend bool operator<(const Vector_column& c1, const Vector_column& c2){
		if (&c1 == &c2) return false;

		auto it1 = c1.column_.begin();
		auto it2 = c2.column_.begin();
		while (it1 != c1.column_.end() && it2 != c2.column_.end()) {
			if constexpr (!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type){
				while (it1 != c1.column_.end() && c1.erasedValues_.find((*it1)->get_row_index()) != c1.erasedValues_.end()) ++it1;
				while (it2 != c2.column_.end() && c2.erasedValues_.find((*it2)->get_row_index()) != c2.erasedValues_.end()) ++it2;
				if (it1 == c1.column_.end() || it2 == c2.column_.end()) break;
			}

			if ((*it1)->get_row_index() != (*it2)->get_row_index())
				return (*it1)->get_row_index() < (*it2)->get_row_index();
			if constexpr (!Master_matrix::Option_list::is_z2){
				if ((*it1)->get_element() != (*it2)->get_element())
					return (*it1)->get_element() < (*it2)->get_element();
			}
			++it1; ++it2;
		}
		if constexpr (!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type){
			while (it1 != c1.column_.end() && c1.erasedValues_.find((*it1)->get_row_index()) != c1.erasedValues_.end()) ++it1;
			while (it2 != c2.column_.end() && c2.erasedValues_.find((*it2)->get_row_index()) != c2.erasedValues_.end()) ++it2;
		}
		return it2 != c2.column_.end();
	}

	// void set_operators(Field_operators* operators){ operators_ = operators; }

	//Disabled with row access.
	Vector_column& operator=(const Vector_column& other);

	friend void swap(Vector_column& col1, Vector_column& col2){
		swap(static_cast<typename Master_matrix::Row_access_option&>(col1),
			 static_cast<typename Master_matrix::Row_access_option&>(col2));
		swap(static_cast<typename Master_matrix::Column_dimension_option&>(col1),
			 static_cast<typename Master_matrix::Column_dimension_option&>(col2));
		swap(static_cast<typename Master_matrix::Chain_column_option&>(col1),
			 static_cast<typename Master_matrix::Chain_column_option&>(col2));
		col1.column_.swap(col2.column_);
		col1.erasedValues_.swap(col2.erasedValues_);
		std::swap(col1.operators_, col2.operators_);
		std::swap(col1.cellPool_, col2.cellPool_);
	}

private:
	using ra_opt = typename Master_matrix::Row_access_option;
	using dim_opt = typename Master_matrix::Column_dimension_option;
	using chain_opt = typename Master_matrix::Chain_column_option;

	Column_type column_;
	std::unordered_set<id_index> erasedValues_;	//TODO: test other containers? Useless when clear(index) is never called, how much is it worth it?
	Field_operators* operators_;
	Cell_constructor* cellPool_;

	// void _clean_values();
	void _delete_cell(Cell* cell);
	void _insert_cell(const Field_element_type& value, id_index rowIndex, Column_type& column);
	void _insert_cell(id_index rowIndex, Column_type& column);
	void _update_cell(const Field_element_type& value, id_index rowIndex, index position);
	void _update_cell(id_index rowIndex, index position);
	template<class Cell_range>
	bool _add(const Cell_range& column);
	template<class Cell_range>
	bool _multiply_and_add(const Field_element_type& val, const Cell_range& column);
	template<class Cell_range>
	bool _multiply_and_add(const Cell_range& column, const Field_element_type& val);

	void _verifyCellConstructor(){
		if (cellPool_ == nullptr){
			if constexpr (std::is_same_v<Cell_constructor, New_cell_constructor<typename Master_matrix::Cell_type> >){
				cellPool_ = &Master_matrix::defaultCellConstructor;
			} else {
				throw std::invalid_argument("Cell constructor pointer cannot be null.");
			}
		}
	}
};

template<class Master_matrix, class Cell_constructor>
inline Vector_column<Master_matrix,Cell_constructor>::Vector_column(Field_operators* operators, Cell_constructor* cellConstructor) 
	: ra_opt(), dim_opt(), chain_opt(), operators_(operators),
	  cellPool_(cellConstructor)
{
	if (operators_ == nullptr && cellPool_ == nullptr) return;	//to allow default constructor which gives a dummy column
	_verifyCellConstructor();
}

template<class Master_matrix, class Cell_constructor>
template<class Container_type>
inline Vector_column<Master_matrix,Cell_constructor>::Vector_column(const Container_type &nonZeroRowIndices, Field_operators* operators, Cell_constructor* cellConstructor)
	: ra_opt(), 
	  dim_opt(nonZeroRowIndices.size() == 0 ? 0 : nonZeroRowIndices.size() - 1), 
	  chain_opt(), 
	  column_(nonZeroRowIndices.size(), nullptr),
	  operators_(operators),
	  cellPool_(cellConstructor)
{
	static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type, 
						"Constructor not available for chain columns, please specify the dimension of the chain.");

	_verifyCellConstructor();

	index i = 0;
	if constexpr (Master_matrix::Option_list::is_z2){
		for (id_index id : nonZeroRowIndices){
			_update_cell(id, i++);
		}
	} else {
		for (const auto& p : nonZeroRowIndices){
			_update_cell(operators_->get_value(p.second), p.first, i++);
		}
	}
}

template<class Master_matrix, class Cell_constructor>
template<class Container_type, class Row_container_type>
inline Vector_column<Master_matrix,Cell_constructor>::Vector_column(
	index columnIndex, const Container_type &nonZeroRowIndices, Row_container_type *rowContainer, Field_operators* operators, Cell_constructor* cellConstructor) 
	: ra_opt(columnIndex, rowContainer), 
	  dim_opt(nonZeroRowIndices.size() == 0 ? 0 : nonZeroRowIndices.size() - 1),
	  chain_opt([&]{
			if constexpr (Master_matrix::Option_list::is_z2){
				return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : *std::prev(nonZeroRowIndices.end());
			} else {
				return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : std::prev(nonZeroRowIndices.end())->first;
			}
		}()), 
	  column_(nonZeroRowIndices.size(), nullptr),
	  operators_(operators),
	  cellPool_(cellConstructor)
{
	static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type, 
						"Constructor not available for chain columns, please specify the dimension of the chain.");

	_verifyCellConstructor();

	index i = 0;
	if constexpr (Master_matrix::Option_list::is_z2){
		for (id_index id : nonZeroRowIndices){
			_update_cell(id, i++);
		}
	} else {
		for (const auto& p : nonZeroRowIndices){
			_update_cell(operators_->get_value(p.second), p.first, i++);
		}
	}
}

template<class Master_matrix, class Cell_constructor>
template<class Container_type>
inline Vector_column<Master_matrix,Cell_constructor>::Vector_column(
	const Container_type &nonZeroRowIndices, dimension_type dimension, Field_operators* operators, Cell_constructor* cellConstructor) 
	: ra_opt(), 
	  dim_opt(dimension),
	  chain_opt([&]{
			if constexpr (Master_matrix::Option_list::is_z2){
				return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : *std::prev(nonZeroRowIndices.end());
			} else {
				return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : std::prev(nonZeroRowIndices.end())->first;
			}
		}()), 
	  column_(nonZeroRowIndices.size(), nullptr),
	  operators_(operators),
	  cellPool_(cellConstructor)
{
	_verifyCellConstructor();

	index i = 0;
	if constexpr (Master_matrix::Option_list::is_z2){
		for (id_index id : nonZeroRowIndices){
			_update_cell(id, i++);
		}
	} else {
		for (const auto& p : nonZeroRowIndices){
			_update_cell(operators_->get_value(p.second), p.first, i++);
		}
	}
}

template<class Master_matrix, class Cell_constructor>
template<class Container_type, class Row_container_type>
inline Vector_column<Master_matrix,Cell_constructor>::Vector_column(
	index columnIndex, const Container_type &nonZeroRowIndices, dimension_type dimension, Row_container_type *rowContainer, Field_operators* operators, Cell_constructor* cellConstructor) 
	: ra_opt(columnIndex, rowContainer), 
	  dim_opt(dimension),
	  chain_opt([&]{
			if constexpr (Master_matrix::Option_list::is_z2){
				return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : *std::prev(nonZeroRowIndices.end());
			} else {
				return nonZeroRowIndices.begin() == nonZeroRowIndices.end() ? -1 : std::prev(nonZeroRowIndices.end())->first;
			}
		}()), 
	  column_(nonZeroRowIndices.size(), nullptr),
	  operators_(operators),
	  cellPool_(cellConstructor)
{
	_verifyCellConstructor();

	index i = 0;
	if constexpr (Master_matrix::Option_list::is_z2){
		for (id_index id : nonZeroRowIndices){
			_update_cell(id, i++);
		}
	} else {
		for (const auto& p : nonZeroRowIndices){
			_update_cell(operators_->get_value(p.second), p.first, i++);
		}
	}
}

template<class Master_matrix, class Cell_constructor>
inline Vector_column<Master_matrix,Cell_constructor>::Vector_column(const Vector_column &column, Field_operators* operators, Cell_constructor* cellConstructor) 
	: ra_opt(), 
	  dim_opt(static_cast<const dim_opt&>(column)), 
	  chain_opt(static_cast<const chain_opt&>(column)),
	  column_(column.column_.size(), nullptr),
	  erasedValues_(column.erasedValues_),
	  operators_(operators == nullptr ? column.operators_ : operators),
	  cellPool_(cellConstructor == nullptr ? column.cellPool_ : cellConstructor)
{
	static_assert(!Master_matrix::Option_list::has_row_access,
			"Simple copy constructor not available when row access option enabled. Please specify the new column index and the row container.");

	index i = 0;
	for (const Cell* cell : column.column_){
		if constexpr (Master_matrix::Option_list::is_z2){
			_update_cell(cell->get_row_index(), i++);
		} else {
			_update_cell(cell->get_element(), cell->get_row_index(), i++);
		}
	}
}

template<class Master_matrix, class Cell_constructor>
template<class Row_container_type>
inline Vector_column<Master_matrix,Cell_constructor>::Vector_column(
	const Vector_column &column, index columnIndex, Row_container_type *rowContainer, Field_operators* operators, Cell_constructor* cellConstructor) 
	: ra_opt(columnIndex, rowContainer), 
	  dim_opt(static_cast<const dim_opt&>(column)), 
	  chain_opt(static_cast<const chain_opt&>(column)),
	  column_(column.column_.size(), nullptr),
	  erasedValues_(column.erasedValues_),
	  operators_(operators == nullptr ? column.operators_ : operators),
	  cellPool_(cellConstructor == nullptr ? column.cellPool_ : cellConstructor)
{
	index i = 0;
	for (const Cell* cell : column.column_){
		if constexpr (Master_matrix::Option_list::is_z2){
			_update_cell(cell->get_row_index(), i++);
		} else {
			_update_cell(cell->get_element(), cell->get_row_index(), i++);
		}
	}
}

template<class Master_matrix, class Cell_constructor>
inline Vector_column<Master_matrix,Cell_constructor>::Vector_column(Vector_column &&column) noexcept 
	: ra_opt(std::move(static_cast<ra_opt&>(column))), 
	  dim_opt(std::move(static_cast<dim_opt&>(column))), 
	  chain_opt(std::move(static_cast<chain_opt&>(column))), 
	  column_(std::move(column.column_)),
	  erasedValues_(std::move(column.erasedValues_)),
	  operators_(std::exchange(column.operators_, nullptr)),
	  cellPool_(std::exchange(column.cellPool_, nullptr))
{}

template<class Master_matrix, class Cell_constructor>
inline Vector_column<Master_matrix,Cell_constructor>::~Vector_column()
{
	for (auto* cell : column_){
		_delete_cell(cell);
	}
}

template<class Master_matrix, class Cell_constructor>
inline std::vector<typename Vector_column<Master_matrix,Cell_constructor>::Field_element_type> 
Vector_column<Master_matrix,Cell_constructor>::get_content(int columnLength) const
{
	if (columnLength < 0 && column_.size() > 0) columnLength = column_.back()->get_row_index() + 1;
	else if (columnLength < 0) return std::vector<Field_element_type>();

	std::vector<Field_element_type> container(columnLength, 0);
	for (auto it = column_.begin(); it != column_.end() && (*it)->get_row_index() < static_cast<id_index>(columnLength); ++it){
		if constexpr (!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type){
			if (erasedValues_.find((*it)->get_row_index()) != erasedValues_.end()) continue;
		}
		if constexpr (Master_matrix::Option_list::is_z2){
			container[(*it)->get_row_index()] = 1;
		} else {
			container[(*it)->get_row_index()] = (*it)->get_element();
		}
	}
	return container;
}

template<class Master_matrix, class Cell_constructor>
inline bool Vector_column<Master_matrix,Cell_constructor>::is_non_zero(id_index rowIndex) const
{
	if constexpr (!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type) 
		if (erasedValues_.find(rowIndex) != erasedValues_.end()) return false;

	return std::binary_search(column_.begin(), column_.end(), 
							  cellPool_->construct(rowIndex), 	//cell gets destroyed with the pool at the end, but I don't know if that's a good solution
							  [](const Cell* a, const Cell* b){
								  return a->get_row_index() < b->get_row_index();
							  });
}

template<class Master_matrix, class Cell_constructor>
inline bool Vector_column<Master_matrix,Cell_constructor>::is_empty() const
{
	if constexpr (!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type){
		return column_.size() == erasedValues_.size();	//assumes that erasedValues is always a subset of column_, which is wrong if someone cleared an non exitsing value...
	} else {
		return column_.empty();
	}
}

template<class Master_matrix, class Cell_constructor>
inline std::size_t Vector_column<Master_matrix,Cell_constructor>::size() const{
	if constexpr (!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type){
		return column_.size() - erasedValues_.size();	//assumes that erasedValues is always a subset of column_, which is wrong if someone cleared an non exitsing value...
	} else {
		return column_.size();
	}
}

template<class Master_matrix, class Cell_constructor>
template<class Map_type>
inline void Vector_column<Master_matrix,Cell_constructor>::reorder(const Map_type &valueMap)
{
	static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type, 
						"Method not available for chain columns.");

	if (erasedValues_.empty()){		//to avoid useless push_backs.
		for (Cell* cell : column_) {
			if constexpr (Master_matrix::Option_list::has_row_access) ra_opt::unlink(cell);
			cell->set_row_index(valueMap.at(cell->get_row_index()));
			if constexpr (Master_matrix::Option_list::has_intrusive_rows && Master_matrix::Option_list::has_row_access) 
				ra_opt::insert_cell(cell->get_row_index(), cell);
		}

		//all cells have to be deleted first, to avoid problem with insertion when row is a set
		if constexpr (!Master_matrix::Option_list::has_intrusive_rows && Master_matrix::Option_list::has_row_access){
			for (Cell* cell : column_) {
				ra_opt::insert_cell(cell->get_row_index(), cell);
			}
		}

		std::sort(column_.begin(), column_.end(), [](const Cell* c1, const Cell* c2){return *c1 < *c2;});
	} else {
		Column_type newColumn;
		for (Cell* cell : column_) {
			if (erasedValues_.find(cell->get_row_index()) == erasedValues_.end()){
				if constexpr (Master_matrix::Option_list::has_row_access) ra_opt::unlink(cell);
				cell->set_row_index(valueMap.at(cell->get_row_index()));
				newColumn.push_back(cell);
				if constexpr (Master_matrix::Option_list::has_intrusive_rows && Master_matrix::Option_list::has_row_access) 
					ra_opt::insert_cell(cell->get_row_index(), cell);
			} else {
				_delete_cell(cell);
			}
		}
		//all cells have to be deleted first, to avoid problem with insertion when row is a set
		if constexpr (!Master_matrix::Option_list::has_intrusive_rows && Master_matrix::Option_list::has_row_access){
			for (Cell* cell : column_) {
				ra_opt::insert_cell(cell->get_row_index(), cell);
			}
		}
		std::sort(newColumn.begin(), newColumn.end(), [](const Cell* c1, const Cell* c2){return *c1 < *c2;});
		erasedValues_.clear();
		column_.swap(newColumn);
	}
}

template<class Master_matrix, class Cell_constructor>
inline void Vector_column<Master_matrix,Cell_constructor>::clear()
{
	static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type, 
						"Method not available for chain columns as a base element should not be empty.");

	for (auto* cell : column_){
		if constexpr (Master_matrix::Option_list::has_row_access) ra_opt::unlink(cell);
		cellPool_->destroy(cell);
	}

	column_.clear();
	erasedValues_.clear();
}

template<class Master_matrix, class Cell_constructor>
inline void Vector_column<Master_matrix,Cell_constructor>::clear(id_index rowIndex)
{
	static_assert(!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type, 
						"Method not available for chain columns.");

	erasedValues_.insert(rowIndex);
}

template<class Master_matrix, class Cell_constructor>
inline typename Vector_column<Master_matrix,Cell_constructor>::id_index Vector_column<Master_matrix,Cell_constructor>::get_pivot()
{
	static_assert(Master_matrix::isNonBasic, "Method not available for base columns.");	//could technically be, but is the notion usefull then?

	if constexpr (Master_matrix::Option_list::is_of_boundary_type){
		if (column_.empty()) return -1;
		if (erasedValues_.empty()) return column_.back()->get_row_index();

		auto it = erasedValues_.find(column_.back()->get_row_index());
		while (!column_.empty() && it != erasedValues_.end()) {
			erasedValues_.erase(it);
			_delete_cell(column_.back());
			column_.pop_back();
			if (!column_.empty()) it = erasedValues_.find(column_.back()->get_row_index());
		}

		if (column_.empty()) return -1;
		return column_.back()->get_row_index();
	} else {
		return chain_opt::get_pivot();
	}
}

template<class Master_matrix, class Cell_constructor>
inline typename Vector_column<Master_matrix,Cell_constructor>::Field_element_type Vector_column<Master_matrix,Cell_constructor>::get_pivot_value()
{
	static_assert(Master_matrix::isNonBasic, "Method not available for base columns.");	//could technically be, but is the notion usefull then?

	if constexpr (Master_matrix::Option_list::is_z2){
		return 1;
	} else {
		if constexpr (Master_matrix::Option_list::is_of_boundary_type){
			if (column_.empty()) return 0;
			if (erasedValues_.empty()) return column_.back()->get_element();

			auto it = erasedValues_.find(column_.back()->get_row_index());
			while (!column_.empty() && it != erasedValues_.end()) {
				erasedValues_.erase(it);
				_delete_cell(column_.back());
				column_.pop_back();
				if (!column_.empty()) it = erasedValues_.find(column_.back()->get_row_index());
			}

			if (column_.empty()) return 0;
			return column_.back()->get_element();
		} else {
			if (chain_opt::get_pivot() == -1) return Field_element_type();
			for (const Cell* cell : column_){
				if (cell->get_row_index() == chain_opt::get_pivot()) return cell->get_element();
			}
			return Field_element_type();	//should never happen if chain column is used properly
		}
	}
}

template<class Master_matrix, class Cell_constructor>
inline typename Vector_column<Master_matrix,Cell_constructor>::iterator
Vector_column<Master_matrix,Cell_constructor>::begin() noexcept
{
	return column_.begin();
}

template<class Master_matrix, class Cell_constructor>
inline typename Vector_column<Master_matrix,Cell_constructor>::const_iterator
Vector_column<Master_matrix,Cell_constructor>::begin() const noexcept
{
	return column_.begin();
}

template<class Master_matrix, class Cell_constructor>
inline typename Vector_column<Master_matrix,Cell_constructor>::iterator
Vector_column<Master_matrix,Cell_constructor>::end() noexcept
{
	return column_.end();
}

template<class Master_matrix, class Cell_constructor>
inline typename Vector_column<Master_matrix,Cell_constructor>::const_iterator
Vector_column<Master_matrix,Cell_constructor>::end() const noexcept
{
	return column_.end();
}

template<class Master_matrix, class Cell_constructor>
inline typename Vector_column<Master_matrix,Cell_constructor>::reverse_iterator
Vector_column<Master_matrix,Cell_constructor>::rbegin() noexcept
{
	return column_.rbegin();
}

template<class Master_matrix, class Cell_constructor>
inline typename Vector_column<Master_matrix,Cell_constructor>::const_reverse_iterator
Vector_column<Master_matrix,Cell_constructor>::rbegin() const noexcept
{
	return column_.rbegin();
}

template<class Master_matrix, class Cell_constructor>
inline typename Vector_column<Master_matrix,Cell_constructor>::reverse_iterator
Vector_column<Master_matrix,Cell_constructor>::rend() noexcept
{
	return column_.rend();
}

template<class Master_matrix, class Cell_constructor>
inline typename Vector_column<Master_matrix,Cell_constructor>::const_reverse_iterator
Vector_column<Master_matrix,Cell_constructor>::rend() const noexcept
{
	return column_.rend();
}

template<class Master_matrix, class Cell_constructor>
template<class Cell_range>
inline Vector_column<Master_matrix,Cell_constructor> &
Vector_column<Master_matrix,Cell_constructor>::operator+=(const Cell_range &column)
{
	static_assert((!Master_matrix::isNonBasic || std::is_same_v<Cell_range, Vector_column>), 
					"For boundary columns, the range has to be a column of same type to help ensure the validity of the base element.");	//could be removed, if we give the responsability to the user.
	static_assert((!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type), 
					"For chain columns, the given column cannot be constant.");

	_add(column);

	return *this;
}

template<class Master_matrix, class Cell_constructor>
inline Vector_column<Master_matrix,Cell_constructor> &
Vector_column<Master_matrix,Cell_constructor>::operator+=(Vector_column &column)
{
	if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type){
		//assumes that the addition never zeros out this column. 
		if (_add(column)){
			chain_opt::swap_pivots(column);
			dim_opt::swap_dimension(column);
		}
	} else {
		_add(column);
	}

	return *this;
}

template<class Master_matrix, class Cell_constructor>
inline Vector_column<Master_matrix,Cell_constructor> &
Vector_column<Master_matrix,Cell_constructor>::operator*=(unsigned int v)
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
		Field_element_type val = operators_->get_value(v);

		if (val == Field_operators::get_additive_identity()) {
			if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type){
				throw std::invalid_argument("A chain column should not be multiplied by 0.");
			} else {
				clear();
			}
			return *this;
		}

		if (val == Field_operators::get_multiplicative_identity()) return *this;

		for (Cell* cell : column_){
			cell->get_element() = operators_->multiply(cell->get_element(), val);
			if constexpr (Master_matrix::Option_list::has_row_access)
				ra_opt::update_cell(*cell);
		}
	}

	return *this;
}

template<class Master_matrix, class Cell_constructor>
template<class Cell_range>
inline Vector_column<Master_matrix,Cell_constructor> &
Vector_column<Master_matrix,Cell_constructor>::multiply_and_add(const Field_element_type& val, const Cell_range& column)
{
	static_assert((!Master_matrix::isNonBasic || std::is_same_v<Cell_range, Vector_column>), 
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

template<class Master_matrix, class Cell_constructor>
inline Vector_column<Master_matrix,Cell_constructor> &
Vector_column<Master_matrix,Cell_constructor>::multiply_and_add(const Field_element_type& val, Vector_column& column)
{
	if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type){
		//assumes that the addition never zeros out this column. 
		if constexpr (Master_matrix::Option_list::is_z2){
			if (val){
				if (_add(column)){
					chain_opt::swap_pivots(column);
					dim_opt::swap_dimension(column);
				}
			} else {
				throw std::invalid_argument("A chain column should not be multiplied by 0.");
			}
		} else {
			if (_multiply_and_add(val, column)){
				chain_opt::swap_pivots(column);
				dim_opt::swap_dimension(column);
			}
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

template<class Master_matrix, class Cell_constructor>
template<class Cell_range>
inline Vector_column<Master_matrix,Cell_constructor> &
Vector_column<Master_matrix,Cell_constructor>::multiply_and_add(const Cell_range& column, const Field_element_type& val)
{
	static_assert((!Master_matrix::isNonBasic || std::is_same_v<Cell_range, Vector_column>), 
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

template<class Master_matrix, class Cell_constructor>
inline Vector_column<Master_matrix,Cell_constructor> &
Vector_column<Master_matrix,Cell_constructor>::multiply_and_add(Vector_column& column, const Field_element_type& val)
{
	if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type){
		//assumes that the addition never zeros out this column. 
		if constexpr (Master_matrix::Option_list::is_z2){
			if (val){
				if (_add(column)){
					chain_opt::swap_pivots(column);
					dim_opt::swap_dimension(column);
				}
			}
		} else {
			if (_multiply_and_add(column, val)){
				chain_opt::swap_pivots(column);
				dim_opt::swap_dimension(column);
			}
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

template<class Master_matrix, class Cell_constructor>
inline Vector_column<Master_matrix,Cell_constructor> &
Vector_column<Master_matrix,Cell_constructor>::operator=(const Vector_column& other)
{
	static_assert(!Master_matrix::Option_list::has_row_access, "= assignement not enabled with row access option.");

	dim_opt::operator=(other);
	chain_opt::operator=(other);

	auto tmpPool = cellPool_;
	cellPool_ = other.cellPool_;

	while (column_.size() > other.column_.size()) {
		if (column_.back() != nullptr){
			if constexpr (Master_matrix::Option_list::has_row_access) ra_opt::unlink(column_.back());
			tmpPool->destroy(column_.back());
		}
		column_.pop_back();
	}
	
	column_.resize(other.column_.size(), nullptr);
	index i = 0;
	for (const Cell* cell : other.column_){
		if (column_[i] != nullptr){
			if constexpr (Master_matrix::Option_list::has_row_access) ra_opt::unlink(column_[i]);
			tmpPool->destroy(column_[i]);
		}
		if constexpr (Master_matrix::Option_list::is_z2){
			_update_cell(cell->get_row_index(), i++);
		} else {
			_update_cell(cell->get_element(), cell->get_row_index(), i++);
		}
	}
	erasedValues_ = other.erasedValues_;
	operators_ = other.operators_;
	
	return *this;
}

// template<class Master_matrix, class Cell_constructor>
// inline void Vector_column<Master_matrix,Cell_constructor>::_clean_values()
// {
// 	if (erasedValues_.empty()) return;

// 	Column_type newColumn;
// 	newColumn.reserve(column_.size());

// 	for (Cell* cell : column_){
// 		if (erasedValues_.find(cell->get_row_index()) == erasedValues_.end())
// 			newColumn.push_back(cell);
// 		else
// 			_delete_cell(cell);
// 	}
// 	erasedValues_.clear();
// 	column_.swap(newColumn);
// }

template<class Master_matrix, class Cell_constructor>
inline void Vector_column<Master_matrix,Cell_constructor>::_delete_cell(Cell* cell)
{
	if constexpr (Master_matrix::Option_list::has_row_access) ra_opt::unlink(cell);
	cellPool_->destroy(cell);
}

template<class Master_matrix, class Cell_constructor>
inline void Vector_column<Master_matrix,Cell_constructor>::_insert_cell(
		const Field_element_type &value, id_index rowIndex, Column_type &column)
{
	if constexpr (Master_matrix::Option_list::has_row_access){
		Cell *new_cell = cellPool_->construct(ra_opt::columnIndex_, rowIndex);
		new_cell->set_element(value);
		column.push_back(new_cell);
		ra_opt::insert_cell(rowIndex, new_cell);
	} else {
		Cell *new_cell = cellPool_->construct(rowIndex);
		new_cell->set_element(value);
		column.push_back(new_cell);
	}
}

template<class Master_matrix, class Cell_constructor>
inline void Vector_column<Master_matrix,Cell_constructor>::_insert_cell(
		id_index rowIndex, Column_type &column)
{
	if constexpr (Master_matrix::Option_list::has_row_access){
		Cell *new_cell = cellPool_->construct(ra_opt::columnIndex_, rowIndex);
		column.push_back(new_cell);
		ra_opt::insert_cell(rowIndex, new_cell);
	} else {
		Cell *new_cell = cellPool_->construct(rowIndex);
		column.push_back(new_cell);
	}
}

template<class Master_matrix, class Cell_constructor>
inline void Vector_column<Master_matrix,Cell_constructor>::_update_cell(
		const Field_element_type &value, id_index rowIndex, index position)
{
	if constexpr (Master_matrix::Option_list::has_row_access){
		Cell *new_cell = cellPool_->construct(ra_opt::columnIndex_, rowIndex);
		new_cell->set_element(value);
		column_[position] = new_cell;
		ra_opt::insert_cell(rowIndex, new_cell);
	} else {
		column_[position] = cellPool_->construct(rowIndex);
		column_[position]->set_element(value);
	}
}

template<class Master_matrix, class Cell_constructor>
inline void Vector_column<Master_matrix,Cell_constructor>::_update_cell(
		id_index rowIndex, index position)
{
	if constexpr (Master_matrix::Option_list::has_row_access){
		Cell *new_cell = cellPool_->construct(ra_opt::columnIndex_, rowIndex);
		column_[position] = new_cell;
		ra_opt::insert_cell(rowIndex, new_cell);
	} else {
		column_[position] = cellPool_->construct(rowIndex);
	}
}

template<class Master_matrix, class Cell_constructor>
template<class Cell_range>
inline bool Vector_column<Master_matrix,Cell_constructor>::_add(const Cell_range &column)
{
	if (column.begin() == column.end()) return false;
	if (column_.empty()){	//chain should never enter here.
		column_.resize(column.size());
		index i = 0;
		for (const Cell& cell : column){
			if constexpr (Master_matrix::Option_list::is_z2){
				_update_cell(cell.get_row_index(), i++);
			} else {
				_update_cell(cell.get_element(), cell.get_row_index(), i++);
			}
		}
		return true;
	}

	auto itTarget = column_.begin();
	auto itSource = column.begin();
	bool pivotIsZeroed = false;
	Column_type newColumn;
	newColumn.reserve(column_.size() + column.size());	//safe upper bound

	while (itTarget != column_.end() && itSource != column.end())
	{
		if constexpr (!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type){
			while (itTarget != column_.end() && erasedValues_.find((*itTarget)->get_row_index()) != erasedValues_.end()){
				_delete_cell(*itTarget);
				++itTarget;
			}
			if (itTarget == column_.end()) break;

			if constexpr (std::is_same_v<Cell_range, Vector_column<Master_matrix,Cell_constructor> >){
				while (itSource != column.end() && column.erasedValues_.find(itSource->get_row_index()) != column.erasedValues_.end()) ++itSource;
				if (itSource == column.end()) break;
			}
		}

		Cell* cellTarget = *itTarget;
		const Cell& cellSource = *itSource;
		if (cellTarget->get_row_index() < cellSource.get_row_index()) {
			newColumn.push_back(cellTarget);
			++itTarget;
		} else if (cellTarget->get_row_index() > cellSource.get_row_index()) {
			if constexpr (Master_matrix::Option_list::is_z2){
				_insert_cell(cellSource.get_row_index(), newColumn);
			} else {
				_insert_cell(cellSource.get_element(), cellSource.get_row_index(), newColumn);
			}
			++itSource;
		} else {
			if constexpr (Master_matrix::Option_list::is_z2){
				if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type){
					if (cellTarget->get_row_index() == chain_opt::get_pivot()) pivotIsZeroed = true;
				}
				_delete_cell(cellTarget);
			} else {
				cellTarget->get_element() = operators_->add(cellTarget->get_element(), cellSource.get_element());
				if (cellTarget->get_element() == Field_operators::get_additive_identity()){
					if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type){
						if (cellTarget->get_row_index() == chain_opt::get_pivot()) pivotIsZeroed = true;
					}
					_delete_cell(cellTarget);
				} else {
					newColumn.push_back(cellTarget);
					if constexpr (Master_matrix::Option_list::has_row_access)
						ra_opt::update_cell(**itTarget);
				}
			}
			++itTarget;
			++itSource;
		}
	}

	while (itSource != column.end()) {
		if constexpr (std::is_same_v<Cell_range, Vector_column<Master_matrix,Cell_constructor> > && (!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type)){
			while (itSource != column.end() && column.erasedValues_.find(itSource->get_row_index()) != column.erasedValues_.end()) ++itSource;
			if (itSource == column.end()) break;
		}
		if constexpr (Master_matrix::Option_list::is_z2){
			_insert_cell(itSource->get_row_index(), newColumn);
		} else {
			_insert_cell(itSource->get_element(), itSource->get_row_index(), newColumn);
		}
		++itSource;
	}

	while (itTarget != column_.end()){
		if constexpr (!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type){
			while (itTarget != column_.end() && erasedValues_.find((*itTarget)->get_row_index()) != erasedValues_.end()){
				_delete_cell(*itTarget);
				++itTarget;
			}
			if (itTarget == column_.end()) break;
		}
		newColumn.push_back(*itTarget);
		itTarget++;
	}

	column_.swap(newColumn);
	erasedValues_.clear();

	return pivotIsZeroed;
}

template<class Master_matrix, class Cell_constructor>
template<class Cell_range>
inline bool Vector_column<Master_matrix,Cell_constructor>::_multiply_and_add(const Field_element_type& val, const Cell_range& column)
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
	if (column_.empty()){	//chain should never enter here.
		column_.resize(column.size());
		index i = 0;
		for (const Cell& cell : column){
			if constexpr (Master_matrix::Option_list::is_z2){
				_update_cell(cell.get_row_index(), i++);
			} else {
				_update_cell(cell.get_element(), cell.get_row_index(), i++);
			}
		}
		return true;
	}

	Column_type newColumn;
	newColumn.reserve(column_.size() + column.size());	//safe upper bound

	auto itTarget = column_.begin();
	auto itSource = column.begin();
	while (itTarget != column_.end() && itSource != column.end())
	{
		if constexpr (!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type){
			while (itTarget != column_.end() && erasedValues_.find((*itTarget)->get_row_index()) != erasedValues_.end()){
				_delete_cell(*itTarget);
				++itTarget;
			}
			if (itTarget == column_.end()) break;

			if constexpr (std::is_same_v<Cell_range, Vector_column<Master_matrix,Cell_constructor> >){
				while (itSource != column.end() && column.erasedValues_.find(itSource->get_row_index()) != column.erasedValues_.end()) ++itSource;
				if (itSource == column.end()) break;
			}
		}

		Cell* cellTarget = *itTarget;
		const Cell& cellSource = *itSource;
		if (cellTarget->get_row_index() < cellSource.get_row_index()) {
			cellTarget->get_element() = operators_->multiply(cellTarget->get_element(), val);
			if constexpr (Master_matrix::Option_list::has_row_access)
				ra_opt::update_cell(**itTarget);
			newColumn.push_back(cellTarget);
			++itTarget;
		} else if (cellTarget->get_row_index() > cellSource.get_row_index()) {
			_insert_cell(cellSource.get_element(), cellSource.get_row_index(), newColumn);
			++itSource;
		} else {
			cellTarget->get_element() = operators_->multiply_and_add(cellTarget->get_element(), val, cellSource.get_element());
			if (cellTarget->get_element() == Field_operators::get_additive_identity()){
				if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type){
					if (cellTarget->get_row_index() == chain_opt::get_pivot()) pivotIsZeroed = true;
				}
				_delete_cell(cellTarget);
			} else {
				if constexpr (Master_matrix::Option_list::has_row_access)
					ra_opt::update_cell(**itTarget);
				newColumn.push_back(cellTarget);
			}
			++itTarget;
			++itSource;
		}
	}

	while (itTarget != column_.end()){
		if constexpr (!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type){
			while (itTarget != column_.end() && erasedValues_.find((*itTarget)->get_row_index()) != erasedValues_.end()){
				_delete_cell(*itTarget);
				++itTarget;
			}
			if (itTarget == column_.end()) break;
		}

		(*itTarget)->get_element() = operators_->multiply((*itTarget)->get_element(), val);
		if constexpr (Master_matrix::Option_list::has_row_access)
			ra_opt::update_cell(**itTarget);
		newColumn.push_back(*itTarget);
		itTarget++;
	}

	while (itSource != column.end()) {
		if constexpr (std::is_same_v<Cell_range, Vector_column<Master_matrix,Cell_constructor> > && (!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type)){
			while (itSource != column.end() && column.erasedValues_.find(itSource->get_row_index()) != column.erasedValues_.end()) ++itSource;
			if (itSource == column.end()) break;
		}

		_insert_cell(itSource->get_element(), itSource->get_row_index(), newColumn);
		++itSource;
	}

	column_.swap(newColumn);
	erasedValues_.clear();

	return pivotIsZeroed;
}

template<class Master_matrix, class Cell_constructor>
template<class Cell_range>
inline bool Vector_column<Master_matrix,Cell_constructor>::_multiply_and_add(const Cell_range& column, const Field_element_type& val)
{
	if (val == 0u || column.begin() == column.end()) {
		return false;
	}

	bool pivotIsZeroed = false;
	Column_type newColumn;
	newColumn.reserve(column_.size() + column.size());	//safe upper bound

	auto itTarget = column_.begin();
	auto itSource = column.begin();
	while (itTarget != column_.end() && itSource != column.end())
	{
		if constexpr (!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type){
			while (itTarget != column_.end() && erasedValues_.find((*itTarget)->get_row_index()) != erasedValues_.end()){
				_delete_cell(*itTarget);
				++itTarget;
			}
			if (itTarget == column_.end()) break;

			if constexpr (std::is_same_v<Cell_range, Vector_column<Master_matrix,Cell_constructor> >){
				while (itSource != column.end() && column.erasedValues_.find(itSource->get_row_index()) != column.erasedValues_.end()) ++itSource;
				if (itSource == column.end()) break;
			}
		}

		Cell* cellTarget = *itTarget;
		const Cell& cellSource = *itSource;
		if (cellTarget->get_row_index() < cellSource.get_row_index()) {
			newColumn.push_back(cellTarget);
			++itTarget;
		} else if (cellTarget->get_row_index() > cellSource.get_row_index()) {
			_insert_cell(operators_->multiply(cellSource.get_element(), val), cellSource.get_row_index(), newColumn);
			++itSource;
		} else {
			cellTarget->get_element() = operators_->multiply_and_add(cellSource.get_element(),  val, cellTarget->get_element());
			if (cellTarget->get_element() == Field_operators::get_additive_identity()){
				if constexpr (Master_matrix::isNonBasic && !Master_matrix::Option_list::is_of_boundary_type){
					if (cellTarget->get_row_index() == chain_opt::get_pivot()) pivotIsZeroed = true;
				}
				_delete_cell(cellTarget);
			} else {
				if constexpr (Master_matrix::Option_list::has_row_access)
						ra_opt::update_cell(**itTarget);
				newColumn.push_back(cellTarget);
			}
			++itTarget;
			++itSource;
		}
	}

	while (itSource != column.end()) {
		if constexpr (std::is_same_v<Cell_range, Vector_column<Master_matrix,Cell_constructor> > && (!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type)){
			while (itSource != column.end() && column.erasedValues_.find(itSource->get_row_index()) != column.erasedValues_.end()) ++itSource;
			if (itSource == column.end()) break;
		}

		_insert_cell(operators_->multiply(itSource->get_element(), val), itSource->get_row_index(), newColumn);
		++itSource;
	}

	while (itTarget != column_.end()){
		if constexpr (!Master_matrix::isNonBasic || Master_matrix::Option_list::is_of_boundary_type){
			while (itTarget != column_.end() && erasedValues_.find((*itTarget)->get_row_index()) != erasedValues_.end()){
				_delete_cell(*itTarget);
				++itTarget;
			}
			if (itTarget == column_.end()) break;
		}

		newColumn.push_back(*itTarget);
		++itTarget;
	}

	column_.swap(newColumn);
	erasedValues_.clear();

	return pivotIsZeroed;
}

} //namespace persistence_matrix
} //namespace Gudhi

template<class Master_matrix, class Cell_constructor>
struct std::hash<Gudhi::persistence_matrix::Vector_column<Master_matrix,Cell_constructor> >
{
	size_t operator()(const Gudhi::persistence_matrix::Vector_column<Master_matrix,Cell_constructor>& column) const
	{
		std::size_t seed = 0;
		for (const auto& cell : column){
			seed ^= std::hash<unsigned int>()(cell.get_row_index() * static_cast<unsigned int>(cell.get_element())) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		}
		return seed;
	}
};

#endif // PM_VECTOR_COLUMN_H

/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_COLUMN_CELL_CONSTRUCTORS_H
#define PM_COLUMN_CELL_CONSTRUCTORS_H

#include <utility>	//std::swap

#include <gudhi/Simple_object_pool.h>

namespace Gudhi {
namespace persistence_matrix {

template<class Cell>
struct New_cell_constructor{
	New_cell_constructor() {}

	template<class...U>
	Cell* construct(U&&...u) const { return new Cell(std::forward<U>(u)...); }

	void destroy(Cell* cell) const { delete cell; }

	friend void swap(New_cell_constructor& col1, New_cell_constructor& col2){}
};

template<class Cell>
struct Pool_cell_constructor{
public:
	Pool_cell_constructor() : cellPool_() {}
	Pool_cell_constructor(const Pool_cell_constructor& col) : cellPool_(col.cellPool_) {}
	Pool_cell_constructor(Pool_cell_constructor&& col) : cellPool_(std::move(col.cellPool_)) {}

	template<class...U>
	Cell* construct(U&&...u) {
		return cellPool_.construct(std::forward<U>(u)...);
	}

	void destroy(Cell* cell) {
		cellPool_.destroy(cell);
	}

	Pool_cell_constructor& operator=(const Pool_cell_constructor& other){ 
		cellPool_ = other.cellPool_;
		return *this;
	}

	friend void swap(Pool_cell_constructor& col1, Pool_cell_constructor& col2){
		std::swap(col1.cellPool_, col2.cellPool_);
	}

private:
	Simple_object_pool<Cell> cellPool_;
};

} //namespace persistence_matrix
} //namespace Gudhi

#endif // PM_COLUMN_CELL_CONSTRUCTORS_H

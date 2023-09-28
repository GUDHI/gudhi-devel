/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022-23 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_COLUMN_DIM_HOLDER_H
#define PM_COLUMN_DIM_HOLDER_H

#include <utility>	//std::swap

namespace Gudhi {
namespace persistence_matrix {

struct Dummy_dimension_holder{
	Dummy_dimension_holder(){}
	template<typename dimension_type>
	Dummy_dimension_holder([[maybe_unused]] dimension_type dim){}
	Dummy_dimension_holder([[maybe_unused]] const Dummy_dimension_holder& col){}
	Dummy_dimension_holder([[maybe_unused]] Dummy_dimension_holder&& col){}

	Dummy_dimension_holder& operator=([[maybe_unused]] const Dummy_dimension_holder& other){ return *this; }

	friend void swap([[maybe_unused]] Dummy_dimension_holder& col1, [[maybe_unused]] Dummy_dimension_holder& col2){}
};

template<class Master_matrix>
struct Column_dimension_holder{
	using dimension_type = typename Master_matrix::dimension_type;

	Column_dimension_holder() : dim_(Master_matrix::Option_list::is_of_boundary_type ? 0 : -1) {}
	Column_dimension_holder(dimension_type dim) : dim_(dim) {}
	Column_dimension_holder(const Column_dimension_holder& col) : dim_(col.dim_) {}
	Column_dimension_holder(Column_dimension_holder&& col) : dim_(std::exchange(col.dim_, -1)) {}

	dimension_type get_dimension() const { return dim_; }

	Column_dimension_holder& operator=(const Column_dimension_holder& other){ 
		dim_ = other.dim_;
		return *this;
	}

	friend void swap(Column_dimension_holder& col1, Column_dimension_holder& col2){
		std::swap(col1.dim_, col2.dim_);
	}

private:
	dimension_type dim_;
};

} //namespace persistence_matrix
} //namespace Gudhi

#endif // PM_COLUMN_DIM_HOLDER_H

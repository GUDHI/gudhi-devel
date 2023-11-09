/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022-23 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_CHAIN_COLUMN_PROP_H
#define PM_CHAIN_COLUMN_PROP_H

#include <utility>	//std::swap

namespace Gudhi {
namespace persistence_matrix {

struct Dummy_chain_properties{
	Dummy_chain_properties(){}
	Dummy_chain_properties([[maybe_unused]] int pivot){}
	Dummy_chain_properties([[maybe_unused]] int pivot, [[maybe_unused]] int pair){}
	Dummy_chain_properties([[maybe_unused]] const Dummy_chain_properties& col){}
	Dummy_chain_properties([[maybe_unused]] Dummy_chain_properties&& col){}

	Dummy_chain_properties& operator=([[maybe_unused]] const Dummy_chain_properties& other){ return *this; }

	friend void swap([[maybe_unused]] Dummy_chain_properties& col1, [[maybe_unused]] Dummy_chain_properties& col2){}
};

template<class Master_matrix>
class Chain_column_extra_properties{
public:
	using index = typename Master_matrix::index;

	Chain_column_extra_properties() : pivot_(-1), pairedColumn_(-1) {}
	Chain_column_extra_properties(int pivot) : pivot_(pivot), pairedColumn_(-1) {}
	Chain_column_extra_properties(int pivot, int pair) : pivot_(pivot), pairedColumn_(pair) {}
	Chain_column_extra_properties(const Chain_column_extra_properties& col) 
		: pivot_(col.pivot_), pairedColumn_(col.pairedColumn_) {}
	Chain_column_extra_properties(Chain_column_extra_properties&& col) 
		: pivot_(std::exchange(col.pivot_, -1)), pairedColumn_(std::exchange(col.pairedColumn_, -1)) {}

	int get_paired_chain_index() const { return pairedColumn_; }
	bool is_paired() const { return pairedColumn_ != -1; }
	void assign_paired_chain(index other_col){ pairedColumn_ = other_col; }
	void unassign_paired_chain() { pairedColumn_ = -1; };

	Chain_column_extra_properties& operator=(const Chain_column_extra_properties& other){ 
		pivot_ = other.pivot_;
		pairedColumn_ = other.pairedColumn_;
		return *this;
	}

	friend void swap(Chain_column_extra_properties& col1, Chain_column_extra_properties& col2){
		std::swap(col1.pivot_, col2.pivot_);
		std::swap(col1.pairedColumn_, col2.pairedColumn_);
	}

protected:
	int get_pivot() const { return pivot_; }
	void swap_pivots(Chain_column_extra_properties& other) { std::swap(pivot_, other.pivot_); }

private:
	int pivot_;			//simplex index associated to the chain
	int pairedColumn_;	//represents the (F, G x H) partition of the columns
};

} //namespace persistence_matrix
} //namespace Gudhi

#endif // PM_CHAIN_COLUMN_PROP_H

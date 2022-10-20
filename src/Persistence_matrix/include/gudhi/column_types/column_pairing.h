/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef COL_PAIRING_H
#define COL_PAIRING_H

#include <utility>

#include "../utilities/utilities.h"

namespace Gudhi {
namespace persistence_matrix {

class Column_pairing
{
public:
	Column_pairing();
	Column_pairing(const Column_pairing &toCopy);
	Column_pairing(Column_pairing&& other) noexcept;

	index get_paired_chain_index() const;
	bool is_paired() const;
	void assign_paired_chain(index other_col);
	void unassign_paired_chain();

	Column_pairing& operator=(Column_pairing other);
	friend void swap(Column_pairing& pairing1,
					 Column_pairing& pairing2){
		std::swap(pairing1.pairedColumn_, pairing2.pairedColumn_);
	}

protected:
	int pairedColumn_;
	static constexpr bool isActive_ = true;
};

inline Column_pairing::Column_pairing() : pairedColumn_(-1)
{}

inline Column_pairing::Column_pairing(const Column_pairing &toCopy)
	: pairedColumn_(toCopy.pairedColumn_)
{}

inline Column_pairing::Column_pairing(Column_pairing &&other) noexcept
	: pairedColumn_(std::exchange(other.pairedColumn_, 0))
{}

inline index Column_pairing::get_paired_chain_index() const
{
	return pairedColumn_;
}

inline bool Column_pairing::is_paired() const
{
	return pairedColumn_ != -1;
}

inline void Column_pairing::assign_paired_chain(index other_col)
{
	pairedColumn_ = other_col;
}

inline void Column_pairing::unassign_paired_chain()
{
	pairedColumn_ = -1;
}

inline Column_pairing& Column_pairing::operator=(Column_pairing other)
{
	std::swap(pairedColumn_, other.pairedColumn_);
	return *this;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // COL_PAIRING_H

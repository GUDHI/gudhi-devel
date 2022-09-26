/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef UTILITIES_H
#define UTILITIES_H

#include <vector>
#include <list>
#include <limits>

#include "Zp_field.h"

namespace Gudhi {
namespace persistence_matrix {

//const double inf = std::numeric_limits<double>::infinity();
//const double negInf = -1 * inf;

using index = unsigned int;
using dimension_type = int;
//using filtration_value_type = double;
//using persistence_pair = std::pair<filtration_value_type, filtration_value_type>;

struct Bar{
	Bar() : dim(-1), birth(-1), death(-1)
	{}

	Bar(dimension_type dim, int birth, int death)
		: dim(dim), birth(birth), death(death)
	{}

	dimension_type dim;
	int birth;
	int death;
};

template<class Field_element_type>
struct CellPairComparator {
	bool operator()(const std::pair<index,Field_element_type>& p1, const std::pair<index,Field_element_type>& p2) const
	{
		return p1.first < p2.first;
	};
};

} //namespace persistence_matrix
} //namespace Gudhi

template<unsigned int characteristic>
struct std::hash<std::pair<unsigned int, Gudhi::persistence_matrix::Zp_field_element<characteristic> > >
{
	size_t operator()(const std::pair<unsigned int, Gudhi::persistence_matrix::Zp_field_element<characteristic> >& p) const
	{
		return std::hash<unsigned int>()(p.first);
	}
};

template<unsigned int characteristic>
bool operator<(
		const std::pair<unsigned int, Gudhi::persistence_matrix::Zp_field_element<characteristic> >& p1,
		const std::pair<unsigned int, Gudhi::persistence_matrix::Zp_field_element<characteristic> >& p2){
	return p1.first < p2.first;
}

#endif // UTILITIES_H

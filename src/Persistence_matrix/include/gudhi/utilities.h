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

namespace Gudhi {
namespace persistence_matrix {

const double inf = std::numeric_limits<double>::infinity();
const double negInf = -1 * inf;

using index = unsigned int;
using filtration_value_type = double;
using dimension_type = int;
using persistence_pair = std::pair<filtration_value_type, filtration_value_type>;
using boundary_type = std::vector<index>;
using boundary_matrix = std::vector<boundary_type>;

//struct Bar{
//    Bar() : dim(-1), birth(-1), death(-1)
//    {}

//    Bar(dimension_type dim, int birth, int death)
//        : dim(dim), birth(birth), death(death)
//    {}

//    dimension_type dim;
//    int birth;
//    int death;
//};

//using barcode_type = std::vector<Bar>;
//using erasable_barcode_type = std::list<Bar>;

} //namespace persistence_matrix
} //namespace Gudhi

#endif // UTILITIES_H

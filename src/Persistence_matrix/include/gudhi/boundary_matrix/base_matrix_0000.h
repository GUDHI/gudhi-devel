/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef BASE_MATRIX_0000_H
#define BASE_MATRIX_0000_H

#include "../utilities.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class Base_matrix : Master_matrix::Base_swap_option, Master_matrix::Base_pairing_option{
public:
	Base_matrix() : Master_matrix::Base_swap_option(matrix_), Master_matrix::Base_pairing_option(){};

private:
	typename Master_matrix::column_container_type matrix_;
};

} //namespace persistence_matrix
} //namespace Gudhi

#endif // BASE_MATRIX_0000_H

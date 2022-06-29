/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef RU_MATRIX_0000_H
#define RU_MATRIX_0000_H

#include "../utilities.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class RU_matrix;

template<class Master_matrix>
class RU_matrix : Master_matrix::RU_pairing_option, Master_matrix::RU_vine_swap_option, Master_matrix::RU_representative_cycles_option{
public:
	RU_matrix()
		: Master_matrix::RU_pairing_option(),
		  Master_matrix::RU_vine_swap_option(reducedMatrixR_, mirrorMatrixU_),
		  Master_matrix::RU_representative_cycles_option(reducedMatrixR_, mirrorMatrixU_){};

private:
	typename Master_matrix::Base_matrix reducedMatrixR_;
	typename Master_matrix::Base_matrix mirrorMatrixU_;
};

} //namespace persistence_matrix
} //namespace Gudhi

#endif // RU_MATRIX_0000_H

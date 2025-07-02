/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CONCEPT_TOPOLOGICAL_DATA_WITH_SCALAR_PRODUCT_H_
#define CONCEPT_TOPOLOGICAL_DATA_WITH_SCALAR_PRODUCT_H_

namespace Gudhi {
namespace Persistence_representations {

/**
 * \brief The concept Topological_data_with_scalar_product describes the requirements
 * for a type to implement a container that allows computations of scalar products.
 */
class Topological_data_with_scalar_product
{
 public:
  double compute_scalar_product(const Topological_data_with_scalar_product& second);
};

}  // namespace Persistence_representations
}  // namespace Gudhi

#endif  // CONCEPT_TOPOLOGICAL_DATA_WITH_SCALAR_PRODUCT_H_

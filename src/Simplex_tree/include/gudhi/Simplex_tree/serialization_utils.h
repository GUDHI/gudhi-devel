/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef SIMPLEX_TREE_SERIALIZATION_UTILS_H_
#define SIMPLEX_TREE_SERIALIZATION_UTILS_H_

#include <cstring>  // for memcpy and std::size_t
#include <iostream>
#include <iomanip>
#include <memory>

namespace Gudhi {

namespace simplex_tree {

template<class SimplicialComplex>
std::size_t get_serialization_size(std::size_t num_simplices) {
  std::size_t vh_byte_size = sizeof(typename SimplicialComplex::Vertex_handle);
  std::size_t fv_byte_size = sizeof(typename SimplicialComplex::Filtration_value);
  return (vh_byte_size + num_simplices * (fv_byte_size + 2 * vh_byte_size));
}

template<class ArgumentType>
std::size_t serialize(char* start, ArgumentType value) {
  std::size_t arg_size = sizeof(ArgumentType);
  memcpy(start, &value, arg_size);
  return arg_size;
}

template<class ArgumentType>
std::size_t deserialize(char* start, ArgumentType& value) {
  std::size_t arg_size = sizeof(ArgumentType);
  memcpy(&value, start, arg_size);
  return arg_size;
}

}  // namespace simplex_tree

}  // namespace Gudhi

#endif  // SIMPLEX_TREE_SERIALIZATION_UTILS_H_

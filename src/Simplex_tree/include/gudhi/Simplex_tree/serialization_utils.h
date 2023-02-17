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

namespace Gudhi {

namespace simplex_tree {

template<class SimplicialComplex>
std::size_t get_serialization_size(std::size_t num_simplices) {
  const std::size_t vh_byte_size = sizeof(typename SimplicialComplex::Vertex_handle);
  const std::size_t fv_byte_size = SimplicialComplex::Options::store_filtration ?
    sizeof(typename SimplicialComplex::Filtration_value) : 0;
  return (vh_byte_size + num_simplices * (fv_byte_size + 2 * vh_byte_size));
}

/** \brief Serialize the given value and insert it at start position.
 * 
 * @param[in] value The value to serialize.
 * @param[in] start Start position where the value is serialized.
 * @return The new position in the array of char for the next serialization.
 * 
 * @warning It is the user's responsibility to provide a pointer to a buffer with enough memory space.
 */
template<class ArgumentType>
char* serialize_trivial(ArgumentType value, char* start) {
  std::size_t arg_size = sizeof(ArgumentType);
  memcpy(start, &value, arg_size);
  return start + arg_size;
}

/** \brief Deserialize at the start position in an array of char and sets the value with it.
 * 
 * @param[in] value The value where to deserialize based on its type.
 * @param[in] start Start position where the value is serialized.
 * @return The new position in the array of char for the next deserialization.
 * 
 * @warning It is the user's responsibility to ensure that the pointer will not go out of bounds.
 */
template<class ArgumentType>
const char* deserialize_trivial(ArgumentType& value, const char* start) {
  std::size_t arg_size = sizeof(ArgumentType);
  memcpy(&value, start, arg_size);
  return (start + arg_size);
}

}  // namespace simplex_tree

}  // namespace Gudhi

#endif  // SIMPLEX_TREE_SERIALIZATION_UTILS_H_

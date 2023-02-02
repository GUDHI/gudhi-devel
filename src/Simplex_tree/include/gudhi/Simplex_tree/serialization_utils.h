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

/** \addtogroup simplex_tree 
 * Serialization and deserialization utils for the Simplex_tree.
 *  @{
 */

/** \brief Computes and return the serialization size in bytes in function of the number of simplices.
 */
template<class SimplicialComplex>
std::size_t get_serialization_size(std::size_t num_simplices) {
  std::size_t vh_byte_size = sizeof(typename SimplicialComplex::Vertex_handle);
  std::size_t fv_byte_size = sizeof(typename SimplicialComplex::Filtration_value);
  return (vh_byte_size + num_simplices * (fv_byte_size + 2 * vh_byte_size));
}

/** \brief Serialize the given value at the start position in an array of char.
 * 
 * @warning It is the user resposibility to ensure that the array of char is wide enough.
 * 
 * @return The new position in the array of char.
 */
template<class ArgumentType>
char* serialize(char* start, ArgumentType value) {
  std::size_t arg_size = sizeof(ArgumentType);
  memcpy(start, &value, arg_size);
  return (start + arg_size);
}

/** \brief Deserialize at the start position in an array of char and sets the value with it.
 * 
 * @return The new position in the array of char.
 * 
 * @warning It is the user resposibility to ensure that the array of char is wide enough.
 */
template<class ArgumentType>
char* deserialize(char* start, ArgumentType& value) {
  std::size_t arg_size = sizeof(ArgumentType);
  memcpy(&value, start, arg_size);
  return (start + arg_size);
}

/** @}*/  // end addtogroup simplex_tree

}  // namespace simplex_tree

}  // namespace Gudhi

#endif  // SIMPLEX_TREE_SERIALIZATION_UTILS_H_

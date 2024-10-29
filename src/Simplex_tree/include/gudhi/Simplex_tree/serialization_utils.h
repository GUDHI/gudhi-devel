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

namespace Gudhi {

namespace simplex_tree {

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

/**
 * @brief Returns the size of the serialization of the given object.
 */
template<class ArgumentType>
constexpr std::size_t get_serialization_size_of([[maybe_unused]] ArgumentType value) {
  return sizeof(ArgumentType);
}

}  // namespace simplex_tree

}  // namespace Gudhi

#endif  // SIMPLEX_TREE_SERIALIZATION_UTILS_H_

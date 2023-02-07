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

/** \brief Serialize the given value and insert it at the end of the buffer.
 */
template<class ArgumentType>
void serialize_trivial(ArgumentType value, std::vector<char>& buffer) {
  buffer.insert(std::end(buffer),
                reinterpret_cast<const char*>(&value),
                reinterpret_cast<const char*>(&value + 1));
}

/** \brief Deserialize at the start position in an array of char and sets the value with it.
 * 
 * @return The new position in the array of char.
 * 
 * @warning It is the user resposibility to ensure that the array of char is wide enough.
 */
template<class ArgumentType>
std::vector<char>::const_iterator deserialize_trivial(std::vector<char>::const_iterator start, ArgumentType& value) {
  std::size_t arg_size = sizeof(ArgumentType);
  // TODO: Not really nice, but I didn't manage to do it with std::copy.
  memcpy(&value, &*start, arg_size);
  return (start + arg_size);
}

/** @}*/  // end addtogroup simplex_tree

}  // namespace simplex_tree

}  // namespace Gudhi

#endif  // SIMPLEX_TREE_SERIALIZATION_UTILS_H_

/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - 2026/08 Hannah Schreiber: Move to common + addition of versions for vectors and pairs
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef GUDHI_SERIALIZATION_UTILS_H_
#define GUDHI_SERIALIZATION_UTILS_H_

#include <cstring>  // for memcpy and std::size_t
#include <vector>
#include <utility>

namespace Gudhi {


/********************************************************/
/***** Signatures to enable recursion on the types ******/
/********************************************************/


// /**
//  * @private
//  * @brief Serialize the given value and insert it at start position.
//  *
//  * @tparam T, T1, T2 Native arithmetic type for the first version and for the others: any of the types
//  * @ref serialize_value_to_char_buffer has a version for.
//  * @param value The value to serialize.
//  * @param start Start position where the value is serialized.
//  * @return The new position in the array of char for the next serialization.
//  *
//  * @warning It is the user's responsibility to provide a pointer to a buffer with enough memory space.
//  */


template <typename T, class>
inline char *serialize_value_to_char_buffer(T value, char *start);
template <typename T>
inline char *serialize_value_to_char_buffer(const std::vector<T> &range, char *start);
template <typename T1, typename T2>
inline char *serialize_value_to_char_buffer(const std::pair<T1, T2> &pair, char *start);


// /**
//  * @private
//  * @brief Deserialize at the start position in an array of char and sets the value with it.
//  *
//  * @tparam T, T1, T2 Native arithmetic type for the first version and for the others: any of the types
//  * @ref deserialize_value_from_char_buffer has a version for.
//  * @param value The value where to deserialize based on its type.
//  * @param start Start position where the value is serialized.
//  * @return The new position in the array of char for the next deserialization.
//  *
//  * @warning It is the user's responsibility to ensure that the pointer will not go out of bounds.
//  */


template <typename T, class>
inline const char *deserialize_value_from_char_buffer(T &value, const char *start);
template <typename T>
inline const char *deserialize_value_from_char_buffer(std::vector<T> &range, const char *start);
template <typename T1, typename T2>
inline const char *deserialize_value_from_char_buffer(std::pair<T1, T2> &pair, const char *start);


// /**
//  * @private
//  * @brief Returns the serialization size of given value needed by the corresponding
//  * @ref serialize_value_to_char_buffer.
//  *
//  * @tparam T, T1, T2 Native arithmetic type for the first version and for the others: any of the types
//  * @ref get_serialization_size_of has a version for.
//  * @param value Value to get the size of.
//  */


template <typename T, class>
inline std::size_t get_serialization_size_of([[maybe_unused]] T value);
template <typename T>
inline std::size_t get_serialization_size_of(const std::vector<T> &range);
template <typename T1, typename T2>
inline std::size_t get_serialization_size_of(const std::pair<T1, T2> &pair);


/****************************************************/
/***** Implementations for the different types ******/
/****************************************************/


template <typename T, class = std::enable_if_t<std::is_arithmetic_v<T>>>
inline char *serialize_value_to_char_buffer(T value, char *start) {
  const std::size_t argSize = sizeof(T);
  memcpy(start, &value, argSize);
  return start + argSize;
}

template <typename T, class = std::enable_if_t<std::is_arithmetic_v<T>>>
inline const char *deserialize_value_from_char_buffer(T &value, const char *start) {
  const std::size_t argSize = sizeof(T);
  memcpy(&value, start, argSize);
  return start + argSize;
}

template <typename T, class = std::enable_if_t<std::is_arithmetic_v<T>>>
inline std::size_t get_serialization_size_of([[maybe_unused]] T value) {
  return sizeof(T);
}

template <typename T>
inline char *serialize_value_to_char_buffer(const std::vector<T> &range, char *start) {
  const std::size_t sizeSize = sizeof(std::size_t);
  const std::size_t length = range.size();
  char *curr = start;
  memcpy(curr, &length, sizeSize);
  curr += sizeSize;
  for (const T &v : range) {
    curr = serialize_value_to_char_buffer(v, curr);
  }
  return curr;
}

template <typename T>
inline const char *deserialize_value_from_char_buffer(std::vector<T> &range, const char *start) {
  const std::size_t sizeSize = sizeof(std::size_t);
  std::size_t length;
  const char *curr = start;
  memcpy(&length, curr, sizeSize);
  curr += sizeSize;
  range.resize(length);
  for (T &v : range) {
    curr = deserialize_value_from_char_buffer(v, curr);
  }
  return curr;
}

template <typename T>
inline std::size_t get_serialization_size_of(const std::vector<T> &range) {
  std::size_t size = sizeof(std::size_t);
  for (const T &v : range) {
    size += get_serialization_size_of(v);
  }
  return size;
}

template <typename T1, typename T2>
inline char *serialize_value_to_char_buffer(const std::pair<T1, T2> &pair, char *start) {
  char *curr = start;
  curr = serialize_value_to_char_buffer(pair.first, curr);
  curr = serialize_value_to_char_buffer(pair.second, curr);
  return curr;
}

template <typename T1, typename T2>
inline const char *deserialize_value_from_char_buffer(std::pair<T1, T2> &pair, const char *start) {
  const char *curr = start;
  curr = deserialize_value_from_char_buffer(pair.first, curr);
  curr = deserialize_value_from_char_buffer(pair.second, curr);
  return curr;
}

template <typename T1, typename T2>
inline std::size_t get_serialization_size_of(const std::pair<T1, T2> &pair) {
  std::size_t size = get_serialization_size_of(pair.first);
  size += get_serialization_size_of(pair.second);
  return size;
}

}  // namespace Gudhi

#endif  // GUDHI_SERIALIZATION_UTILS_H_

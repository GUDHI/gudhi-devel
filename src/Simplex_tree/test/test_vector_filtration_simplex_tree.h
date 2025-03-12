/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef SIMPLEX_TREE_TEST_CUSTOM_SIMPLEX_TREE_H_
#define SIMPLEX_TREE_TEST_CUSTOM_SIMPLEX_TREE_H_

#include <cstddef>
#include <stdexcept>
#include <vector>
#include <limits>
#include <cstdint>

#include <gudhi/Simplex_tree.h>

namespace Gudhi {

class Vector_filtration_value : public std::vector<int>
{
  using Base = std::vector<int>;

 public:
  using const_iterator = Base::const_iterator;

  Vector_filtration_value() : Base() {}
  Vector_filtration_value(std::size_t count) : Base(count) {}
  Vector_filtration_value(std::initializer_list<int> init) : Base(init) {}
  Vector_filtration_value(const_iterator start, const_iterator end) : Base(start, end) {}
  // Vector_filtration_value(const Gudhi::simplex_tree::empty_filtration_value_t& e) : Base(0) {}

  friend bool unify_lifetimes(Vector_filtration_value& f1, const Vector_filtration_value& f2) {
    int max = std::numeric_limits<int>::max();
    bool f1_is_inf = f1.size() == 1 && f1[0] == max;
    bool f2_is_inf = f2.size() == 1 && f2[0] == max;
    if (f1_is_inf){
      if (f2_is_inf){
        return false;
      }
      f1 = f2;
      return true;
    }
    if (f2_is_inf){
      return false;
    }
    // same size otherwise
    unsigned int i = 0;
    bool modified = false;
    for (int v : f2){
      if (v < f1[i]){
        f1[i] = v;
        modified = true;
      }
      ++i;
    }
    
    return modified;
  }

  friend bool intersect_lifetimes(Vector_filtration_value& f1, const Vector_filtration_value& f2) {
    if (f1 < f2) {
      f1 = f2;
      return true;
    }
    return false;
  }

  friend char* serialize_value_to_char_buffer(const Vector_filtration_value& value, char* start)
  {
    const auto length = value.size();
    const std::size_t arg_size = sizeof(int) * length;
    const std::size_t type_size = sizeof(Vector_filtration_value::size_type);
    memcpy(start, &length, type_size);
    memcpy(start + type_size, value.data(), arg_size);
    return start + arg_size + type_size;
  }

  friend const char* deserialize_value_from_char_buffer(Vector_filtration_value& value, const char* start)
  {
    const std::size_t type_size = sizeof(Vector_filtration_value::size_type);
    Vector_filtration_value::size_type length;
    memcpy(&length, start, type_size);
    std::size_t arg_size = sizeof(int) * length;
    value.resize(length);
    memcpy(value.data(), start + type_size, arg_size);
    return start + arg_size + type_size;
  }

  friend std::size_t get_serialization_size_of(const Vector_filtration_value& value) {
    return sizeof(Vector_filtration_value::size_type) + sizeof(int) * value.size();
  }

  friend std::ostream &operator<<(std::ostream &stream, const Vector_filtration_value &f) {
    if (f.empty()) {
      stream << "[]";
      return stream;
    }
    stream << "[";
    for (std::size_t i = 0; i < f.size() - 1; i++) {
      stream << f[i] << ", ";
    }
    stream << f.back();
    stream << "]";
    return stream;
  }
};

struct Simplex_tree_options_custom_fil_values_default {
  typedef linear_indexing_tag Indexing_tag;
  typedef std::int16_t Vertex_handle;
  typedef Vector_filtration_value Filtration_value;
  typedef std::int32_t Simplex_key;
  static const bool store_key = false;
  static const bool store_filtration = true;
  static const bool contiguous_vertices = false;
  static const bool link_nodes_by_label = false;
  static const bool stable_simplex_handles = false;
};

struct Simplex_tree_options_custom_fil_values_fast_persistence {
  typedef linear_indexing_tag Indexing_tag;
  typedef std::int16_t Vertex_handle;
  typedef Vector_filtration_value Filtration_value;
  typedef std::int32_t Simplex_key;
  static const bool store_key = true;
  static const bool store_filtration = true;
  static const bool contiguous_vertices = true;
  static const bool link_nodes_by_label = false;
  static const bool stable_simplex_handles = false;
};

struct Simplex_tree_options_custom_fil_values_full_featured {
  typedef linear_indexing_tag Indexing_tag;
  typedef std::int16_t Vertex_handle;
  typedef Vector_filtration_value Filtration_value;
  typedef std::int32_t Simplex_key;
  static const bool store_key = true;
  static const bool store_filtration = true;
  static const bool contiguous_vertices = false;
  static const bool link_nodes_by_label = true;
  static const bool stable_simplex_handles = true;
};

}  // namespace Gudhi

namespace std {

template<>
class numeric_limits<Gudhi::Vector_filtration_value> {
 public:
  static constexpr bool has_infinity = true;

  static Gudhi::Vector_filtration_value infinity() noexcept {
    return {std::numeric_limits<int>::max()};
  };

  static Gudhi::Vector_filtration_value max() noexcept(false) {
    throw std::logic_error(
        "The maximal value cannot be represented with no finite numbers of parameters.");
  };
};

}  // namespace std

#endif  // SIMPLEX_TREE_TEST_CUSTOM_SIMPLEX_TREE_H_

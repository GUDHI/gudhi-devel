/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef INCLUDE_FILTRATION_VALUE_TYPES_SUPPORT_H_
#define INCLUDE_FILTRATION_VALUE_TYPES_SUPPORT_H_

#include <string>
#include <type_traits>  // for std::is_same_v

namespace Gudhi {

  // Developper notice:
  // If you need to add a new type:
  // 1. you have to add a condition to get_class_name_for_filtration_type() for the python class name to be suffixed
  //    with the type as the class names must be different.
  // 2. you need to add the new type in Supported::List
  // 
  // Supported_for_data_structure::List is automatically updated (nothing to be done here)

  template<class Filtration_value>
  constexpr std::string get_class_name_for_filtration_type(const std::string& root_class_name)
  {
    if constexpr (std::is_same_v<Filtration_value, float>)
      return root_class_name + std::string("_float32");
    if constexpr (std::is_same_v<Filtration_value, double>)
      return root_class_name + std::string("_float64");
  }

  // To get the identity of a type
  template<class T>
  struct identity { using type = T; };

  // A list of Types
  template<class... Ts>
  struct type_list {};

  struct Supported {
    using List = type_list<float, double>;
  };

  // Data_structure must accept Filtration_value as a template parameter (Data_structure<Filtration_value> shall compile)
  // Then for_each_filtration_value_type(Supported_for_data_structure<Data_structure>::List>,my_lambda);
  // will call my_lambda for all Data_structure<Filtration_value>, where Filtration_value is one of
  // Supported
  // 
  // Equivalent to:
  //
  // template<template<typename> class Data_structure>
  // struct Supported_for_data_structure {
  //   using List = type_list<Data_structure<float>, Data_structure<double>>;
  // };
  template<template<typename> class Data_structure>
  struct Supported_for_data_structure {
   private:
    // metafunction (partial specialization): type_list<Ts...> -> type_list<Data_structure<Ts>...>
    template<class L> struct wrap;
    template<class... Ts>
    struct wrap<type_list<Ts...>> {
      using type = type_list<Data_structure<Ts>...>;
    };

   public:
    using List = typename wrap<typename Supported::List>::type;
  };

  // Calls f(identity<T>{}) for every T in the list
  template<class... Ts, class F>
  void for_each_filtration_value_type(type_list<Ts...>, F&& f) {
    (f(identity<Ts>{}), ...);
  }

}  // namespace Gudhi

#endif  // INCLUDE_FILTRATION_VALUE_TYPES_SUPPORT_H_

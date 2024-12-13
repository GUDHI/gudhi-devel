/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_COMMON_BOOST_TYPE_LISTS_H
#define PM_COMMON_BOOST_TYPE_LISTS_H

#include <boost/mp11.hpp>

#include <gudhi/persistence_matrix_options.h>

using Gudhi::persistence_matrix::Column_types;

struct ct_intrusive_list {
  static constexpr const Column_types t = Column_types::INTRUSIVE_LIST;
};

struct ct_intrusive_set {
  static constexpr const Column_types t = Column_types::INTRUSIVE_SET;
};

struct ct_list {
  static constexpr const Column_types t = Column_types::LIST;
};

struct ct_set {
  static constexpr const Column_types t = Column_types::SET;
};

struct ct_heap {
  static constexpr const Column_types t = Column_types::HEAP;
};

struct ct_unordered_set {
  static constexpr const Column_types t = Column_types::UNORDERED_SET;
};

struct ct_vector {
  static constexpr const Column_types t = Column_types::VECTOR;
};

struct ct_naive_vector {
  static constexpr const Column_types t = Column_types::NAIVE_VECTOR;
};

struct ct_small_vector {
  static constexpr const Column_types t = Column_types::SMALL_VECTOR;
};

struct true_value {
  static constexpr const bool t = true;
};

struct false_value {
  static constexpr const bool t = false;
};

template <class T>
using get_type = typename T::type;

template <template <typename...> class... F>
using mp_list_q = boost::mp11::mp_list<boost::mp11::mp_quote<F>...>;

using bool_value_list = boost::mp11::mp_list<true_value, false_value>;
using true_value_list = boost::mp11::mp_list<true_value>;
using false_value_list = boost::mp11::mp_list<false_value>;

#endif  // PM_COMMON_BOOST_TYPE_LISTS_H

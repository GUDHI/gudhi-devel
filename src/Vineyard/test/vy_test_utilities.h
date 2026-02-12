/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef VY_TEST_UTILITIES_H
#define VY_TEST_UTILITIES_H

#include <tuple>
#include <vector>
#include <limits>

#include <gudhi/vineyard_base.h>

using BC = std::vector<std::vector<int>>;
using DC = std::vector<int>;
template<typename T>
using FC = std::vector<T>;

struct Chain_vineyard_options : Gudhi::vineyard::Default_vineyard_options {
  static constexpr bool is_RU = false;
};

struct RU_vineyard_options : Gudhi::vineyard::Default_vineyard_options {
  static constexpr bool is_RU = true;
};

template<typename T = double>
inline std::tuple<BC, DC, FC<T> > build_simple_input_complex()
{
  T inf = std::numeric_limits<T>::has_infinity ? std::numeric_limits<T>::infinity() : std::numeric_limits<T>::max();

  BC bc = {{}, {}, {}, {0, 1}, {1, 2}, {0, 2}, {3, 4, 5}, {}, {1, 7}};
  DC dc = {0, 0, 0, 1, 1, 1, 2, 0, 1};
  FC<T> fc = {1, 2, 1, 6, 4, 3, inf, inf, inf};

  return std::tuple<BC, DC, FC<T> >(bc, dc, fc);
}

template <class Bar, class Barcode>
std::vector<Bar> get_barcode(const Barcode& bc)
{
  std::vector<Bar> barcode(bc.begin(), bc.end());
  std::sort(barcode.begin(), barcode.end(), [](const Bar& b1, const Bar& b2) {
    if (b1.dim == b2.dim) return b1.birth < b2.birth;
    return b1.dim < b2.dim;
  });
  return barcode;
}

template <class Cycle>
std::vector<int> get_cycle(const Cycle& c)
{
  std::vector<int> cycle(c.begin(), c.end());
  std::sort(cycle.begin(), cycle.end());
  return cycle;
}

template <class Cycles>
std::vector<std::vector<int>> get_all_cycles(const Cycles& cs)
{
  std::vector<std::vector<int>> cycles(cs.size());
  unsigned int i = 0;
  for (const auto& c : cs) {
    cycles[i] = get_cycle(c);
    ++i;
  }
  std::sort(cycles.begin(), cycles.end());
  return cycles;
}

template <class V>
std::vector<std::vector<int>> get_all_cycles_individually(V& vy)
{
  std::vector<std::vector<int>> cycles(vy.get_current_barcode().size());
  for (unsigned int i = 0; i < cycles.size(); ++i) {
    cycles[i] = get_cycle(vy.get_current_representative_cycle(i));
  }
  std::sort(cycles.begin(), cycles.end());
  return cycles;
}

#endif  // VY_TEST_UTILITIES_H

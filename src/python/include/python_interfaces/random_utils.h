/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef INCLUDE_RANDOM_UTILS_PYTHON_H_
#define INCLUDE_RANDOM_UTILS_PYTHON_H_

#include <iostream>
#include <cstdint>

#include <nanobind/nanobind.h>  // for nb::capsule

#include <numpy/random/bitgen.h>  // for bitgen_t

#include <CGAL/Random.h>  // for CGAL::get_default_random()
#include <thread>

namespace nb = nanobind;

namespace Gudhi {
namespace random {

class RandomGenerator {
 public:
  RandomGenerator() { std::cout << "CGAL::get_default_random: " << &(CGAL::get_default_random()) << std::endl; }
  RandomGenerator(uint64_t seed) { CGAL::get_default_random() = CGAL::Random(seed); std::cout << "seed = " << seed << " - thread id " << std::this_thread::get_id() << " - CGAL::get_default_random: " << &(CGAL::get_default_random()) << std::endl; }
  RandomGenerator(RandomGenerator& other) = delete;
  RandomGenerator(RandomGenerator&& other) = delete;
  RandomGenerator& operator=(const RandomGenerator& other) = delete;
  RandomGenerator& operator=(RandomGenerator&& other) = delete;
  
  CGAL::Random* get_default_random() { return &(CGAL::get_default_random()); }
  uint32_t next_uint32() { return CGAL::get_default_random().get_int(0, UINT32_MAX); }
  uint64_t next_uint64() { return CGAL::get_default_random().get_int(0, UINT64_MAX); }
  double   next_double() { return CGAL::get_default_random().get_double(); }
};

uint32_t next_uint32(void* st) {
  return (uint32_t) static_cast<RandomGenerator*>(st)->next_uint32();
}

uint64_t next_uint64(void* st) {
  return (uint64_t) static_cast<RandomGenerator*>(st)->next_uint64();
}

double next_double(void* st) {
  return (double) static_cast<RandomGenerator*>(st)->next_double();
}

void setup_bitgen(RandomGenerator* rng, nb::capsule capsule) {
  bitgen_t* bg = static_cast<bitgen_t*>(capsule.data());
  bg->state = rng;
  bg->next_uint64 = next_uint64;
  bg->next_uint32 = next_uint32;
  bg->next_double = next_double;
  bg->next_raw    = next_uint64;
};

}  // namespace random
}  // namespace Gudhi

#endif  // INCLUDE_RANDOM_UTILS_PYTHON_H_

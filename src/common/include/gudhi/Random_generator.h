/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef RANDOM_GENERATOR_H_
#define RANDOM_GENERATOR_H_

#include <cstdint>
#include <random>

// Must be done before #include <CGAL/Random.h> - Requires also to compile without -fvisibility=hidden
namespace CGAL {
  class Random;
#ifdef _WIN32
  #if defined RANDOM_DLL_IMPORT
    extern __declspec(dllimport) Random& get_default_random();
  #endif
  #if defined RANDOM_DLL_EXPORT
    __declspec(dllexport) Random& get_default_random();
  #endif
#else
  Random& get_default_random() __attribute__((visibility("default")));
#endif
}

#include <CGAL/Random.h>  // for CGAL::get_default_random()

namespace Gudhi {
namespace random {

// Let's hide CGAL::Random
using Random = CGAL::Random;

class Random_generator {
 public:
  // CGAL::Random() default constructor sets the seed with std::time which is not convenient.
  // Sets it with std::random_device is more advised
  Random_generator() {
    unsigned int random_seed = rand_dev();
    CGAL::get_default_random() = Random(random_seed);
  }
  Random_generator(unsigned int seed) { CGAL::get_default_random() = Random(seed); }
  Random_generator(Random_generator& other) = delete;
  Random_generator(Random_generator&& other) = delete;
  Random_generator& operator=(const Random_generator& other) = delete;
  Random_generator& operator=(Random_generator&& other) = delete;
  
  Random* get_default_random() { return &(CGAL::get_default_random()); }
  uint32_t next_uint32() { return CGAL::get_default_random().uniform_int<uint32_t>(0, UINT32_MAX); }
  uint64_t next_uint64() { return CGAL::get_default_random().uniform_int<uint64_t>(0, UINT64_MAX); }
  double   next_double() { return CGAL::get_default_random().uniform_real<double>(); }

 private:
  std::random_device rand_dev;
};

uint32_t next_uint32(void* st) {
  return (uint32_t) static_cast<Random_generator*>(st)->next_uint32();
}

uint64_t next_uint64(void* st) {
  return (uint64_t) static_cast<Random_generator*>(st)->next_uint64();
}

double next_double(void* st) {
  return (double) static_cast<Random_generator*>(st)->next_double();
}

}  // namespace random
}  // namespace Gudhi

#endif  // RANDOM_GENERATOR_H_

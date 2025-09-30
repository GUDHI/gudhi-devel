/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <cstdint>
#include <nanobind/nanobind.h>

#include <numpy/random/bitgen.h>

#include <CGAL/Random.h>

namespace nb = nanobind;

namespace Gudhi {
namespace random {

class RandomGenerator {
 public:
  RandomGenerator() : engine() {}
  RandomGenerator(uint64_t seed) : engine(seed) { CGAL::get_default_random() = CGAL::Random(seed); }

  uint32_t next_uint32() { return engine.get_int(0, UINT32_MAX); }
  uint64_t next_uint64() { return engine.get_int(0, UINT64_MAX); }
  double   next_double() { return engine.get_double(); }

 private:
  CGAL::Random engine;
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


NB_MODULE(_random_ext, m) {
  // Based on https://doc.cgal.org/latest/Generator/index.html
  m.attr("__license__") = "LGPL v3";
  
  nb::class_<Gudhi::random::RandomGenerator> (m, "RandomGenerator")
      .def(nb::init<>())
      .def(nb::init<long>())
      .def("setup_bitgen", &Gudhi::random::setup_bitgen);
}

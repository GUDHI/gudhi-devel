/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <nanobind/nanobind.h>
#include <numpy/random/bitgen.h>  // for bitgen_t

// For Windows, where Random_generator is a client of random here (random.dll is the provider)
#define RANDOM_DLL_IMPORT

#include <gudhi/Random_generator.h>

namespace nb = nanobind;

namespace Gudhi {
namespace random {

void setup_bitgen(Random_generator* rng, nb::capsule capsule) {
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
  
  nb::class_<Gudhi::random::Random_generator> (m, "GudhiRandomGenerator")
      .def(nb::init<>())
      .def(nb::init<unsigned int>())
      .def("setup_bitgen", &Gudhi::random::setup_bitgen);
}

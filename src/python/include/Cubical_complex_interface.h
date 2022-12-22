/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef INCLUDE_CUBICAL_COMPLEX_INTERFACE_H_
#define INCLUDE_CUBICAL_COMPLEX_INTERFACE_H_

#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Bitmap_cubical_complex_base.h>
#include <gudhi/Bitmap_cubical_complex_periodic_boundary_conditions_base.h>

#include <iostream>
#include <vector>
#include <string>

namespace Gudhi {

namespace cubical_complex {

/* These classes are only needed to circumvent the protected access. */

class Cubical_complex_interface : public Bitmap_cubical_complex<Bitmap_cubical_complex_base<double>> {
  typedef Bitmap_cubical_complex<Bitmap_cubical_complex_base<double>> Base;
 public:
  using Base::Base; // inheriting constructors
  using Base::data;

  // not const because cython does not handle const very well
  std::vector<unsigned>& shape() { return this->sizes; };
};

class Periodic_cubical_complex_interface : public Bitmap_cubical_complex<Bitmap_cubical_complex_periodic_boundary_conditions_base<double>> {
  typedef Bitmap_cubical_complex<Bitmap_cubical_complex_periodic_boundary_conditions_base<double>> Base;
 public:
  using Base::Base; // inheriting constructors
  using Base::data;

  // not const because cython does not handle const very well
  std::vector<unsigned>& shape() { return this->sizes; };
  std::vector<bool>& periodicities() { return this->directions_in_which_periodic_b_cond_are_to_be_imposed; }
};

}  // namespace cubical_complex

}  // namespace Gudhi

#endif  // INCLUDE_CUBICAL_COMPLEX_INTERFACE_H_


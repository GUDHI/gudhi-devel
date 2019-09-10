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

template<typename CubicalComplexOptions = Bitmap_cubical_complex_base<double>>
class Cubical_complex_interface : public Bitmap_cubical_complex<CubicalComplexOptions> {
 public:
  Cubical_complex_interface(const std::vector<unsigned>& dimensions,
                            const std::vector<double>& top_dimensional_cells)
  : Bitmap_cubical_complex<CubicalComplexOptions>(dimensions, top_dimensional_cells) {
  }

  Cubical_complex_interface(const std::vector<unsigned>& dimensions,
                            const std::vector<double>& top_dimensional_cells,
                            const std::vector<bool>& periodic_dimensions)
  : Bitmap_cubical_complex<CubicalComplexOptions>(dimensions, top_dimensional_cells, periodic_dimensions) {
  }

  Cubical_complex_interface(const std::string& perseus_file)
  : Bitmap_cubical_complex<CubicalComplexOptions>(perseus_file.c_str()) {
  }
};

}  // namespace cubical_complex

}  // namespace Gudhi

#endif  // INCLUDE_CUBICAL_COMPLEX_INTERFACE_H_


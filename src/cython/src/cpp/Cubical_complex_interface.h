/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016  INRIA Saclay (France)
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CUBICAL_COMPLEX_INTERFACE_H
#define	CUBICAL_COMPLEX_INTERFACE_H

#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Bitmap_cubical_complex_base.h>
#include <gudhi/Bitmap_cubical_complex_periodic_boundary_conditions_base.h>

#include <vector>
#include <utility>  // std::pair
#include <iostream>

namespace Gudhi {

namespace Cubical_complex {

template<typename CubicalComplexOptions = Bitmap_cubical_complex_base<double>>
class Cubical_complex_interface : public Bitmap_cubical_complex<CubicalComplexOptions> {
 public:

  Cubical_complex_interface(const std::vector<unsigned>& dimensions,
                            const std::vector<double>& top_dimensional_cells)
  : Bitmap_cubical_complex<CubicalComplexOptions>(dimensions, top_dimensional_cells) {
  }

};

}  // namespace Cubical_complex

}  // namespace Gudhi

#endif  // CUBICAL_COMPLEX_INTERFACE_H


/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Mathieu Carriere
 *
 *    Copyright (C) 2018 INRIA
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

#ifndef INCLUDE_KERNELS_INTERFACE_H_
#define INCLUDE_KERNELS_INTERFACE_H_

#include <gudhi/Sliced_Wasserstein.h>

#include <iostream>
#include <vector>
#include <utility>  // for std::pair

namespace Gudhi {

namespace persistence_diagram {

  double sw(const std::vector<std::pair<double, double>>& diag1,
            const std::vector<std::pair<double, double>>& diag2,
            double sigma, int N) {
    Gudhi::Persistence_representations::Sliced_Wasserstein sw1(diag1, sigma, N);
    Gudhi::Persistence_representations::Sliced_Wasserstein sw2(diag2, sigma, N);
    return sw1.compute_scalar_product(sw2);
  }

}  // namespace persistence_diagram

}  // namespace Gudhi


#endif  // INCLUDE_KERNELS_INTERFACE_H_

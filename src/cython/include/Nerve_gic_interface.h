/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2018 Inria
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

#ifndef INCLUDE_NERVE_GIC_INTERFACE_H_
#define INCLUDE_NERVE_GIC_INTERFACE_H_

#include <gudhi/Simplex_tree.h>
#include <gudhi/distance_functions.h>
#include <gudhi/GIC.h>

#include "Simplex_tree_interface.h"

#include <iostream>
#include <vector>
#include <string>

namespace Gudhi {

namespace cover_complex {

class Nerve_gic_interface : public Cover_complex<std::vector<double>> {
 public:
  void create_simplex_tree(Simplex_tree_interface<>* simplex_tree) {
    create_complex(*simplex_tree);
    simplex_tree->initialize_filtration();
  }
  void set_cover_from_Euclidean_Voronoi(int m) {
    set_cover_from_Voronoi(Gudhi::Euclidean_distance(), m);
  }
};

}  // namespace cover_complex

}  // namespace Gudhi

#endif  // INCLUDE_NERVE_GIC_INTERFACE_H_

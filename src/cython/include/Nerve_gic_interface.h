/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2018 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
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
  double set_graph_from_automatic_euclidean_rips(int N) {
    return set_graph_from_automatic_rips(Gudhi::Euclidean_distance(), N);
  }
  void set_graph_from_euclidean_rips(double threshold) {
    set_graph_from_rips(threshold, Gudhi::Euclidean_distance());
  }
};

}  // namespace cover_complex

}  // namespace Gudhi

#endif  // INCLUDE_NERVE_GIC_INTERFACE_H_

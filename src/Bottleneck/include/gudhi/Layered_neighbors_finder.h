/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Francois Godi
 *
 *    Copyright (C) 2015  INRIA Saclay (France)
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

#ifndef SRC_BOTTLENECK_INCLUDE_GUDHI_LAYERED_NEIGHBORS_FINDER_H_
#define SRC_BOTTLENECK_INCLUDE_GUDHI_LAYERED_NEIGHBORS_FINDER_H_

#include <vector>

#include "Neighbors_finder.h"

// Layered_neighbors_finder is a data structure used to find if a query point from U has neighbors in V in a given
// vlayer of the vlayered persistence diagrams graph. V's points have to be added manually using their index.
// A neighbor returned is automatically removed.

namespace Gudhi {

namespace bottleneck {

class Layered_neighbors_finder {
 public:
  Layered_neighbors_finder(const Persistence_diagrams_graph& g, double r);
  void add(int v_point_index, int vlayer);
  int pull_near(int u_point_index, int vlayer);
  int vlayers_number() const;

 private:
  const Persistence_diagrams_graph& g;
  const double r;
  std::vector<Neighbors_finder> neighbors_finder;
};

Layered_neighbors_finder::Layered_neighbors_finder(const Persistence_diagrams_graph& g, double r) :
    g(g), r(r), neighbors_finder() { }

inline void Layered_neighbors_finder::add(int v_point_index, int vlayer) {
  for (int l = neighbors_finder.size(); l <= vlayer; l++)
    neighbors_finder.emplace_back(Neighbors_finder(g, r));
  neighbors_finder.at(vlayer).add(v_point_index);
}

inline int Layered_neighbors_finder::pull_near(int u_point_index, int vlayer) {
  if (static_cast<int> (neighbors_finder.size()) <= vlayer)
    return null_point_index();
  return neighbors_finder.at(vlayer).pull_near(u_point_index);
}

inline int Layered_neighbors_finder::vlayers_number() const {
  return neighbors_finder.size();
}

}  // namespace bottleneck

}  // namespace Gudhi

#endif  // SRC_BOTTLENECK_INCLUDE_GUDHI_LAYERED_NEIGHBORS_FINDER_H_

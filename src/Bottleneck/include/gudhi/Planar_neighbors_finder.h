/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Francois Godi
 *
 *    Copyright (C) 2015  INRIA Sophia-Antipolis (France)
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

#ifndef SRC_BOTTLENECK_INCLUDE_GUDHI_PLANAR_NEIGHBORS_FINDER_H_
#define SRC_BOTTLENECK_INCLUDE_GUDHI_PLANAR_NEIGHBORS_FINDER_H_

#include <list>
#include <set>

#include "Persistence_diagrams_graph.h"

namespace Gudhi {

namespace bottleneck {

// Planar_neighbors_finder is a data structure used to find if a query point from U has planar neighbors in V with the
// planar distance.
// V's points have to be added manually using their index. A neighbor returned is automatically removed but we can also
// remove points manually using their index.

class Abstract_planar_neighbors_finder {
 public:
  Abstract_planar_neighbors_finder(const Persistence_diagrams_graph& g, double r);
  virtual ~Abstract_planar_neighbors_finder() = 0;
  virtual void add(int v_point_index) = 0;
  virtual void remove(int v_point_index) = 0;
  virtual bool contains(int v_point_index) const = 0;
  virtual int pull_near(int u_point_index) = 0;
  virtual std::unique_ptr< std::list<int> > pull_all_near(int u_point_index);

 protected:
  const Persistence_diagrams_graph& g;
  const double r;
};


// Naive_pnf is a nave implementation of Abstract_planar_neighbors_finder

class Naive_pnf : public Abstract_planar_neighbors_finder {
 public:
  Naive_pnf(const Persistence_diagrams_graph& g, double r);
  void add(int v_point_index);
  void remove(int v_point_index);
  bool contains(int v_point_index) const;
  int pull_near(int u_point_index);

 private:
  std::set<int> candidates;
};


// Planar_neighbors_finder is the used Abstract_planar_neighbors_finder's implementation
typedef Naive_pnf Planar_neighbors_finder;

Abstract_planar_neighbors_finder::Abstract_planar_neighbors_finder(const Persistence_diagrams_graph& g, double r) :
    g(g), r(r) { }

/* inline */ Abstract_planar_neighbors_finder::~Abstract_planar_neighbors_finder() { }

/* inline */ std::unique_ptr< std::list<int> > Abstract_planar_neighbors_finder::pull_all_near(int u_point_index) {
  std::unique_ptr< std::list<int> > all_pull(new std::list<int>);
  int last_pull = pull_near(u_point_index);
  while (last_pull != null_point_index()) {
    all_pull->emplace_back(last_pull);
    last_pull = pull_near(u_point_index);
  }
  return all_pull;
}

Naive_pnf::Naive_pnf(const Persistence_diagrams_graph& g, double r) :
    Abstract_planar_neighbors_finder(g, r), candidates() { }

/* inline */ void Naive_pnf::add(int v_point_index) {
  candidates.emplace(v_point_index);
}

/* inline */ void Naive_pnf::remove(int v_point_index) {
  candidates.erase(v_point_index);
}

/* inline */ bool Naive_pnf::contains(int v_point_index) const {
  return (candidates.count(v_point_index) > 0);
}

/* inline */ int Naive_pnf::pull_near(int u_point_index) {
  for (auto it = candidates.begin(); it != candidates.end(); ++it)
    if (g.distance(u_point_index, *it) <= r) {
      int tmp = *it;
      candidates.erase(it);
      return tmp;
    }
  return null_point_index();
}

}  // namespace bottleneck

}  // namespace Gudhi

#endif  // SRC_BOTTLENECK_INCLUDE_GUDHI_PLANAR_NEIGHBORS_FINDER_H_

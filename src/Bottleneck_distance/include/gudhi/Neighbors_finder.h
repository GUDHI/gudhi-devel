/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author:       Francois Godi
 *
 *    Copyright (C) 2015 Inria
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

#ifndef NEIGHBORS_FINDER_H_
#define NEIGHBORS_FINDER_H_

// Inclusion order is important for CGAL patch
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits.h>

#include <gudhi/Persistence_graph.h>
#include <gudhi/Internal_point.h>

// Specific for Windows
#define NOMINMAX

#include <unordered_set>
#include <vector>
#include <algorithm>  // for std::max
#include <cmath>  // for std::abs

namespace Gudhi {

namespace persistence_diagram {

/** \internal \brief Variant of CGAL::Fuzzy_iso_box to ensure that the box ic closed. It isn't clear how necessary that is.
 */
struct Square_query {
  typedef CGAL::Dimension_tag<2> D;
  typedef Internal_point Point_d;
  typedef double FT;
  bool contains(Point_d p) const {
    return std::max(std::abs(p.x()-c.x()), std::abs(p.y()-c.y())) <= size;
  }
  bool inner_range_intersects(CGAL::Kd_tree_rectangle<FT, D> const&r) const {
    return
      r.max_coord(0) >= c.x() - size &&
      r.min_coord(0) <= c.x() + size &&
      r.max_coord(1) >= c.y() - size &&
      r.min_coord(1) <= c.y() + size;
  }
  bool outer_range_contains(CGAL::Kd_tree_rectangle<FT, D> const&r) const {
    return
      r.min_coord(0) >= c.x() - size &&
      r.max_coord(0) <= c.x() + size &&
      r.min_coord(1) >= c.y() - size &&
      r.max_coord(1) <= c.y() + size;
  }
  Point_d c;
  FT size;
};

/** \internal \brief data structure used to find any point (including projections) in V near to a query point from U
 * (which can be a projection).
 *
 * V points have to be added manually using their index and before the first pull. A neighbor pulled is automatically
 * removed.
 *
 * \ingroup bottleneck_distance
 */
class Neighbors_finder {
  typedef CGAL::Dimension_tag<2> D;
  typedef CGAL::Search_traits<double, Internal_point, const double*, Construct_coord_iterator, D> Traits;
  typedef CGAL::Kd_tree<Traits> Kd_tree;

 public:
  /** \internal \brief Constructor taking the near distance definition as parameter. */
  Neighbors_finder(const Persistence_graph& g, double r);
  /** \internal \brief A point added will be possibly pulled. */
  void add(int v_point_index);
  /** \internal \brief Returns and remove a V point near to the U point given as parameter, null_point_index() if
   * there isn't such a point. */
  int pull_near(int u_point_index);
  /** \internal \brief Returns and remove all the V points near to the U point given as parameter. */
  std::vector<int> pull_all_near(int u_point_index);

 private:
  const Persistence_graph& g;
  const double r;
  Kd_tree kd_t;
  std::unordered_set<int> projections_f;
};

/** \internal \brief data structure used to find any point (including projections) in V near to a query point from U
 * (which can be a projection) in a layered graph layer given as parmeter.
 *
 * V points have to be added manually using their index and before the first pull. A neighbor pulled is automatically
 * removed.
 *
 * \ingroup bottleneck_distance
 */
class Layered_neighbors_finder {
 public:
  /** \internal \brief Constructor taking the near distance definition as parameter. */
  Layered_neighbors_finder(const Persistence_graph& g, double r);
  /** \internal \brief A point added will be possibly pulled. */
  void add(int v_point_index, int vlayer);
  /** \internal \brief Returns and remove a V point near to the U point given as parameter, null_point_index() if
   * there isn't such a point. */
  int pull_near(int u_point_index, int vlayer);
  /** \internal \brief Returns the number of layers. */
  int vlayers_number() const;

 private:
  const Persistence_graph& g;
  const double r;
  std::vector<std::unique_ptr<Neighbors_finder>> neighbors_finder;
};

inline Neighbors_finder::Neighbors_finder(const Persistence_graph& g, double r) :
    g(g), r(r), kd_t(), projections_f() { }

inline void Neighbors_finder::add(int v_point_index) {
  if (g.on_the_v_diagonal(v_point_index))
    projections_f.emplace(v_point_index);
  else
    kd_t.insert(g.get_v_point(v_point_index));
}

inline int Neighbors_finder::pull_near(int u_point_index) {
  int tmp;
  int c = g.corresponding_point_in_v(u_point_index);
  if (g.on_the_u_diagonal(u_point_index) && !projections_f.empty()) {
    // Any pair of projection is at distance 0
    tmp = *projections_f.cbegin();
    projections_f.erase(tmp);
  } else if (projections_f.count(c) && (g.distance(u_point_index, c) <= r)) {
    // Is the query point near to its projection ?
    tmp = c;
    projections_f.erase(tmp);
  } else {
    // Is the query point near to a V point in the plane ?
    Internal_point u_point = g.get_u_point(u_point_index);
    auto neighbor = kd_t.search_any_point(Square_query{u_point, r});
    if (!neighbor)
      return null_point_index();
    tmp = neighbor->point_index;
    auto point = g.get_v_point(tmp);
    int idx = point.point_index;
    kd_t.remove(point, [idx](Internal_point const&p){return p.point_index == idx;});
  }
  return tmp;
}

inline std::vector<int> Neighbors_finder::pull_all_near(int u_point_index) {
  std::vector<int> all_pull;
  int last_pull = pull_near(u_point_index);
  while (last_pull != null_point_index()) {
    all_pull.emplace_back(last_pull);
    last_pull = pull_near(u_point_index);
  }
  return all_pull;
}

inline Layered_neighbors_finder::Layered_neighbors_finder(const Persistence_graph& g, double r) :
    g(g), r(r), neighbors_finder() { }

inline void Layered_neighbors_finder::add(int v_point_index, int vlayer) {
  for (int l = neighbors_finder.size(); l <= vlayer; l++)
    neighbors_finder.emplace_back(std::unique_ptr<Neighbors_finder>(new Neighbors_finder(g, r)));
  neighbors_finder.at(vlayer)->add(v_point_index);
}

inline int Layered_neighbors_finder::pull_near(int u_point_index, int vlayer) {
  if (static_cast<int> (neighbors_finder.size()) <= vlayer)
    return null_point_index();
  return neighbors_finder.at(vlayer)->pull_near(u_point_index);
}

inline int Layered_neighbors_finder::vlayers_number() const {
  return static_cast<int> (neighbors_finder.size());
}

}  // namespace persistence_diagram

}  // namespace Gudhi

#endif  // NEIGHBORS_FINDER_H_

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

#ifndef SRC_BOTTLENECK_INCLUDE_GUDHI_NEIGHBORS_FINDER_H_
#define SRC_BOTTLENECK_INCLUDE_GUDHI_NEIGHBORS_FINDER_H_

#include <unordered_set>
#include <list>

#include "Planar_neighbors_finder.h"

namespace Gudhi {

namespace bottleneck {

/** \internal \brief data structure used to find any point (including projections) in V near to a query point from U
 * (which can be a projection).
 *
 * V points have to be added manually using their index and before the first pull. A neighbor pulled is automatically removed.
 *
 * \ingroup bottleneck_distance
 */
class Neighbors_finder {
public:
    /** \internal \brief Constructor taking the near distance definition as parameter. */
    Neighbors_finder(double r);
    /** \internal \brief A point added will be possibly pulled. */
    void add(int v_point_index);
    /** \internal \brief Returns and remove a V point near to the U point given as parameter, null_point_index() if there isn't such a point. */
    int pull_near(int u_point_index);
    /** \internal \brief Returns and remove all the V points near to the U point given as parameter. */
    std::unique_ptr< std::list<int> > pull_all_near(int u_point_index);

private:
    const double r;
    Planar_neighbors_finder planar_neighbors_f;
    std::unordered_set<int> projections_f;
    void remove(int v_point_index);
    bool contains(int v_point_index);
};

inline Neighbors_finder::Neighbors_finder(double r) :
    r(r), planar_neighbors_f(r), projections_f() { }

inline void Neighbors_finder::add(int v_point_index) {
    if (G::on_the_v_diagonal(v_point_index))
        projections_f.emplace(v_point_index);
    else
        planar_neighbors_f.add(v_point_index);
}

inline void Neighbors_finder::remove(int v_point_index) {
    if(v_point_index == null_point_index())
        return;
    if (G::on_the_v_diagonal(v_point_index))
        projections_f.erase(v_point_index);
    else
        planar_neighbors_f.remove(v_point_index);
}

inline bool Neighbors_finder::contains(int v_point_index) {
    return planar_neighbors_f.contains(v_point_index) || (projections_f.count(v_point_index)>0);
}

inline int Neighbors_finder::pull_near(int u_point_index) {
    int tmp;
    int c = G::corresponding_point_in_v(u_point_index);
    if (G::on_the_u_diagonal(u_point_index) && !projections_f.empty())
        //All projections are at distance 0
        tmp = *projections_f.cbegin();
    else if (contains(c) && (G::distance(u_point_index, c) <= r))
        //Is the query point near to its projection ?
        tmp = c;
    else
        //Is the query point near to a V point in the plane ?
        tmp = planar_neighbors_f.pull_near(u_point_index);
    remove(tmp);
    return tmp;
}

inline std::unique_ptr< std::list<int> > Neighbors_finder::pull_all_near(int u_point_index) {
    std::unique_ptr< std::list<int> > all_pull = std::move(planar_neighbors_f.pull_all_near(u_point_index));
    int last_pull = pull_near(u_point_index);
    while (last_pull != null_point_index()) {
        all_pull->emplace_back(last_pull);
        last_pull = pull_near(u_point_index);
    }
    return all_pull;
}

}  // namespace bottleneck

}  // namespace Gudhi

#endif  // SRC_BOTTLENECK_INCLUDE_GUDHI_NEIGHBORS_FINDER_H_

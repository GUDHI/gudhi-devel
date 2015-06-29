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

#ifndef SRC_BOTTLENECK_INCLUDE_GUDHI_LAYERED_NEIGHBORS_FINDER_H_
#define SRC_BOTTLENECK_INCLUDE_GUDHI_LAYERED_NEIGHBORS_FINDER_H_

#include <vector>

#include "Neighbors_finder.h"

namespace Gudhi {

namespace bottleneck {

/** \internal \brief data structure used to find any point (including projections) in V near to a query point from U
 * (which can be a projection) in a layered graph layer given as parmeter.
 *
 * V points have to be added manually using their index and before the first pull. A neighbor pulled is automatically removed.
 *
 * \ingroup bottleneck_distance
 */
class Layered_neighbors_finder {
public:
    /** \internal \brief Constructor taking the near distance definition as parameter. */
    Layered_neighbors_finder(double r);
    /** \internal \brief A point added will be possibly pulled. */
    void add(int v_point_index, int vlayer);
    /** \internal \brief Returns and remove a V point near to the U point given as parameter, null_point_index() if there isn't such a point. */
    int pull_near(int u_point_index, int vlayer);
    /** \internal \brief Returns the number of layers. */
    int vlayers_number() const;

private:
    const double r;
    std::vector<Neighbors_finder> neighbors_finder;
};

Layered_neighbors_finder::Layered_neighbors_finder(double r) :
    r(r), neighbors_finder() { }

inline void Layered_neighbors_finder::add(int v_point_index, int vlayer) {
    for (int l = neighbors_finder.size(); l <= vlayer; l++)
        neighbors_finder.emplace_back(Neighbors_finder(r));
    neighbors_finder.at(vlayer).add(v_point_index);
}

inline int Layered_neighbors_finder::pull_near(int u_point_index, int vlayer) {
    if (static_cast<int> (neighbors_finder.size()) <= vlayer)
        return null_point_index();
    return neighbors_finder.at(vlayer).pull_near(u_point_index);
}

inline int Layered_neighbors_finder::vlayers_number() const {
    return static_cast<int>(neighbors_finder.size());
}

}  // namespace bottleneck

}  // namespace Gudhi

#endif  // SRC_BOTTLENECK_INCLUDE_GUDHI_LAYERED_NEIGHBORS_FINDER_H_

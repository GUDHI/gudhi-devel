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

#ifndef SRC_BOTTLENECK_INCLUDE_GUDHI_PERSISTENCE_DIAGRAMS_GRAPH_H_
#define SRC_BOTTLENECK_INCLUDE_GUDHI_PERSISTENCE_DIAGRAMS_GRAPH_H_

#include <vector>
#include <set>
#include <cmath>
#include <utility>
#include <algorithm>
#include <math.h>
#include <memory>
#include "Internal_point.h"

namespace Gudhi {

namespace bipartite_graph_matching {


/** \internal \brief Structure representing an euclidean bipartite graph containing
 *  the points from the two persistence diagrams (including the projections).
 *
 * \ingroup bottleneck_distance
 */
class Persistence_diagrams_graph {
public:
    /** \internal \brief Initializer taking 2 Point (concept) ranges as parameters. */
    template<typename Persistence_diagram1, typename Persistence_diagram2>
    static void initialize(const Persistence_diagram1& diag1, const Persistence_diagram2& diag2, double e);
    /** \internal \brief Is the given point from U the projection of a point in V ? */
    static bool on_the_u_diagonal(int u_point_index);
    /** \internal \brief Is the given point from V the projection of a point in U ? */
    static bool on_the_v_diagonal(int v_point_index);
    /** \internal \brief Given a point from V, returns the corresponding (projection or projector) point in U. */
    static int corresponding_point_in_u(int v_point_index);
    /** \internal \brief Given a point from U, returns the corresponding (projection or projector) point in V. */
    static int corresponding_point_in_v(int u_point_index);
    /** \internal \brief Given a point from U and a point from V, returns the distance between those points. */
    static double distance(int u_point_index, int v_point_index);
    /** \internal \brief Returns size = |U| = |V|. */
    static int size();
    /** \internal \brief Returns the O(n^2) sorted distances between the points. */
    static std::shared_ptr< std::vector<double> > sorted_distances();

private:
    static std::vector<Internal_point> u;
    static std::vector<Internal_point> v;
    static Internal_point get_u_point(int u_point_index);
    static Internal_point get_v_point(int v_point_index);

    friend class Naive_pnf;
    friend class Cgal_pnf;
};

/** \internal \typedef \brief Shorter alias */
typedef Persistence_diagrams_graph G;

// static initialization
std::vector<Internal_point> G::u = [] {return std::vector<Internal_point>();}();
std::vector<Internal_point> G::v = [] {return std::vector<Internal_point>();}();

template<typename Persistence_diagram1, typename Persistence_diagram2>
inline void G::initialize(const Persistence_diagram1 &diag1,
                          const Persistence_diagram2 &diag2, double e){
    u.clear();
    v.clear();
    for (auto it = diag1.cbegin(); it != diag1.cend(); ++it)
        if (it->y() - it->x() > e)
            u.emplace_back(*it);
    for (auto it = diag2.cbegin(); it != diag2.cend(); ++it)
        if (it->y() - it->x() > e)
            v.emplace_back(*it);
    if (u.size() < v.size())
        swap(u, v);
}

inline bool G::on_the_u_diagonal(int u_point_index) {
    return u_point_index >= static_cast<int> (u.size());
}

inline bool G::on_the_v_diagonal(int v_point_index) {
    return v_point_index >= static_cast<int> (v.size());
}

inline int G::corresponding_point_in_u(int v_point_index) {
    return on_the_v_diagonal(v_point_index) ?
                v_point_index - static_cast<int> (v.size()) : v_point_index + static_cast<int> (u.size());
}

inline int G::corresponding_point_in_v(int u_point_index) {
    return on_the_u_diagonal(u_point_index) ?
                u_point_index - static_cast<int> (u.size()) : u_point_index + static_cast<int> (v.size());
}

inline double G::distance(int u_point_index, int v_point_index) {
    if (on_the_u_diagonal(u_point_index) && on_the_v_diagonal(v_point_index))
        return 0;
    Internal_point p_u = get_u_point(u_point_index);
    Internal_point p_v = get_v_point(v_point_index);
    return std::max(std::fabs(p_u.x() - p_v.x()), std::fabs(p_u.y() - p_v.y()));
}

inline int G::size() {
    return static_cast<int> (u.size() + v.size());
}

inline std::shared_ptr< std::vector<double> > G::sorted_distances() {
    // could be optimized
    std::set<double> sorted_distances;
    for (int u_point_index = 0; u_point_index < size(); ++u_point_index)
        for (int v_point_index = 0; v_point_index < size(); ++v_point_index)
            sorted_distances.emplace(distance(u_point_index, v_point_index));
    std::shared_ptr< std::vector<double> > sd_up(new std::vector<double>(sorted_distances.cbegin(), sorted_distances.cend()));
    return sd_up;
}

inline Internal_point G::get_u_point(int u_point_index) {
    if (!on_the_u_diagonal(u_point_index))
        return u.at(u_point_index);
    Internal_point projector = v.at(corresponding_point_in_v(u_point_index));
    double m = (projector.x() + projector.y()) / 2;
    return Internal_point(m, m);
}

inline Internal_point G::get_v_point(int v_point_index) {
    if (!on_the_v_diagonal(v_point_index))
        return v.at(v_point_index);
    Internal_point projector = u.at(corresponding_point_in_u(v_point_index));
    double m = (projector.x() + projector.y()) / 2;
    return Internal_point(m, m);
}

}  // namespace bipartite_graph_matching

}  // namespace Gudhi

#endif  // SRC_BOTTLENECK_INCLUDE_GUDHI_PERSISTENCE_DIAGRAMS_GRAPH_H_

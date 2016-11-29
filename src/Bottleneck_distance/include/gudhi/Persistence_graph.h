/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author:       Francois Godi
 *
 *    Copyright (C) 2015  INRIA (France)
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

#ifndef PERSISTENCE_GRAPH_H_
#define PERSISTENCE_GRAPH_H_

#include <vector>
#include <set>
#include <gudhi/Internal_point.h>

namespace Gudhi {

namespace bottleneck_distance {


/** \internal \brief Structure representing an euclidean bipartite graph containing
 *  the points from the two persistence diagrams (including the projections).
 *
 * \ingroup bottleneck_distance
 */
class Persistence_graph {
public:
    /** \internal \brief Initializer taking 2 Point (concept) ranges as parameters. */
    template<typename Persistence_diagram1, typename Persistence_diagram2>
    Persistence_graph(const Persistence_diagram1& diag1, const Persistence_diagram2& diag2, double e);
    /** \internal \brief Is the given point from U the projection of a point in V ? */
    bool on_the_u_diagonal(int u_point_index) const;
    /** \internal \brief Is the given point from V the projection of a point in U ? */
    bool on_the_v_diagonal(int v_point_index) const;
    /** \internal \brief Given a point from V, returns the corresponding (projection or projector) point in U. */
    int corresponding_point_in_u(int v_point_index) const;
    /** \internal \brief Given a point from U, returns the corresponding (projection or projector) point in V. */
    int corresponding_point_in_v(int u_point_index) const;
    /** \internal \brief Given a point from U and a point from V, returns the distance between those points. */
    double distance(int u_point_index, int v_point_index) const;
    /** \internal \brief Returns size = |U| = |V|. */
    int size() const;
    /** \internal \brief Returns the O(n^2) sorted distances between the points. */
    std::vector<double> sorted_distances() const;
    /** \internal \brief Returns an upper bound of the diameter of the convex hull */
    double diameter_bound() const;
    /** \internal \brief Returns the corresponding internal point */
    Internal_point get_u_point(int u_point_index) const;
    /** \internal \brief Returns the corresponding internal point */
    Internal_point get_v_point(int v_point_index) const;

private:
    std::vector<Internal_point> u;
    std::vector<Internal_point> v;
};

template<typename Persistence_diagram1, typename Persistence_diagram2>
Persistence_graph::Persistence_graph(const Persistence_diagram1 &diag1,
                                     const Persistence_diagram2 &diag2, double e)
    : u(), v()
{
    for (auto it = diag1.cbegin(); it != diag1.cend(); ++it)
        if (it->second != std::numeric_limits<double>::infinity() && it->second - it->first > e)
            u.push_back(Internal_point(std::get<0>(*it), std::get<1>(*it), u.size()));
    for (auto it = diag2.cbegin(); it != diag2.cend(); ++it)
        if (it->second != std::numeric_limits<double>::infinity() && it->second - it->first > e)
            v.push_back(Internal_point(std::get<0>(*it), std::get<1>(*it), v.size()));
    if (u.size() < v.size())
        swap(u, v);
}

inline bool Persistence_graph::on_the_u_diagonal(int u_point_index) const {
    return u_point_index >= static_cast<int> (u.size());
}

inline bool Persistence_graph::on_the_v_diagonal(int v_point_index) const {
    return v_point_index >= static_cast<int> (v.size());
}

inline int Persistence_graph::corresponding_point_in_u(int v_point_index) const {
    return on_the_v_diagonal(v_point_index) ?
                v_point_index - static_cast<int> (v.size()) : v_point_index + static_cast<int> (u.size());
}

inline int Persistence_graph::corresponding_point_in_v(int u_point_index) const {
    return on_the_u_diagonal(u_point_index) ?
                u_point_index - static_cast<int> (u.size()) : u_point_index + static_cast<int> (v.size());
}

inline double Persistence_graph::distance(int u_point_index, int v_point_index) const {
    if (on_the_u_diagonal(u_point_index) && on_the_v_diagonal(v_point_index))
        return 0.;
    Internal_point p_u = get_u_point(u_point_index);
    Internal_point p_v = get_v_point(v_point_index);
    return std::max(std::fabs(p_u.x() - p_v.x()), std::fabs(p_u.y() - p_v.y()));
}

inline int Persistence_graph::size() const {
    return static_cast<int> (u.size() + v.size());
}

inline std::vector<double> Persistence_graph::sorted_distances() const {
    // could be optimized
    std::set<double> sorted_distances;
    sorted_distances.emplace(0.);
    for (int u_point_index = 0; u_point_index < size(); ++u_point_index)
        for (int v_point_index = 0; v_point_index < size(); ++v_point_index)
            sorted_distances.emplace(distance(u_point_index, v_point_index));
    sorted_distances.emplace(std::numeric_limits<double>::infinity());
    return std::vector<double>(sorted_distances.begin(),sorted_distances.end());
}

inline Internal_point Persistence_graph::get_u_point(int u_point_index) const {
    if (!on_the_u_diagonal(u_point_index))
        return u.at(u_point_index);
    Internal_point projector = v.at(corresponding_point_in_v(u_point_index));
    double m = (projector.x() + projector.y()) / 2;
    return Internal_point(m,m,u_point_index);
}

inline Internal_point Persistence_graph::get_v_point(int v_point_index) const {
    if (!on_the_v_diagonal(v_point_index))
        return v.at(v_point_index);
    Internal_point projector = u.at(corresponding_point_in_u(v_point_index));
    double m = (projector.x() + projector.y()) / 2;
    return Internal_point(m,m,v_point_index);
}

inline double Persistence_graph::diameter_bound() const {
    double max = 0.;
    for(auto it = u.cbegin(); it != u.cend(); it++)
        max = std::max(max,it->y());
    for(auto it = v.cbegin(); it != v.cend(); it++)
        max = std::max(max,it->y());
    return max;
}

}  // namespace bottleneck_distance

}  // namespace Gudhi

#endif  // PERSISTENCE_GRAPH_H_

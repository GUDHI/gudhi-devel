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
#include <utility>  // for pair<>
#include <algorithm>  // for max

namespace Gudhi {

namespace bottleneck {

// Diagram_point is the type of the persistence diagram's points
typedef std::pair<double, double> Diagram_point;

// Return the used index for encoding none of the points
int null_point_index();

// Persistence_diagrams_graph is the interface beetwen any external representation of the two persistence diagrams and
// the bottleneck distance computation. An interface is necessary to ensure basic functions complexity.

class Persistence_diagrams_graph {
 public:
    // Persistence_diagram1 and 2 are the types of any externals representations of persistence diagrams.
    // They have to have an iterator over points, which have to have fields first (for birth) and second (for death).
    template<typename Persistence_diagram1, typename Persistence_diagram2>
    Persistence_diagrams_graph(const Persistence_diagram1& diag1, const Persistence_diagram2& diag2, double e);
    Persistence_diagrams_graph();
    bool on_the_u_diagonal(int u_point_index) const;
    bool on_the_v_diagonal(int v_point_index) const;
    int corresponding_point_in_u(int v_point_index) const;
    int corresponding_point_in_v(int u_point_index) const;
    double distance(int u_point_index, int v_point_index) const;
    int size() const;
    std::unique_ptr< std::vector<double> > sorted_distances();

 private:
    std::vector<Diagram_point> u;
    std::vector<Diagram_point> v;
    Diagram_point get_u_point(int u_point_index) const;
    Diagram_point get_v_point(int v_point_index) const;
};

/* inline */ int null_point_index() {
    return -1;
}

template<typename Persistence_diagram1, typename Persistence_diagram2>
Persistence_diagrams_graph::Persistence_diagrams_graph(const Persistence_diagram1 &diag1,
                                                       const Persistence_diagram2 &diag2, double e)
    : u(), v() {
    for (auto it = diag1.cbegin(); it != diag1.cend(); ++it)
        if (it->second - it->first > e)
            u.emplace_back(*it);
    for (auto it = diag2.cbegin(); it != diag2.cend(); ++it)
        if (it->second - it->first > e)
            v.emplace_back(*it);
    if (u.size() < v.size())
        swap(u, v);
}

Persistence_diagrams_graph::Persistence_diagrams_graph()
    : u(), v() { }

/* inline */ bool Persistence_diagrams_graph::on_the_u_diagonal(int u_point_index) const {
    return u_point_index >= static_cast<int> (u.size());
}

/* inline */ bool Persistence_diagrams_graph::on_the_v_diagonal(int v_point_index) const {
    return v_point_index >= static_cast<int> (v.size());
}

/* inline */ int Persistence_diagrams_graph::corresponding_point_in_u(int v_point_index) const {
    return on_the_v_diagonal(v_point_index) ?
                v_point_index - static_cast<int> (v.size()) : v_point_index + static_cast<int> (u.size());
}

/* inline */ int Persistence_diagrams_graph::corresponding_point_in_v(int u_point_index) const {
    return on_the_u_diagonal(u_point_index) ?
                u_point_index - static_cast<int> (u.size()) : u_point_index + static_cast<int> (v.size());
}

/* inline */ double Persistence_diagrams_graph::distance(int u_point_index, int v_point_index) const {
  // could be optimized for the case where one point is the projection of the other
  if (on_the_u_diagonal(u_point_index) && on_the_v_diagonal(v_point_index))
    return 0;
  Diagram_point p_u = get_u_point(u_point_index);
  Diagram_point p_v = get_v_point(v_point_index);
  return std::max(std::fabs(p_u.first - p_v.first), std::fabs(p_u.second - p_v.second));
}

/* inline */ int Persistence_diagrams_graph::size() const {
    return static_cast<int> (u.size() + v.size());
}

/* inline */ std::unique_ptr< std::vector<double> > Persistence_diagrams_graph::sorted_distances() {
    // could be optimized
    std::set<double> sorted_distances;
    for (int u_point_index = 0; u_point_index < size(); ++u_point_index)
        for (int v_point_index = 0; v_point_index < size(); ++v_point_index)
            sorted_distances.emplace(distance(u_point_index, v_point_index));
    std::unique_ptr< std::vector<double> > sd_up(new std::vector<double>(sorted_distances.cbegin(), sorted_distances.cend()));
    return sd_up;
}

/* inline */ Diagram_point Persistence_diagrams_graph::get_u_point(int u_point_index) const {
    if (!on_the_u_diagonal(u_point_index))
        return u.at(u_point_index);
    Diagram_point projector = v.at(corresponding_point_in_v(u_point_index));
    double x = (projector.first + projector.second) / 2;
    return Diagram_point(x, x);
}

/* inline */ Diagram_point Persistence_diagrams_graph::get_v_point(int v_point_index) const {
    if (!on_the_v_diagonal(v_point_index))
        return v.at(v_point_index);
    Diagram_point projector = u.at(corresponding_point_in_u(v_point_index));
    double x = (projector.first + projector.second) / 2;
    return Diagram_point(x, x);
}

}  // namespace bottleneck

}  // namespace Gudhi

#endif  // SRC_BOTTLENECK_INCLUDE_GUDHI_PERSISTENCE_DIAGRAMS_GRAPH_H_

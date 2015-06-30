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

#ifndef SRC_BOTTLENECK_INCLUDE_GUDHI_GRID_CELL_H_
#define SRC_BOTTLENECK_INCLUDE_GUDHI_GRID_CELL_H_

#include <list>
#include <set>
#include <map>

#include "Persistence_diagrams_graph.h"

namespace Gudhi {

namespace bipartite_graph_matching {

/** \internal \brief TODO
 *
 * \ingroup bottleneck_distance
 */
class Grid_cell {
public:
    Grid_cell(double r);
    void add(int v_point_index);
    void remove(int v_point_index);
    bool contains(int v_point_index) const;
    int pull_center();
    int pull_xi(int u_point_index);
    int pull_xd(int u_point_index);
    int pull_yi(int u_point_index);
    int pull_yd(int u_point_index);
    int pull_xi_yi(int u_point_index);
    int pull_xi_yd(int u_point_index);
    int pull_xd_yi(int u_point_index);
    int pull_xd_yd(int u_point_index);

private:
    double r;
    std::set<int, G::Compare_x> xi_order;
    std::set<int, G::Compare_y> yi_order;
    struct Hidden_points_tree{int point; std::list<Hidden_points_tree> hidden;};
    typedef std::map<int, std::list<Hidden_points_tree>, G::Compare_x> Corner_tree;
    Corner_tree xi_yi_order;
    Corner_tree xi_yd_order;
    Corner_tree xd_yi_order;
    Corner_tree xd_yd_order;
    void remove_aux(Corner_tree t, int v_point_index);
    void build_xi_yi();
    void build_xi_yd();
    void build_xd_yi();
    void build_xd_yd();
};


inline Grid_cell::Grid_cell(double r)
    : r(r), xi_order(G::Compare_x(r)), yi_order(G::Compare_y(r)), xi_yi_order(G::Compare_x(r)),
      xi_yd_order(G::Compare_x(r)), xd_yi_order(G::Compare_x(r)), xd_yd_order(G::Compare_x(r)) {}

inline void Grid_cell::add(int v_point_index){
    xi_order.emplace(v_point_index);
}

inline bool Grid_cell::contains(int v_point_index) const{
    return (xi_order.count(v_point_index) > 0);
}

inline void Grid_cell::remove_aux(Corner_tree t, int v_point_index){
    if(t.empty())
        return;
    std::list<Hidden_points_tree> hidden_points = t.at(v_point_index);
    t.erase(v_point_index);
    for(auto it = hidden_points.begin(); it != hidden_points.end(); ++it)
        t.emplace(it->point,it->hidden);

}

inline void Grid_cell::remove(int v_point_index){
    xi_order.erase(v_point_index);
    yi_order.erase(v_point_index);
    remove_aux(xi_yi_order,v_point_index);
    remove_aux(xi_yd_order,v_point_index);
    remove_aux(xd_yi_order,v_point_index);
    remove_aux(xd_yd_order,v_point_index);
}

//factorization needed \/ \/ \/

inline int Grid_cell::pull_center(){
    if(xi_order.empty())
        return null_point_index();
    int v_point_index = *xi_order.begin();
    remove(v_point_index);
    return v_point_index;
}

inline int Grid_cell::pull_xi(int u_point_index){
    if(xi_order.empty())
        return null_point_index();
    int v_point_index = *xi_order.begin(); //!
    if(G::distance(u_point_index,v_point_index)<=r){
        remove(v_point_index);
        return v_point_index;
    }
    return null_point_index();
}

inline int Grid_cell::pull_xd(int u_point_index){
    if(xi_order.empty())
        return null_point_index();
    int v_point_index = *xi_order.rbegin(); //!
    if(G::distance(u_point_index,v_point_index)<=r){
        remove(v_point_index);
        return v_point_index;
    }
    return null_point_index();
}

inline int Grid_cell::pull_yi(int u_point_index){
    if(xi_order.empty())
        return null_point_index();
    if(yi_order.empty())
        for(auto it = xi_order.begin(); it!= xi_order.end(); ++it)
            yi_order.emplace(*it);
    int v_point_index = *yi_order.begin(); //!
    if(G::distance(u_point_index,v_point_index)<=r){
        remove(v_point_index);
        return v_point_index;
    }
    return null_point_index();
}

inline int Grid_cell::pull_yd(int u_point_index){
    if(xi_order.empty())
        return null_point_index();
    if(yi_order.empty())
        for(auto it = xi_order.begin(); it!= xi_order.end(); ++it)
            yi_order.emplace(*it);
    int v_point_index = *yi_order.rbegin(); //!
    if(G::distance(u_point_index,v_point_index)<=r){
        remove(v_point_index);
        return v_point_index;
    }
    return null_point_index();
}

inline int Grid_cell::pull_xi_yi(int u_point_index){
    if(xi_order.empty())
        return null_point_index();
    if(xi_yi_order.empty())
        build_xi_yi();
    auto it = xi_yi_order.upper_bound(u_point_index+G::size());
    if(it==xi_yi_order.cbegin()) //!
        return null_point_index();
    it--; //!
    int v_point_index = it->first;
    for(auto it2=it->second.begin();it2!=it->second.end();it2++)
        xi_yi_order.emplace(it2->point,it2->hidden);
    if(G::distance(u_point_index,v_point_index)<=r){
        remove(v_point_index);
        return v_point_index;
    }
    return null_point_index();
}

}  // namespace bipartite_graph_matching

}  // namespace Gudhi

#endif  // SRC_BOTTLENECK_INCLUDE_GUDHI_GRID_CELL_H_

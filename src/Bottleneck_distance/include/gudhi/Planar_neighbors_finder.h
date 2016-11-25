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

#ifndef PLANAR_NEIGHBORS_FINDER_H_
#define PLANAR_NEIGHBORS_FINDER_H_

// Inclusion order is important for CGAL patch
#include "CGAL/Kd_tree_node.h"
#include "CGAL/Kd_tree.h"
#include "CGAL/Orthogonal_incremental_neighbor_search.h"

#include <CGAL/Weighted_Minkowski_distance.h>
#include <CGAL/Search_traits.h>

#include <gudhi/Persistence_diagrams_graph.h>
#include <gudhi/Construct_coord_iterator.h>


namespace Gudhi {

namespace bottleneck_distance {

/** \internal \brief Structure used to find any point in V near (according to the planar distance) to a query point from U.
 *
 * V points have to be added manually using their index and before the first remove/pull. A neighbor pulled is automatically removed. but we can also
 * remove points manually using their index.
 *
 * \ingroup bottleneck_distance
 */
class Naive_pnf {
public:
    /** \internal \brief Constructor taking the near distance definition as parameter. */
    Naive_pnf(double r_);
    /** \internal \brief A point added will be possibly pulled. */
    void add(int v_point_index);
    /** \internal \brief A point manually removed will no longer be possibly pulled. */
    void remove(int v_point_index);
    /** \internal \brief Can the point given as parameter be returned ? */
    bool contains(int v_point_index) const;
    /** \internal \brief Provide and remove a V point near to the U point given as parameter, null_point_index() if there isn't such a point. */
    int pull_near(int u_point_index);
    /** \internal \brief Provide and remove all the V points near to the U point given as parameter. */
    std::vector<int> pull_all_near(int u_point_index);

private:
    double r;
    std::pair<int,int> get_v_key(int v_point_index) const;
    std::multimap<std::pair<int,int>,int> grid;
};

class Planar_neighbors_finder {

    typedef CGAL::Dimension_tag<2> D;
    typedef CGAL::Search_traits<double, Internal_point, const double*, CGAL::Construct_coord_iterator, D> Traits;
    typedef CGAL::Weighted_Minkowski_distance<Traits> Distance;
    typedef CGAL::Orthogonal_incremental_neighbor_search<Traits, Distance> K_neighbor_search;
    typedef K_neighbor_search::Tree Kd_tree;


public:
    /** \internal \brief Constructor taking the near distance definition as parameter. */
    Planar_neighbors_finder(double r_);
    /** \internal \brief A point added will be possibly pulled. */
    void add(int v_point_index);
    /** \internal \brief A point manually removed will no longer be possibly pulled. */
    void remove(int v_point_index);
    /** \internal \brief Can the point given as parameter be returned ? */
    bool contains(int v_point_index) const;
    /** \internal \brief Provide a V point near to the U point given as parameter, null_point_index() if there isn't such a point. */
    int pull_near(int u_point_index);
    /** \internal \brief Provide and remove all the V points near to the U point given as parameter. */
    virtual std::vector<int> pull_all_near(int u_point_index);

private:
    double r;
    std::set<int> contents;
    Kd_tree kd_t;
};


/** \internal \brief Constructor taking the near distance definition as parameter. */
inline Planar_neighbors_finder::Planar_neighbors_finder(double r_)
    : r(r_), contents(), kd_t() {}


/** \internal \brief A point added will be possibly pulled. */
inline void Planar_neighbors_finder::add(int v_point_index){
    if(v_point_index == null_point_index())
        return;
    contents.insert(v_point_index);
    kd_t.insert(G::get_v_point(v_point_index));
}

/** \internal \brief A point manually removed will no longer be possibly pulled. */
inline void Planar_neighbors_finder::remove(int v_point_index){
        contents.erase(v_point_index);
        kd_t.remove(G::get_v_point(v_point_index));
}

/** \internal \brief Can the point given as parameter be returned ? */
inline bool Planar_neighbors_finder::contains(int v_point_index) const{
    if(v_point_index == null_point_index())
        return false;
    return contents.count(v_point_index)>0;
}

/** \internal \brief Provide and remove a V point near to the U point given as parameter, null_point_index() if there isn't such a point. */
inline int Planar_neighbors_finder::pull_near(int u_point_index){
    Internal_point u_point = G::get_u_point(u_point_index);
    std::vector<double> w = {1., 1.};
    K_neighbor_search search(kd_t, u_point, 0., true, Distance(0, 2, w));
    auto it = search.begin();
    if(it==search.end() || G::distance(u_point_index, it->first.point_index) > r)
        return null_point_index();
    int tmp = it->first.point_index;
    if(!contains(tmp))
        std::cout << "!! A kd_tree returns a point (Point_index:" << tmp << ") previously removed !!" << std::endl;
    remove(tmp);
    return tmp;
}

inline std::vector<int> Planar_neighbors_finder::pull_all_near(int u_point_index) {
    std::vector<int> all_pull;
    int last_pull = pull_near(u_point_index);
    while (last_pull != null_point_index()) {
        all_pull.emplace_back(last_pull);
        last_pull = pull_near(u_point_index);
    }
    return all_pull;
}


}  // namespace bottleneck_distance

}  // namespace Gudhi

#endif  // PLANAR_NEIGHBORS_FINDER_H_

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
#include <map>
#include <CGAL/Search_traits.h>
#include <CGAL/Orthogonal_incremental_neighbor_search.h>
#include <CGAL/Weighted_Minkowski_distance.h>
#include "../CGAL/Miscellaneous.h"
#include "Persistence_diagrams_graph.h"


namespace Gudhi {

namespace bipartite_graph_matching {

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
    virtual std::shared_ptr< std::list<int> > pull_all_near(int u_point_index);

private:
    double r;
    std::pair<int,int> get_v_key(int v_point_index) const;
    std::multimap<std::pair<int,int>,int> grid;
};

class Cgal_pnf {

    typedef CGAL::Dimension_tag<2> D;
    typedef CGAL::Search_traits<double, Internal_point, const double*, CGAL::Construct_coord_iterator, D> Traits;
    typedef CGAL::Weighted_Minkowski_distance<Traits> Distance;
    typedef CGAL::Orthogonal_incremental_neighbor_search<Traits, Distance> K_neighbor_search;
    typedef K_neighbor_search::Tree Kd_tree;


public:
    /** \internal \brief Constructor taking the near distance definition as parameter. */
    Cgal_pnf(double r_);
    /** \internal \brief A point added will be possibly pulled. */
    void add(int v_point_index);
    /** \internal \brief A point manually removed will no longer be possibly pulled. */
    void remove(int v_point_index);
    /** \internal \brief Can the point given as parameter be returned ? */
    bool contains(int v_point_index) const;
    /** \internal \brief Provide a V point near to the U point given as parameter, null_point_index() if there isn't such a point. */
    int pull_near(int u_point_index);
    /** \internal \brief Provide and remove all the V points near to the U point given as parameter. */
    virtual std::shared_ptr< std::list<int> > pull_all_near(int u_point_index);

private:
    double r;
    std::set<int> contents;
    Kd_tree kd_t;
};

/** \internal \typedef \brief Planar_neighbors_finder is the used implementation. */
typedef Cgal_pnf Planar_neighbors_finder;

inline Naive_pnf::Naive_pnf(double r_) :
    r(r_), grid() { }


inline std::pair<int,int> Naive_pnf::get_v_key(int v_point_index) const{
    Internal_point v_point = G::get_v_point(v_point_index);
    return std::make_pair(static_cast<int>(v_point.x()/r), static_cast<int>(v_point.y()/r));
}

inline void Naive_pnf::add(int v_point_index) {
    grid.emplace(get_v_key(v_point_index),v_point_index);
}

inline void Naive_pnf::remove(int v_point_index) {
    if(!contains(v_point_index))
    for(auto it = grid.find(get_v_key(v_point_index)); it!=grid.end(); it++)
        if(it->second==v_point_index){
            grid.erase(it);
            return;
        }
}

inline bool Naive_pnf::contains(int v_point_index) const {
    if(v_point_index == null_point_index())
        return false;
    for(auto it = grid.find(get_v_key(v_point_index)); it!=grid.end(); it++)
        if(it->second==v_point_index)
            return true;
    return false;
}

inline int Naive_pnf::pull_near(int u_point_index) {
    Internal_point u_point = G::get_u_point(u_point_index);
    int i0 = static_cast<int>(u_point.x()/r);
    int j0 = static_cast<int>(u_point.y()/r);
    for(int i = 1; i<= 3; i++)
        for(int j = 1; j<= 3; j++)
            for(auto it = grid.find(std::make_pair(i0 +(i%3)-1, j0+(j%3)-1)); it!=grid.end(); it++)
                if (G::distance(u_point_index, it->second) <= r) {
                    int tmp = it->second;
                    return tmp;
                }
    return null_point_index();
}

inline std::shared_ptr< std::list<int> > Naive_pnf::pull_all_near(int u_point_index) {
    std::shared_ptr< std::list<int> > all_pull(new std::list<int>);
    Internal_point u_point = G::get_u_point(u_point_index);
    int i0 = static_cast<int>(u_point.x()/r);
    int j0 = static_cast<int>(u_point.y()/r);
    for(int i = 1; i<= 3; i++)
        for(int j = 1; j<= 3; j++)
            for(auto it = grid.find(std::make_pair(i0 +(i%3)-1, j0+(j%3)-1)); it!=grid.end(); it++)
                if (G::distance(u_point_index, it->second) <= r) {
                    int tmp = it->second;
                    grid.erase(it);
                    all_pull->emplace_back(tmp);
                }
    return all_pull;
}


/** \internal \brief Constructor taking the near distance definition as parameter. */
inline Cgal_pnf::Cgal_pnf(double r_)
    : r(r_), contents(), kd_t() {}


/** \internal \brief A point added will be possibly pulled. */
inline void Cgal_pnf::add(int v_point_index){
    contents.insert(v_point_index);
    kd_t.insert(G::get_v_point(v_point_index));
}

/** \internal \brief A point manually removed will no longer be possibly pulled. */
inline void Cgal_pnf::remove(int v_point_index){
    contents.erase(v_point_index);
    kd_t.remove(G::get_v_point(v_point_index));
}

/** \internal \brief Can the point given as parameter be returned ? */
inline bool Cgal_pnf::contains(int v_point_index) const{
    return contents.count(v_point_index)>0;
}

/** \internal \brief Provide and remove a V point near to the U point given as parameter, null_point_index() if there isn't such a point. */
inline int Cgal_pnf::pull_near(int u_point_index){
    Internal_point u_point = G::get_u_point(u_point_index);
    std::vector<double> w = {1., 1.};
    K_neighbor_search search(kd_t, u_point, 0., true, Distance(0, 2, w));
    auto it = search.begin();
    if(it==search.end() || G::distance(u_point_index, it->first.point_index) > r)
        return null_point_index();
    return it->first.point_index;
}

inline std::shared_ptr< std::list<int> > Cgal_pnf::pull_all_near(int u_point_index) {
    std::shared_ptr< std::list<int> > all_pull(new std::list<int>);
    int last_pull = pull_near(u_point_index);
    while (last_pull != null_point_index()) {
        all_pull->emplace_back(last_pull);
        remove(last_pull);
        last_pull = pull_near(u_point_index);
    }
    return all_pull;
}


}  // namespace bipartite_graph_matching

}  // namespace Gudhi

#endif  // SRC_BOTTLENECK_INCLUDE_GUDHI_PLANAR_NEIGHBORS_FINDER_H_

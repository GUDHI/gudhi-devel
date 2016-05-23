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

#ifndef SRC_BOTTLENECK_INCLUDE_GUDHI_INTERNAL_POINT_H_
#define SRC_BOTTLENECK_INCLUDE_GUDHI_INTERNAL_POINT_H_

//namespace Gudhi {

//namespace bipartite_graph_matching {

/** \internal \brief Returns the used index for encoding none of the points */
int null_point_index();

/** \internal \typedef \brief Internal_point is the internal points representation, indexes used outside. */
struct Internal_point {
    double vec[2];
    int point_index;
    Internal_point() { vec[0]= vec[1] = 0.; point_index = null_point_index(); }
    Internal_point(double x, double y) { vec[0]=x; vec[1]=y; point_index = null_point_index(); }
    double x() const { return vec[ 0 ]; }
    double y() const { return vec[ 1 ]; }
    double& x() { return vec[ 0 ]; }
    double& y() { return vec[ 1 ]; }
    bool operator==(const Internal_point& p) const
    {
        return (x() == p.x()) && (y() == p.y());
    }
    bool  operator!=(const Internal_point& p) const { return ! (*this == p); }
};

namespace CGAL {
template <>
struct Kernel_traits<Internal_point> {
    struct Kernel {
        typedef double FT;
        typedef double RT;
    };
};
}

struct Construct_coord_iterator {
    typedef  const double* result_type;
    const double* operator()(const Internal_point& p) const
    { return static_cast<const double*>(p.vec); }
    const double* operator()(const Internal_point& p, int)  const
    { return static_cast<const double*>(p.vec+2); }
};
/*
struct Distance {
    typedef Internal_point Query_item;
    typedef Internal_point Point_d;
    typedef double FT;
    typedef CGAL::Dimension_tag<2> D;

    FT transformed_distance(const Query_item& q, const Point_d& p) const {
        FT d0= std::abs(q.x()-p.x());
        FT d1= std::abs(q.y()-p.y());
        return std::max(d0,d1);
    }
   FT min_distance_to_rectangle(const Query_item& q,
                                     const CGAL::Kd_tree_rectangle<FT,D>& b) const {
        FT d0(0.), d1(0.);
        if (q.x() < b.min_coord(0))
            d0 += (b.min_coord(0)-q.x());
        if (q.x() > b.max_coord(0))
            d0 += (q.x()-b.max_coord(0));
        if (q.y() < b.min_coord(1))
            d1 += (b.min_coord(1)-q.y());
        if (q.y() > b.max_coord(1))
            d1 += (q.y()-b.max_coord(1));
        return std::max(d0,d1);
    }
    FT min_distance_to_rectangle(const Query_item& q, const CGAL::Kd_tree_rectangle<FT,D>& b,std::vector<FT>& dists){
        dists[0] = dists[1] = 0.;
        if (q.x() < b.min_coord(0))
            dists[0] = (b.min_coord(0)- q.x());
        if (q.x() > b.max_coord(0))
            dists[0] = (q.x()-b.max_coord(0));
        if (q.y() < b.min_coord(1))
            dists[1] = (b.min_coord(1)-q.y());
        if (q.y() > b.max_coord(1))
            dists[1] = (q.y()-b.max_coord(1));
        return std::max(dists[0],dists[1]);
    }
    FT max_distance_to_rectangle(const Query_item& q, const CGAL::Kd_tree_rectangle<FT,D>& b) const {
        FT d0 = (q.x() >= (b.min_coord(0)+b.max_coord(0))/2.) ?
                    (q.x()-b.min_coord(0)) : (b.max_coord(0)-q.x());
        FT d1 = (q.y() >= (b.min_coord(1)+b.max_coord(1))/2.) ?
                    (q.y()-b.min_coord(1)) : (b.max_coord(1)-q.y());
        return std::max(d0, d1);
    }
    FT max_distance_to_rectangle(const Query_item& q, const CGAL::Kd_tree_rectangle<FT,D>& b,std::vector<FT>& dists){
        dists[0] = (q.x() >= (b.min_coord(0)+b.max_coord(0))/2.0) ?
                    (q.x()-b.min_coord(0)) : (b.max_coord(0)-q.x());
        dists[1] = (q.y() >= (b.min_coord(1)+b.max_coord(1))/2.0) ?
                    (q.y()-b.min_coord(1)) : (b.max_coord(1)-q.y());
        return std::max(dists[0], dists[1]);
    }
    FT new_distance(FT& dist, FT old_off, FT new_off,
                        int)  const {
        return dist + new_off - old_off; //works ?
    }
    FT transformed_distance(FT d) const { return d; }
    FT inverse_of_transformed_distance(FT d) { return d; }
}; // end of struct Distance
*/
inline int null_point_index() {
    return -1;
}

//}  // namespace bipartite_graph_matching

//}  // namespace Gudhi

#endif  // SRC_BOTTLENECK_INCLUDE_GUDHI_PERSISTENCE_DIAGRAMS_GRAPH_H_

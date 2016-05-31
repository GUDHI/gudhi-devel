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

namespace Gudhi {

namespace bipartite_graph_matching {

/** \internal \brief Returns the used index for encoding none of the points */
int null_point_index();

/** \internal \typedef \brief Internal_point is the internal points representation, indexes used outside. */
struct Internal_point {
    double vec[2];
    int point_index;
    Internal_point() {}
    Internal_point(double x, double y, int p_i = null_point_index()) { vec[0]=x; vec[1]=y; point_index = p_i; }
    double x() const { return vec[ 0 ]; }
    double y() const { return vec[ 1 ]; }
    double& x() { return vec[ 0 ]; }
    double& y() { return vec[ 1 ]; }
    bool operator==(const Internal_point& p) const
    {
        return point_index==p.point_index;
    }
    bool  operator!=(const Internal_point& p) const { return ! (*this == p); }
/* 
Useless
    friend void swap(Internal_point& a, Internal_point& b)
    {
        using std::swap;
        double x_tmp = a.vec[0];
        double y_tmp = a.vec[1];
        int pi_tmp = a.point_index;
        a.vec[0] = b.vec[0];
        a.vec[1] = b.vec[1];
        a.point_index = b.point_index;
        b.vec[0] =  x_tmp;
        b.vec[1] = y_tmp;
        b.point_index = pi_tmp;
    }
*/
};

inline int null_point_index() {
    return -1;
}

}  // namespace bipartite_graph_matching

}  // namespace Gudhi

#endif  // SRC_BOTTLENECK_INCLUDE_GUDHI_INTERNAL_POINT_H_

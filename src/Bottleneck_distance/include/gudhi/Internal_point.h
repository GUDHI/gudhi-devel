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

#ifndef INTERNAL_POINT_H_
#define INTERNAL_POINT_H_

namespace Gudhi {

namespace persistence_diagram {

/** \internal \brief Returns the used index for encoding none of the points */
int null_point_index();

/** \internal \typedef \brief Internal_point is the internal points representation, indexes used outside. */
struct Internal_point {
    double vec[2];
    int point_index;
    Internal_point() {}
    Internal_point(double x, double y, int p_i) { vec[0]=x; vec[1]=y; point_index = p_i; }
    double x() const { return vec[ 0 ]; }
    double y() const { return vec[ 1 ]; }
    double& x() { return vec[ 0 ]; }
    double& y() { return vec[ 1 ]; }
    bool operator==(const Internal_point& p) const
    {
        return point_index==p.point_index;
    }
    bool  operator!=(const Internal_point& p) const { return !(*this == p); }
};

inline int null_point_index() {
    return -1;
}

struct Construct_coord_iterator {
    typedef  const double* result_type;
    const double* operator()(const Internal_point& p) const
    { return p.vec; }
    const double* operator()(const Internal_point& p, int)  const
    { return p.vec+2; }
};

}  // namespace persistence_diagram

}  // namespace Gudhi





#endif  // INTERNAL_POINT_H_

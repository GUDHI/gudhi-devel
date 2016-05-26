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

#ifndef SRC_BOTTLENECK_INCLUDE_CGAL_MISCELLANEOUS_H_
#define SRC_BOTTLENECK_INCLUDE_CGAL_MISCELLANEOUS_H_

#include <gudhi/Internal_point.h>

namespace CGAL {

typedef Gudhi::bipartite_graph_matching::Internal_point Internal_point;

template <>
struct Kernel_traits<Internal_point> {
    struct Kernel {
        typedef double FT;
        typedef double RT;
    };
};


struct Construct_coord_iterator {
    typedef  const double* result_type;
    const double* operator()(const Internal_point& p) const
    { return static_cast<const double*>(p.vec); }
    const double* operator()(const Internal_point& p, int)  const
    { return static_cast<const double*>(p.vec+2); }
};

} //namespace CGAL

#endif  // SRC_BOTTLENECK_INCLUDE_CGAL_MISCELLANEOUS_H_

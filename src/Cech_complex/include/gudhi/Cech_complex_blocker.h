/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2018 Inria
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

#ifndef CECH_COMPLEX_BLOCKER_H_
#define CECH_COMPLEX_BLOCKER_H_

#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/Cech_complex.h>

#include <Miniball/Miniball.hpp>

#include <iostream>
#include <vector>
#include <cmath>  // for std::sqrt

namespace Gudhi {

namespace cech_complex {

// Just declaring Cech_complex class because used and not yet defined.
template<typename SimplicialComplexForCechComplex, typename ForwardPointRange>
class Cech_complex;

template <typename SimplicialComplexForCech, typename ForwardPointRange>
class Cech_blocker {
 private:
  using Point = std::vector<double>;
  using Point_cloud = std::vector<Point>;
  using Point_iterator = Point_cloud::const_iterator;
  using Coordinate_iterator = Point::const_iterator;
  using Min_sphere = Miniball::Miniball<Miniball::CoordAccessor<Point_iterator, Coordinate_iterator>>;
  using Simplex_handle = typename SimplicialComplexForCech::Simplex_handle;
  using Filtration_value = typename SimplicialComplexForCech::Filtration_value;
  using Cech_complex = Gudhi::cech_complex::Cech_complex<SimplicialComplexForCech, ForwardPointRange>;

 public:
  bool operator()(Simplex_handle sh) {
    Point_cloud points;
    for (auto vertex : simplicial_complex_.simplex_vertex_range(sh)) {
      points.push_back(cc_ptr_->point(vertex));
#ifdef DEBUG_TRACES
      std::cout << "#(" << vertex << ")#";
#endif  // DEBUG_TRACES
    }
    Min_sphere ms(cc_ptr_->dimension(), points.begin(),points.end());
    Filtration_value radius = std::sqrt(ms.squared_radius());
#ifdef DEBUG_TRACES
    std::cout << "radius = " << radius << " - " << (radius > cc_ptr_->threshold()) << std::endl;
#endif  // DEBUG_TRACES
    simplicial_complex_.assign_filtration(sh, radius);
    return (radius > cc_ptr_->threshold());
  }
  Cech_blocker(SimplicialComplexForCech& simplicial_complex, Cech_complex* cc_ptr)
    : simplicial_complex_(simplicial_complex),
      cc_ptr_(cc_ptr) {
  }
 private:
  SimplicialComplexForCech simplicial_complex_;
  Cech_complex* cc_ptr_;
};

}  // namespace cech_complex

}  // namespace Gudhi

#endif  // CECH_COMPLEX_BLOCKER_H_

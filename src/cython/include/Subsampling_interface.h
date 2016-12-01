/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016 INRIA
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

#ifndef SUBSAMPLING_INTERFACE_H
#define	SUBSAMPLING_INTERFACE_H

#include <gudhi/choose_n_farthest_points.h>
#include <gudhi/Points_off_io.h>
#include <CGAL/Epick_d.h>

#include <vector>
#include <iostream>

namespace Gudhi {

namespace subsampling {

using Subsampling_dynamic_kernel = CGAL::Epick_d< CGAL::Dynamic_dimension_tag >;
using Subsampling_point_d = Subsampling_dynamic_kernel::Point_d;
using Subsampling_ft = Subsampling_dynamic_kernel::FT;

std::vector<std::vector<double>> subsampling_n_farthest_points(std::vector<std::vector<double>>& points, unsigned nb_points) {
  std::vector<Subsampling_point_d> input, output;
  for (auto point : points)
      input.push_back(Subsampling_point_d(point.size(), point.begin(), point.end()));
  std::vector<std::vector<double>> landmarks;
  Subsampling_dynamic_kernel k;
  choose_n_farthest_points(k, points, nb_points, std::back_inserter(landmarks));
  std::cout << "output " << landmarks.size() << std::endl;


  return landmarks;
}

}  // namespace subsampling

}  // namespace Gudhi

#endif  // SUBSAMPLING_INTERFACE_H


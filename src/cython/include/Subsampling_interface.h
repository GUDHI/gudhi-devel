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

#ifndef INCLUDE_SUBSAMPLING_INTERFACE_H_
#define INCLUDE_SUBSAMPLING_INTERFACE_H_

#include <gudhi/choose_n_farthest_points.h>
#include <gudhi/pick_n_random_points.h>
#include <gudhi/sparsify_point_set.h>
#include <gudhi/Points_off_io.h>
#include <CGAL/Epick_d.h>

#include <iostream>
#include <vector>
#include <string>

namespace Gudhi {

namespace subsampling {

using Subsampling_dynamic_kernel = CGAL::Epick_d< CGAL::Dynamic_dimension_tag >;
using Subsampling_point_d = Subsampling_dynamic_kernel::Point_d;
using Subsampling_ft = Subsampling_dynamic_kernel::FT;

// ------ choose_n_farthest_points ------
std::vector<std::vector<double>> subsampling_n_farthest_points(const std::vector<std::vector<double>>& points,
                                                               unsigned nb_points) {
  std::vector<std::vector<double>> landmarks;
  Subsampling_dynamic_kernel k;
  choose_n_farthest_points(k, points, nb_points, std::back_inserter(landmarks));

  return landmarks;
}

std::vector<std::vector<double>> subsampling_n_farthest_points(const std::vector<std::vector<double>>& points,
                                                               unsigned nb_points, unsigned starting_point) {
  std::vector<std::vector<double>> landmarks;
  Subsampling_dynamic_kernel k;
  choose_n_farthest_points(k, points, nb_points, starting_point, std::back_inserter(landmarks));

  return landmarks;
}

std::vector<std::vector<double>> subsampling_n_farthest_points_from_file(const std::string& off_file,
                                                                         unsigned nb_points) {
  Gudhi::Points_off_reader<std::vector<double>> off_reader(off_file);
  std::vector<std::vector<double>> points = off_reader.get_point_cloud();
  return subsampling_n_farthest_points(points, nb_points);
}

std::vector<std::vector<double>> subsampling_n_farthest_points_from_file(const std::string& off_file,
                                                                         unsigned nb_points, unsigned starting_point) {
    Gudhi::Points_off_reader<std::vector<double>> off_reader(off_file);
    std::vector<std::vector<double>> points = off_reader.get_point_cloud();
    return subsampling_n_farthest_points(points, nb_points, starting_point);
}

// ------ pick_n_random_points ------
std::vector<std::vector<double>> subsampling_n_random_points(const std::vector<std::vector<double>>& points,
                                                             unsigned nb_points) {
  std::vector<std::vector<double>> landmarks;
  pick_n_random_points(points, nb_points, std::back_inserter(landmarks));

  return landmarks;
}

std::vector<std::vector<double>> subsampling_n_random_points_from_file(const std::string& off_file,
                                                                       unsigned nb_points) {
  Gudhi::Points_off_reader<std::vector<double>> off_reader(off_file);
  std::vector<std::vector<double>> points = off_reader.get_point_cloud();
  return subsampling_n_random_points(points, nb_points);
}

// ------ sparsify_point_set ------
std::vector<std::vector<double>> subsampling_sparsify_points(const std::vector<std::vector<double>>& points,
                                                             double min_squared_dist) {
  std::vector<Subsampling_point_d> input, output;
  for (auto point : points)
      input.push_back(Subsampling_point_d(point.size(), point.begin(), point.end()));
  Subsampling_dynamic_kernel k;
  sparsify_point_set(k, input, min_squared_dist, std::back_inserter(output));

  std::vector<std::vector<double>> landmarks;
  for (auto point : output)
      landmarks.push_back(std::vector<double>(point.cartesian_begin(), point.cartesian_end()));
  return landmarks;
}

std::vector<std::vector<double>> subsampling_sparsify_points_from_file(const std::string& off_file,
                                                                       double min_squared_dist) {
  Gudhi::Points_off_reader<std::vector<double>> off_reader(off_file);
  std::vector<std::vector<double>> points = off_reader.get_point_cloud();
  return subsampling_sparsify_points(points, min_squared_dist);
}

}  // namespace subsampling

}  // namespace Gudhi

#endif  // INCLUDE_SUBSAMPLING_INTERFACE_H_

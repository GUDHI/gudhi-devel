/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef INCLUDE_SUBSAMPLING_INTERFACE_H_
#define INCLUDE_SUBSAMPLING_INTERFACE_H_

#include <gudhi/distance_functions.h>
#include <gudhi/choose_n_farthest_points.h>
#include <gudhi/pick_n_random_points.h>
#include <gudhi/Points_off_io.h>

// Default value is undefined
#define _GUDHI_SUBSAMPLING_USE_CGAL 0

#if __has_include(<CGAL/version.h>) && __has_include(<Eigen/Jacobi>)
#include <CGAL/version.h>  // for CGAL_VERSION_NR
#if CGAL_VERSION_NR >= 1041101000
#undef _GUDHI_SUBSAMPLING_USE_CGAL  // To remove warning
#define _GUDHI_SUBSAMPLING_USE_CGAL 1
#include <CGAL/Epick_d.h>
#include <gudhi/sparsify_point_set.h>
#endif
#endif

#include <iostream>
#include <vector>
#include <string>

namespace Gudhi {

namespace subsampling {

// ------ choose_n_farthest_points ------
std::vector<std::vector<double>> subsampling_n_farthest_points(bool metric, const std::vector<std::vector<double>>& points,
                                                               std::size_t nb_points, std::size_t starting_point = random_starting_point) {
  std::vector<std::vector<double>> landmarks;
  if (metric)
    choose_n_farthest_points_metric(Euclidean_distance(), points, nb_points, starting_point, std::back_inserter(landmarks));
  else
    choose_n_farthest_points(Euclidean_distance(), points, nb_points, starting_point, std::back_inserter(landmarks));

  return landmarks;
}

std::vector<std::vector<double>> subsampling_n_farthest_points_from_file(bool metric, const std::string& off_file,
                                                                         std::size_t nb_points, std::size_t starting_point = random_starting_point) {
    Gudhi::Points_off_reader<std::vector<double>> off_reader(off_file);
    std::vector<std::vector<double>> points = off_reader.get_point_cloud();
    return subsampling_n_farthest_points(metric, points, nb_points, starting_point);
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

#if _GUDHI_SUBSAMPLING_USE_CGAL == 1
using Subsampling_dynamic_kernel = CGAL::Epick_d< CGAL::Dynamic_dimension_tag >;
using Subsampling_point_d = Subsampling_dynamic_kernel::Point_d;

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
#else
// ------ sparsify_point_set ------
std::vector<std::vector<double>> subsampling_sparsify_points(const std::vector<std::vector<double>>&, double) {
  return std::vector<std::vector<double>>{};
}

std::vector<std::vector<double>> subsampling_sparsify_points_from_file(const std::string&, double) {
    return std::vector<std::vector<double>>{};
}
#endif

}  // namespace subsampling

}  // namespace Gudhi

#endif  // INCLUDE_SUBSAMPLING_INTERFACE_H_

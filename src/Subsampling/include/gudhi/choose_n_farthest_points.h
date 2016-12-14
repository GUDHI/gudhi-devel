/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
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

#ifndef CHOOSE_N_FARTHEST_POINTS_H_
#define CHOOSE_N_FARTHEST_POINTS_H_

#include <boost/range.hpp>

#include <gudhi/Kd_tree_search.h>

#include <gudhi/Clock.h>

#include <CGAL/Search_traits.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Fuzzy_sphere.h>

#include <iterator>
#include <algorithm>  // for sort
#include <vector>
#include <random>
#include <limits>  // for numeric_limits<>

namespace Gudhi {

namespace subsampling {

/** 
 *  \ingroup subsampling
 *  \brief Subsample by a greedy strategy of iteratively adding the farthest point from the
 *  current chosen point set to the subsampling. 
 *  The iteration starts with the landmark `starting point`.
 *  \tparam Kernel must provide a type Kernel::Squared_distance_d which is a model of the 
 *          concept <a target="_blank"
 *   href="http://doc.cgal.org/latest/Kernel_d/classKernel__d_1_1Squared__distance__d.html">Kernel_d::Squared_distance_d</a>
 *   concept.
 *  It must also contain a public member 'squared_distance_d_object' of this type.
 *  \tparam Point_range Range whose value type is Kernel::Point_d.  It must provide random-access 
 *         via `operator[]` and the points should be stored contiguously in memory.
 *  \tparam OutputIterator Output iterator whose value type is Kernel::Point_d.
 *  \details It chooses `final_size` points from a random access range `input_pts` and
 *  outputs it in the output iterator `output_it`.
 * @param[in] k A kernel object.
 * @param[in] input_pts Const reference to the input points.
 * @param[in] final_size The size of the subsample to compute.
 * @param[in] starting_point The seed in the farthest point algorithm.
 * @param[out] output_it The output iterator.
 *  
 */
template < typename Kernel,
typename Point_range,
typename OutputIterator>
void choose_n_farthest_points(Kernel const &k,
                              Point_range const &input_pts,
                              std::size_t final_size,
                              std::size_t starting_point,
                              OutputIterator output_it) {
  typename Kernel::Squared_distance_d sqdist = k.squared_distance_d_object();

  std::size_t nb_points = boost::size(input_pts);
  assert(nb_points >= final_size);

  std::size_t current_number_of_landmarks = 0;  // counter for landmarks
  const double infty = std::numeric_limits<double>::infinity();  // infinity (see next entry)
  std::vector< double > dist_to_L(nb_points, infty);  // vector of current distances to L from input_pts

  std::size_t curr_max_w = starting_point;

  for (current_number_of_landmarks = 0; current_number_of_landmarks != final_size; current_number_of_landmarks++) {
    // curr_max_w at this point is the next landmark
    *output_it++ = input_pts[curr_max_w];
    std::size_t i = 0;
    for (auto& p : input_pts) {
      double curr_dist = sqdist(p, *(std::begin(input_pts) + curr_max_w));
      if (curr_dist < dist_to_L[i])
        dist_to_L[i] = curr_dist;
      ++i;
    }
    // choose the next curr_max_w
    double curr_max_dist = 0;  // used for defining the furhest point from L
    for (i = 0; i < dist_to_L.size(); i++)
      if (dist_to_L[i] > curr_max_dist) {
        curr_max_dist = dist_to_L[i];
        curr_max_w = i;
      }
  }
}

/** 
 *  \ingroup subsampling
 *  \brief Subsample by a greedy strategy of iteratively adding the farthest point from the
 *  current chosen point set to the subsampling. 
 *  The iteration starts with a random landmark.
 *  \tparam Kernel must provide a type Kernel::Squared_distance_d which is a model of the 
 *          concept <a target="_blank"
 *   href="http://doc.cgal.org/latest/Kernel_d/classKernel__d_1_1Squared__distance__d.html">Kernel_d::Squared_distance_d</a>
 *   concept.
 *  It must also contain a public member 'squared_distance_d_object' of this type.
 *  \tparam Point_range Range whose value type is Kernel::Point_d.  It must provide random-access 
 *         via `operator[]` and the points should be stored contiguously in memory.
 *  \tparam OutputIterator Output iterator whose value type is Kernel::Point_d.
 *  \details It chooses `final_size` points from a random access range `input_pts` and
 *  outputs it in the output iterator `output_it`.
 * @param[in] k A kernel object.
 * @param[in] input_pts Const reference to the input points.
 * @param[in] final_size The size of the subsample to compute.
 * @param[out] output_it The output iterator.
 *  
 */
template < typename Kernel,
typename Point_range,
typename OutputIterator>
void choose_n_farthest_points(Kernel const& k,
                              Point_range const &input_pts,
                              unsigned final_size,
                              OutputIterator output_it) {
  // Choose randomly the first landmark
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, final_size);
  int starting_point = dis(gen);
  choose_n_farthest_points(k, input_pts, final_size, starting_point, output_it);
}

}  // namespace subsampling

}  // namespace Gudhi

#endif  // CHOOSE_N_FARTHEST_POINTS_H_

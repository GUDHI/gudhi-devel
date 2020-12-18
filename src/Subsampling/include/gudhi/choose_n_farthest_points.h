/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CHOOSE_N_FARTHEST_POINTS_H_
#define CHOOSE_N_FARTHEST_POINTS_H_

#include <boost/range.hpp>

#include <gudhi/Null_output_iterator.h>

#include <iterator>
#include <vector>
#include <random>
#include <limits>  // for numeric_limits<>

namespace Gudhi {

namespace subsampling {

/**
 *  \ingroup subsampling
 */
enum : std::size_t {
/**
 *  Argument for `choose_n_farthest_points` to indicate that the starting point should be picked randomly.
 */
  random_starting_point = std::size_t(-1)
};

/** 
 *  \ingroup subsampling
 *  \brief Subsample by a greedy strategy of iteratively adding the farthest point from the
 *  current chosen point set to the subsampling. 
 *  \details
 *  The iteration starts with the landmark `starting point` or, if `starting point==random_starting_point`,
 *  with a random landmark.
 *  It chooses `final_size` points from a random access range
 *  `input_pts` (or the number of distinct points if `final_size` is larger)
 *  and outputs them in the output iterator `output_it`. It also
 *  outputs the distance from each of those points to the set of previous
 *  points in `dist_it`.
 *  \tparam Distance must provide an operator() that takes 2 points (value type of the range)
 *  and returns their distance (or some more general proximity measure) as a `double`.
 *  \tparam Point_range Random access range of points.
 *  \tparam PointOutputIterator Output iterator whose value type is the point type.
 *  \tparam DistanceOutputIterator Output iterator for distances.
 * @param[in] dist A distance function.
 * @param[in] input_pts The input points.
 * @param[in] final_size The size of the subsample to compute.
 * @param[in] starting_point The seed in the farthest point algorithm.
 * @param[out] output_it The output iterator for points.
 * @param[out] dist_it The optional output iterator for distances.
 *
 * \warning Older versions of this function took a CGAL kernel as argument. Users need to replace `k` with
 * `k.squared_distance_d_object()` in the first argument of every call to `choose_n_farthest_points`.
 *  
 */
template < typename Distance,
typename Point_range,
typename PointOutputIterator,
typename DistanceOutputIterator = Null_output_iterator>
void choose_n_farthest_points(Distance dist,
                              Point_range const &input_pts,
                              std::size_t final_size,
                              std::size_t starting_point,
                              PointOutputIterator output_it,
                              DistanceOutputIterator dist_it = {}) {
  std::size_t nb_points = boost::size(input_pts);
  if (final_size > nb_points)
    final_size = nb_points;

  // Tests to the limit
  if (final_size < 1)
    return;

  if (starting_point == random_starting_point) {
    // Choose randomly the first landmark
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<std::size_t> dis(0, nb_points - 1);
    starting_point = dis(gen);
  }

  std::size_t current_number_of_landmarks = 0;  // counter for landmarks
  static_assert(std::numeric_limits<double>::has_infinity, "the number type needs to support infinity()");
  // FIXME: don't hard-code the type as double. For Epeck_d, we also want to handle types that do not have an infinity.
  const double infty = std::numeric_limits<double>::infinity();  // infinity (see next entry)
  std::vector< double > dist_to_L(nb_points, infty);  // vector of current distances to L from input_pts

  std::size_t curr_max_w = starting_point;

  for (current_number_of_landmarks = 0; current_number_of_landmarks != final_size; current_number_of_landmarks++) {
    // curr_max_w at this point is the next landmark
    *output_it++ = input_pts[curr_max_w];
    *dist_it++ = dist_to_L[curr_max_w];
    std::size_t i = 0;
    for (auto&& p : input_pts) {
      double curr_dist = dist(p, input_pts[curr_max_w]);
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
    // If all that remains are duplicates of points already taken, stop.
    if (curr_max_dist == 0) break;
  }
}

}  // namespace subsampling

}  // namespace Gudhi

#endif  // CHOOSE_N_FARTHEST_POINTS_H_

/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siargey Kachanovich, Marc Glisse
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - 2022/11 Glisse: New *_metric variant
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CHOOSE_N_FARTHEST_POINTS_H_
#define CHOOSE_N_FARTHEST_POINTS_H_

#include <boost/range.hpp>
#include <boost/heap/d_ary_heap.hpp>
#include <boost/unordered_set.hpp> // preferably with boost 1.79+ for Fibonacci hashing

#include <gudhi/Null_output_iterator.h>

#include <iterator>
#include <vector>
#include <utility>
#include <random>
#include <limits>

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
 *  \brief Subsample by an iterative, greedy strategy.
 *  \details
 *  The algorithm starts with the landmark `starting point` or, if `starting point==random_starting_point`,
 *  with a landmark chosen randomly from the point set.
 *  At each iteration, it finds the point farthest from the set of current landmarks,
 *  outputs this distance, and promotes the point to landmark.
 *  It stops after finding `final_size` landmarks.
 *  \tparam Distance must provide an operator() that takes 2 points (value type of the range)
 *  and returns their distance (or some more general proximity measure) as a `double`.
 *  \tparam Point_range Random access range of points.
 *  \tparam PointOutputIterator Output iterator whose value type is the point type.
 *  \tparam DistanceOutputIterator Output iterator for distances.
 * @param[in] dist Distance function.
 * @param[in] input_pts The input points.
 * @param[in] final_size The size of the subsample to compute (reduced to the number of input points if `final_size` is larger).
 * @param[in] starting_point The seed in the farthest point algorithm.
 * @param[out] output_it The output iterator where landmarks are written.
 * @param[out] dist_it The optional output iterator where the distance from a landmark to the previous landmarks is written.
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

  // FIXME: don't hard-code the type as double. For Epeck_d, we also want to handle types that do not have an infinity.
  typedef double FT;
  static_assert(std::numeric_limits<FT>::has_infinity, "the number type needs to support infinity()");

  *output_it++ = input_pts[starting_point];
  *dist_it++ = std::numeric_limits<FT>::infinity();
  if (final_size == 1) return;

  std::vector<std::size_t> points(nb_points);  // map from remaining points to indexes in input_pts
  std::vector< FT > dist_to_L(nb_points);  // vector of current distances to L from points
  for(std::size_t i = 0; i < nb_points; ++i) {
    points[i] = i;
    dist_to_L[i] = dist(input_pts[i], input_pts[starting_point]);
  }
  // The indirection through points makes the program a bit slower. Some alternatives:
  // - the original code never removed points and counted on them not
  //   reappearing because of a self-distance of 0. This causes unnecessary
  //   computations when final_size is large. It also causes trouble if there are
  //   input points at distance 0 from each other.
  // - copy input_pts and update the local copy when removing points.

  std::size_t curr_max_w = starting_point;

  for (std::size_t current_number_of_landmarks = 1; current_number_of_landmarks != final_size; current_number_of_landmarks++) {
    std::size_t latest_landmark = points[curr_max_w];
    // To remove the latest landmark at index curr_max_w, replace it
    // with the last point and reduce the length of the vector.
    std::size_t last = points.size() - 1;
    if (curr_max_w != last) {
      points[curr_max_w] = points[last];
      dist_to_L[curr_max_w] = dist_to_L[last];
    }
    points.pop_back();

    // Update distances to L.
    std::size_t i = 0;
    for (auto p : points) {
      FT curr_dist = dist(input_pts[p], input_pts[latest_landmark]);
      if (curr_dist < dist_to_L[i])
        dist_to_L[i] = curr_dist;
      ++i;
    }
    // choose the next landmark
    curr_max_w = 0;
    FT curr_max_dist = dist_to_L[curr_max_w];  // used for defining the furthest point from L
    for (i = 1; i < points.size(); i++)
      if (dist_to_L[i] > curr_max_dist) {
        curr_max_dist = dist_to_L[i];
        curr_max_w = i;
      }
    *output_it++ = input_pts[points[curr_max_w]];
    *dist_it++ = dist_to_L[curr_max_w];
  }
}


// How bad is it to use the triangle inequality with inexact double computations?
// Hopefully we still get a net-tree.

template<class FT> struct Landmark_info;
template<class FT>
struct Compare_landmark_radius {
  std::vector<Landmark_info<FT>>* landmarks_p;
  Compare_landmark_radius(std::vector<Landmark_info<FT>>* p): landmarks_p(p) {}
  bool operator()(std::size_t, std::size_t) const;
};
// I compared all the heaps in boost. Fibonacci is not bad, but d_ary is by far the fastest.
// Arity of 3 is clearly faster than 2, and speed doesn't change much when increasing the arity even more.
template<class FT>
using radius_priority_ds =
  boost::heap::d_ary_heap<std::size_t, boost::heap::arity<7>, boost::heap::compare<Compare_landmark_radius<FT>>,
                          boost::heap::mutable_<true>, boost::heap::constant_time_size<false>>;
template<class FT>
struct Landmark_info {
  std::size_t far; FT radius;
  // The points that are closer to this landmark than to other landmarks
  std::vector<std::pair<std::size_t, FT>> voronoi;
  // For a landmark A, the list of landmarks B such that picking a Voronoi
  // point of A as a new landmark might steal a Voronoi point from B, or vice versa.
  std::vector<std::pair<std::size_t, FT>> neighbors;
  // Note that above we cache the distances. This is always good for Voronoi, and for neighbors it is neutral in 2D and helps in 4D.
  typename radius_priority_ds<FT>::handle_type position_in_queue;
};
template<class FT>
bool Compare_landmark_radius<FT>::operator()(std::size_t a, std::size_t b)const{ return (*landmarks_p)[a].radius < (*landmarks_p)[b].radius; }

/**
 *  \ingroup subsampling
 *  \brief Subsample by an iterative, greedy strategy.
 *  \details
 *  This computes the same thing as `choose_n_farthest_points()`, but relies on the triangle
 *  inequality to reduce the amount of computation when the doubling dimension and spread are small.
 *  In the worst case, this can be much slower than `choose_n_farthest_points()` though.
 *  See \cite sheehy20onehop and its references for details about this algorithm.
 *  \tparam Distance must provide an operator() that takes 2 points (value type of the range)
 *  and returns their distance as a `double`. It must be a true metric (\a not squared Euclidean),
 *  the algorithm relies on the triangle inequality.
 *  \tparam Point_range Random access range of points.
 *  \tparam PointOutputIterator Output iterator whose value type is the point type.
 *  \tparam DistanceOutputIterator Output iterator for distances.
 * @param[in] dist_ Distance function.
 * @param[in] input_pts The input points.
 * @param[in] final_size The size of the subsample to compute (reduced to the number of input points if `final_size` is larger).
 * @param[in] starting_point The seed in the farthest point algorithm.
 * @param[out] output_it The output iterator where landmarks are written.
 * @param[out] dist_it The optional output iterator where the distance from a landmark to the previous landmarks is written.
 */
template < typename Distance,
typename Point_range,
typename PointOutputIterator,
typename DistanceOutputIterator = Null_output_iterator>
void choose_n_farthest_points_metric(Distance dist_,
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

  // FIXME: don't hard-code the type as double. For Epeck_d, we also want to handle types that do not have an infinity.
  typedef double FT;
  static_assert(std::numeric_limits<FT>::has_infinity, "the number type needs to support infinity()");

  *output_it++ = input_pts[starting_point];
  *dist_it++ = std::numeric_limits<FT>::infinity();
  if (final_size == 1) return;

  auto dist = [&](std::size_t a, std::size_t b){ return dist_(input_pts[a], input_pts[b]); };

  std::vector<Landmark_info<FT>> landmarks(nb_points);
  radius_priority_ds<FT> radius_priority(&landmarks);

  auto compute_radius = [&](std::size_t i)
  {
    FT r = -std::numeric_limits<FT>::infinity(); // -2 * diameter should suffice
    std::size_t jmax = -1;
    for(auto [ j, d ] : landmarks[i].voronoi) {
      if (d > r) {
        r = d;
        jmax = j;
      }
    }
    landmarks[i].radius = r;
    landmarks[i].far = jmax;
  };
  auto update_radius = [&](std::size_t i)
  {
    compute_radius(i);
    radius_priority.decrease(landmarks[i].position_in_queue);
  };

  {
    // Initialize everything with starting_point
    auto& ini = landmarks[starting_point];
    ini.voronoi.reserve(nb_points - 1);
    for (std::size_t i = 0; i < nb_points; ++i)
      if (i != starting_point)
        ini.voronoi.emplace_back(i, dist(starting_point, i));
    compute_radius(starting_point);
    ini.position_in_queue = radius_priority.push(starting_point);
  }
  // outside the loop to recycle the allocation
  std::vector<std::size_t> modified_neighbors;
  boost::unordered_set<std::size_t> l_neighbors; // Should we use an allocator?
  for (std::size_t current_number_of_landmarks = 1; current_number_of_landmarks != final_size; current_number_of_landmarks++) {
    std::size_t l_parent = radius_priority.top();
    auto& parent_info = landmarks[l_parent];
    std::size_t l = parent_info.far;
    FT radius = parent_info.radius;
    auto& info = landmarks[l];
    *output_it++ = input_pts[l];
    *dist_it++ = radius;
    l_neighbors.clear();
    modified_neighbors.clear();
    // If a Voronoi point X of A can steal a Voronoi point Y from B, then
    // BY > XY >= AB - AX - BY, so AB < AX + 2 * BY. Symmetrized.
    auto max_dist = [](FT a, FT b){ return a + b + std::max(a, b); }; // tighter than 3 * radius
    // Check if any Voronoi points from ngb need to move to l
    auto handle_neighbor_voronoi = [&](std::size_t ngb)
    {
      auto& ngb_info = landmarks[ngb];
      auto it = std::remove_if(ngb_info.voronoi.begin(), ngb_info.voronoi.end(), [&](auto wd)
          {
            std::size_t w = wd.first;
            FT d = wd.second;
            FT newd = dist(l, w);
            if (newd < d) {
              if (w != l) // w==l can only happen for ngb==l_parent
                info.voronoi.emplace_back(w, newd);
              return true;
            }
            return false;
          });
      if (it != ngb_info.voronoi.end()) { // modified, always true for ngb==l_parent
        ngb_info.voronoi.erase(it, ngb_info.voronoi.end());
        modified_neighbors.push_back(ngb);
        // We only need to recompute the radius if far was removed, which we can test here with
        //   if (dist(l, ngb_info.far) < ngb_info.radius)
        // to avoid a costly test for each w in the loop above, but it does not seem to help.
        update_radius(ngb);
        // if (ngb_info.voronoi.empty()) radius_priority.erase(ngb_info.position_in_queue);
      }
      // We could easily return true/false here to say whether anything was modified, if useful.
    };
    auto handle_neighbor_neighbors = [&](std::size_t ngb)
    {
        auto& ngb_info = landmarks[ngb];
        std::remove_if(ngb_info.neighbors.begin(), ngb_info.neighbors.end(), [&](auto near_){
            std::size_t near = near_.first;
            FT d = near_.second;
            // Conservative 3 * radius: we could use the old radii of ngb and near, but not the new ones.
            if (d <= 3 * radius) { l_neighbors.insert(near); }
            // Here it is safe to use the new radii.
            return d >= max_dist(ngb_info.radius, landmarks[near].radius);
            });
    };
    // First update the Voronoi diagram, so we can compute all the updated
    // radii before pruning neighbor lists. The main drawback is that we have
    // to store modified_neighbors, and we don't have access to the old radii
    // in handle_neighbor_neighbors.
    handle_neighbor_voronoi(l_parent);
    // Should we make this loop a remove_if? We already remove in the next loop.
    for (auto ngb_ : parent_info.neighbors) {
      std::size_t ngb = ngb_.first;
      //if(ngb_.second <= max_dist(radius, landmarks[ngb].radius)) // radius from before update_radius(l_parent)
      //if(ngb_.second <= radius + 2 * landmarks[ngb].radius) // no need to symmetrize
      // If X can steal a Voronoi point Y from B, then BY > XY >= BX - BY, so BX < 2 * BY.
      if(dist(l, ngb) < 2 * landmarks[ngb].radius)
        handle_neighbor_voronoi(ngb);
    }
    // If there are too many neighbors (of neighbors), this could be quadratic (?).
    // Testing every landmark would then be faster, linear.
    // TODO: find a good heuristic to switch to the linear case.
    for (std::size_t ngb : modified_neighbors)
      handle_neighbor_neighbors(ngb);
    // Finalize the new Voronoi cell.
    compute_radius(l);
    info.position_in_queue = radius_priority.push(l);
    // The parent is an obvious candidate to be a neighbor.
    l_neighbors.insert(l_parent);
    for (std::size_t ngb : l_neighbors) {
      FT d = dist(l, ngb);
      if (d < max_dist(info.radius, landmarks[ngb].radius)) {
        info.neighbors.emplace_back(ngb, d);
        landmarks[ngb].neighbors.emplace_back(l, d);
      }
    }
  }
}

}  // namespace subsampling

}  // namespace Gudhi

#endif  // CHOOSE_N_FARTHEST_POINTS_H_

/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clement Jamin
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

#ifndef SPARSIFY_POINT_SET_H_
#define SPARSIFY_POINT_SET_H_

#include <gudhi/Kd_tree_search.h>
#ifdef GUDHI_SUBSAMPLING_PROFILING
#include <gudhi/Clock.h>
#endif

#include <cstddef>
#include <vector>

namespace Gudhi {

namespace subsampling {

/**
 *  \ingroup subsampling
 *  \brief Outputs a subset of the input points so that the 
 *         squared distance between any two points
 *         is greater than or equal to `min_squared_dist`.
 *
 * \tparam Kernel must be a model of the <a target="_blank"
 *   href="http://doc.cgal.org/latest/Spatial_searching/classSearchTraits.html">SearchTraits</a>
 *   concept, such as the <a target="_blank"
 *   href="http://doc.cgal.org/latest/Kernel_d/classCGAL_1_1Epick__d.html">CGAL::Epick_d</a> class, which
 *   can be static if you know the ambiant dimension at compile-time, or dynamic if you don't.
 * \tparam Point_range Range whose value type is Kernel::Point_d.  It must provide random-access 
 *         via `operator[]` and the points should be stored contiguously in memory.
 * \tparam OutputIterator Output iterator whose value type is Kernel::Point_d.
 *
 * @param[in] k A kernel object.
 * @param[in] input_pts Const reference to the input points.
 * @param[in] min_squared_dist Minimum squared distance separating the output points.
 * @param[out] output_it The output iterator.
 */
template <typename Kernel, typename Point_range, typename OutputIterator>
void
sparsify_point_set(
                   const Kernel &k, Point_range const& input_pts,
                   typename Kernel::FT min_squared_dist,
                   OutputIterator output_it) {
  typedef typename Gudhi::spatial_searching::Kd_tree_search<
      Kernel, Point_range> Points_ds;

  typename Kernel::Squared_distance_d sqdist = k.squared_distance_d_object();

#ifdef GUDHI_SUBSAMPLING_PROFILING
  Gudhi::Clock t;
#endif

  Points_ds points_ds(input_pts);

  std::vector<bool> dropped_points(input_pts.size(), false);

  // Parse the input points, and add them if they are not too close to
  // the other points
  std::size_t pt_idx = 0;
  for (typename Point_range::const_iterator it_pt = input_pts.begin();
       it_pt != input_pts.end();
       ++it_pt, ++pt_idx) {
    if (dropped_points[pt_idx])
      continue;

    *output_it++ = *it_pt;

    auto ins_range = points_ds.query_incremental_nearest_neighbors(*it_pt);

    // If another point Q is closer that min_squared_dist, mark Q to be dropped
    for (auto const& neighbor : ins_range) {
      std::size_t neighbor_point_idx = neighbor.first;
      // If the neighbor is too close, we drop the neighbor
      if (neighbor.second < min_squared_dist) {
        // N.B.: If neighbor_point_idx < pt_idx,
        // dropped_points[neighbor_point_idx] is already true but adding a
        // test doesn't make things faster, so why bother?
        dropped_points[neighbor_point_idx] = true;
      } else {
        break;
      }
    }
  }

#ifdef GUDHI_SUBSAMPLING_PROFILING
  t.end();
  std::cerr << "Point set sparsified in " << t.num_seconds()
      << " seconds." << std::endl;
#endif
}

}  // namespace subsampling
}  // namespace Gudhi

#endif  // SPARSIFY_POINT_SET_H_

/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clement Jamin
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef SPARSIFY_POINT_SET_H_
#define SPARSIFY_POINT_SET_H_

#include <boost/version.hpp>
#if BOOST_VERSION < 106600
# include <boost/function_output_iterator.hpp>
#else
# include <boost/iterator/function_output_iterator.hpp>
#endif

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
 *         is greater than `min_squared_dist`.
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

#ifdef GUDHI_SUBSAMPLING_PROFILING
  Gudhi::Clock t;
#endif

  Points_ds points_ds(input_pts);

  std::vector<bool> dropped_points(input_pts.size(), false);

  // Parse the input points, and add them if they are not too close to
  // the other points
  std::size_t pt_idx = 0;
  for (auto const& pt : input_pts) {
    if (dropped_points[pt_idx++])
      continue;

    *output_it++ = pt;

    // If another point Q is closer that min_squared_dist, mark Q to be dropped
    auto drop = [&dropped_points] (std::ptrdiff_t neighbor_point_idx) { dropped_points[neighbor_point_idx] = true; };
    points_ds.all_near_neighbors2(pt, min_squared_dist, min_squared_dist, boost::make_function_output_iterator(std::ref(drop)));
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

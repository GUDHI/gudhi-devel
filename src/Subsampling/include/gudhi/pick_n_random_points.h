/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PICK_N_RANDOM_POINTS_H_
#define PICK_N_RANDOM_POINTS_H_

#ifdef GUDHI_SUBSAMPLING_PROFILING
# include <gudhi/Clock.h>
#endif

#include <boost/range/size.hpp>

#include <cstddef>
#include <random>     // random_device, mt19937
#include <algorithm>  // shuffle
#include <numeric>    // iota
#include <iterator>
#include <vector>


namespace Gudhi {

namespace subsampling {

/**
 *  \ingroup subsampling
 * \brief Subsample a point set by picking random vertices.
 *
 *  \details It chooses `final_size` distinct points from a random access range `points`
 *  and outputs them to the output iterator `output_it`.
 *  Point_container::iterator should be ValueSwappable and RandomAccessIterator.
 */
template <typename Point_container,
typename OutputIterator>
void pick_n_random_points(Point_container const &points,
                          std::size_t final_size,
                          OutputIterator output_it) {
#ifdef GUDHI_SUBSAMPLING_PROFILING
  Gudhi::Clock t;
#endif

  std::random_device rd;
  std::mt19937 g(rd());

#if __cplusplus >= 201703L
  std::sample(std::begin(points), std::end(points), output_it, final_size, g);
#else
  std::size_t nbP = boost::size(points);
  if (final_size > nbP)
      final_size = nbP;

  std::vector<int> landmarks(nbP);
  std::iota(landmarks.begin(), landmarks.end(), 0);

  std::shuffle(landmarks.begin(), landmarks.end(), g);
  landmarks.resize(final_size);

  for (int l : landmarks)
    *output_it++ = points[l];
#endif

#ifdef GUDHI_SUBSAMPLING_PROFILING
  t.end();
  std::cerr << "Random landmark choice took " << t.num_seconds()
      << " seconds." << std::endl;
#endif
}

}  // namespace subsampling

}  // namespace Gudhi

#endif  // PICK_N_RANDOM_POINTS_H_

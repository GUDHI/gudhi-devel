/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2016 Inria
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

#ifndef PICK_N_RANDOM_POINTS_H_
#define PICK_N_RANDOM_POINTS_H_

#include <gudhi/Clock.h>

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

  std::size_t nbP = boost::size(points);
  if (final_size > nbP)
      final_size = nbP;

  std::vector<int> landmarks(nbP);
  std::iota(landmarks.begin(), landmarks.end(), 0);

  std::random_device rd;
  std::mt19937 g(rd());

  std::shuffle(landmarks.begin(), landmarks.end(), g);
  landmarks.resize(final_size);

  for (int l : landmarks)
    *output_it++ = points[l];

#ifdef GUDHI_SUBSAMPLING_PROFILING
  t.end();
  std::cerr << "Random landmark choice took " << t.num_seconds()
      << " seconds." << std::endl;
#endif
}

}  // namespace subsampling

}  // namespace Gudhi

#endif  // PICK_N_RANDOM_POINTS_H_

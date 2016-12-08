/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2015  INRIA
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

#ifndef CONSTRUCT_CLOSEST_LANDMARK_TABLE_H_
#define CONSTRUCT_CLOSEST_LANDMARK_TABLE_H_

#include <boost/range/size.hpp>

#include <queue>  // for priority_queue<>
#include <utility>  // for pair<>
#include <iterator>
#include <vector>
#include <set>

namespace Gudhi {

namespace witness_complex {

  /**
   *  \ingroup witness_complex
   * \brief Construct the closest landmark tables for all witnesses.
   *  \details Output a table 'knn', each line of which represents a witness and
   *   consists of landmarks sorted by
   *   euclidean distance from the corresponding witness.
   *
   *  The type WitnessContainer is a random access range and
   *  the type LandmarkContainer is a range.
   *  The type KNearestNeighbors can be seen as 
   *  Witness_range<Closest_landmark_range<Vertex_handle>>, where
   *  Witness_range and Closest_landmark_range are random access ranges and
   *  Vertex_handle is the label type of a vertex in a simplicial complex.
   *  Closest_landmark_range needs to have push_back operation.
   */

  template <typename FiltrationValue,
            typename WitnessContainer,
            typename LandmarkContainer,
            typename KNearestNeighbours>
  void construct_closest_landmark_table(WitnessContainer const &points,
                                        LandmarkContainer const &landmarks,
                                        KNearestNeighbours &knn) {
    int nbP = boost::size(points);
    assert(nbP >= boost::size(landmarks));

    int dim = boost::size(*std::begin(points));
    typedef std::pair<double, int> dist_i;
    typedef bool (*comp)(dist_i, dist_i);
    knn = KNearestNeighbours(nbP);
    for (int points_i = 0; points_i < nbP; points_i++) {
      std::priority_queue<dist_i, std::vector<dist_i>, comp> l_heap([](dist_i j1, dist_i j2) {
          return j1.first > j2.first;
        });
      typename LandmarkContainer::const_iterator landmarks_it;
      int landmarks_i = 0;
      for (landmarks_it = landmarks.begin(), landmarks_i = 0; landmarks_it != landmarks.end();
           ++landmarks_it, landmarks_i++) {
        dist_i dist = std::make_pair(Euclidean_distance()(points[points_i], *landmarks_it),
                                     landmarks_i);
        l_heap.push(dist);
      }
      for (int i = 0; i < dim + 1; i++) {
        dist_i dist = l_heap.top();
        knn[points_i].push_back(dist.second);
        l_heap.pop();
      }
    }
  }

}  // namespace witness_complex

}  // namespace Gudhi

#endif  // CONSTRUCT_CLOSEST_LANDMARK_TABLE_H_

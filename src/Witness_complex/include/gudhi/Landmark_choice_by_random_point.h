/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2015  INRIA Sophia Antipolis-Méditerranée (France)
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

#ifndef LANDMARK_CHOICE_BY_RANDOM_POINT_H_
#define LANDMARK_CHOICE_BY_RANDOM_POINT_H_

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
   * \brief Landmark choice strategy by taking random vertices for landmarks.
   *  \details It chooses nbL distinct landmarks from a random access range `points`
   *  and outputs a matrix {witness}*{closest landmarks} in knn.
   *
   *  The type KNearestNeighbors can be seen as 
   *  Witness_range<Closest_landmark_range<Vertex_handle>>, where
   *  Witness_range and Closest_landmark_range are random access ranges and
   *  Vertex_handle is the label type of a vertex in a simplicial complex.
   *  Closest_landmark_range needs to have push_back operation.
   */

  template <typename KNearestNeighbours,
            typename Point_random_access_range>
  void landmark_choice_by_random_point(Point_random_access_range const &points,
                                       int nbL,
                                       KNearestNeighbours &knn) {
    int nbP = boost::size(points);
    assert(nbP >= nbL);
    std::set<int> landmarks;
    int current_number_of_landmarks = 0;  // counter for landmarks

    // TODO(SK) Consider using rand_r(...) instead of rand(...) for improved thread safety
    int chosen_landmark = rand() % nbP;
    for (current_number_of_landmarks = 0; current_number_of_landmarks != nbL; current_number_of_landmarks++) {
      while (landmarks.find(chosen_landmark) != landmarks.end())
        chosen_landmark = rand() % nbP;
      landmarks.insert(chosen_landmark);
    }

    int dim = boost::size(*std::begin(points));
    typedef std::pair<double, int> dist_i;
    typedef bool (*comp)(dist_i, dist_i);
    knn = KNearestNeighbours(nbP);
    for (int points_i = 0; points_i < nbP; points_i++) {
      std::priority_queue<dist_i, std::vector<dist_i>, comp> l_heap([](dist_i j1, dist_i j2) {
          return j1.first > j2.first;
        });
      std::set<int>::iterator landmarks_it;
      int landmarks_i = 0;
      for (landmarks_it = landmarks.begin(), landmarks_i = 0; landmarks_it != landmarks.end();
           ++landmarks_it, landmarks_i++) {
        dist_i dist = std::make_pair(euclidean_distance(points[points_i], points[*landmarks_it]), landmarks_i);
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

#endif  // LANDMARK_CHOICE_BY_RANDOM_POINT_H_

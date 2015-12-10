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

#ifndef GUDHI_LANDMARK_CHOICE_BY_RANDOM_POINT_H_
#define GUDHI_LANDMARK_CHOICE_BY_RANDOM_POINT_H_

/** 
 *  \class Landmark_choice_by_random_point
 *  \brief The class `Landmark_choice_by_random_point` allows to construct the matrix
 *  of closest landmarks per witness by iteratively choosing a random non-chosen witness
 *  as a new landmark. 
 *  \ingroup witness_complex
 */

class Landmark_choice_by_random_point {

public:
  
    /** \brief Landmark choice strategy by taking random vertices for landmarks.
     *  \details It chooses nbL distinct landmarks from a random access range `points`
     *  and outputs a matrix {witness}*{closest landmarks} in knn.
     */
    
    template <typename KNearestNeighbours,
              typename Point_random_access_range>                 
    Landmark_choice_by_random_point(Point_random_access_range &points, int nbL, KNearestNeighbours &knn)
    {
      int nbP = points.end() - points.begin();
      std::set<int> landmarks;
      int current_number_of_landmarks=0;                        // counter for landmarks 

      int chosen_landmark = rand()%nbP;
      for (current_number_of_landmarks = 0; current_number_of_landmarks != nbL; current_number_of_landmarks++)
        {
          while (landmarks.find(chosen_landmark) != landmarks.end())
            chosen_landmark = rand()% nbP;
          landmarks.insert(chosen_landmark);
        }

      int dim = points.begin()->size();
      typedef std::pair<double,int> dist_i;
      typedef bool (*comp)(dist_i,dist_i);
      knn = KNearestNeighbours(nbP);
      for (int points_i = 0; points_i < nbP; points_i++)
        {
          std::priority_queue<dist_i, std::vector<dist_i>, comp> l_heap([&](dist_i j1, dist_i j2){return j1.first > j2.first;});
          std::set<int>::iterator landmarks_it;
          int landmarks_i = 0;
          for (landmarks_it = landmarks.begin(), landmarks_i=0; landmarks_it != landmarks.end(); landmarks_it++, landmarks_i++)
            {
              dist_i dist = std::make_pair(euclidean_distance(points[points_i],points[*landmarks_it]), landmarks_i);
              l_heap.push(dist);
            }
          for (int i = 0; i < dim+1; i++)
            {
              dist_i dist = l_heap.top();
              knn[points_i].push_back(dist.second);
              l_heap.pop();
            }
        }
    }

};

#endif

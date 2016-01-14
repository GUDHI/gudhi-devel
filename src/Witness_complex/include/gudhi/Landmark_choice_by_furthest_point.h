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

#ifndef GUDHI_LANDMARK_CHOICE_BY_FURTHEST_POINT_H_
#define GUDHI_LANDMARK_CHOICE_BY_FURTHEST_POINT_H_

/**
 *  \class Landmark_choice_by_furthest_point
 *  \brief The class `Landmark_choice_by_furthest_point` allows to construct the matrix
 *  of closest landmarks per witness by iteratively choosing the furthest witness
 *  from the set of already chosen landmarks as the new landmark. 
 *  \ingroup witness_complex
 */

class Landmark_choice_by_furthest_point {

private:
  typedef std::vector<int> typeVectorVertex;
  
public:
  
/** 
 *  \brief Landmark choice strategy by iteratively adding the furthest witness from the
 *  current landmark set as the new landmark. 
 *  \details It chooses nbL landmarks from a random access range `points` and
 *  writes {witness}*{closest landmarks} matrix in `knn`.
 */

    template <typename KNearestNeighbours,
              typename Point_random_access_range>
    Landmark_choice_by_furthest_point(Point_random_access_range &points,
                                      int nbL,
                                      KNearestNeighbours &knn)
    {
      int nb_points = points.end() - points.begin();
      std::vector<std::vector<double>> wit_land_dist(nb_points, std::vector<double>());    // distance matrix witness x landmarks
      typeVectorVertex  chosen_landmarks;                       // landmark list
      
      knn = KNearestNeighbours(nb_points, std::vector<int>());                             
      int current_number_of_landmarks=0;                        // counter for landmarks 
      double curr_max_dist = 0;                                 // used for defining the furhest point from L
      const double infty = std::numeric_limits<double>::infinity(); // infinity (see next entry)
      std::vector< double > dist_to_L(nb_points,infty);         // vector of current distances to L from points
      //int dim = points.begin()->size();
      
      int rand_int = rand() % nb_points;
      int curr_max_w = rand_int; //For testing purposes a pseudo-random number is used here
      
      for (current_number_of_landmarks = 0; current_number_of_landmarks != nbL; current_number_of_landmarks++)
        {
          //curr_max_w at this point is the next landmark
          chosen_landmarks.push_back(curr_max_w);
          unsigned i = 0;
          for (auto& p: points)
            {
              double curr_dist = euclidean_distance(p, *(points.begin() + chosen_landmarks[current_number_of_landmarks]));
              wit_land_dist[i].push_back(curr_dist);
              knn[i].push_back(current_number_of_landmarks);
              if (curr_dist < dist_to_L[i])
                dist_to_L[i] = curr_dist;
              ++i;
            }
          curr_max_dist = 0;
          for (i = 0; i < dist_to_L.size(); i++)
            if (dist_to_L[i] > curr_max_dist)
              {
                curr_max_dist = dist_to_L[i];
                curr_max_w = i;
              }
        }
      for (unsigned i = 0; i < points.size(); ++i)
        std::sort(knn[i].begin(),
                  knn[i].end(),
                  [&wit_land_dist, i](int a, int b)
                  { return wit_land_dist[i][a] < wit_land_dist[i][b]; });
    }

};

#endif

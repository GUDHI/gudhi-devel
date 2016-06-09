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

#ifndef LANDMARK_CHOICE_BY_RANDOM_KNN_H_
#define LANDMARK_CHOICE_BY_RANDOM_KNN_H_

#include <utility>  // for pair<>
#include <vector>
#include <cstddef> // for ptrdiff_t type

//#include <CGAL/Cartesian_d.h>
#include <CGAL/Search_traits.h>
#include <CGAL/Search_traits_adapter.h>
//#include <CGAL/property_map.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Orthogonal_incremental_neighbor_search.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Euclidean_distance.h>
//#include <CGAL/Kernel_d/Vector_d.h>
#include <CGAL/Random.h>


namespace Gudhi {

namespace witness_complex {

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> K;
typedef K::FT FT;
typedef K::Point_d Point_d;
typedef CGAL::Search_traits< FT,
                             Point_d,
                             typename K::Cartesian_const_iterator_d,
                             typename K::Construct_cartesian_const_iterator_d > Traits_base;
typedef CGAL::Euclidean_distance<Traits_base> Euclidean_distance;
typedef CGAL::Search_traits_adapter< std::ptrdiff_t,
                                     Point_d*,
                                     Traits_base> STraits;
typedef CGAL::Distance_adapter< std::ptrdiff_t,
                                Point_d*,
                                Euclidean_distance > Distance_adapter;
typedef CGAL::Orthogonal_incremental_neighbor_search< STraits,
                                                      Distance_adapter > Neighbor_search;
  typedef CGAL::Orthogonal_k_neighbor_search< STraits > Neighbor_search2;
typedef Neighbor_search::Tree Tree;

  
  /** \brief Landmark choice strategy by taking random vertices for landmarks.
   *  \details It chooses nbL distinct landmarks from a random access range `points`
   *  and outputs a matrix {witness}*{closest landmarks} in knn.
   */
  template <typename KNearestNeighbours,
            typename Point_random_access_range,
            typename Distance_matrix>
  void landmark_choice_by_random_knn(Point_random_access_range const & points,
                                     int nbL,
                                     FT alpha,
                                     unsigned limD,
                                     KNearestNeighbours & knn,
                                     Distance_matrix & distances) {
    int nbP = points.end() - points.begin();
    assert(nbP >= nbL);
    std::vector<Point_d> landmarks;
    std::vector<int> landmarks_ind;
    Point_d p;
    int chosen_landmark;
    CGAL::Random rand;
    // TODO(SK) Consider using rand_r(...) instead of rand(...) for improved thread safety
    int current_number_of_landmarks = 0;  // counter for landmarks
    for (; current_number_of_landmarks != nbL; current_number_of_landmarks++) {
      do chosen_landmark = rand.get_int(0,nbP);
      while (std::find(landmarks_ind.begin(), landmarks_ind.end(), chosen_landmark) != landmarks_ind.end());
      p = points[chosen_landmark];
      landmarks.push_back(p);
      landmarks_ind.push_back(chosen_landmark);
    }
    // std::cout << "Choice finished!" << std::endl;
    
    //int dim = points.begin()->size();
    knn = KNearestNeighbours(nbP);
    distances = Distance_matrix(nbP);
    STraits traits(&(landmarks[0]));
    CGAL::Distance_adapter<std::ptrdiff_t,Point_d*,Euclidean_distance> adapter(&(landmarks[0]));
    Euclidean_distance ed;
    Tree landmark_tree(boost::counting_iterator<std::ptrdiff_t>(0),
                       boost::counting_iterator<std::ptrdiff_t>(nbL),
                       typename Tree::Splitter(),
                       traits);
    for (int points_i = 0; points_i < nbP; points_i++) {
      Point_d const & w = points[points_i];
      Neighbor_search search(landmark_tree,
                             w,
                             FT(0),
                             true,
                             adapter);
      Neighbor_search::iterator search_it = search.begin();
      // Neighbor_search2 search(landmark_tree,
      //                         w, limD+1,
      //                         FT(0),
      //                         true,
      //                         adapter);
      // Neighbor_search2::iterator search_it = search.begin();
      
      while (knn[points_i].size() < limD) {
        distances[points_i].push_back(sqrt(search_it->second));
        knn[points_i].push_back((search_it++)->first);
      }
      FT dtow = distances[points_i][limD-1];
      
      if (alpha != 0)
        while (search_it != search.end() && search_it->second < dtow + alpha) {
          distances[points_i].push_back(sqrt(search_it->second));
          knn[points_i].push_back((search_it++)->first);
        }
      std::cout << "k = " << knn[points_i].size() << std::endl;
    }
  }

}  // namespace witness_complex

}  // namespace Gudhi

#endif  // LANDMARK_CHOICE_BY_RANDOM_POINT_H_

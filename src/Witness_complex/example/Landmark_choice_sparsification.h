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

#ifndef LANDMARK_CHOICE_BY_SPARSIFICATION_H_
#define LANDMARK_CHOICE_BY_SPARSIFICATION_H_

#include <utility>  // for pair<>
#include <vector>
#include <cstddef> // for ptrdiff_t type
#include <algorithm>

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
#include <CGAL/Fuzzy_sphere.h>

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
  typedef CGAL::Fuzzy_sphere<STraits> Fuzzy_sphere;

  /** \brief Landmark selection function by as a sub-epsilon-net of the
   *  given set of points.
   */
  template <typename Point_random_access_range>
  void landmark_choice_by_sparsification(Point_random_access_range & points,
                                         unsigned nbL,
                                         FT mu_epsilon,
                                         Point_random_access_range & landmarks)
  {
    unsigned nbP = points.end() - points.begin();
    assert(nbP >= nbL);
    CGAL::Random rand;
    // TODO(SK) Consider using rand_r(...) instead of rand(...) for improved thread safety
    STraits points_traits(&(points[0]));
    CGAL::Distance_adapter<std::ptrdiff_t,Point_d*,Euclidean_distance> points_adapter(&(points[0]));
    std::vector<bool> dropped_points(nbP, false);
    
    Tree witness_tree(boost::counting_iterator<std::ptrdiff_t>(0),
                      boost::counting_iterator<std::ptrdiff_t>(nbP),
                      typename Tree::Splitter(),
                      points_traits);
    
    for (unsigned points_i = 0; points_i < nbP; points_i++) {
      if (dropped_points[points_i])
        continue;
      Point_d & w = points[points_i];
      Fuzzy_sphere fs(w, mu_epsilon, 0, points_traits);
      std::vector<int> close_neighbors;
      witness_tree.search(std::insert_iterator<std::vector<int>>(close_neighbors,close_neighbors.begin()),fs);
      for (int i: close_neighbors)
        dropped_points[i] = true;
    }

    for (unsigned points_i = 0; points_i < nbP; points_i++) {
      if (dropped_points[points_i])
        landmarks.push_back(points[points_i]);
    }

    if (nbL < landmarks.size()) {
      std::random_shuffle(landmarks.begin(), landmarks.end());
      landmarks.resize(nbL);
    }
  }
    


  
  /** \brief Landmark choice strategy by taking random vertices for landmarks.
   *  \details It chooses nbL distinct landmarks from a random access range `points`
   *  and outputs a matrix {witness}*{closest landmarks} in knn.
   */
  template <typename KNearestNeighbours,
            typename Point_random_access_range,
            typename Distance_matrix>
  void build_distance_matrix(Point_random_access_range const & points,
                             Point_random_access_range & landmarks,
                             FT alpha,
                             unsigned limD,
                             KNearestNeighbours & knn,
                             Distance_matrix & distances)
  {
    int nbP = points.end() - points.begin();
    knn = KNearestNeighbours(nbP);
    distances = Distance_matrix(nbP);
    STraits traits(&(landmarks[0]));
    CGAL::Distance_adapter<std::ptrdiff_t,Point_d*,Euclidean_distance> adapter(&(landmarks[0]));
    Euclidean_distance ed;
    Tree landmark_tree(boost::counting_iterator<std::ptrdiff_t>(0),
                       boost::counting_iterator<std::ptrdiff_t>(landmarks.size()),
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
      
      while (knn[points_i].size() <= limD) {
        distances[points_i].push_back(search_it->second); //!sq_dist
        knn[points_i].push_back((search_it++)->first);
      }
      FT dtow = distances[points_i][limD];
      
      while (search_it != search.end() && search_it->second < dtow + alpha) {
        distances[points_i].push_back(search_it->second);
        knn[points_i].push_back((search_it++)->first);
      }
      //std::cout << "k = " << knn[points_i].size() << std::endl;
    }
  }
 
  /*
  template <typename Kernel, typename Point_container>
  std::vector<typename Point_container::value_type>
  sparsify_point_set(const Kernel &k,
                     Point_container const& input_pts,
                     typename Kernel::FT min_squared_dist)
  {
    typedef typename CGAL::Tangential_complex_::Point_cloud_data_structure<Kernel, Point_container> Points_ds;
    typedef typename Points_ds::INS_iterator      INS_iterator;
    typedef typename Points_ds::INS_range         INS_range;

  typename Kernel::Squared_distance_d sqdist = k.squared_distance_d_object();

  // Create the output container
  std::vector<typename Point_container::value_type> output;

  Points_ds points_ds(input_pts);

  std::vector<bool> dropped_points(input_pts.size(), false);

  // Parse the input points, and add them if they are not too close to
  // the other points
  std::size_t pt_idx = 0;
  for (typename Point_container::const_iterator it_pt = input_pts.begin() ;
       it_pt != input_pts.end();
       ++it_pt, ++pt_idx)
  {
    if (dropped_points[pt_idx])
      continue;

    output.push_back(*it_pt);

    INS_range ins_range = points_ds.query_incremental_ANN(*it_pt);

    // If another point Q is closer that min_squared_dist, mark Q to be dropped
    for (INS_iterator nn_it = ins_range.begin() ;
        nn_it != ins_range.end() ;
        ++nn_it)
    {
      std::size_t neighbor_point_idx = nn_it->first;
      // If the neighbor is too close, we drop the neighbor
      if (nn_it->second < min_squared_dist)
      {
        // N.B.: If neighbor_point_idx < pt_idx, 
        // dropped_points[neighbor_point_idx] is already true but adding a
        // test doesn't make things faster, so why bother?
        dropped_points[neighbor_point_idx] = true;
      }
      else
        break;
    }
  }

  return output;
}
  */  


}  // namespace witness_complex

}  // namespace Gudhi

#endif  // LANDMARK_CHOICE_BY_RANDOM_POINT_H_

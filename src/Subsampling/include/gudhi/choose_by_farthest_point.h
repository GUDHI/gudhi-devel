/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
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

#ifndef CHOOSE_BY_FARTHEST_POINT_H_
#define CHOOSE_BY_FARTHEST_POINT_H_

#include <boost/range.hpp>

#include <gudhi/Spatial_tree_data_structure.h>

#include <gudhi/Clock.h>

#include <CGAL/Search_traits.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Fuzzy_sphere.h>

#include <iterator>
#include <algorithm>  // for sort
#include <vector>
#include <random>

namespace Gudhi {

namespace subsampling {
  /** 
   *  \ingroup subsampling
   *  \brief Subsample by a greedy strategy of iteratively adding the farthest point from the
   *  current chosen point set to the subsampling. 
   *  \details It chooses `final_size` points from a random access range `points` and
   *  outputs it in the output iterator `output_it`.
   *  
   */

  template < typename Kernel,
             typename Point_container,
             typename OutputIterator>
  void choose_by_farthest_point_old( Kernel& k,
                                     Point_container const &points,
                                     int final_size,
                                     int starting_point,
                                     OutputIterator output_it)
  {
    typename Kernel::Squared_distance_d sqdist = k.squared_distance_d_object();
    
    int nb_points = boost::size(points);
    assert(nb_points >= final_size);

    int current_number_of_landmarks = 0;  // counter for landmarks
    double curr_max_dist = 0;  // used for defining the furhest point from L
    const double infty = std::numeric_limits<double>::infinity();  // infinity (see next entry)
    std::vector< double > dist_to_L(nb_points, infty);  // vector of current distances to L from points

    int curr_max_w = starting_point;

    for (current_number_of_landmarks = 0; current_number_of_landmarks != final_size; current_number_of_landmarks++) {
      // curr_max_w at this point is the next landmark
      *output_it++ = points[curr_max_w];
      // std::cout << curr_max_w << "\n";
      unsigned i = 0;
      for (auto& p : points) {
        double curr_dist = sqdist(p, *(std::begin(points) + curr_max_w));
        if (curr_dist < dist_to_L[i])
          dist_to_L[i] = curr_dist;
        ++i;
      }
      // choose the next curr_max_w
      curr_max_dist = 0;
      for (i = 0; i < dist_to_L.size(); i++)
        if (dist_to_L[i] > curr_max_dist) {
          curr_max_dist = dist_to_L[i];
          curr_max_w = i;
        }
    }
  }
  
  template < typename Kernel,
             typename Point_container,
             typename OutputIterator>
  void choose_by_farthest_point_old( Kernel& k,
                                     Point_container const &points,
                                     int final_size,
                                     OutputIterator output_it)
  {
    // Choose randomly the first landmark 
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(1, 6);
    int starting_point = dis(gen);
    choose_by_farthest_point_old(k, points, final_size, starting_point, output_it);
  }

  template < typename Kernel,
             typename Point_container,
             typename OutputIterator>
  void choose_by_farthest_point( Kernel& k,
                                 Point_container const &points,
                                 int final_size,
                                 int starting_point,
                                 OutputIterator output_it)
  {
    // typedef typename Kernel::Point_d Point_d;
    // typedef typename Kernel::FT FT;
    // typedef CGAL::Search_traits<
    //   FT, Point_d,
    //   typename Kernel::Cartesian_const_iterator_d,
    //   typename Kernel::Construct_cartesian_const_iterator_d> Traits_base;

    // typedef CGAL::Search_traits_adapter< std::ptrdiff_t, Point_d*, Traits_base > STraits;
    // typedef CGAL::Fuzzy_sphere< STraits > Fuzzy_sphere;
    
    typename Kernel::Squared_distance_d sqdist = k.squared_distance_d_object();
    
    int nb_points = boost::size(points);
    assert(nb_points >= final_size);

    Clock t;
    Gudhi::spatial_searching::Spatial_tree_data_structure< Kernel, Point_container> tree(points);
    t.end();
    //std::cout << "Constructed the Kd tree: " << t.num_seconds()  << " s." << std::endl;
  
    //CGAL::Fuzzy_sphere< CGAL::Search_trai>
    
    int current_number_of_landmarks = 0;  // counter for landmarks
    const double infty = std::numeric_limits<double>::infinity();  // infinity (see next entry)
    double curr_max_dist = infty;  // used for defining the furhest point from L
    std::vector< double > dist_to_L(nb_points, infty);  // vector of current distances to L from points
    
    // Choose randomly the first landmark 
    int curr_max_w = starting_point;

    for (current_number_of_landmarks = 0; current_number_of_landmarks != final_size; current_number_of_landmarks++) {
      // curr_max_w at this point is the next landmark
      *output_it++ = points[curr_max_w];
      // std::cout << curr_max_w << "\n";
      //for (auto& p : points) {
      auto search = tree.query_incremental_ANN(points[curr_max_w]);
      auto search_it = search.begin();
      while (search_it != search.end() && search_it->second <= curr_max_dist ) {
        //std::cout << search_it->second << " " << curr_max_dist << "\n";
        if (dist_to_L[search_it->first] > search_it->second)
          dist_to_L[search_it->first] = search_it->second;
        search_it++;
      }
      // choose the next curr_max_w
      curr_max_dist = 0;
      for (unsigned i = 0; i < dist_to_L.size(); i++)
        if (dist_to_L[i] > curr_max_dist) {
          curr_max_dist = dist_to_L[i];
          curr_max_w = i;
        }
    }
  }
  
  template < typename Kernel,
             typename Point_container,
             typename OutputIterator>
  void choose_by_farthest_point( Kernel& k,
                                 Point_container const &points,
                                 int final_size,
                                 OutputIterator output_it)
  {    
    // Choose randomly the first landmark 
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(1, 6);
    int starting_point = dis(gen);

    choose_by_farthest_point_old(k, points, final_size, starting_point, output_it);
  }

  
} // namespace subsampling
  
}  // namespace Gudhi

#endif  // CHOOSE_BY_FARTHEST_POINT_H_

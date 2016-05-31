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

#ifndef LANDMARK_CHOICE_BY_FARTHEST_POINT_H_
#define LANDMARK_CHOICE_BY_FARTHEST_POINT_H_

#include <gudhi/Spatial_tree_data_structure.h>

#include <iterator>
#include <algorithm>  // for sort
#include <vector>
#include <random>
#include <boost/heap/fibonacci_heap.hpp>

namespace Gudhi {


  template < typename Point_d,
             typename Heap,
             typename Tree,
             typename Presence_table >
  void update_heap( Point_d &l,
                    unsigned nbL,
                    Heap &heap,
                    Tree &tree,
                    Presence_table &table)
  {
    auto search = tree.query_incremental_ANN(l);
    for (auto w: search) {
      if (table[w.first].first)
        if (w.second < table[w.first].second->second) {
          heap.update(table[w.first].second, w);
        }
    }
  }
  
  /** 
   *  \ingroup witness_complex
   *  \brief Landmark choice strategy by iteratively adding the farthest witness from the
   *  current landmark set as the new landmark. 
   *  \details It chooses nbL landmarks from a random access range `points` and
   *  writes {witness}*{closest landmarks} matrix in `knn`.
   *
   *  The type KNearestNeighbors can be seen as 
   *  Witness_range<Closest_landmark_range<Vertex_handle>>, where
   *  Witness_range and Closest_landmark_range are random access ranges 
   *  
   */

  template < typename Kernel,
             typename Point_container,
             typename OutputIterator>
  void landmark_choice_by_farthest_point( Kernel& k,
                                          Point_container const &points,
                                          int nbL,
                                          OutputIterator output_it)
  {

    // typedef typename Kernel::FT FT;
    // typedef std::pair<unsigned, FT> Heap_node;
    
    // struct R_max_compare
    // {
    //   bool operator()(const Heap_node &rmh1, const Heap_node &rmh2) const
    //   {
    //     return rmh1.second < rmh2.second;
    //   }
    // };
    
    // typedef boost::heap::fibonacci_heap<Heap_node, boost::heap::compare<R_max_compare>> Heap;
    // typedef Spatial_tree_data_structure<Kernel, Point_container> Tree;
    // typedef std::vector< std::pair<bool, Heap_node*> > Presence_table;

    typename Kernel::Squared_distance_d sqdist = k.squared_distance_d_object();
    
  //   Tree tree(points);
  //   Heap heap;
  //   Presence_table table(points.size());
  //   for (auto p: table)
  //     std::cout << p.first << "\n";
  //   int number_landmarks = 0; // number of treated landmarks

  //   double curr_max_dist = 0;                                      // used for defining the furhest point from L
  //   const double infty = std::numeric_limits<double>::infinity();  // infinity (see next entry)
  //   std::vector< double > dist_to_L(points.size(), infty);         // vector of current distances to L from points
    
  //   // Choose randomly the first landmark 
  //   std::random_device rd;
  //   std::mt19937 gen(rd());
  //   std::uniform_int_distribution<> dis(1, 6);
  //   int curr_landmark = dis(gen);
    
  //   do {
  //     *output_landmarks++ = points[curr_landmark];
  //     std::cout << curr_landmark << "\n";
  //     number_landmarks++;
  //   }
  //   while (number_landmarks < nbL);
  // }

    int nb_points = boost::size(points);
    assert(nb_points >= nbL);

    int current_number_of_landmarks = 0;  // counter for landmarks
    double curr_max_dist = 0;  // used for defining the furhest point from L
    const double infty = std::numeric_limits<double>::infinity();  // infinity (see next entry)
    std::vector< double > dist_to_L(nb_points, infty);  // vector of current distances to L from points

    // Choose randomly the first landmark 
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(1, 6);
    int curr_max_w = dis(gen);

    
    for (current_number_of_landmarks = 0; current_number_of_landmarks != nbL; current_number_of_landmarks++) {
      // curr_max_w at this point is the next landmark
      *output_it++ = points[curr_max_w];
      std::cout << curr_max_w << "\n";
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
  
}  // namespace Gudhi

#endif  // LANDMARK_CHOICE_BY_FARTHEST_POINT_H_

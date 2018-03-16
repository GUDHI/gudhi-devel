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

#ifndef SIMPLEX_WITH_COFACES_H_
#define SIMPLEX_WITH_COFACES_H_

#include <forward_list>
#include <map>
#include <vector>
#include <algorithm>
// #include <utility>
// #include "gudhi/reader_utils.h"
// #include "gudhi/distance_functions.h"
// #include <gudhi/Dim_list_iterator.h>
// #include <set>
// #include <queue>
// #include <limits>
// #include <math.h>
// #include <ctime>
// #include <iostream>

template <class Simplex_with_cofaces_>
struct Coface_compare {
  typedef typename std::map<std::vector<std::size_t>, std::pair<Simplex_with_cofaces_, double> > Simplex_map;
  typedef typename Simplex_map::iterator Map_iterator;
  bool operator()(const Map_iterator& lhs, const Map_iterator& rhs) const { 
    return lhs->second.second > rhs->second.second ||
      (lhs->second.second == rhs->second.second && std::lexicographical_compare(rhs->first.begin(),
                                                                                rhs->first.end(),
                                                                                lhs->first.begin(),
                                                                                lhs->first.end()));
  }
};

template <class Set_of_set_iterators_>
struct Filtered_range_compare {
  typedef typename Set_of_set_iterators_::iterator Set_iterator;
  bool operator()(const Set_iterator& lhs, const Set_iterator& rhs) const { 
    return (*lhs)->second.second > (*rhs)->second.second ||
      ((*lhs)->second.second == (*rhs)->second.second && std::lexicographical_compare((*rhs)->first.begin(),
                                                                                (*rhs)->first.end(),
                                                                                (*lhs)->first.begin(),
                                                                                (*lhs)->first.end()));
  }
};

class Simplex_with_cofaces {
  
public:
  typedef std::map<std::vector<std::size_t>, std::pair<Simplex_with_cofaces, double> > Simplex_map;
  typedef typename Simplex_map::iterator Map_iterator;
  typedef std::vector<Map_iterator> Filtered_range;
  typedef typename Filtered_range::iterator Filtered_range_iterator;
  // typedef std::set<Filtered_range_iterator, Filtered_range_compare<Filtered_range> > Set_of_facets;
  typedef std::forward_list<Map_iterator> List_of_facets;
  typedef typename std::map<Map_iterator, List_of_facets::iterator, Coface_compare<Simplex_with_cofaces> > Map_of_cofaces;
  
private:
  // std::vector<std::size_t> vertices_;
  List_of_facets facets_;
  Map_of_cofaces cofaces_;
  
public:

  Simplex_with_cofaces() {}
  
  // Simplex_with_cofaces(std::vector<std::size_t>& vertices, typename Set_of_facets::iterator face_it)
  //   : vertices_(vertices) {
  //   cofaces_.emplace(face_it);
  // }

  void push_front_facet(Map_iterator facet_it) {
    return facets_.push_front(facet_it);
  }

  void emplace_coface(Map_iterator& coface_it, typename List_of_facets::iterator set_it) {
    cofaces_.emplace(std::make_pair(coface_it, set_it));
  }

  std::size_t number_of_cofaces() const {
    return cofaces_.size();
  }

  void remove_coface(Map_iterator& coface_it) {
    cofaces_.erase(coface_it);
  }
  
  // std::vector<std::size_t> const& vertices() const {
  //   return vertices_;
  // }

  List_of_facets& facets() {
    return facets_;
  }

  typename std::pair<Map_iterator, List_of_facets::iterator> const coface() const {
    return *(cofaces_.begin());
  }
  
}; //class Simplex_with_cofaces

#endif

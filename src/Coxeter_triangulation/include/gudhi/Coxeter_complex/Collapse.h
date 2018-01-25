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

#ifndef COLLAPSE_H_
#define COLLAPSE_H_

#include <map>
#include "Simplex_with_cofaces.h"

/* The simplices are sorted by dimension in reverse order */
class Collapse {

private:

  typedef std::map<std::vector<std::size_t>, Simplex_with_cofaces> Simplex_map;
  typedef typename Simplex_map::iterator Map_iterator; 
  Simplex_map* cofaces_;
  Simplex_map* facets_;  
  std::map<int, int> label_map; // Used to output contiguous labels for vertices
  int max_index = 0;
  
public:

  template <class Range_of_simplices,
            class Simplex_tree>
  Collapse(Range_of_simplices& simplices, Simplex_tree& output)
    : max_index(0) {
    // std::cout << "Started collapses.\n";
    Simplex_tree collapsed_tree;
    typename Range_of_simplices::iterator current_it = simplices.begin();
    if (current_it == simplices.end())
      return;
    unsigned d = current_it->size() - 1;
    cofaces_ = new Simplex_map();
    for (unsigned curr_dim = d; curr_dim > 0; curr_dim--) {
      while (current_it != simplices.end() && current_it->size() == curr_dim + 1)
        cofaces_->emplace(std::make_pair(*current_it++, Simplex_with_cofaces()));
      facets_  = new Simplex_map();
      auto coface_it = cofaces_->begin();
      for (; coface_it != cofaces_->end(); ++coface_it)
        for (unsigned i = 0; i <= curr_dim; i++) {
          std::vector<std::size_t> facet;
          for (unsigned j = 0; j <= curr_dim; j++)
            if (i != j)
              facet.push_back(coface_it->first[j]);
          auto fac_res = facets_->emplace(std::make_pair(facet, Simplex_with_cofaces()));
          coface_it->second.push_front_facet(fac_res.first);
          fac_res.first->second.emplace_coface(coface_it, coface_it->second.facets().begin());          
        }
      /* ERASE PART: TO FIX */
      // coface_it = cofaces_->begin();
      // while (coface_it != cofaces_->end())
      //   if (coface_it->second.number_of_cofaces() != 0) {
      //     // std::cout << "Coface " << coface_it->first << " erased\n";
      //     auto list_it = coface_it->second.facets().begin();
      //     while (list_it != coface_it->second.facets().end()) {
      //       (*list_it)->second.remove_coface(coface_it);
      //       if ((*list_it)->second.number_of_cofaces() == 0)
      //         facets_->erase(*(list_it++));
      //       else
      //         list_it++;
      //     }
      //     cofaces_->erase(coface_it++);
      //   }
      //   else
      //     coface_it++;

      // std::cout << "Coface map is ready, size = " << cofaces_->size() << "\n";
      // for (auto cf: *cofaces_) {
      //   std::cout << cf.first << ": " << cf.second.number_of_cofaces() << " cofaces\n";
      // }
      // std::cout << "Facet map is ready, size = " << facets_->size() << "\n";
      // for (auto cf: *facets_) {
      //   std::cout << cf.first << ": " << cf.second.number_of_cofaces() << " cofaces\n";
      // }
      auto facet_it = facets_->begin();
      while (facet_it != facets_->end())
        if (facet_it->second.number_of_cofaces() == 1)
          elementary_collapse(facet_it++);
        else
          facet_it++;
      for (auto cf: *cofaces_) {
        // std::cout << "Coface " << cf.first << " inserted\n";
        std::vector<std::size_t> vertices;
        for (std::size_t v: cf.first) {
          auto m_it = label_map.find(v);
          if (m_it == label_map.end()) {
            label_map.emplace(std::make_pair(v, max_index));
            vertices.push_back(max_index++);
          }
          else
            vertices.push_back(m_it->second);
        }
        output.insert_simplex_and_subfaces(vertices);
      }
      delete cofaces_;
      cofaces_ = facets_;
    }
    // The simplex tree has vertices with non-contiguous labels
    for (auto cf: *cofaces_) {
      // std::cout << "Coface " << cf.first << " inserted\n"; 
      std::vector<std::size_t> vertices;
      for (std::size_t v: cf.first) {
        auto m_it = label_map.find(v);
        if (m_it == label_map.end()) {
          label_map.emplace(std::make_pair(v, max_index));
          vertices.push_back(max_index++);
        }
        else
          vertices.push_back(m_it->second);
      }
      output.insert_simplex_and_subfaces(vertices);
    }
    delete cofaces_;
  }

  void elementary_collapse(Map_iterator facet_it) {
    typedef typename Simplex_with_cofaces::List_of_facets List_of_facets;

    Map_iterator coface_it = facet_it->second.coface().first;
    // std::cout << "EColapse: facet = " << facet_it->first << ", coface = " << coface_it->first << "\n";
    typename List_of_facets::iterator col_it = facet_it->second.coface().second;
    List_of_facets& facet_list = coface_it->second.facets();
    auto list_it = facet_list.begin();
    while (list_it != col_it) {
      (*list_it)->second.remove_coface(coface_it);
      if ((*list_it)->second.number_of_cofaces() == 1)
        elementary_collapse(*(list_it++));
      else
        list_it++;
    }
    list_it++; // skip col_it
    while (list_it != facet_list.end())
      (*(list_it++))->second.remove_coface(coface_it);
    facets_->erase(facet_it);
    cofaces_->erase(coface_it);
  }
  
}; //class Collapse


#endif

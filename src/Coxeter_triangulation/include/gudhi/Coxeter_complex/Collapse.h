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

/* The simplices are sorted by decreasing dimension and decreasing filtration */
/* The filtration values of the input simplices should be coherent */
class Collapse {

private:

  // typedef std::pair<std::vector<std::size_t>, double> Key_pair;
  // struct Key_compare {
  //   bool operator() (const Key_pair& lhs, const Key_pair& rhs) {
  //     return lhs.second < rhs.second || (lhs.second == rhs.second && lhs.first < rhs.first);
  //   }
  // }
  // typedef std::map< Key_pair, Simplex_with_cofaces, Key_compare> Simplex_map;
  typedef std::map<std::vector<std::size_t>, std::pair<Simplex_with_cofaces, double> > Simplex_map;
  typedef typename Simplex_map::iterator Map_iterator; 
  typedef std::vector<Map_iterator> Filtered_range;
  typedef typename Filtered_range::iterator Filtered_range_iterator;
  Simplex_map* cofaces_;
  Simplex_map* facets_;
  std::map<int, int> label_map; // Used to output contiguous labels for vertices
  int max_index = 0;
  
public:
  // debug
  int max_simplices_out = 0;

  template <class Range_of_simplices,
            class Simplex_tree>
  Collapse(Range_of_simplices& simplices, Simplex_tree& output)
    : max_index(0) {
    // std::cout << "Started collapses.\n";
    Simplex_tree collapsed_tree;
    typename Range_of_simplices::iterator current_it = simplices.begin();
    if (current_it == simplices.end())
      return;
    unsigned d = current_it->first.size() - 1;
    cofaces_ = new Simplex_map();
    for (unsigned curr_dim = d; curr_dim > 0; curr_dim--) {
      while (current_it != simplices.end() && current_it->first.size() == curr_dim + 1) {
        // std::cout << "Adding coface " << current_it->first << " filtration " << current_it->second << "\n";
        cofaces_->emplace(std::make_pair(current_it->first, std::make_pair(Simplex_with_cofaces(), current_it->second)));
        // if (ins_res.second)
        //   coface_filtered_range_.push_back(ins_res.first);
        current_it++;
      }
      Filtered_range coface_filtered_range_;
      for (auto cf_it = cofaces_->begin(); cf_it != cofaces_->end(); ++cf_it) {
        // std::cout << "Coface " << cf_it->first << " added to cofaces_\n";
        coface_filtered_range_.push_back(cf_it);
      }
      
      std::sort(coface_filtered_range_.begin(), coface_filtered_range_.end(), Coface_compare<Simplex_with_cofaces>());
      // for (auto cf: coface_filtered_range_) 
      //   std::cout << cf->first << ": " << cf->second.second << std::endl;    
      // std::cout << "\n";

      facets_  = new Simplex_map();
      Filtered_range facet_filtered_range_;
      for (auto coface_it: coface_filtered_range_) {
        // std::cout << "\nAdding facets of " << coface_it->first << " filtration " << coface_it->second.second << "\n";
        for (unsigned i = 0; i <= curr_dim; i++) {
          std::vector<std::size_t> facet;
          for (unsigned j = 0; j <= curr_dim; j++)
            if (i != j)
              facet.push_back(coface_it->first[j]);
          // std::cout << "Facet " << facet << "\n";
          auto fac_res = facets_->emplace(std::make_pair(facet, std::make_pair(Simplex_with_cofaces(), 0)));
          // if (fac_res.second)
          //   std::cout << "Successfully added " << facet << " to facets_\n";
          fac_res.first->second.second = coface_it->second.second;
          // std::cout << "Filtration value of " << facet << " is updated to " << fac_res.first->second.second << "\n";
          if (fac_res.second)
            facet_filtered_range_.push_back(fac_res.first);
          coface_it->second.first.push_front_facet(fac_res.first);
          fac_res.first->second.first.emplace_coface(coface_it, coface_it->second.first.facets().begin());
        }
      }
      std::sort(facet_filtered_range_.begin(), facet_filtered_range_.end(), Coface_compare<Simplex_with_cofaces>());
      
      auto facet_it = facet_filtered_range_.begin();
      while (facet_it != facet_filtered_range_.end()) {
        // std::cout << "Attempting to collapse " << (*facet_it)->first << " filtration " << (*facet_it)->second.second << "\n";
        if ((*facet_it)->second.first.number_of_cofaces() == 1 && (*facet_it)->second.second == (*facet_it)->second.first.coface().first->second.second)
          elementary_collapse(*facet_it++);
        else
          facet_it++;
      }
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
        if (!cf.second.first.number_of_cofaces())
          max_simplices_out++;
        // std::cout << "Coface after relabeling: " << vertices << "\n"; 
        output.insert_simplex_and_subfaces(vertices, cf.second.second);
      }
      delete cofaces_;
      cofaces_ = facets_;
    }
    // The simplex tree has vertices with non-contiguous labels
    for (auto cf: *cofaces_) {
      // std::cout << "Coface " << cf.first << " to be inserted\n"; 
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
      if (!cf.second.first.number_of_cofaces())
        max_simplices_out++;
      // std::cout << "Coface after relabeling: " << vertices << "\n"; 
      output.insert_simplex_and_subfaces(vertices, cf.second.second);
    }
    delete cofaces_;
    // std::cout << "Number of max simplices after collapse: " << max_simplices_out << "\n";
  }

  void elementary_collapse(Map_iterator facet_it) {
    typedef typename Simplex_with_cofaces::List_of_facets List_of_facets;

    Map_iterator coface_it = facet_it->second.first.coface().first;
    // std::cout << "EColapse: facet = " << facet_it->first << ", coface = " << coface_it->first << "\n";
    typename List_of_facets::iterator col_it = facet_it->second.first.coface().second;
    List_of_facets& facet_list = coface_it->second.first.facets();
    auto list_it = facet_list.begin();
    while (list_it != facet_list.end())
      if (list_it != col_it) {
        (*list_it)->second.first.remove_coface(coface_it);
        if (Coface_compare<Simplex_with_cofaces>()(*list_it, *col_it))
          if ((*list_it)->second.first.number_of_cofaces() == 1 && (*list_it)->second.second == (*list_it)->second.first.coface().first->second.second)
            elementary_collapse(*list_it++);
          else 
            list_it++;
        else
          list_it++;
      }
      else
        list_it++;
    facets_->erase(facet_it);
    cofaces_->erase(coface_it);
  }
  
}; //class Collapse


#endif

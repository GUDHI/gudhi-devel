/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author:       Mathieu Carriere
 *
 *    Copyright (C) 2017  INRIA
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

#ifndef GIC_H_
#define GIC_H_

#include <gudhi/Debug_utils.h>
#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/reader_utils.h>
#include <gudhi/Simplex_tree.h>

#include <boost/graph/adjacency_list.hpp>

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <limits>  // for numeric_limits
#include <utility>  // for pair<>


namespace Gudhi {

namespace graph_induced_complex {


/**
 * \class Graph_induced_complex
 * \brief Graph induced complex data structure.
 *
 * \ingroup graph_induced_complex
 *
 * \details
 *
 *
 */

class Graph_induced_complex {

 private:
  typedef int Cover_t;

 private:
   std::vector<std::vector<Cover_t> > cliques;

 public:
   template<typename SimplicialComplexForGIC>
   void create_complex(SimplicialComplexForGIC & complex) {
     size_t sz = cliques.size();
     for(int i = 0; i < sz; i++)  complex.insert_simplex_and_subfaces(cliques[i]);
   }

 public:
   void find_all_simplices(std::vector<std::vector<Cover_t> > & cliques, const std::vector<std::vector<Cover_t> > & cover_elts, int & token, std::vector<Cover_t> & simplex_tmp){
     int num_nodes = cover_elts.size();
     if(token == num_nodes-1){
       int num_clus = cover_elts[token].size();
       for(int i = 0; i < num_clus; i++){
         std::vector<Cover_t> simplex = simplex_tmp; simplex.push_back(cover_elts[token][i]);
         cliques.push_back(simplex);
       }
     }
     else{
       int num_clus = cover_elts[token].size();
       for(int i = 0; i < num_clus; i++){
         std::vector<Cover_t> simplex = simplex_tmp; simplex.push_back(cover_elts[token][i]);
         find_all_simplices(cliques, cover_elts, ++token, simplex);
       }
     }
   }

 public:
   /** \brief Graph_induced_complex constructor from a graph and a cover.
    *
    * @param[in] graph built on point cloud.
    * @param[in] cover of points.
    *
    * \tparam Cover must be a range for which `std::begin` and `std::end` return input iterators on a
    * `Cover_value`.
    *
    */
   template<typename Cover>
   Graph_induced_complex(Simplex_tree & st, const Cover& C, const int& max_dim) {

     // Construct the Simplex Tree corresponding to the graph
     st.expansion(max_dim);

     // Find complexes of GIC
     cliques.clear();
     for (auto simplex : st.complex_simplex_range()) {
       std::vector<std::vector<Cover_t> > cover_elts;
       for (auto vertex : st.simplex_vertex_range(simplex)) {
         cover_elts.push_back(C[vertex]);
         find_all_simplices(cliques,cover_elts);
       }
     }
   }

};

}

}

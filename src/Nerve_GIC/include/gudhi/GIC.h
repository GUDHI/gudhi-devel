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
#include <gudhi/Rips_complex.h>
#include <gudhi/Points_off_io.h>
#include <gudhi/distance_functions.h>

#include <boost/graph/adjacency_list.hpp>

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <limits>  // for numeric_limits
#include <utility>  // for pair<>
#include <functional>  // for greater<>
#include <stdexcept>
#include <initializer_list>
#include <algorithm>  // for std::max
#include <cstdint>  // for std::uint32_t

using Simplex_tree = Gudhi::Simplex_tree<>;
using Filtration_value = Simplex_tree::Filtration_value;
using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;
using Point = std::vector<float>;

bool simplex_comp(const std::vector<int>& s1, const std::vector<int>& s2){
  if(s1.size() == s2.size()){
    return s1[0] < s2[0];
  }
  else  return s1.size() < s2.size();
}

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

 private:
   std::map<int, std::vector<Cover_t> > Cov;

 private:
   int maximal_dim;

 private:
   Simplex_tree<> st;

 public:
   void set_cover(const std::string& cover_file_name){
     int vertex_id = 0; Cover_t cov; std::vector<Cover_t> cov_elts, cov_number;
     std::ifstream input(cover_file_name); std::string line;
     while(std::getline(input,line)){
       cov_elts.clear(); std::stringstream stream(line);
       while(stream >> cov){cov_elts.push_back(cov); cov_number.push_back(cov);}
       Cov.insert(std::pair<int, std::vector<Cover_t> >(vertex_id, cov_elts)); vertex_id++;
     }
     std::vector<Cover_t>::iterator it;
     std::sort(cov_number.begin(),cov_number.end()); it = std::unique(cov_number.begin(),cov_number.end());
     cov_number.resize(std::distance(cov_number.begin(),it)); maximal_dim = cov_number.size();
     return;
   }

 public:
   void set_graph_simplex_tree(const double& threshold, const std::string& off_file_name){
     Points_off_reader<Point> off_reader(off_file_name);
     Rips_complex rips_complex_from_points(off_reader.get_point_cloud(), threshold, Euclidean_distance());
     rips_complex_from_points.create_complex(st, 1);
     return;}

 public:
   template<typename SimplicialComplexForGIC>
   void create_complex(SimplicialComplexForGIC & complex) {
     size_t sz = cliques.size(); int dimension = 0;
     for(int i = 0; i < sz; i++){
       complex.insert_simplex_and_subfaces(cliques[i]);
       if(dimension < cliques[i].size()-1)  dimension = cliques[i].size()-1;
     }
     complex.set_dimension(dimension);
   }

 public:
   void find_all_simplices(std::vector<std::vector<Cover_t> > & cliques, const std::vector<std::vector<Cover_t> > & cover_elts,\
                           int & token, std::vector<Cover_t> & simplex_tmp){
     int num_nodes = cover_elts.size();
     if(token == num_nodes-1){
       int num_clus = cover_elts[token].size();
       for(int i = 0; i < num_clus; i++){
         std::vector<Cover_t> simplex = simplex_tmp; simplex.push_back(cover_elts[token][i]);
         std::vector<Cover_t>::iterator it;
         std::sort(simplex.begin(),simplex.end()); it = std::unique(simplex.begin(),simplex.end());
         simplex.resize(std::distance(simplex.begin(),it));
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
   void find_simplices() {

     // Find IDs of edges to remove
     std::vector<int> simplex_to_remove; int simplex_id = 0;
     for (auto simplex : st.complex_simplex_range()) {
       if(st.dimension(simplex) == 1){
         std::vector<std::vector<Cover_t> > comp;
         for(auto vertex : st.simplex_vertex_range(simplex))  comp.push_back(Cov[vertex]);
         if(comp[0].size() == 1 && comp[1].size() == 1 && comp[0][0] == comp[1][0])  simplex_to_remove.push_back(simplex_id);

       }
       simplex_id++;
     }

     // Remove edges
     int current_id = 0; auto simplex_tmp = st.complex_simplex_range().begin();
     if(simplex_to_remove[current_id] == 0){st.remove_maximal_simplex(*simplex_tmp); current_id++;}
     auto simplex = st.complex_simplex_range().begin();
     for(int i = 1; i < --simplex_id; i++){
         int j = i+1; auto simplex_tmp = simplex; simplex_tmp++;
         if(j == simplex_to_remove[current_id]){st.remove_maximal_simplex(*simplex_tmp); current_id++; i++;}
         simplex++;
     }

     // Build the Simplex Tree corresponding to the graph
     st.expansion(maximal_dim);

     // Find simplices of GIC
     cliques.clear();
     for (auto simplex : st.complex_simplex_range()) {
       std::vector<std::vector<Cover_t> > cover_elts; int token = 0; std::vector<Cover_t> sim;
       for (auto vertex : st.simplex_vertex_range(simplex))  cover_elts.push_back(Cov[vertex]);
       find_all_simplices(cliques,cover_elts,token,sim);
     }
     std::vector<std::vector<Cover_t> >::iterator it;
     std::sort(cliques.begin(),cliques.end()); it = std::unique(cliques.begin(),cliques.end());
     cliques.resize(std::distance(cliques.begin(),it));
   }

};

} // namespace graph_induced_complex

} // namespace Gudhi

#endif  // GIC_H_

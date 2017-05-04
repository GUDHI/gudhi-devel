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

std::map<int, double> func;

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
   std::vector<std::vector<Cover_t> > simplices;

 private:
   std::map<int, std::vector<Cover_t> > cover;

 private:
   int maximal_dim;

 private:
   std::map<Cover_t,int> cover_fct;

 private:
   Simplex_tree<> st;
   std::map<int,std::vector<int> > adjacency_matrix;

 // Simplex comparator
 private:
   bool simplex_comp(const std::vector<Cover_t>& s1, const std::vector<Cover_t>& s2){
     if(s1.size() == s2.size()){
       return s1[0] < s2[0];
     }
     else  return s1.size() < s2.size();
   }

 // Point comparator
 private:
   static bool functional_comp(const int& a, const int& b){
     if(func[a] == func[b])  return a < b;
     else  return func[a] < func[b];
   }

 // DFS
 private:
   void dfs(std::map<int,std::vector<int> >& G, const int& p, std::vector<int>& cc, std::map<int,bool>& visit){
     cc.push_back(p);
     visit[p] = true; int neighb = G[p].size();
     for (int i = 0; i < neighb; i++)
       if (  visit.find(G[p][i]) != visit.end() )
         if(  !(visit[G[p][i]])  )
           dfs(G,G[p][i],cc,visit);
   }

 // *******************************************************************************************************************
 // Graphs.
 // *******************************************************************************************************************

 public: // Set graph from file.
  void set_graph_from_file(const std::string& graph_file_name){
    int neighb; int vid; std::ifstream input(graph_file_name); std::string line; std::vector<int> edge(2);
    while(std::getline(input,line)){
      std::stringstream stream(line); stream >> vid; edge[0] = vid;
      while(stream >> neighb){edge[1] = neighb; st.insert_simplex_and_subfaces(edge);}
    }
  }

 public:
  void set_graph_from_OFF(const std::string& off_file_name){
    int numpts, numedges, numfaces, i; std::vector<int> edge(2); double x; int num; std::vector<int> simplex;
    std::ifstream input(off_file_name); std::string line; getline(input, line);
    input >> numpts; input >> numfaces; input >> numedges;
    i = 0;  while(i < numpts){input >> x; input >> x; input >> x; i++;}
    i = 0;  while(i < numfaces){
      simplex.clear(); input >> num;
      for(int j = 0; j < num; j++){int k; input >> k; simplex.push_back(k);}
      for(int j = 0; j < num; j++){
        for(int k = j+1; k < num; k++){
          edge[0] = simplex[j]; edge[1] = simplex[k];
          st.insert_simplex_and_subfaces(edge);
        }
      }
      i++;
    }

  }

 public: // Set graph from Rips complex.
   void set_graph_from_rips(const double& threshold, const std::string& off_file_name){
     Points_off_reader<Point> off_reader(off_file_name);
     Rips_complex rips_complex_from_points(off_reader.get_point_cloud(), threshold, Euclidean_distance());
     rips_complex_from_points.create_complex(st, 1);
   }


 // *******************************************************************************************************************
 // Functions.
 // *******************************************************************************************************************

 public: // Set function from file.
   void set_function_from_file(const std::string& func_file_name){
     int vertex_id = 0; std::ifstream input(func_file_name); std::string line; double f;
     while(std::getline(input,line)){
       std::stringstream stream(line); stream >> f;
       func.insert(std::pair<int,double>(vertex_id, f)); vertex_id++;
     }
   }

 public: // Set function from kth coordinate
   void set_function_from_coordinate(const int& k, const std::string& off_file_name){
     Points_off_reader<Point> off_reader(off_file_name);
     //std::vector<Point> cloud = off_reader.get_point_cloud();
     int numpts = off_reader.get_point_cloud().size(); //cloud.size();
     for(int i = 0; i < numpts; i++)  func.insert(std::pair<int,double>(i,off_reader.get_point_cloud()[i][k])); //std::cout << cloud[i][k] << std::endl;
   }

 public: // Set function from vector.
   void set_function_from_vector(const std::vector<double>& function){
     for(int i = 0; i < function.size(); i++)  func.insert(std::pair<int,double>(i, function[i]));
   }

 // *******************************************************************************************************************
 // Covers.
 // *******************************************************************************************************************

 public: // Set cover from file.
   void set_cover_from_file(const std::string& cover_file_name){
     int vertex_id = 0; Cover_t cov; std::vector<Cover_t> cov_elts, cov_number;
     std::ifstream input(cover_file_name); std::string line;
     while(std::getline(input,line)){
       cov_elts.clear(); std::stringstream stream(line);
       while(stream >> cov){cov_elts.push_back(cov); cov_number.push_back(cov);}
       cover.insert(std::pair<int, std::vector<Cover_t> >(vertex_id, cov_elts)); vertex_id++;
     }
     std::vector<Cover_t>::iterator it;
     std::sort(cov_number.begin(),cov_number.end()); it = std::unique(cov_number.begin(),cov_number.end());
     cov_number.resize(std::distance(cov_number.begin(),it)); maximal_dim = cov_number.size();
   }

 public: // Set cover with preimages of function.
   void set_cover_from_function(const double& resolution, const double& gain, const bool& token){

     // Read function values and compute min and max
     std::map<int, double>::iterator it; double maxf, minf; minf = std::numeric_limits<float>::max(); maxf = std::numeric_limits<float>::min();
     for(it = func.begin(); it != func.end(); it++){minf = std::min(minf, it->second); maxf = std::max(maxf, it->second);}
     int num_pts = func.size(); //std::cout << minf << " " << maxf << std::endl;

     // Compute cover of im(f)
     std::vector<std::pair<double,double> > intervals; int res;
     if(!token){
       double incr = (maxf-minf)/resolution; double x = minf; double alpha = (incr*gain)/(2-2*gain);
       double y = minf + incr + alpha; std::pair<double, double> interm(x,y); intervals.push_back(interm);
       for(int i = 1; i < resolution-1; i++){
         x = minf + i*incr - alpha;
         y = minf + (i+1)*incr + alpha;
         std::pair<double, double> inter(x,y); intervals.push_back(inter);
       }
       x = minf + (resolution-1)*incr - alpha; y = maxf;
       std::pair<double, double> interM(x,y); intervals.push_back(interM); res = intervals.size();
       //for(int i = 0; i < res; i++)  std::cout << intervals[i].first << " " << intervals[i].second << std::endl;
     }
     else{
       double x = minf; double y = x + resolution;
       while(y <= maxf && maxf - (y-gain*resolution) >= resolution){
         std::pair<double, double> inter(x,y); intervals.push_back(inter);
         x = y - gain*resolution;
         y = x + resolution;
       }
       std::pair<double, double> interM(x,maxf); intervals.push_back(interM); res = intervals.size();
     }

     // Sort points according to function values
     std::vector<int> points(num_pts); for(int i = 0; i < num_pts; i++)  points[i] = i;
     std::sort(points.begin(),points.end(),functional_comp);

     // Build adjacency matrix
     std::vector<int> empty;
     for(int i = 0; i < num_pts; i++)  adjacency_matrix.insert(std::pair<int,std::vector<int> >(points[i],empty));
     for (auto simplex : st.complex_simplex_range()) {
       if(st.dimension(simplex) == 1){
         std::vector<int> vertices;
         for(auto vertex : st.simplex_vertex_range(simplex))  vertices.push_back(vertex);
         adjacency_matrix[vertices[0]].push_back(vertices[1]); adjacency_matrix[vertices[1]].push_back(vertices[0]);
       }
     }

     int id = 0; int pos = 0; double min_prop_int; double max_prop_int;

     for(int i = 0; i < res; i++){       

       // Find points in the preimage
       std::map<int,std::vector<int> > prop; prop.clear();
       std::pair<double, double> inter1 = intervals[i]; min_prop_int = inter1.first;
       int tmp = pos;

       if(i != res-1){

         if(i != 0){
           std::pair<double, double> inter3 = intervals[i-1]; min_prop_int = inter3.second;
           while(func[points[tmp]] < inter3.second && tmp != num_pts){
             prop.insert(std::pair<int,std::vector<int> >(points[tmp],adjacency_matrix[points[tmp]]));
             tmp++;
           }
         }

         std::pair<double, double> inter2 = intervals[i+1]; max_prop_int = inter2.first;
         while(func[points[tmp]] < inter2.first && tmp != num_pts){
           prop.insert(std::pair<int,std::vector<int> >(points[tmp],adjacency_matrix[points[tmp]]));
           tmp++;
         }

         pos = tmp;
         while(func[points[tmp]] < inter1.second && tmp != num_pts){
           prop.insert(std::pair<int,std::vector<int> >(points[tmp],adjacency_matrix[points[tmp]]));
           tmp++;
         }

       }

       else{

         std::pair<double, double> inter3 = intervals[i-1];
         min_prop_int = inter3.second; max_prop_int = inter1.second;
         while(func[points[tmp]] < inter3.second && tmp != num_pts){
           prop.insert(std::pair<int,std::vector<int> >(points[tmp],adjacency_matrix[points[tmp]]));
           tmp++;
         }

         while(tmp != num_pts){
           prop.insert(std::pair<int,std::vector<int> >(points[tmp],adjacency_matrix[points[tmp]]));
           tmp++;
         }

       }

       // Compute the connected components with DFS
       std::map<int,bool> visit; //std::cout << i << std::endl;
       for(std::map<int, std::vector<int> >::iterator it = prop.begin(); it != prop.end(); it++)
         visit.insert(std::pair<int,bool>(it->first, false));
       if (!(prop.empty())){
         for(std::map<int, std::vector<int> >::iterator it = prop.begin(); it != prop.end(); it++){
           if (  !(visit[it->first])  ){
             std::vector<int> cc; cc.clear();
             dfs(prop,it->first,cc,visit); int cci = cc.size(); //std::cout << cci << " ";
             for(int j = 0; j < cci; j++)  cover[cc[j]].push_back(id);
             cover_fct[id] = i;
             id++;
           }
         }
       }
       //std::cout << std::endl;
     }

     maximal_dim = id;

   }


 // *******************************************************************************************************************
 // *******************************************************************************************************************


 public:
   template<typename SimplicialComplexForGIC>
   void create_complex(SimplicialComplexForGIC & complex) {
     size_t sz = simplices.size(); int dimension = 0;
     for(int i = 0; i < sz; i++){
       complex.insert_simplex_and_subfaces(simplices[i]);
       if(dimension < simplices[i].size()-1)  dimension = simplices[i].size()-1;
     }
     complex.set_dimension(dimension);
   }

 public:
   void find_all_simplices(const std::vector<std::vector<Cover_t> > & cover_elts){
     int num_nodes = cover_elts.size();
     std::vector<Cover_t> simplex;
     for(int i = 0; i < num_nodes; i++)
       for(int j = 0; j < cover_elts[i].size(); j++)
         simplex.push_back(cover_elts[i][j]);
     std::sort(simplex.begin(),simplex.end()); std::vector<Cover_t>::iterator it = std::unique(simplex.begin(),simplex.end());
     simplex.resize(std::distance(simplex.begin(),it));
     simplices.push_back(simplex); //for(int i = 0; i < simplex.size(); i++)  std::cout << simplex[i] << " "; std::cout << std::endl;
   }

 public:
   void find_Nerve_simplices(){
     simplices.clear();
     for(std::map<int,std::vector<Cover_t> >::iterator it = cover.begin(); it!=cover.end(); it++){
       simplices.push_back(it->second); //std::cout << it->second[0] << std::endl;
     }
     std::vector<std::vector<Cover_t> >::iterator it;
     std::sort(simplices.begin(),simplices.end()); it = std::unique(simplices.begin(),simplices.end());
     simplices.resize(std::distance(simplices.begin(),it));
   }

 public:
   void find_GIC_simplices() {

     // Find IDs of edges to remove
     std::vector<int> simplex_to_remove; int simplex_id = 0;
     for (auto simplex : st.complex_simplex_range()) {
       if(st.dimension(simplex) == 1){
         std::vector<std::vector<Cover_t> > comp;
         for(auto vertex : st.simplex_vertex_range(simplex))  comp.push_back(cover[vertex]);
         if(comp[0].size() == 1 && comp[1].size() == 1 && comp[0] == comp[1])  simplex_to_remove.push_back(simplex_id);
       }
       simplex_id++;
     }

     // Remove edges
     int current_id = 1;
     auto simplex = st.complex_simplex_range().begin(); int num_rem = 0;
     for(int i = 0; i < simplex_id-1; i++){
         int j = i+1; auto simplex_tmp = simplex; simplex_tmp++;
         if(j == simplex_to_remove[current_id]){st.remove_maximal_simplex(*simplex_tmp); current_id++; num_rem++;}
         else  simplex++;
     } simplex = st.complex_simplex_range().begin();
     for(int i = 0; i < simplex_to_remove[0]; i++)  simplex++;  st.remove_maximal_simplex(*simplex);

     // Build the Simplex Tree corresponding to the graph
     st.expansion(maximal_dim);

     // Find simplices of GIC
     simplices.clear();
     for (auto simplex : st.complex_simplex_range()) {
       if(!st.has_children(simplex)){
         std::vector<std::vector<Cover_t> > cover_elts;
         for (auto vertex : st.simplex_vertex_range(simplex))  cover_elts.push_back(cover[vertex]);
         find_all_simplices(cover_elts);
       }
     }
     std::vector<std::vector<Cover_t> >::iterator it;
     std::sort(simplices.begin(),simplices.end()); it = std::unique(simplices.begin(),simplices.end());
     simplices.resize(std::distance(simplices.begin(),it));
   }

 public:
   void find_GIC_simplices_with_functional_minimal_cover(){

     int v1, v2;

     // Loop on all points.
     for(std::map<int,std::vector<Cover_t> >::iterator it = cover.begin(); it != cover.end(); it++){

       int vid = it->first; std::vector<int> neighbors = adjacency_matrix[vid]; int num_neighb = neighbors.size();
       //std::cout << vid << " " << num_neighb << std::endl;

       // Find cover of current point (vid).
       if(cover[vid].size() == 2)  v1 = std::min(cover[vid][0],cover[vid][1]); else  v1 = cover[vid][0];

       // Loop on neighbors.
       for(int i = 0; i < num_neighb; i++){

         int neighb = neighbors[i];

         // Find cover of neighbor (neighb).
         if(cover[neighb].size() == 2)  v2 = std::max(cover[neighb][0],cover[neighb][1]); else  v2 = cover[neighb][0];

         // If neighbor is in next interval, add edge.
         if(cover_fct[v2] == cover_fct[v1] + 1){
           std::vector<int> edge(2); edge[0] = v1; edge[1] = v2; //std::cout << v1 << " " << v2 << std::endl;
           simplices.push_back(edge);
         }
       }
     }
     std::vector<std::vector<Cover_t> >::iterator it;
     std::sort(simplices.begin(),simplices.end()); it = std::unique(simplices.begin(),simplices.end());
     simplices.resize(std::distance(simplices.begin(),it));
   }

};

} // namespace graph_induced_complex

} // namespace Gudhi

#endif  // GIC_H_

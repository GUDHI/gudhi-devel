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
#include <random>

#define CONSTANT 10
#define ETA 0.001
#define MASK 0

using Simplex_tree = Gudhi::Simplex_tree<>;
using Filtration_value = Simplex_tree::Filtration_value;
using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;
using Point = std::vector<float>;

std::map<int, double> func;
std::map<int, double> func_color;

namespace Gudhi {

namespace graph_induced_complex {


/**
 * \class Graph_induced_complex
 * \brief Graph induced complex data structure.
 *
 * \ingroup graph_induced_complex
 *
 * \details
 * The data structure is a simplicial complex, representing a
 * Graph Induced simplicial Complex (GIC) or a Nerve,
 * and whose simplices are computed with a cover C of a point
 * cloud P, which often comes from the preimages of intervals
 * covering the image of a function f defined on P.
 * These intervals are parameterized by their resolution
 * (either their length or their number)
 * and their gain (percentage of overlap).
 * To compute a GIC, one also needs a graph G built on top of P,
 * whose cliques with vertices belonging to different elements of C
 * correspond to the simplices of the GIC.
 *
 */

class Graph_induced_complex {

 private:
   bool verbose;
   typedef int Cover_t;
   std::vector<std::vector<Cover_t> > simplices;
   std::map<int, std::vector<Cover_t> > cover;
   int maximal_dim;
   int data_dimension;
   std::map<Cover_t,int> cover_fct;
   std::map<Cover_t,std::pair<int,double> > cover_color;
   Simplex_tree<> st;
   std::map<int,std::vector<int> > adjacency_matrix;
   int resolution_int;
   double resolution_double;
   double gain;

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

 // Find random number in [0,1].
 double GetUniform(){
   static std::default_random_engine re;
   static std::uniform_real_distribution<double> Dist(0,1);
   return Dist(re);
 }

 // Subsample points.
 void SampleWithoutReplacement(int populationSize, int sampleSize, std::vector<int> & samples){
   int& n = sampleSize; int& N = populationSize;
   int t = 0; int m = 0; double u;
   while (m < n){
     u = GetUniform();
     if ( (N - t)*u >= n - m )
       t++;
      else{samples[m] = t; t++; m++;}
   }
 }

 public:
   void set_verbose(bool verb = 0){verbose = verb;}

 // *******************************************************************************************************************
 // Graphs.
 // *******************************************************************************************************************

 public: // Set graph from file.
   /** \brief Creates the graph G from a file containing the edges.
    *
    * @param[in] graph_file_name name of the input graph file.
    *
    */
   void set_graph_from_file(const std::string& graph_file_name){
     int neighb; int vid; std::ifstream input(graph_file_name); std::string line; std::vector<int> edge(2);
     while(std::getline(input,line)){
       std::stringstream stream(line); stream >> vid; edge[0] = vid;
       while(stream >> neighb){edge[1] = neighb; st.insert_simplex_and_subfaces(edge);}
     }
   }

 public: // Set graph from OFF file.
   /** \brief Creates the graph G from the triangulation given by an .OFF file.
    *
    * @param[in] off_file_name name of the input .OFF file.
    *
    */
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
   /** \brief Creates the graph G from a Rips complex.
    *
    * @param[in] threshold threshold value for the Rips complex.
    * @param[in] off_file_name name of the input .OFF file.
    *
    */
   void set_graph_from_rips(const double& threshold, const std::string& off_file_name){
     Points_off_reader<Point> off_reader(off_file_name);
     Rips_complex rips_complex_from_points(off_reader.get_point_cloud(), threshold, Euclidean_distance());
     rips_complex_from_points.create_complex(st, 1); data_dimension = off_reader.get_point_cloud()[0].size();
   }

 public: // Automatic tuning of Rips complex.
   /** \brief Creates the graph G from a Rips complex whose threshold value is automatically tuned with subsampling.
    *
    * @param[in] off_file_name name of the input .OFF file.
    * @param[in] N number of subsampling iteration (default value 100).
    *
    */
   void set_graph_from_automatic_rips(const std::string& off_file_name, int N = 100){

     Points_off_reader<Point> off_reader(off_file_name);
     int n = off_reader.get_point_cloud().size();
     int m = floor(n/pow(log(n)/log(CONSTANT),1+ETA)); m = std::min(m,n-1);
     std::vector<int> samples(m); double delta = 0; int dim = off_reader.get_point_cloud()[0].size(); data_dimension = dim;

     if(verbose)  std::cout << n << " points in R^" << dim << std::endl;
     if(verbose)  std::cout << "Subsampling " << m << " points" << std::endl;

     std::vector<std::vector<double> > dist(n); std::vector<double> dumb(n);
     for(int i = 0; i < n; i++)  dist[i] = dumb;
     double d;

     char distances[100]; sprintf(distances,"%s_dist",(char*) off_file_name.c_str());
     std::ifstream input(distances, std::ios::out | std::ios::binary);

     if(input.good()){
       if(verbose)  std::cout << "Reading distances..." << std::endl;
       for(int i = 0; i < n; i++){
         for (int j = 0; j < n; j++){
           input.read((char*) &d,8); dist[i][j] = d;
         }
       }
       input.close();
     }
     else{
       if(verbose)  std::cout << "Computing distances..." << std::endl;
       input.close(); std::ofstream output(distances, std::ios::out | std::ios::binary);
       for(int i = 0; i < n; i++){
         if( (int) floor( 100*(i*1.0+1)/(n*1.0) ) %10 == 0  )
           if(verbose)  std::cout << "\r" << floor( 100*(i*1.0+1)/(n*1.0) ) << "%" << std::flush;
         for (int j = 0; j < n; j++){
           double dis = 0; for(int k = 0; k < dim; k++)  dis += pow(off_reader.get_point_cloud()[i][k]-off_reader.get_point_cloud()[j][k],2);
           dis = std::sqrt(dis); dist[i][j] = dis; output.write((char*) &dis,8);
         }
       }
       output.close(); if(verbose)  std::cout << std::endl;
     }

     for(int i = 0; i < n; i++){
       for(int j = i; j < n; j++){
         double dis = 0; for(int k = 0; k < dim; k++)  dis += pow(off_reader.get_point_cloud()[i][k]-off_reader.get_point_cloud()[j][k],2);
         dist[i][j] = std::sqrt(dis); dist[j][i] = std::sqrt(dis);
       }
     }

     #pragma omp parallel for
     for (int i = 0; i < N; i++){

       SampleWithoutReplacement(n,m,samples);
       double hausdorff_dist = 0;
       for (int j = 0; j < n; j++){
         double mj = dist[j][samples[0]]; for (int k = 1; k < m; k++)  mj = std::min(mj, dist[j][samples[k]]);
         hausdorff_dist = std::max(hausdorff_dist, mj);
       }
       delta += hausdorff_dist/N;

     }

     if(verbose)  std::cout << "delta = " << delta << std::endl;
     Rips_complex rips_complex_from_points(off_reader.get_point_cloud(), delta, Euclidean_distance());
     rips_complex_from_points.create_complex(st, 1);

   }


 // *******************************************************************************************************************
 // Functions.
 // *******************************************************************************************************************

 public: // Set function from file.
   /** \brief Creates the function f from a file containing the function values.
    *
    * @param[in] func_file_name name of the input function file.
    *
    */
   void set_function_from_file(const std::string& func_file_name){
     int vertex_id = 0; std::ifstream input(func_file_name); std::string line; double f;
     while(std::getline(input,line)){
       std::stringstream stream(line); stream >> f;
       func.insert(std::pair<int,double>(vertex_id, f)); vertex_id++;
     }
   }

 public: // Set function from kth coordinate
   /** \brief Creates the function f from the k-th coordinate of the point cloud P.
    *
    * @param[in] k coordinate to use (start at 0).
    * @param[in] off_file_name name of the input .OFF file.
    *
    */
   void set_function_from_coordinate(const int& k, const std::string& off_file_name){
     Points_off_reader<Point> off_reader(off_file_name);
     int numpts = off_reader.get_point_cloud().size();
     for(int i = 0; i < numpts; i++)  func.insert(std::pair<int,double>(i,off_reader.get_point_cloud()[i][k]));
   }

 public: // Set function from vector.
   /** \brief Creates the function f from a vector stored in memory.
    *
    * @param[in] function input vector of values.
    *
    */
   void set_function_from_vector(const std::vector<double>& function){
     for(int i = 0; i < function.size(); i++)  func.insert(std::pair<int,double>(i, function[i]));
   }

 // *******************************************************************************************************************
 // Covers.
 // *******************************************************************************************************************

 public: // Set cover from file.
   /** \brief Creates the cover C from a file containing the cover elements of each point (the order has to be the same
    * as in the input file!).
    *
    * @param[in] cover_file_name name of the input cover file.
    *
    */
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

 public: // Automatic tuning of resolution for Mapper Delta.
   /** \brief Computes the optimal length of intervals for a Mapper Delta.
    */
   void set_automatic_resolution_for_GICMAP(){
     double reso = 0;
     for (auto simplex : st.complex_simplex_range()) {
       if(st.dimension(simplex) == 1){
         std::vector<int> vertices;
         for(auto vertex : st.simplex_vertex_range(simplex))  vertices.push_back(vertex);
         reso = std::max(reso, std::abs(func[vertices[0]] - func[vertices[1]]));
       }
     }
     if(verbose)  std::cout << "resolution = " << reso << std::endl;
     resolution_double = reso;
   }

 public:
   /** \brief Sets a length of intervals from a value stored in memory.
    *
    * @param[in] reso length of intervals.
    *
    */
   void set_resolution_double(double reso){resolution_double = reso;}
   /** \brief Sets a number of intervals from a value stored in memory.
    *
    * @param[in] reso number of intervals.
    *
    */
   void set_resolution_int(int reso){resolution_int = reso;}
   /** \brief Sets a gain from a value stored in memory (default value 0.3).
    *
    * @param[in] g gain.
    *
    */
   void set_gain(double g = 0.3){gain = g;}

 public: // Automatic tuning of resolution for Mapper Point.
   /** \brief Computes the optimal length of intervals for a standard Mapper.
    */
   void set_automatic_resolution_for_MAP(const double& gain){
     double reso = 0;
     for (auto simplex : st.complex_simplex_range()) {
       if(st.dimension(simplex) == 1){
         std::vector<int> vertices;
         for(auto vertex : st.simplex_vertex_range(simplex))  vertices.push_back(vertex);
         reso = std::max(reso, (std::abs(func[vertices[0]] - func[vertices[1]]))/gain);
       }
     }
     if(verbose)  std::cout << "resolution = " << reso << std::endl;
     resolution_double = reso;
   }

 public: // Set cover with preimages of function.
   /** \brief Creates a cover C from the preimages of the function f.
    *
    * @param[in] token boolean specifying whether we use the length or the number of intervals for the cover of im(f).
    *
    */
   void set_cover_from_function(const bool& token){

     // Read function values and compute min and max
     std::map<int, double>::iterator it; double maxf, minf; minf = std::numeric_limits<float>::max(); maxf = std::numeric_limits<float>::min();
     for(it = func.begin(); it != func.end(); it++){minf = std::min(minf, it->second); maxf = std::max(maxf, it->second);}
     int num_pts = func.size(); if(verbose)  std::cout << "Min function value = " << minf << " and Max function value = " << maxf << std::endl;

     // Compute cover of im(f)
     std::vector<std::pair<double,double> > intervals; int res;

     if(!token){ // Case we use an integer for the number of intervals.
       double incr = (maxf-minf)/resolution_int; double x = minf; double alpha = (incr*gain)/(2-2*gain);
       double y = minf + incr + alpha; std::pair<double, double> interm(x,y); intervals.push_back(interm);
       for(int i = 1; i < resolution_int-1; i++){
         x = minf + i*incr - alpha;
         y = minf + (i+1)*incr + alpha;
         std::pair<double, double> inter(x,y); intervals.push_back(inter);
       }
       x = minf + (resolution_int-1)*incr - alpha; y = maxf;
       std::pair<double, double> interM(x,y); intervals.push_back(interM); res = intervals.size();
       for(int i = 0; i < res; i++)
         if(verbose)  std::cout << "Interval " << i << " = [" << intervals[i].first << ", " << intervals[i].second << "]" << std::endl;
     }

     else{ // Case we use a double for the length of the intervals.
       double x = minf; double y = x + resolution_double;
       while(y <= maxf && maxf - (y-gain*resolution_double) >= resolution_double){
         std::pair<double, double> inter(x,y); intervals.push_back(inter);
         x = y - gain*resolution_double;
         y = x + resolution_double;
       }
       std::pair<double, double> interM(x,maxf); intervals.push_back(interM); res = intervals.size();
       for(int i = 0; i < res; i++)  if(verbose)  std::cout << "Interval " << i << " = [" << intervals[i].first << ", " << intervals[i].second << "]" << std::endl;
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
       std::map<int,bool> visit; if(verbose)  std::cout << "Preimage of interval " << i << std::endl;
       for(std::map<int, std::vector<int> >::iterator it = prop.begin(); it != prop.end(); it++)
         visit.insert(std::pair<int,bool>(it->first, false));
       if (!(prop.empty())){
         for(std::map<int, std::vector<int> >::iterator it = prop.begin(); it != prop.end(); it++){
           if (  !(visit[it->first])  ){
             std::vector<int> cc; cc.clear();
             dfs(prop,it->first,cc,visit); int cci = cc.size(); if(verbose)  std::cout << "one CC with " << cci << " points, ";
             double average_col = 0;
             for(int j = 0; j < cci; j++){cover[cc[j]].push_back(id); average_col += func_color[cc[j]]/cci;}
             cover_fct[id] = i; cover_color[id] = std::pair<int,double>(cci,average_col);
             id++;
           }
         }
       }
       if(verbose)  std::cout << std::endl;
     }

     maximal_dim = id;

   }

 // *******************************************************************************************************************
 // Visualization.
 // *******************************************************************************************************************

 public: // Set color from file.
   /** \brief Computes the function used to color the nodes of the simplicial complex from a file containing the function values.
    *
    * @param[in] color_file_name name of the input color file.
    *
    */
   void set_color_from_file(const std::string& color_file_name){
     int vertex_id = 0; std::ifstream input(color_file_name); std::string line; double f;
     while(std::getline(input,line)){
       std::stringstream stream(line); stream >> f;
       func_color.insert(std::pair<int,double>(vertex_id, f)); vertex_id++;
     }
   }

 public: // Set color from kth coordinate
   /** \brief Computes the function used to color the nodes of the simplicial complex from the k-th coordinate.
    *
    * @param[in] k coordinate to use (start at 0).
    * @param[in] off_file_name name of the input .OFF file.
    *
    */
   void set_color_from_coordinate(int k = 0, const std::string& off_file_name){
     Points_off_reader<Point> off_reader(off_file_name);
     int numpts = off_reader.get_point_cloud().size();
     for(int i = 0; i < numpts; i++)  func_color.insert(std::pair<int,double>(i,off_reader.get_point_cloud()[i][k]));
   }

 public: // Set color from vector.
   /** \brief Computes the function used to color the nodes of the simplicial complex from a vector stored in memory.
    *
    * @param[in] color input vector of values.
    *
    */
   void set_color_from_vector(const std::vector<double>& color){
     for(int i = 0; i < color.size(); i++)  func_color.insert(std::pair<int,double>(i, color[i]));
   }

 public: // Create a .dot file that can be compiled with neato to produce a .pdf file
   /** \brief Creates a .dot file for neato once the simplicial complex is computed to get a .pdf output.
    */
   void plot_with_neato(){
     char mapp[11] = "mapper.dot";   std::ofstream graphic(mapp); graphic << "graph Mapper {" << std::endl;
     double maxv, minv; maxv = std::numeric_limits<double>::min(); minv = std::numeric_limits<double>::max();
     for (std::map<Cover_t,std::pair<int,double> >::iterator iit = cover_color.begin(); iit != cover_color.end(); iit++){
       maxv = std::max(maxv, iit->second.second);  minv = std::min(minv, iit->second.second);
     }
     int k = 0; std::vector<int> nodes; nodes.clear();
     for (std::map<Cover_t,std::pair<int,double> >::iterator iit = cover_color.begin(); iit != cover_color.end(); iit++){
       if(iit->second.first > MASK){
         nodes.push_back(iit->first);
         graphic << iit->first << "[shape=circle fontcolor=black color=black label=\"" \
           << iit->first << ":" << iit->second.first << "\" style=filled fillcolor=\"" \
           << (1-(maxv-iit->second.second)/(maxv-minv))*0.6 << ", 1, 1\"]" << std::endl;
         k++;
       }
     }
     int ke = 0; int num_simplices = simplices.size();
     for (int i = 0; i < num_simplices; i++)
       if (simplices[i].size() == 2)
         if (cover_color[simplices[i][0]].first > MASK && cover_color[simplices[i][1]].first > MASK){
           graphic << "  " << simplices[i][0] << " -- " << simplices[i][1] << " [weight=15];" << std::endl; ke++;}
     graphic << "}"; graphic.close();
   }

 public: // Create a .txt file that can be compiled with KeplerMapper to produce a .html file
   /** \brief Creates a .html file for KeplerMapper once the simplicial complex is computed to get a nice visualization in browser.
    */
   void plot_with_KeplerMapper(){

     int num_simplices = simplices.size(); int num_edges = 0;
     char mapp[11] = "mapper.txt";  std::ofstream graphic(mapp);
     for (int i = 0; i < num_simplices; i++)
       if (simplices[i].size() == 2)
         if (cover_color[simplices[i][0]].first > MASK && cover_color[simplices[i][1]].first > MASK)
           num_edges++;

     graphic << "Cloud" << std::endl;
     graphic << "Function" << std::endl;
     graphic << 0 << " " << 0 << std::endl;
     graphic << cover_color.size() << " " << num_edges << std::endl;

     for (std::map<Cover_t,std::pair<int,double> >::iterator iit = cover_color.begin(); iit != cover_color.end(); iit++)
       graphic << iit->first << " " << iit->second.second << " " << iit->second.first << std::endl;

     for (int i = 0; i < num_simplices; i++)
       if (simplices[i].size() == 2)
         if (cover_color[simplices[i][0]].first > MASK && cover_color[simplices[i][1]].first > MASK)
           graphic << simplices[i][0] << " " << simplices[i][1] << std::endl;
     graphic.close();
   }

 // *******************************************************************************************************************
 // *******************************************************************************************************************


 public:
   /** \brief Creates the simplicial complex.
    *
    * \tparam SimplicialComplexForGIC must meet `SimplicialComplexForRips` concept.
    * @param[in] complex SimplicialComplexForGIC to be created.
    *
    */
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
   /** \brief Finds the maximal clique formed by different elements of the cover in a set of points.
    *
    * @param[in] cover_elts vector of points represented by vectors of cover elements (the ones to which they belong).
    *
    */
   void find_all_simplices(const std::vector<std::vector<Cover_t> > & cover_elts){
     int num_nodes = cover_elts.size();
     std::vector<Cover_t> simplex;
     for(int i = 0; i < num_nodes; i++)
       for(int j = 0; j < cover_elts[i].size(); j++)
         simplex.push_back(cover_elts[i][j]);
     std::sort(simplex.begin(),simplex.end()); std::vector<Cover_t>::iterator it = std::unique(simplex.begin(),simplex.end());
     simplex.resize(std::distance(simplex.begin(),it));
     simplices.push_back(simplex);
   }

 public:
   /** \brief Computes the simplices in the Nerve of the cover C.
    */
   void find_Nerve_simplices(){
     for(std::map<int,std::vector<Cover_t> >::iterator it = cover.begin(); it!=cover.end(); it++)  simplices.push_back(it->second);
     std::vector<std::vector<Cover_t> >::iterator it;
     std::sort(simplices.begin(),simplices.end()); it = std::unique(simplices.begin(),simplices.end());
     simplices.resize(std::distance(simplices.begin(),it));
   }

 public:
   /** \brief Computes the simplices in the GIC of the graph G and the cover C.
    */
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
   /** \brief Computes the simplices in the Mapper Delta.
    */
   void find_GICMAP_simplices_with_functional_minimal_cover(){

     int v1, v2;

     // Loop on all points.
     for(std::map<int,std::vector<Cover_t> >::iterator it = cover.begin(); it != cover.end(); it++){

       int vid = it->first; std::vector<int> neighbors = adjacency_matrix[vid]; int num_neighb = neighbors.size();

       // Find cover of current point (vid).
       if(cover[vid].size() == 2)  v1 = std::min(cover[vid][0],cover[vid][1]); else  v1 = cover[vid][0];
       std::vector<int> node(1); node[0] = v1; simplices.push_back(node);

       // Loop on neighbors.
       for(int i = 0; i < num_neighb; i++){

         int neighb = neighbors[i];

         // Find cover of neighbor (neighb).
         if(cover[neighb].size() == 2)  v2 = std::max(cover[neighb][0],cover[neighb][1]); else  v2 = cover[neighb][0];

         // If neighbor is in next interval, add edge.
         if(cover_fct[v2] == cover_fct[v1] + 1){
           std::vector<int> edge(2); edge[0] = v1; edge[1] = v2;
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

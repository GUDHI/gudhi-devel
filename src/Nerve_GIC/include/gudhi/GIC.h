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

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <limits>  // for numeric_limits
#include <utility>  // for pair<>
#include <algorithm>  // for std::max
#include <random>
#include <cassert>

namespace Gudhi {

namespace cover_complex {

using Simplex_tree = Gudhi::Simplex_tree<>;
using Filtration_value = Simplex_tree::Filtration_value;
using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;


/**
 * \class Cover_complex
 * \brief Cover complex data structure.
 *
 * \ingroup cover_complex
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
template<typename Point>
class Cover_complex {

 private:
   //Graph_induced_complex(std::map<int, double> fun){func = fun;}
   bool verbose = false; // whether to display information.
   std::vector<Point> point_cloud;
   std::vector<std::vector<int> > one_skeleton;
   typedef int Cover_t; // elements of cover C are indexed by integers.
   std::vector<std::vector<Cover_t> > simplices;
   std::map<int, std::vector<Cover_t> > cover;
   std::map<Cover_t, std::vector<int> > cover_back;
   int maximal_dim; // maximal dimension of output simplicial complex.
   int data_dimension; // dimension of input data.
   int n; // number of points.
   std::map<Cover_t,int> cover_fct; // integer-valued function that allows to state if two elements of the cover are consecutive or not.
   std::map<Cover_t,std::pair<int,double> > cover_color; // size and coloring of the vertices of the output simplicial complex.
   Simplex_tree st;

   std::map<int,std::vector<int> > adjacency_matrix;
   std::vector<std::vector<double> > distances;

   int resolution_int = -1;
   double resolution_double = -1;
   double gain = -1;
   double rate_constant = 10; // Constant in the subsampling.
   double rate_power = 0.001; // Power in the subsampling.
   int mask = 0; // Ignore nodes containing less than mask points.

   std::map<int, double> func;
   std::map<int, double> func_color;
   std::vector<int> voronoi_subsamples;
   std::string cover_name;
   std::string point_cloud_name;
   std::string color_name;
   std::string type; // Nerve or GIC
   bool functional_cover = false; // whether we use a cover with preimages of a function or not

   // Point comparator
   struct Less{
     Less(std::map<int, double> func){Fct = func;}
     std::map<int, double> Fct;
     bool operator()(int a, int b){if(Fct[a] == Fct[b])  return a < b;  else  return Fct[a] < Fct[b];}
   };

 // DFS
 private:
   void dfs(std::map<int,std::vector<int> >& G, int p, std::vector<int>& cc, std::map<int,bool>& visit){
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
   int t = 0; int m = 0; double u;
   while (m < sampleSize){
     u = GetUniform();
     if ( (populationSize - t)*u >= sampleSize - m )
       t++;
      else{samples[m] = t; t++; m++;}
   }
 }

 private:
   void fill_adjacency_matrix_from_st(){
     std::vector<int> empty;
     for(int i = 0; i < n; i++)  adjacency_matrix[i] = empty;
     for (auto simplex : st.complex_simplex_range()) {
       if(st.dimension(simplex) == 1){
         std::vector<int> vertices;
         for(auto vertex : st.simplex_vertex_range(simplex))  vertices.push_back(vertex);
         adjacency_matrix[vertices[0]].push_back(vertices[1]); adjacency_matrix[vertices[1]].push_back(vertices[0]);
       }
     }
   }

  public:
    /** \brief Specifies whether the type of the output simplicial complex.
     *
     * @param[in] t string (either "GIC" or "Nerve").
     *
     */
    void set_type(const std::string & t){type = t;}

 public:
   /** \brief Specifies whether the program should display information or not.
    *
    * @param[in] verb boolean (true = display info, false = do not display info).
    *
    */
   void set_verbose(bool verb = false){verbose = verb;}

 public:
   /** \brief Sets the constants used to subsample the data set. These constants are
    * explained in \cite Carriere17c.
    *
    * @param[in] constant double.
    * @param[in] power double.
    *
    */
   void set_subsampling(double constant, double power){rate_constant = constant; rate_power = power;}

 public:
   /** \brief Sets the mask, which is a threshold integer such that nodes in the complex that contain a number of data points which is less than or equal to
    * this threshold are not displayed.
    *
    * @param[in] nodemask integer.
    *
    */
   void set_mask(int nodemask){mask = nodemask;}

 public:
   /** \brief Reads and stores the input point cloud.
    *
    * @param[in] off_file_name name of the input .OFF or .nOFF file.
    *
    */
   bool read_point_cloud(const std::string & off_file_name){

     point_cloud_name = off_file_name; std::ifstream input(off_file_name); std::string line;

     char comment = '#';
     while(comment == '#'){getline(input, line); if(!line.empty() && !std::all_of(line.begin(),line.end(),isspace)) comment = line[line.find_first_not_of(' ')];}
     if(std::strcmp((char*) line.c_str(),"nOFF")==0){
       comment = '#'; while(comment == '#'){getline(input, line); if(!line.empty() && !std::all_of(line.begin(),line.end(),isspace)) comment = line[line.find_first_not_of(' ')];}
       std::stringstream stream(line); stream >> data_dimension;
     }
     else data_dimension = 3;

     comment = '#'; int numedges, numfaces, i, num;
     while(comment == '#'){getline(input, line); if(!line.empty() && !std::all_of(line.begin(),line.end(),isspace)) comment = line[line.find_first_not_of(' ')];}
     std::stringstream stream(line); stream >> n; stream >> numfaces; stream >> numedges;

     i = 0; while(i < n){
       getline(input, line);
       if(!line.empty() && line[line.find_first_not_of(' ')] != '#' && !std::all_of(line.begin(),line.end(),isspace)){
         std::vector<double> point; std::istringstream iss(line);
         point.assign(std::istream_iterator<double>(iss), std::istream_iterator<double>());
         point_cloud.emplace_back(point.begin(),point.begin()+data_dimension); i++;
       }
     }

     i = 0;  while(i < numfaces){
       getline(input, line);
       if(!line.empty() && line[line.find_first_not_of(' ')] != '#' && !std::all_of(line.begin(),line.end(),isspace)){
         std::vector<int> simplex; std::istringstream iss(line);
         simplex.assign(std::istream_iterator<int>(iss), std::istream_iterator<int>());
         num = simplex[0]; std::vector<int> edge(2);
         for(int j = 1; j <= num; j++){ for(int k = j+1; k <= num; k++){ edge[0] = simplex[j]; edge[1] = simplex[k]; one_skeleton.push_back(edge); } }
         i++;
       }
     }

     return input.is_open();
   }

 // *******************************************************************************************************************
 // Graphs.
 // *******************************************************************************************************************

 public: // Set graph from file.
   /** \brief Creates a graph G from a file containing the edges.
    *
    * @param[in] graph_file_name name of the input graph file.
    * The graph file contains one edge per line,
    * each edge being represented by the IDs of its two nodes.
    *
    */
   void set_graph_from_file(const std::string & graph_file_name){
     int neighb; std::ifstream input(graph_file_name);
     std::string line; int edge[2]; int n = 0;
     while(std::getline(input,line)){
       std::stringstream stream(line); stream >> edge[0];
       while(stream >> neighb){edge[1] = neighb; st.insert_simplex_and_subfaces(edge);}
       n++;
     }

     fill_adjacency_matrix_from_st();

   }


 public: // Set graph from OFF file.
   /** \brief Creates a graph G from the triangulation given by the input .OFF file.
    *
    */
   void set_graph_from_OFF(){

     int num_edges = one_skeleton.size();
     if(num_edges > 0){
       for(int i = 0; i < num_edges; i++)  st.insert_simplex_and_subfaces(one_skeleton[i]);
       fill_adjacency_matrix_from_st();
     }
     else  std::cout << "No triangulation read in OFF file!" << std::endl;

   }

 public: // Set graph from Rips complex.
   /** \brief Creates a graph G from a Rips complex.
    *
    * @param[in] threshold threshold value for the Rips complex.
    * @param[in] distance distance used to compute the Rips complex.
    *
    */
   template<typename Distance> void set_graph_from_rips(double threshold, Distance distance){

     Rips_complex rips_complex_from_points(point_cloud, threshold, distance);
     rips_complex_from_points.create_complex(st, 1);
     fill_adjacency_matrix_from_st();

   }



 public: // Pairwise distances.
   /** \private \brief Computes all pairwise distances.
    */
   template<typename Distance> void compute_pairwise_distances(Distance ref_distance){

     double d; std::vector<double> zeros(n); for(int i = 0; i < n; i++)  distances.push_back(zeros);
     std::string distance = point_cloud_name; distance.append("_dist");
     std::ifstream input(distance.c_str(), std::ios::out | std::ios::binary);

     if(input.good()){
       if(verbose)  std::cout << "Reading distances..." << std::endl;
       for(int i = 0; i < n; i++){
         for (int j = i; j < n; j++){
           input.read((char*) &d,8); distances[i][j] = d; distances[j][i] = d;
         }
       }
       input.close();
     }

     else{
       if(verbose)  std::cout << "Computing distances..." << std::endl;
       input.close(); std::ofstream output(distance, std::ios::out | std::ios::binary);
       for(int i = 0; i < n; i++){
         int state = (int) floor( 100*(i*1.0+1)/n ) %10;
         if( state == 0  && verbose)  std::cout << "\r" << state << "%" << std::flush;
         for (int j = i; j < n; j++){
           double dis = ref_distance(point_cloud[i],point_cloud[j]);
           distances[i][j] = dis; distances[j][i] = dis;
           output.write((char*) &dis,8);
         }
       }
       output.close(); if(verbose)  std::cout << std::endl;
     }

   }

 public: // Automatic tuning of Rips complex.
   /** \brief Creates a graph G from a Rips complex whose threshold value is automatically tuned with subsampling---see \cite Carriere17c.
    *
    * @param[in] distance distance between data points.
    * @param[in]  N number of subsampling iteration (the default reasonable value is 100, but there is no guarantee on how to choose it).
    * @result delta threshold used for computing the Rips complex.
    *
    */
   template<typename Distance> double set_graph_from_automatic_rips(Distance distance, int N = 100){

     int m = floor(n/   std::exp((1+rate_power)*std::log(std::log(n)/std::log(rate_constant)))   );
     m = std::min(m,n-1);
     std::vector<int> samples(m); double delta = 0;

     if(verbose)  std::cout << n << " points in R^" << data_dimension << std::endl;
     if(verbose)  std::cout << "Subsampling " << m << " points" << std::endl;

     if(distances.size() == 0)  compute_pairwise_distances(distance);

     //#pragma omp parallel for
     for (int i = 0; i < N; i++){

       SampleWithoutReplacement(n,m,samples);
       double hausdorff_dist = 0;
       for (int j = 0; j < n; j++){
         double mj = distances[j][samples[0]]; for (int k = 1; k < m; k++)  mj = std::min(mj, distances[j][samples[k]]);
         hausdorff_dist = std::max(hausdorff_dist, mj);
       }
       delta += hausdorff_dist/N;

     }

     if(verbose)  std::cout << "delta = " << delta << std::endl;
     Rips_complex rips_complex_from_points(point_cloud, delta, distance);
     rips_complex_from_points.create_complex(st, 1);
     fill_adjacency_matrix_from_st();

     return delta;

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
   void set_function_from_file(const std::string & func_file_name){
     int vertex_id = 0; std::ifstream input(func_file_name); std::string line; double f;
     while(std::getline(input,line)){
       std::stringstream stream(line); stream >> f;
       func.emplace(vertex_id, f); vertex_id++;
     }
     functional_cover = true;
     cover_name = func_file_name;
   }

 public: // Set function from kth coordinate
   /** \brief Creates the function f from the k-th coordinate of the point cloud P.
    *
    * @param[in] k coordinate to use (start at 0).
    *
    */
   void set_function_from_coordinate(int k){
     for(int i = 0; i < n; i++)  func.emplace(i,point_cloud[i][k]);
     char coordinate[100]; sprintf(coordinate, "coordinate %d", k);
     functional_cover = true;
     cover_name = coordinate;
   }

 public: // Set function from vector.
   /** \brief Creates the function f from a vector stored in memory.
    *
    * @param[in] function input vector of values.
    *
    */
   template<class InputRange> void set_function_from_range(InputRange const & function){
     functional_cover = true;
     int index = 0; for(auto v : function){func.emplace(index, v); index++;}
   }

 // *******************************************************************************************************************
 // Covers.
 // *******************************************************************************************************************

 public: // Automatic tuning of resolution.
   /** \brief Computes the optimal length of intervals
    * (i.e. the smallest interval length avoiding discretization artifacts---see \cite Carriere17c) for a functional cover.
    *
    * @result reso interval length used to compute the cover.
    *
    */
   double set_automatic_resolution(){

     if(!functional_cover){std::cout << "Cover needs to come from the preimages of a function." << std::endl; return 0;}
     if(type != "Nerve" && type != "GIC"){std::cout << "Type of complex needs to be specified." << std::endl; return 0;}

     double reso = 0;

     if(type == "GIC"){
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

     if(type == "Nerve"){
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

     return reso;

   }

 public:
   /** \brief Sets a length of intervals from a value stored in memory.
    *
    * @param[in] reso length of intervals.
    *
    */
   void set_resolution_with_interval_length(double reso){resolution_double = reso;}
   /** \brief Sets a number of intervals from a value stored in memory.
    *
    * @param[in] reso number of intervals.
    *
    */
   void set_resolution_with_interval_number(int reso){resolution_int = reso;}
   /** \brief Sets a gain from a value stored in memory (default value 0.3).
    *
    * @param[in] g gain.
    *
    */
   void set_gain(double g = 0.3){gain = g;}


 public: // Set cover with preimages of function.
   /** \brief Creates a cover C from the preimages of the function f.
    *
    */
   void set_cover_from_function(){

     if(resolution_double == -1 && resolution_int == -1){ std::cout << "Number and/or length of intervals not specified" << std::endl; return; }
     if(gain == -1){ std::cout << "Gain not specified" << std::endl; return; }

     // Read function values and compute min and max
     std::map<int, double>::iterator it;
     double maxf, minf; minf = std::numeric_limits<float>::max(); maxf = std::numeric_limits<float>::min();
     for(it = func.begin(); it != func.end(); it++){minf = std::min(minf, it->second); maxf = std::max(maxf, it->second);}
     int n = func.size(); if(verbose)  std::cout << "Min function value = " << minf << " and Max function value = " << maxf << std::endl;

     // Compute cover of im(f)
     std::vector<std::pair<double,double> > intervals; int res;

     if(resolution_double == -1){ // Case we use an integer for the number of intervals.
       double incr = (maxf-minf)/resolution_int; double x = minf; double alpha = (incr*gain)/(2-2*gain);
       double y = minf + incr + alpha; std::pair<double, double> interm(x,y); intervals.push_back(interm);
       for(int i = 1; i < resolution_int-1; i++){
         x = minf + i*incr - alpha;
         y = minf + (i+1)*incr + alpha;
         std::pair<double, double> inter(x,y); intervals.push_back(inter);
       }
       x = minf + (resolution_int-1)*incr - alpha; y = maxf;
       std::pair<double, double> interM(x,y); intervals.push_back(interM); res = intervals.size();
       if(verbose)
         for(int i = 0; i < res; i++)  std::cout << "Interval " << i << " = [" << intervals[i].first << ", " << intervals[i].second << "]" << std::endl;
     }

     else{
       if(resolution_int == -1){ // Case we use a double for the length of the intervals.
         double x = minf; double y = x + resolution_double;
         while(y <= maxf && maxf - (y-gain*resolution_double) >= resolution_double){
           std::pair<double, double> inter(x,y); intervals.push_back(inter);
           x = y - gain*resolution_double;
           y = x + resolution_double;
         }
         std::pair<double, double> interM(x,maxf); intervals.push_back(interM); res = intervals.size();
         if(verbose)
           for(int i = 0; i < res; i++)  std::cout << "Interval " << i << " = [" << intervals[i].first << ", " << intervals[i].second << "]" << std::endl;
       }

       else{ // Case we use an integer and a double for the length of the intervals.
         double x = minf; double y = x + resolution_double; int count = 0;
         while(count < resolution_int && y <= maxf && maxf - (y-gain*resolution_double) >= resolution_double){
             std::pair<double, double> inter(x,y); intervals.push_back(inter);  count++;
             x = y - gain*resolution_double;
             y = x + resolution_double;
           }
           res = intervals.size();
           if(verbose)
             for(int i = 0; i < res; i++)  std::cout << "Interval " << i << " = [" << intervals[i].first << ", " << intervals[i].second << "]" << std::endl;
       }
     }

     // Sort points according to function values
     std::vector<int> points(n); for(int i = 0; i < n; i++)  points[i] = i;
     std::sort(points.begin(),points.end(),Less(this->func));
     int id = 0; int pos = 0;

     for(int i = 0; i < res; i++){       

       // Find points in the preimage
       std::map<int,std::vector<int> > prop;
       std::pair<double, double> inter1 = intervals[i];
       int tmp = pos;

       if(i != res-1){

         if(i != 0){
           std::pair<double, double> inter3 = intervals[i-1];
           while(func[points[tmp]] < inter3.second && tmp != n){
             prop[points[tmp]] = adjacency_matrix[points[tmp]];
             tmp++;
           }
         }

         std::pair<double, double> inter2 = intervals[i+1];
         while(func[points[tmp]] < inter2.first && tmp != n){
           prop[points[tmp]] = adjacency_matrix[points[tmp]];
           tmp++;
         }

         pos = tmp;
         while(func[points[tmp]] < inter1.second && tmp != n){
           prop[points[tmp]] = adjacency_matrix[points[tmp]];
           tmp++;
         }

       }

       else{

         std::pair<double, double> inter3 = intervals[i-1];
         while(func[points[tmp]] < inter3.second && tmp != n){
           prop[points[tmp]] = adjacency_matrix[points[tmp]];
           tmp++;
         }

         while(tmp != n){
           prop[points[tmp]] = adjacency_matrix[points[tmp]];
           tmp++;
         }

       }

       // Compute the connected components with DFS
       std::map<int,bool> visit; if(verbose)  std::cout << "Preimage of interval " << i << std::endl;
       for(std::map<int, std::vector<int> >::iterator it = prop.begin(); it != prop.end(); it++)
         visit[it->first] = false;
       if (!(prop.empty())){
         for(std::map<int, std::vector<int> >::iterator it = prop.begin(); it != prop.end(); it++){
           if (  !(visit[it->first])  ){
             std::vector<int> cc; cc.clear();
             dfs(prop,it->first,cc,visit); int cci = cc.size(); if(verbose)  std::cout << "one CC with " << cci << " points, ";
             double average_col = 0;
             for(int j = 0; j < cci; j++){cover[cc[j]].push_back(id); cover_back[id].push_back(cc[j]); average_col += func_color[cc[j]]/cci;}
             cover_fct[id] = i; cover_color[id] = std::pair<int,double>(cci,average_col);
             id++;
           }
         }
       }
       if(verbose)  std::cout << std::endl;
     }

     maximal_dim = id-1;

   }

 public: // Set cover from file.
   /** \brief Creates the cover C from a file containing the cover elements of each point (the order has to be the same
   * as in the input file!).
   *
   * @param[in] cover_file_name name of the input cover file.
   *
   */
   void set_cover_from_file(const std::string & cover_file_name){
     int vertex_id = 0; Cover_t cov; std::vector<Cover_t> cov_elts, cov_number;
     std::ifstream input(cover_file_name); std::string line;
     while(std::getline(input,line)){
       cov_elts.clear(); std::stringstream stream(line);
       while(stream >> cov){
         cov_elts.push_back(cov); cov_number.push_back(cov);
         cover_fct[cov] = cov; cover_color[cov].second += func_color[vertex_id]; cover_color[cov].first++; cover_back[cov].push_back(vertex_id);
       }
       cover[vertex_id] = cov_elts; vertex_id++;
     }
     std::vector<Cover_t>::iterator it;
     std::sort(cov_number.begin(),cov_number.end()); it = std::unique(cov_number.begin(),cov_number.end());
     cov_number.resize(std::distance(cov_number.begin(),it)); maximal_dim = cov_number.size()-1;
     for(int i = 0; i <= maximal_dim; i++)  cover_color[i].second /= cover_color[i].first;
     cover_name = cover_file_name;
   }

 public: // Set cover from Voronoi
   /** \brief Creates the cover C from the Voronoï cells of a subsampling of the point cloud.
   *
   * @param[in] distance distance between the points.
   * @param[in] m number of points in the subsample.
   *
   */
   template<typename Distance> void set_cover_from_Voronoi(Distance distance, int m = 100){

     voronoi_subsamples.resize(m); SampleWithoutReplacement(n,m,voronoi_subsamples);
     if(distances.size() == 0)  compute_pairwise_distances(distance);
     std::vector<double> mindist(n); for(int j = 0; j < n; j++)  mindist[j] = std::numeric_limits<double>::max();

     // Compute the geodesic distances to subsamples with Dijkstra
     for(int i = 0; i < m; i++){
       if(verbose)  std::cout << "Computing geodesic distances to seed " << i << "..." << std::endl;
       int seed = voronoi_subsamples[i];
       std::vector<double> dist(n); std::vector<int> process(n);
       for(int j = 0; j < n; j++){ dist[j] = std::numeric_limits<double>::max(); process[j] = j; }
       dist[seed] = 0; int curr_size = process.size(); int min_point, min_index; double min_dist;
       std::vector<int> neighbors; int num_neighbors;

       while(curr_size > 0){
         min_dist = std::numeric_limits<double>::max(); min_index = -1; min_point = -1;
         for(int j = 0; j < curr_size; j++){
           if(dist[process[j]] < min_dist){
             min_point = process[j]; min_dist = dist[process[j]]; min_index = j;
           }
         }
         assert(min_index != -1); process.erase(process.begin() + min_index);
         assert(min_point != -1); neighbors = adjacency_matrix[min_point]; num_neighbors = neighbors.size();
         for(int j = 0; j < num_neighbors; j++){
           double d = dist[min_point] + distances[min_point][neighbors[j]];
           dist[neighbors[j]] = std::min(dist[neighbors[j]], d);
         }
         curr_size = process.size();
       }

       for(int j = 0; j < n; j++)
         if(mindist[j] > dist[j]){
           mindist[j] = dist[j];
           if(cover[j].size() == 0)  cover[j].push_back(i);
           else cover[j][0] = i;
         }
     }

     for(int i = 0; i < n; i++){ cover_back[cover[i][0]].push_back(i); cover_color[cover[i][0]].second += func_color[i]; cover_color[cover[i][0]].first++; }
     for(int i = 0; i < m; i++)  cover_color[i].second /= cover_color[i].first;
     maximal_dim = m-1; cover_name = "Voronoi";

   }

 public: // return subset of data corresponding to a node
   /** \brief Returns the data subset corresponding to a specific node of the created complex.
   *
   * @param[in] c ID of the node.
   * @result cover_back(c) vector of IDs of data points.
   *
   */
   const std::vector<int> & subpopulation(Cover_t c){
     return cover_back[c];
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
   void set_color_from_file(const std::string & color_file_name){
     int vertex_id = 0; std::ifstream input(color_file_name); std::string line; double f;
     while(std::getline(input,line)){
       std::stringstream stream(line); stream >> f;
       func_color.emplace(vertex_id, f); vertex_id++;
     }
     color_name = color_file_name;
   }

 public: // Set color from kth coordinate
   /** \brief Computes the function used to color the nodes of the simplicial complex from the k-th coordinate.
    *
    * @param[in] k coordinate to use (start at 0).
    *
    */
   void set_color_from_coordinate(int k = 0){
     for(int i = 0; i < n; i++)  func_color.emplace(i, point_cloud[i][k]);
     color_name = "coordinate "; color_name.append(std::to_string(k));
   }

 public: // Set color from vector.
   /** \brief Computes the function used to color the nodes of the simplicial complex from a vector stored in memory.
    *
    * @param[in] color input vector of values.
    *
    */
   void set_color_from_vector(std::vector<double> color){
     for(unsigned int i = 0; i < color.size(); i++)  func_color.emplace(i, color[i]);
   }

 public: // Create a .dot file that can be compiled with neato to produce a .pdf file.
   /** \brief Creates a .dot file called SC.dot for neato (part of the graphviz package) once the simplicial complex is computed to get a visualization
    * of its 1-skeleton in a .pdf file.
    */
   void plot_DOT(){
     char mapp[11] = "SC.dot";   std::ofstream graphic(mapp); graphic << "graph GIC {" << std::endl;
     double maxv, minv; maxv = std::numeric_limits<double>::min(); minv = std::numeric_limits<double>::max();
     for (std::map<Cover_t,std::pair<int,double> >::iterator iit = cover_color.begin(); iit != cover_color.end(); iit++){
       maxv = std::max(maxv, iit->second.second);  minv = std::min(minv, iit->second.second);
     }
     int k = 0; std::vector<int> nodes; nodes.clear();
     for (std::map<Cover_t,std::pair<int,double> >::iterator iit = cover_color.begin(); iit != cover_color.end(); iit++){
       if(iit->second.first > mask){
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
         if (cover_color[simplices[i][0]].first > mask && cover_color[simplices[i][1]].first > mask){
           graphic << "  " << simplices[i][0] << " -- " << simplices[i][1] << " [weight=15];" << std::endl; ke++;}
     graphic << "}"; graphic.close();
     std::cout << "SC.dot generated. It can be visualized with e.g. neato." << std::endl;
   }

 public: // Create a .txt file that can be compiled with KeplerMapper.
   /** \brief Creates a .txt file called SC.txt describing the 1-skeleton, which can then be plotted with e.g. KeplerMapper.
    */
   void write_info(){

     int num_simplices = simplices.size(); int num_edges = 0;
     char mapp[11] = "SC.txt";  std::ofstream graphic(mapp);
     for (int i = 0; i < num_simplices; i++)
       if (simplices[i].size() == 2)
         if (cover_color[simplices[i][0]].first > mask && cover_color[simplices[i][1]].first > mask)
           num_edges++;

     graphic << point_cloud_name << std::endl;
     graphic << cover_name << std::endl;
     graphic << color_name << std::endl;
     graphic << resolution_double << " " << gain << std::endl;
     graphic << cover_color.size() << " " << num_edges << std::endl;

     for (std::map<Cover_t,std::pair<int,double> >::iterator iit = cover_color.begin(); iit != cover_color.end(); iit++)
       graphic << iit->first << " " << iit->second.second << " " << iit->second.first << std::endl;

     for (int i = 0; i < num_simplices; i++)
       if (simplices[i].size() == 2)
         if (cover_color[simplices[i][0]].first > mask && cover_color[simplices[i][1]].first > mask)
           graphic << simplices[i][0] << " " << simplices[i][1] << std::endl;
     graphic.close();
     std::cout << "SC.txt generated. It can be visualized with e.g. python KeplerMapperVisuFromTxtFile.py and firefox." << std::endl;
   }


 public: // Create a .off file that can be visualized (e.g. with Geomview).
   /** \brief Creates a .off file called SC.off for 3D visualization, which contains the 2-skeleton of the GIC.
    * This function assumes that the cover has been computed with Voronoi. If data points are in 1D or 2D,
    * the remaining coordinates of the points embedded in 3D are set to 0.
    */
   void plot_OFF(){

     assert(cover_name == "Voronoi");
     char gic[11] = "SC.off";  std::ofstream graphic(gic);
     graphic << "OFF" << std::endl; int m = voronoi_subsamples.size(); int numedges = 0; int numfaces = 0;
     std::vector<std::vector<int> > edges, faces; int numsimplices = simplices.size();
     for (int i = 0; i < numsimplices; i++) {
       if(simplices[i].size() == 2){ numedges++; edges.push_back(simplices[i]); }
       if(simplices[i].size() == 3){ numfaces++; faces.push_back(simplices[i]); }
     }
     graphic << m << " " << numedges + numfaces << std::endl;
     for(int i = 0; i < m; i++){
       if(data_dimension <= 3){
         for(int j = 0; j < data_dimension; j++)  graphic << point_cloud[voronoi_subsamples[i]][j] << " ";
         for(int j = data_dimension; j < 3; j++)  graphic << 0 << " ";
         graphic << std::endl;
       }
       else  for(int j = 0; j < 3; j++)  graphic << point_cloud[voronoi_subsamples[i]][j] << " ";
     }
     for(int i = 0; i < numedges; i++)  graphic << 2 << " " << edges[i][0] << " " << edges[i][1] << std::endl;
     for(int i = 0; i < numfaces; i++)  graphic << 3 << " " << faces[i][0] << " " << faces[i][1] << " " << faces[i][2] << std::endl;
     graphic.close();
     std::cout << "SC.off generated. It can be visualized with e.g. geomview." << std::endl;

   }

 // *******************************************************************************************************************
 // *******************************************************************************************************************


 public:
   /** \brief Creates the simplicial complex.
    *
    * @param[in] complex SimplicialComplexForRips to be created.
    *
    */
   template<typename SimplicialComplexForRips>
   void create_complex(SimplicialComplexForRips & complex) {
     unsigned int dimension = 0;
     for(auto const& simplex : simplices){
       complex.insert_simplex_and_subfaces(simplex);
       if(dimension < simplex.size()-1)  dimension = simplex.size()-1;
     }
     complex.set_dimension(dimension);
   }

 public:
   /** \brief Computes the simplices of the simplicial complex.
    */
   void find_simplices(){

     if(type != "Nerve" && type != "GIC"){std::cout << "Type of complex needs to be specified." << std::endl; return;}

     if(type == "Nerve"){
       for(std::map<int,std::vector<Cover_t> >::iterator it = cover.begin(); it!=cover.end(); it++)  simplices.push_back(it->second);
       std::vector<std::vector<Cover_t> >::iterator it;
       std::sort(simplices.begin(),simplices.end()); it = std::unique(simplices.begin(),simplices.end());
       simplices.resize(std::distance(simplices.begin(),it));
     }

     if(type == "GIC"){

       if(functional_cover){

         // Computes the simplices in the GIC by looking at all the edges of the graph and adding the
         // corresponding edges in the GIC if the images of the endpoints belong to consecutive intervals.

         if (gain >= 0.5)
           throw std::invalid_argument("the output of this function is correct ONLY if the cover is minimal, i.e. the gain is less than 0.5.");

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

       else{

         // Find IDs of edges to remove
         std::vector<int> simplex_to_remove; int simplex_id = 0;
         for (auto simplex : st.complex_simplex_range()) {
           if(st.dimension(simplex) == 1){
             std::vector<std::vector<Cover_t> > comp;
             for(auto vertex : st.simplex_vertex_range(simplex))  comp.push_back(cover[vertex]);
             if(comp[0].size() == 1 && comp[0] == comp[1])  simplex_to_remove.push_back(simplex_id);
           }
           simplex_id++;
         }

         // Remove edges
         if(simplex_to_remove.size() > 1){
           int current_id = 1;
           auto simplex = st.complex_simplex_range().begin(); int num_rem = 0;
           for(int i = 0; i < simplex_id-1; i++){
             int j = i+1; auto simplex_tmp = simplex; simplex_tmp++;
             if(j == simplex_to_remove[current_id]){st.remove_maximal_simplex(*simplex_tmp); current_id++; num_rem++;}
             else  simplex++;
           } simplex = st.complex_simplex_range().begin();
           for(int i = 0; i < simplex_to_remove[0]; i++)  simplex++;  st.remove_maximal_simplex(*simplex);
         }

         // Build the Simplex Tree corresponding to the graph
         st.expansion(maximal_dim);

         // Find simplices of GIC
         simplices.clear();
         for (auto simplex : st.complex_simplex_range()) {
           if(!st.has_children(simplex)){

             std::vector<Cover_t> simplx;
             for (auto vertex : st.simplex_vertex_range(simplex)){
               unsigned int sz = cover[vertex].size();
               for(unsigned int i = 0; i < sz; i++){
                 simplx.push_back(cover[vertex][i]);
               }
             }

             std::sort(simplx.begin(),simplx.end()); std::vector<Cover_t>::iterator it = std::unique(simplx.begin(),simplx.end());
             simplx.resize(std::distance(simplx.begin(),it));
             simplices.push_back(simplx);

           }
         }
         std::vector<std::vector<Cover_t> >::iterator it;
         std::sort(simplices.begin(),simplices.end()); it = std::unique(simplices.begin(),simplices.end());
         simplices.resize(std::distance(simplices.begin(),it));

       }
     }

   }

};

} // namespace cover_complex

} // namespace Gudhi

#endif  // GIC_H_

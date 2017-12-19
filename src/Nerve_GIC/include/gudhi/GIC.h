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
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Bottleneck.h>

#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/graph/graph_utility.hpp>

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <limits>     // for numeric_limits
#include <utility>    // for pair<>
#include <algorithm>  // for std::max
#include <random>
#include <cassert>

using namespace boost;
using namespace std;

namespace Gudhi {

namespace cover_complex {

using Simplex_tree        = Gudhi::Simplex_tree<>;
using Filtration_value    = Simplex_tree::Filtration_value;
using Rips_complex        = Gudhi::rips_complex::Rips_complex<Filtration_value>;
using PersistenceDiagram  = vector<pair<double,double> >;
using Graph               = subgraph<adjacency_list<setS, vecS, undirectedS, no_property, property<edge_index_t, int, property<edge_weight_t, double> > > >;
using vertex_t            = graph_traits<Graph>::vertex_descriptor;
using IndexMap            = property_map<Graph, vertex_index_t>::type;
using WeightMap           = property_map<Graph, edge_weight_t>::type;

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
template <typename Point>
class Cover_complex {
 private:

  bool verbose = false;  // whether to display information.

  vector<Point> point_cloud;
  int maximal_dim;     // maximal dimension of output simplicial complex.
  int data_dimension;  // dimension of input data.
  int n;               // number of points.

  vector<vector<double> > distances;

  map<int, double> func;          // function used to compute the output simplicial complex.
  map<int, double> func_color;    // function used to compute the colors of the nodes of the output simplicial complex.
  bool functional_cover = false;  // whether we use a cover with preimages of a function or not.

  Graph one_skeleton_OFF;   // one-skeleton given by the input OFF file (if it exists).
  Graph one_skeleton;       // one-skeleton used to compute the connected components.
  vector<vertex_t> vertices;
  vector<vector<int> > simplices;

  vector<int> voronoi_subsamples;

  PersistenceDiagram PD;
  vector<double> distribution;

  map<int, vector<int> > cover;
  map<int, vector<int> > cover_back;
  map<int, double> cover_std;                // standard function (induced by func) used to compute the extended persistence diagram of the output simplicial complex.
  map<int, int> cover_fct;                   // integer-valued function that allows to state if two elements of the cover are consecutive or not.
  map<int, pair<int, double> > cover_color;  // size and coloring (induced by func_color) of the vertices of the output simplicial complex.

  int resolution_int = -1;
  double resolution_double = -1;
  double gain = -1;
  double rate_constant = 10;  // Constant in the subsampling.
  double rate_power = 0.001;  // Power in the subsampling.
  int mask = 0;               // Ignore nodes containing less than mask points.

  string cover_name;
  string point_cloud_name;
  string color_name;
  string type;               // Nerve or GIC

  // Point comparator
  struct Less {
    Less(map<int, double> func) { Fct = func; }
    map<int, double> Fct;
    bool operator()(int a, int b) {
      if (Fct[a] == Fct[b])
        return a < b;
      else
        return Fct[a] < Fct[b];
    }
  };

  // Find random number in [0,1].
  double GetUniform() {
    thread_local default_random_engine re;
    thread_local uniform_real_distribution<double> Dist(0, 1);
    return Dist(re);
  }

  // Subsample points.
  void SampleWithoutReplacement(int populationSize, int sampleSize, vector<int> & samples) {
    int t = 0; int m = 0; double u;
    while (m < sampleSize){
      u = GetUniform();
      if ((populationSize - t) * u >= sampleSize - m)  t++;
      else{ samples[m] = t; t++; m++;}
    }
  }


 public:
  /** \brief Specifies whether the type of the output simplicial complex.
   *
   * @param[in] t string (either "GIC" or "Nerve").
   *
   */
  void set_type(const std::string& t) { type = t; }

 public:
  /** \brief Specifies whether the program should display information or not.
   *
   * @param[in] verb boolean (true = display info, false = do not display info).
   *
   */
  void set_verbose(bool verb = false) { verbose = verb; }

 public:
  /** \brief Sets the constants used to subsample the data set. These constants are
   * explained in \cite Carriere17c.
   *
   * @param[in] constant double.
   * @param[in] power double.
   *
   */
  void set_subsampling(double constant, double power) {
    rate_constant = constant;
    rate_power = power;
  }

 public:
  /** \brief Sets the mask, which is a threshold integer such that nodes in the complex that contain a number of data
   * points which is less than or equal to
   * this threshold are not displayed.
   *
   * @param[in] nodemask integer.
   *
   */
  void set_mask(int nodemask) { mask = nodemask; }

 public:
  /** \brief Reads and stores the input point cloud.
   *
   * @param[in] off_file_name name of the input .OFF or .nOFF file.
   *
   */
  bool read_point_cloud(const string & off_file_name) {
    point_cloud_name = off_file_name;
    ifstream input(off_file_name);
    string line;

    char comment = '#';
    while (comment == '#') {
      getline(input, line);
      if (!line.empty() && !all_of(line.begin(), line.end(), (int(*)(int))isspace)) comment = line[line.find_first_not_of(' ')];
    }
    if (strcmp((char*)line.c_str(), "nOFF") == 0) {
      comment = '#';
      while (comment == '#') {
        getline(input, line);
        if (!line.empty() && !all_of(line.begin(), line.end(), (int(*)(int))isspace))
          comment = line[line.find_first_not_of(' ')];
      }
      stringstream stream(line);
      stream >> data_dimension;
    } else {
      data_dimension = 3;
    }

    comment = '#';
    int numedges, numfaces, i, dim;
    while (comment == '#') {
      getline(input, line);
      if (!line.empty() && !all_of(line.begin(), line.end(), (int(*)(int))isspace)) comment = line[line.find_first_not_of(' ')];
    }
    stringstream stream(line);
    stream >> n;
    stream >> numfaces;
    stream >> numedges;

    i = 0;
    while (i < n) {
      getline(input, line);
      if (!line.empty() && line[line.find_first_not_of(' ')] != '#' && !all_of(line.begin(), line.end(), (int(*)(int))isspace)) {
        istringstream iss(line); vector<double> point; point.assign(istream_iterator<double>(iss), istream_iterator<double>());
        point_cloud.emplace_back(point.begin(), point.begin() + data_dimension);
        add_vertex(one_skeleton_OFF); vertices.push_back(add_vertex(one_skeleton));
        i++;
      }
    }

    i = 0;
    while (i < numfaces) {
      getline(input, line);
      if (!line.empty() && line[line.find_first_not_of(' ')] != '#' && !all_of(line.begin(), line.end(), (int(*)(int))isspace)) {
        vector<int> simplex; istringstream iss(line);
        simplex.assign(istream_iterator<int>(iss), istream_iterator<int>()); dim = simplex[0];
        for (int j = 1; j <= dim; j++)
          for (int k = j + 1; k <= dim; k++)
            add_edge(vertices[simplex[j]], vertices[simplex[k]], one_skeleton_OFF);
        i++;
      }
    }

    return input.is_open();
  }

  // *******************************************************************************************************************
  // Graphs.
  // *******************************************************************************************************************

 public:  // Set graph from file.
          /** \brief Creates a graph G from a file containing the edges.
           *
           * @param[in] graph_file_name name of the input graph file.
           * The graph file contains one edge per line,
           * each edge being represented by the IDs of its two nodes.
           *
           */
  void set_graph_from_file(const string & graph_file_name){
    int neighb; ifstream input(graph_file_name); string line; int source;
    while (getline(input, line)){
      stringstream stream(line); stream >> source;
      while (stream >> neighb)  add_edge(vertices[source], vertices[neighb], one_skeleton);
    }
  }

 public:  // Set graph from OFF file.
          /** \brief Creates a graph G from the triangulation given by the input .OFF file.
           *
           */
  void set_graph_from_OFF() {
    if(num_edges(one_skeleton_OFF))  one_skeleton = one_skeleton_OFF;
    else  cout << "No triangulation read in OFF file!" << endl;
  }

 public:  // Set graph from Rips complex.
          /** \brief Creates a graph G from a Rips complex.
           *
           * @param[in] threshold threshold value for the Rips complex.
           * @param[in] distance distance used to compute the Rips complex.
           *
           */
  template <typename Distance>
  void set_graph_from_rips(double threshold, Distance distance) {
    if(distances.size() == 0)  compute_pairwise_distances(distance);
    for(int i = 0; i < n; i++){
      for(int j = i+1; j < n; j++){
        if(distances[i][j] <= threshold){
          add_edge(vertices[i], vertices[j], one_skeleton);
          put(edge_weight, one_skeleton, edge(vertices[i], vertices[j], one_skeleton).first, distances[i][j]);
        }
      }
    }
  }

 public:
  void set_graph_weights(){
    IndexMap index = get(vertex_index, one_skeleton); WeightMap weight = get(edge_weight, one_skeleton);
    graph_traits<Graph>::edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(one_skeleton); ei != ei_end; ++ei)
      put(weight, *ei, distances[index[source(*ei, one_skeleton)]][index[target(*ei, one_skeleton)]]);
  }

 public:  // Pairwise distances.
          /** \private \brief Computes all pairwise distances.
           */
  template <typename Distance>
  void compute_pairwise_distances(Distance ref_distance) {
    double d; vector<double> zeros(n); for (int i = 0; i < n; i++) distances.push_back(zeros);
    string distance = point_cloud_name;
    distance.append("_dist");
    ifstream input(distance.c_str(), ios::out | ios::binary);

    if (input.good()) {
      if (verbose) std::cout << "Reading distances..." << std::endl;
      for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
          input.read((char*)&d, 8);
          distances[i][j] = d; distances[j][i] = d;
        }
      }
      input.close();
    } else {
      if (verbose) cout << "Computing distances..." << endl;
      input.close(); ofstream output(distance, ios::out | ios::binary);
      for (int i = 0; i < n; i++) {
        int state = (int)floor(100 * (i * 1.0 + 1) / n) % 10;
        if (state == 0 && verbose) cout << "\r" << state << "%" << flush;
        for (int j = i; j < n; j++) {
          double dis = ref_distance(point_cloud[i], point_cloud[j]);
          distances[i][j] = dis; distances[j][i] = dis;
          output.write((char*)&dis, 8);
        }
      }
      output.close();
      if (verbose) cout << endl;
    }
  }

 public:  // Automatic tuning of Rips complex.
  /** \brief Creates a graph G from a Rips complex whose threshold value is automatically tuned with subsampling---see
   * \cite Carriere17c.
   *
   * @param[in] distance distance between data points.
   * @param[in]  N number of subsampling iteration (the default reasonable value is 100, but there is no guarantee on
   * how to choose it).
   * @result delta threshold used for computing the Rips complex.
   *
   */
  template <typename Distance>
  double set_graph_from_automatic_rips(Distance distance, int N = 100) {
    int m = floor(n / exp((1 + rate_power) * log(log(n) / log(rate_constant))));
    m = min(m, n - 1);
    vector<int> samples(m);
    double delta = 0;

    if (verbose) cout << n << " points in R^" << data_dimension << endl;
    if (verbose) cout << "Subsampling " << m << " points" << endl;

    if (distances.size() == 0) compute_pairwise_distances(distance);

    // #pragma omp parallel for
    for (int i = 0; i < N; i++) {
      SampleWithoutReplacement(n, m, samples);
      double hausdorff_dist = 0;
      for (int j = 0; j < n; j++) {
        double mj = distances[j][samples[0]];
        for (int k = 1; k < m; k++) mj = min(mj, distances[j][samples[k]]);
        hausdorff_dist = max(hausdorff_dist, mj);
      }
      delta += hausdorff_dist / N;
    }

    if (verbose) cout << "delta = " << delta << endl;
    set_graph_from_rips(delta, distance);
    return delta;
  }

  // *******************************************************************************************************************
  // Functions.
  // *******************************************************************************************************************

 public:  // Set function from file.
          /** \brief Creates the function f from a file containing the function values.
           *
           * @param[in] func_file_name name of the input function file.
           *
           */
  void set_function_from_file(const std::string& func_file_name) {
    int i = 0; ifstream input(func_file_name); string line; double f;
    while (getline(input, line)) {
      stringstream stream(line); stream >> f; func.emplace(i, f); i++;
    }
    functional_cover = true;
    cover_name = func_file_name;
  }

 public:  // Set function from kth coordinate
          /** \brief Creates the function f from the k-th coordinate of the point cloud P.
           *
           * @param[in] k coordinate to use (start at 0).
           *
           */
  void set_function_from_coordinate(int k) {
    for (int i = 0; i < n; i++) func.emplace(i, point_cloud[i][k]);
    char coordinate[100];
    sprintf(coordinate, "coordinate %d", k);
    functional_cover = true;
    cover_name = coordinate;
  }

 public:  // Set function from vector.
          /** \brief Creates the function f from a vector stored in memory.
           *
           * @param[in] function input vector of values.
           *
           */
  template <class InputRange>
  void set_function_from_range(InputRange const& f) {
    for (int i = 0; i < n; i++) func.emplace(i, f[i]);
    functional_cover = true;
  }

  // *******************************************************************************************************************
  // Covers.
  // *******************************************************************************************************************

 public:  // Automatic tuning of resolution.
          /** \brief Computes the optimal length of intervals
           * (i.e. the smallest interval length avoiding discretization artifacts---see \cite Carriere17c) for a functional
           * cover.
           *
           * @result reso interval length used to compute the cover.
           *
           */
  double set_automatic_resolution() {
    if (!functional_cover) {
      std::cout << "Cover needs to come from the preimages of a function." << std::endl;
      return 0;
    }
    if (type != "Nerve" && type != "GIC") {
      std::cout << "Type of complex needs to be specified." << std::endl;
      return 0;
    }

    double reso = 0; IndexMap index = get(vertex_index, one_skeleton);

    if (type == "GIC") {
      graph_traits<Graph>::edge_iterator ei, ei_end;
      for (tie(ei, ei_end) = edges(one_skeleton); ei != ei_end; ++ei)
        reso = max(reso, abs(func[index[source(*ei, one_skeleton)]] - func[index[target(*ei, one_skeleton)]]));
      if (verbose) cout << "resolution = " << reso << endl;
      resolution_double = reso;
    }

    if (type == "Nerve") {
      graph_traits<Graph>::edge_iterator ei, ei_end;
      for (tie(ei, ei_end) = edges(one_skeleton); ei != ei_end; ++ei)
        reso = max(reso, abs(func[index[source(*ei, one_skeleton)]] - func[index[target(*ei, one_skeleton)]]) / gain);
      if (verbose) cout << "resolution = " << reso << endl;
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
  void set_resolution_with_interval_length(double reso) { resolution_double = reso; }
  /** \brief Sets a number of intervals from a value stored in memory.
   *
   * @param[in] reso number of intervals.
   *
   */
  void set_resolution_with_interval_number(int reso) { resolution_int = reso; }
  /** \brief Sets a gain from a value stored in memory (default value 0.3).
   *
   * @param[in] g gain.
   *
   */
  void set_gain(double g = 0.3) { gain = g; }

 public:  // Set cover with preimages of function.
          /** \brief Creates a cover C from the preimages of the function f.
           *
           */
  void set_cover_from_function() {
    if (resolution_double == -1 && resolution_int == -1) {
      cout << "Number and/or length of intervals not specified" << endl;
      return;
    }
    if (gain == -1) {
      cout << "Gain not specified" << endl;
      return;
    }

    // Read function values and compute min and max
    double minf = numeric_limits<float>::max(); double maxf = numeric_limits<float>::lowest();
    for (int i = 0; i < n; i++) {
      minf = min(minf, func[i]); maxf = max(maxf, func[i]);
    }
    if (verbose) cout << "Min function value = " << minf << " and Max function value = " << maxf << endl;

    // Compute cover of im(f)
    vector<pair<double, double> > intervals; int res;

    if (resolution_double == -1) {  // Case we use an integer for the number of intervals.
      double incr = (maxf - minf) / resolution_int;
      double x = minf;
      double alpha = (incr * gain) / (2 - 2 * gain);
      double y = minf + incr + alpha;
      pair<double, double> interm(x, y);
      intervals.push_back(interm);
      for (int i = 1; i < resolution_int - 1; i++) {
        x = minf + i * incr - alpha;
        y = minf + (i + 1) * incr + alpha;
        pair<double, double> inter(x, y);
        intervals.push_back(inter);
      }
      x = minf + (resolution_int - 1) * incr - alpha;
      y = maxf;
      pair<double, double> interM(x, y);
      intervals.push_back(interM);
      res = intervals.size();
      if (verbose) {
        for (int i = 0; i < res; i++)
          cout << "Interval " << i << " = [" << intervals[i].first << ", " << intervals[i].second << "]" << endl;
      }
    } else {
      if (resolution_int == -1) {  // Case we use a double for the length of the intervals.
        double x = minf;
        double y = x + resolution_double;
        while (y <= maxf && maxf - (y - gain * resolution_double) >= resolution_double) {
          pair<double, double> inter(x, y);
          intervals.push_back(inter);
          x = y - gain * resolution_double;
          y = x + resolution_double;
        }
        pair<double, double> interM(x, maxf);
        intervals.push_back(interM);
        res = intervals.size();
        if (verbose) {
          for (int i = 0; i < res; i++)
            cout << "Interval " << i << " = [" << intervals[i].first << ", " << intervals[i].second << "]" << endl;
        }
      } else {  // Case we use an integer and a double for the length of the intervals.
        double x = minf;
        double y = x + resolution_double;
        int count = 0;
        while (count < resolution_int && y <= maxf && maxf - (y - gain * resolution_double) >= resolution_double) {
          pair<double, double> inter(x, y);
          intervals.push_back(inter);
          count++;
          x = y - gain * resolution_double;
          y = x + resolution_double;
        }
        res = intervals.size();
        if (verbose) {
          for (int i = 0; i < res; i++)
            cout << "Interval " << i << " = [" << intervals[i].first << ", " << intervals[i].second << "]" << endl;
        }
      }
    }

    // Sort points according to function values
    vector<int> points(n); for (int i = 0; i < n; i++) points[i] = i;
    sort(points.begin(), points.end(), Less(this->func));

    int id = 0; int pos = 0; int maxc = -1; IndexMap index = get(vertex_index, one_skeleton);

    for (int i = 0; i < res; i++) {

      // Find points in the preimage
      vector<int> indices; pair<double, double> inter1 = intervals[i];
      int tmp = pos; double u, v; Graph G = one_skeleton.create_subgraph();

      if (i != res - 1) {

        if (i != 0) {
          pair<double, double> inter3 = intervals[i - 1];
          while (func[points[tmp]] < inter3.second && tmp != n){
            add_vertex(index[vertices[points[tmp]]], G); indices.push_back(points[tmp]); tmp++;
          }
          u = inter3.second;
        }
        else  u = inter1.first;

        pair<double, double> inter2 = intervals[i + 1];
        while (func[points[tmp]] < inter2.first && tmp != n){
          add_vertex(index[vertices[points[tmp]]], G); indices.push_back(points[tmp]); tmp++;
        }

        v = inter2.first;

        pos = tmp;
        while (func[points[tmp]] < inter1.second && tmp != n){
          add_vertex(index[vertices[points[tmp]]], G); indices.push_back(points[tmp]); tmp++;
        }

      } else {
        pair<double, double> inter3 = intervals[i - 1];
        while (func[points[tmp]] < inter3.second && tmp != n){
          add_vertex(index[vertices[points[tmp]]], G); indices.push_back(points[tmp]); tmp++;
        }

        while (tmp != n){
          add_vertex(index[vertices[points[tmp]]], G); indices.push_back(points[tmp]); tmp++;
        }

        u = inter3.second; v = inter1.second;

      }

      int num = num_vertices(G); vector<int> component(num);

      // Compute connected components
      connected_components(G, &component[0]); int maxct = maxc + 1;

      // Update covers
      for(int j = 0; j < num; j++){
        maxc = max(maxc, maxct + component[j]);
        cover       [indices[j]]               .push_back(maxct + component[j]);
        cover_back  [maxct + component[j]]     .push_back(indices[j]);
        cover_fct   [maxct + component[j]]     = i;
        cover_std   [maxct + component[j]]     = 0.5*(u+v);
        cover_color [maxct + component[j]]     .second += func_color[indices[j]]; //= pair<int, double>(cci, average_col);
        cover_color [maxct + component[j]]     .first += 1;
      }
    }

    maximal_dim = id - 1;
    for (map<int, pair<int, double> >::iterator iit = cover_color.begin(); iit != cover_color.end(); iit++)
      iit->second.second /= iit->second.first;
  }

 public:  // Set cover from file.
  /** \brief Creates the cover C from a file containing the cover elements of each point (the order has to be the same
  * as in the input file!).
  *
  * @param[in] cover_file_name name of the input cover file.
  *
  */
  void set_cover_from_file(const string & cover_file_name) {
    int i = 0; int cov; vector<int> cov_elts, cov_number;
    ifstream input(cover_file_name); string line;
    while (getline(input, line)) {
      cov_elts.clear();
      stringstream stream(line);
      while (stream >> cov) {
        cov_elts.push_back(cov);
        cov_number.push_back(cov);
        cover_fct   [cov]  = cov;
        cover_color [cov]  .second += func_color[i];
        cover_color [cov]  .first++;
        cover_back  [cov]  .push_back(i);
      }
      cover[i] = cov_elts; i++;
    }

    sort(cov_number.begin(), cov_number.end());
    vector<int>::iterator it = unique(cov_number.begin(), cov_number.end());
    cov_number.resize(distance(cov_number.begin(), it));

    maximal_dim = cov_number.size() - 1;
    for (int i = 0; i <= maximal_dim; i++) cover_color[i].second /= cover_color[i].first;
    cover_name = cover_file_name;
  }

 public:  // Set cover from Voronoi
          /** \brief Creates the cover C from the Voronoï cells of a subsampling of the point cloud.
          *
          * @param[in] distance distance between the points.
          * @param[in] m number of points in the subsample.
          *
          */
  template <typename Distance>
  void set_cover_from_Voronoi(Distance distance, int m = 100){

    voronoi_subsamples.resize(m); SampleWithoutReplacement(n, m, voronoi_subsamples);
    if (distances.size() == 0) compute_pairwise_distances(distance); set_graph_weights();
    WeightMap weight = get(edge_weight, one_skeleton); IndexMap index = get(vertex_index, one_skeleton);
    vector<double> mindist(n); for (int j = 0; j < n; j++) mindist[j] = numeric_limits<double>::max();

    // Compute the geodesic distances to subsamples with Dijkstra
    for (int i = 0; i < m; i++) {

      if (verbose) cout << "Computing geodesic distances to seed " << i << "..." << endl;
      int seed = voronoi_subsamples[i]; vector<double> dmap(n);
      dijkstra_shortest_paths(one_skeleton, vertices[seed], weight_map(weight).distance_map(make_iterator_property_map(dmap.begin(), index)));

      for (int j = 0; j < n; j++)
        if (mindist[j] > dmap[j]) {
          mindist[j] = dmap[j];
          if (cover[j].size() == 0)  cover[j].push_back(i);  else  cover[j][0] = i;
        }
    }

    for (int i = 0; i < n; i++) {
      cover_back[cover[i][0]].push_back(i);
      cover_color[cover[i][0]].second += func_color[i];
      cover_color[cover[i][0]].first++;
    }
    for (int i = 0; i < m; i++) cover_color[i].second /= cover_color[i].first;
    maximal_dim = m - 1;
    cover_name = "Voronoi";

  }

 public:  // return subset of data corresponding to a node
          /** \brief Returns the data subset corresponding to a specific node of the created complex.
          *
          * @param[in] c ID of the node.
          * @result cover_back(c) vector of IDs of data points.
          *
          */
  const vector<int> & subpopulation(int c) { return cover_back[c]; }

  // *******************************************************************************************************************
  // Visualization.
  // *******************************************************************************************************************

 public:  // Set color from file.
  /** \brief Computes the function used to color the nodes of the simplicial complex from a file containing the function
   * values.
   *
   * @param[in] color_file_name name of the input color file.
   *
   */
  void set_color_from_file(const string & color_file_name) {
    int i = 0;
    ifstream input(color_file_name);
    string line;
    double f;
    while (getline(input, line)) {
      stringstream stream(line);
      //stream >> one_skeleton[vertices[i]].color;
      stream >> f;
      func_color.emplace(i, f);
      i++;
    }
    color_name = color_file_name;
  }

 public:  // Set color from kth coordinate
          /** \brief Computes the function used to color the nodes of the simplicial complex from the k-th coordinate.
           *
           * @param[in] k coordinate to use (start at 0).
           *
           */
  void set_color_from_coordinate(int k = 0) {
    for (int i = 0; i < n; i++) func_color[i] = point_cloud[i][k];
    color_name = "coordinate ";
    color_name.append(to_string(k));
  }

 public:  // Set color from vector.
  /** \brief Computes the function used to color the nodes of the simplicial complex from a vector stored in memory.
   *
   * @param[in] color input vector of values.
   *
   */
  void set_color_from_vector(vector<double> c) {
    for (unsigned int i = 0; i < c.size(); i++) func_color[i] = c[i];
  }

 public:  // Create a .dot file that can be compiled with neato to produce a .pdf file.
  /** \brief Creates a .dot file called SC.dot for neato (part of the graphviz package) once the simplicial complex is
   * computed to get a visualization
   * of its 1-skeleton in a .pdf file.
   */
  void plot_DOT() {

    char mapp[100]; sprintf(mapp, "%s_sc.dot",point_cloud_name.c_str()); ofstream graphic(mapp);

    double maxv = numeric_limits<double>::lowest(); double minv = numeric_limits<double>::max();
    for (map<int, pair<int, double> >::iterator iit = cover_color.begin(); iit != cover_color.end(); iit++) {
      maxv = max(maxv, iit->second.second); minv = min(minv, iit->second.second);
    }

    int k = 0; vector<int> nodes; nodes.clear();

    graphic << "graph GIC {" << endl;
    for (map<int, pair<int, double> >::iterator iit = cover_color.begin(); iit != cover_color.end(); iit++) {
      if (iit->second.first > mask) {
        nodes.push_back(iit->first);
        graphic << iit->first << "[shape=circle fontcolor=black color=black label=\"" << iit->first << ":"
                << iit->second.first << "\" style=filled fillcolor=\""
                << (1 - (maxv - iit->second.second) / (maxv - minv)) * 0.6 << ", 1, 1\"]" << std::endl;
        k++;
      }
    }
    int ke = 0;
    int num_simplices = simplices.size();
    for (int i = 0; i < num_simplices; i++)
      if (simplices[i].size() == 2) {
        if (cover_color[simplices[i][0]].first > mask && cover_color[simplices[i][1]].first > mask) {
          graphic << "  " << simplices[i][0] << " -- " << simplices[i][1] << " [weight=15];" << std::endl;
          ke++;
        }
      }
    graphic << "}";
    graphic.close();
    cout << ".dot file generated. It can be visualized with e.g. neato." << endl;
  }

 public:  // Create a .txt file that can be compiled with KeplerMapper.
          /** \brief Creates a .txt file called SC.txt describing the 1-skeleton, which can then be plotted with e.g.
           * KeplerMapper.
           */
  void write_info() {

    int num_simplices = simplices.size(); int num_edges = 0;
    char mapp[100]; sprintf(mapp, "%s_sc.txt",point_cloud_name.c_str()); ofstream graphic(mapp);

    for (int i = 0; i < num_simplices; i++)
      if (simplices[i].size() == 2)
        if (cover_color[simplices[i][0]].first > mask && cover_color[simplices[i][1]].first > mask) num_edges++;

    graphic << point_cloud_name << endl;
    graphic << cover_name << endl;
    graphic << color_name << endl;
    graphic << resolution_double << " " << gain << endl;
    graphic << cover_color.size() << " " << num_edges << endl;

    for (map<int, pair<int, double> >::iterator iit = cover_color.begin(); iit != cover_color.end(); iit++)
      graphic << iit->first << " " << iit->second.second << " " << iit->second.first << std::endl;

    for (int i = 0; i < num_simplices; i++)
      if (simplices[i].size() == 2)
        if (cover_color[simplices[i][0]].first > mask && cover_color[simplices[i][1]].first > mask)
          graphic << simplices[i][0] << " " << simplices[i][1] << std::endl;
    graphic.close();
    cout << ".txt generated. It can be visualized with e.g. python KeplerMapperVisuFromTxtFile.py and firefox." << endl;

  }

 public:  // Create a .off file that can be visualized (e.g. with Geomview).
          /** \brief Creates a .off file called SC.off for 3D visualization, which contains the 2-skeleton of the GIC.
           * This function assumes that the cover has been computed with Voronoi. If data points are in 1D or 2D,
           * the remaining coordinates of the points embedded in 3D are set to 0.
           */
  void plot_OFF() {

    assert(cover_name == "Voronoi");

    int m = voronoi_subsamples.size(); int numedges = 0; int numfaces = 0; vector<vector<int> > edges, faces;
    int numsimplices = simplices.size();

    char gic[100]; sprintf(gic, "%s_sc.off",point_cloud_name.c_str()); ofstream graphic(gic);

    graphic << "OFF" << std::endl;
    for (int i = 0; i < numsimplices; i++) {
      if (simplices[i].size() == 2) {
        numedges++;
        edges.push_back(simplices[i]);
      }
      if (simplices[i].size() == 3) {
        numfaces++;
        faces.push_back(simplices[i]);
      }
    }
    graphic << m << " " << numedges + numfaces << std::endl;
    for (int i = 0; i < m; i++) {
      if (data_dimension <= 3) {
        for (int j = 0; j < data_dimension; j++) graphic << point_cloud[voronoi_subsamples[i]][j] << " ";
        for (int j = data_dimension; j < 3; j++) graphic << 0 << " ";
        graphic << std::endl;
      } else {
        for (int j = 0; j < 3; j++) graphic << point_cloud[voronoi_subsamples[i]][j] << " ";
      }
    }
    for (int i = 0; i < numedges; i++) graphic << 2 << " " << edges[i][0] << " " << edges[i][1] << std::endl;
    for (int i = 0; i < numfaces; i++)
      graphic << 3 << " " << faces[i][0] << " " << faces[i][1] << " " << faces[i][2] << std::endl;
    graphic.close();
    cout << ".off generated. It can be visualized with e.g. geomview." << endl;
  }

  // *******************************************************************************************************************
  // *******************************************************************************************************************

 public:
  /** \brief Creates the simplicial complex.
   *
   * @param[in] complex SimplicialComplex to be created.
   *
   */
  template <typename SimplicialComplex>
  void create_complex(SimplicialComplex& complex){
    unsigned int dimension = 0;
    for (auto const& simplex : simplices) {
      int numvert = simplex.size(); double filt = std::numeric_limits<double>::lowest();
      for(int i = 0; i < numvert; i++)  filt = std::max(cover_color[simplex[i]].second, filt);
      complex.insert_simplex_and_subfaces(simplex, filt);
      if (dimension < simplex.size() - 1) dimension = simplex.size() - 1;
    }
    complex.set_dimension(dimension);
  }

 public:
  /** \brief Computes the extended persistence diagram of the complex.
   *
   */
  template <typename SimplicialComplex>
  void compute_PD(){

    SimplicialComplex streef, streeb; unsigned int dimension = 0;
    for (auto const & simplex : simplices) {
      int numvert = simplex.size(); double filtM = numeric_limits<double>::lowest(); double filtm = filtM;
      for(int i = 0; i < numvert; i++){filtM = max(cover_std[simplex[i]], filtM); filtm = max(-cover_std[simplex[i]], filtm);}
      streef.insert_simplex_and_subfaces(simplex, filtM); streeb.insert_simplex_and_subfaces(simplex, filtm);
      if (dimension < simplex.size() - 1) dimension = simplex.size() - 1;
    } streef.set_dimension(dimension); streeb.set_dimension(dimension);

    streef.initialize_filtration();
    Gudhi::persistent_cohomology::Persistent_cohomology<SimplicialComplex, Gudhi::persistent_cohomology::Field_Zp> pcohf(streef);
    pcohf.init_coefficients(2); pcohf.compute_persistent_cohomology();
    pcohf.output_diagram();

    streeb.initialize_filtration();
    Gudhi::persistent_cohomology::Persistent_cohomology<SimplicialComplex, Gudhi::persistent_cohomology::Field_Zp> pcohb(streeb);
    pcohb.init_coefficients(2); pcohb.compute_persistent_cohomology();
    pcohb.output_diagram();

    //PD = pcohf.get_persistent_pairs();

  }

 public:
  /** \brief Computes bootstrapped distances distribution.
   *
   * @param[in] N number of bootstrap iterations.
   *
   */
  template <typename SimplicialComplex>
  void compute_distribution(int N = 100){

    if(distribution.size() >= N) std::cout << "Already done!" << std::endl;
    else{
      for(int i = 0; i < N-distribution.size(); i++){

        Cover_complex Cboot; Cboot.n = this->n; std::vector<int> boot(this->n);
        for(int j = 0; j < this->n; j++){
          double u = GetUniform(); int id = std::floor(u*(this->n)); boot[j] = id;
          Cboot.point_cloud[j] = this->point_cloud[id]; Cboot.func.emplace(j,this->func[id]);
        }
        for(int j = 0; j < n; j++){
          vector<double> dist(n);
          for(int k = 0; k < n; k++)
            dist[k] = distances[boot[j]][boot[k]];
          Cboot.distances.push_back(dist);
        }

        Cboot.set_graph_from_automatic_rips(Gudhi::Euclidean_distance());
        Cboot.set_automatic_resolution(); Cboot.set_gain(); Cboot.set_cover_from_function();
        Cboot.find_simplices(); Cboot.compute_PD<Gudhi::Simplex_tree<> >();

        distribution.push_back(Gudhi::persistence_diagram::bottleneck_distance(this->PD,Cboot.PD));

      }

      std::sort(distribution.begin(), distribution.end());

    }
  }

 public:
  /** \brief Computes the bottleneck distance corresponding to a specific confidence level.
   *
   * @param[in] alpha Confidence level.
   *
   */
  double compute_distance_from_confidence_level(double alpha){
    int N = distribution.size();
    return distribution[std::floor(alpha*N)];
  }

 public:
  /** \brief Computes the confidence level of a specific bottleneck distance.
   *
   * @param[in] d Bottleneck distance.
   *
   */
  double compute_confidence_level_from_distance(double d){
    int N = distribution.size();
    for(int i = 0; i < N; i++)  if(distribution[i] > d)  return i*1.0/N;
  }

 public:
  /** \brief Computes the p-value, i.e. the opposite of the confidence level of the largest bottleneck
   * distance preserving the points in the persistence diagram of the output simplicial complex.
   *
   */
  double compute_p_value(){
    double distancemin = -std::numeric_limits<double>::lowest();
    int N = PD.size(); for(int i = 0; i < N; i++)  distancemin = std::min(distancemin, 0.5*(PD[i].second - PD[i].first));
    return 1-compute_confidence_level_from_distance(distancemin);
  }

 public:
  /** \brief Computes the simplices of the simplicial complex.
   */
  void find_simplices() {
    if (type != "Nerve" && type != "GIC") {
      cout << "Type of complex needs to be specified." << endl;
      return;
    }

    if (type == "Nerve") {
      for (map<int, vector<int> >::iterator it = cover.begin(); it != cover.end(); it++)
        simplices.push_back(it->second);
      sort(simplices.begin(), simplices.end());
      vector<vector<int> >::iterator it = unique(simplices.begin(), simplices.end());
      simplices.resize(distance(simplices.begin(), it));
    }

    if (type == "GIC") {

      IndexMap index = get(vertex_index, one_skeleton);

      if (functional_cover) {
        // Computes the simplices in the GIC by looking at all the edges of the graph and adding the
        // corresponding edges in the GIC if the images of the endpoints belong to consecutive intervals.

        if (gain >= 0.5)
          throw invalid_argument("the output of this function is correct ONLY if the cover is minimal, i.e. the gain is less than 0.5.");

        // Loop on all edges.
        graph_traits<Graph>::edge_iterator ei, ei_end;
        for (tie(ei, ei_end) = edges(one_skeleton); ei != ei_end; ++ei){
          int nums = cover[index[source(*ei, one_skeleton)]].size();
          for(int i = 0; i < nums; i++){
            int vs = cover[index[source(*ei, one_skeleton)]][i];
            int numt = cover[index[target(*ei, one_skeleton)]].size();
            for(int j = 0; j < numt; j++){
              int vt = cover[index[target(*ei, one_skeleton)]][j];
              if(cover_fct[vs] == cover_fct[vt] + 1 || cover_fct[vt] == cover_fct[vs] + 1){
                vector<int> edge(2); edge[0] = vs; edge[1] = vt; simplices.push_back(edge);
              }
            }
          }
        }
        sort(simplices.begin(), simplices.end());
        vector<vector<int> >::iterator it = unique(simplices.begin(), simplices.end());
        simplices.resize(distance(simplices.begin(), it));

      } else {

        // Find edges to keep
        Simplex_tree st; graph_traits<Graph>::edge_iterator ei, ei_end;
        for (tie(ei, ei_end) = edges(one_skeleton); ei != ei_end; ++ei)
          if( !( cover[index[target(*ei, one_skeleton)]].size() == 1 &&
                 cover[index[target(*ei, one_skeleton)]] == cover[index[source(*ei, one_skeleton)]]) ){
            vector<int> edge(2); edge[0] = index[source(*ei, one_skeleton)]; edge[1] = index[target(*ei, one_skeleton)];
            st.insert_simplex_and_subfaces(edge);
          }

        //st.insert_graph(one_skeleton);

        // Build the Simplex Tree corresponding to the graph
        st.expansion(maximal_dim);

        // Find simplices of GIC
        simplices.clear();
        for (auto simplex : st.complex_simplex_range()) {
          if (!st.has_children(simplex)) {
            vector<int> simplx;
            for (auto vertex : st.simplex_vertex_range(simplex)) {
              unsigned int sz = cover[vertex].size();
              for (unsigned int i = 0; i < sz; i++) {
                simplx.push_back(cover[vertex][i]);
              }
            }
            sort(simplx.begin(), simplx.end());
            vector<int>::iterator it = unique(simplx.begin(), simplx.end());
            simplx.resize(distance(simplx.begin(), it));
            simplices.push_back(simplx);
          }
        }
        sort(simplices.begin(), simplices.end());
        vector<vector<int> >::iterator it = unique(simplices.begin(), simplices.end());
        simplices.resize(distance(simplices.begin(), it));

      }
    }
  }
};

}  // namespace cover_complex

}  // namespace Gudhi

#endif  // GIC_H_

















/*Old code.

 private:
  void fill_adjacency_matrix_from_st() {
    std::vector<int> empty;
    for (int i = 0; i < n; i++) adjacency_matrix[i] = empty;
    for (auto simplex : st.complex_simplex_range()) {
      if (st.dimension(simplex) == 1) {
        std::vector<int> vertices;
        for (auto vertex : st.simplex_vertex_range(simplex)) vertices.push_back(vertex);
        adjacency_matrix[vertices[0]].push_back(vertices[1]);
        adjacency_matrix[vertices[1]].push_back(vertices[0]);
      }
    }
  }

std::vector<int> simplex_to_remove;
int simplex_id = 0;
for (auto simplex : st.complex_simplex_range()) {
  if (st.dimension(simplex) == 1) {
    std::vector<std::vector<int> > comp;
    for (auto vertex : st.simplex_vertex_range(simplex)) comp.push_back(cover[vertex]);
    if (comp[0].size() == 1 && comp[0] == comp[1]) simplex_to_remove.push_back(simplex_id);
  }
  simplex_id++;
}

// Remove edges
if (simplex_to_remove.size() > 1) {
  int current_id = 1;
  auto simplex = st.complex_simplex_range().begin();
  int num_rem = 0;
  for (int i = 0; i < simplex_id - 1; i++) {
    int j = i + 1;
    auto simplex_tmp = simplex;
    simplex_tmp++;
    if (j == simplex_to_remove[current_id]) {
      st.remove_maximal_simplex(*simplex_tmp);
      current_id++;
      num_rem++;
    } else {
      simplex++;
    }
  }
  simplex = st.complex_simplex_range().begin();
  for (int i = 0; i < simplex_to_remove[0]; i++) simplex++;
  st.remove_maximal_simplex(*simplex);
}


if(cover[index[source(*ei, one_skeleton)]].size() == 1){
  vs = cover[index[source(*ei, one_skeleton)]][0];
  vm = cover_fct[vs];
}
else{
  vs0 = cover[index[source(*ei, one_skeleton)]][0];
  vs1 = cover[index[source(*ei, one_skeleton)]][1];
  vm = min(cover_fct[vs0], cover_fct[vs1]);
  if(vm == cover_fct[vs0])  vs = vs0;  else  vs = vs1;
}

if(cover[index[target(*ei, one_skeleton)]].size() == 1){
  vt = cover[index[target(*ei, one_skeleton)]][0];
  vM = cover_fct[vt];
}
else{
  vt0 = cover[index[target(*ei, one_skeleton)]][0];
  vt1 = cover[index[target(*ei, one_skeleton)]][1];
  vM = max(cover_fct[vt0], cover_fct[vt1]);
  if(vM == cover_fct[vt0])  vt = vt0;  else  vt = vt1;
}

if(vM == vm + 1){
  //if(max(cover_fct[cover[index[target(*ei, one_skeleton)]][0]], cover_fct[cover[index[target(*ei, one_skeleton)]][1]])== min(cover_fct[index[source(*ei, one_skeleton)]][0], cover_fct[index[source(*ei, one_skeleton)]][1]) + 1){
  vector<int> edge(2); edge[0] = vs; edge[1] = vt;
  simplices.push_back(edge);
}


for (map<int, vector<int> >::iterator it = cover.begin(); it != cover.end(); it++) {
int vid = it->first;
vector<int> neighbors = adjacency_matrix[vid];
int num_neighb = neighbors.size();

// Find cover of current point (vid).
if (cover[vid].size() == 2)  v1 = std::min(cover[vid][0], cover[vid][1]);
else  v1 = cover[vid][0];
vector<int> node(1);  node[0] = v1;
simplices.push_back(node);

// Loop on neighbors.
for (int i = 0; i < num_neighb; i++) {
  int neighb = neighbors[i];

  // Find cover of neighbor (neighb).
  if (cover[neighb].size() == 2)  v2 = std::max(cover[neighb][0], cover[neighb][1]);
  else  v2 = cover[neighb][0];

  // If neighbor is in next interval, add edge.
  if (cover_fct[v2] == cover_fct[v1] + 1) {
    vector<int> edge(2); edge[0] = v1; edge[1] = v2;
    simplices.push_back(edge); break;
  }
}
}


  std::vector<double> dist(n);
  std::vector<int> process(n);
  for (int j = 0; j < n; j++) {
    dist[j] = std::numeric_limits<double>::max();
    process[j] = j;
  }
  dist[seed] = 0;
  int curr_size = process.size();
  int min_point, min_index;
  double min_dist;
  std::vector<int> neighbors;
  int num_neighbors;

  while (curr_size > 0) {
    min_dist = std::numeric_limits<double>::max();
    min_index = -1;
    min_point = -1;
    for (int j = 0; j < curr_size; j++) {
      if (dist[process[j]] < min_dist) {
        min_point = process[j];
        min_dist = dist[process[j]];
        min_index = j;
      }
    }
    assert(min_index != -1);
    process.erase(process.begin() + min_index);
    assert(min_point != -1);
    neighbors = adjacency_matrix[min_point];
    num_neighbors = neighbors.size();
    for (int j = 0; j < num_neighbors; j++) {
      double d = dist[min_point] + distances[min_point][neighbors[j]];
      dist[neighbors[j]] = std::min(dist[neighbors[j]], d);
    }
    curr_size = process.size();
  }


  // Compute the connected components with DFS
  std::map<int, bool> visit;
  if (verbose) std::cout << "Preimage of interval " << i << std::endl;
  for (std::map<int, std::vector<int> >::iterator it = prop.begin(); it != prop.end(); it++)
    visit[it->first] = false;
  if (!(prop.empty())) {
    for (std::map<int, std::vector<int> >::iterator it = prop.begin(); it != prop.end(); it++) {
      if (!(visit[it->first])) {
        std::vector<int> cc;
        cc.clear();
        dfs(prop, it->first, cc, visit);
        int cci = cc.size();
        if (verbose) std::cout << "one CC with " << cci << " points, ";
        double average_col = 0;
        for (int j = 0; j < cci; j++) {
          cover[cc[j]].push_back(id);
          cover_back[id].push_back(cc[j]);
          average_col += func_color[cc[j]] / cci;
        }
        cover_fct[id]   = i;
        cover_std[id]   = std::pair<int, double>(cci, 0.5*(u+v));
        cover_color[id] = std::pair<int, double>(cci, average_col);
        id++;
      }
    }
  }


  // DFS
  private:
   void dfs(std::map<int, std::vector<int> >& G, int p, std::vector<int>& cc, std::map<int, bool>& visit) {
     cc.push_back(p);
     visit[p] = true;
     int neighb = G[p].size();
     for (int i = 0; i < neighb; i++)
       if (visit.find(G[p][i]) != visit.end())
         if (!(visit[G[p][i]])) dfs(G, G[p][i], cc, visit);
   }
*/

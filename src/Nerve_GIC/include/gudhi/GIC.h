/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author:       Mathieu Carriere
 *
 *    Copyright (C) 2017 Inria
 *
 *    Modification(s):
 *      - 2019/08 Vincent Rouvreau: Fix issue #10 for CGAL
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef GIC_H_
#define GIC_H_

#ifdef GUDHI_USE_TBB
#include <tbb/parallel_for.h>
#include <mutex>
#endif

#if __has_include(<CGAL/version.h>)
# define GUDHI_GIC_USE_CGAL 1
# include <gudhi/Bottleneck.h>
#elif __has_include(<hera/bottleneck.h>)
# define GUDHI_GIC_USE_HERA 1
# ifdef _MSC_VER
// https://github.com/grey-narn/hera/issues/3
// ssize_t is a non-standard type (well, posix)
#  include <type_traits>
#  include <cstdlib>
using ssize_t = std::make_signed_t<std::size_t>;
# endif
# include <hera/bottleneck.h>
#endif

#include <gudhi/Debug_utils.h>
#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/reader_utils.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Rips_complex.h>
#include <gudhi/Points_off_io.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Persistent_cohomology.h>

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
#include <utility>    // for std::pair<>
#include <algorithm>  // for (std::max)
#include <random>
#include <cassert>
#include <cmath>

namespace Gudhi {

namespace cover_complex {

using Simplex_tree = Gudhi::Simplex_tree<>;
using Filtration_value = Simplex_tree::Filtration_value;
using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;
using Persistence_diagram = std::vector<std::pair<double, double> >;
using Graph = boost::subgraph<
    boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, boost::no_property,
                          boost::property<boost::edge_index_t, int, boost::property<boost::edge_weight_t, double> > > >;
using Vertex_t = boost::graph_traits<Graph>::vertex_descriptor;
using Index_map = boost::property_map<Graph, boost::vertex_index_t>::type;
using Weight_map = boost::property_map<Graph, boost::edge_weight_t>::type;

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
  std::string type;      // Nerve or GIC

  std::vector<Point> point_cloud;               // input point cloud.
  std::vector<std::vector<double> > distances;  // all pairwise distances.
  int maximal_dim;                              // maximal dimension of output simplicial complex.
  int data_dimension;                           // dimension of input data.
  int num_points;                                        // number of points.

  std::vector<double> func;       // function used to compute the output simplicial complex.
  std::vector<double> func_color; // function used to compute the colors of the nodes of the output simplicial complex.
  bool functional_cover = false;  // whether we use a cover with preimages of a function or not.

  Graph one_skeleton_OFF;          // one-skeleton given by the input OFF file (if it exists).
  Graph one_skeleton;              // one-skeleton used to compute the connected components.
  std::vector<Vertex_t> vertices;  // vertices of one_skeleton.

  std::vector<std::vector<int> > simplices;  // simplices of output simplicial complex.
  std::vector<int> voronoi_subsamples;       // Voronoi germs (in case of Voronoi cover).

  Persistence_diagram PD;
  std::vector<double> distribution;

  std::vector<std::vector<int> >
      cover;  // function associating to each data point the vector of cover elements to which it belongs.
  std::map<int, std::vector<int> >
      cover_back;  // inverse of cover, in order to get the data points associated to a specific cover element.
  std::map<int, double> cover_std;  // standard function (induced by func) used to compute the extended persistence
                                    // diagram of the output simplicial complex.
  std::map<int, int>
      cover_fct;  // integer-valued function that allows to state if two elements of the cover are consecutive or not.
  std::map<int, std::pair<int, double> >
      cover_color;  // size and coloring (induced by func_color) of the vertices of the output simplicial complex.

  int resolution_int = -1;
  double resolution_double = -1;
  double gain = -1;
  double rate_constant = 10;  // Constant in the subsampling.
  double rate_power = 0.001;  // Power in the subsampling.
  int mask = 0;               // Ignore nodes containing less than mask points.

  std::map<int, int> name2id, name2idinv;

  std::string cover_name;
  std::string point_cloud_name;
  std::string color_name;

  // Remove all edges of a graph.
  void remove_edges(Graph& G) {
    boost::graph_traits<Graph>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(G); ei != ei_end; ++ei) boost::remove_edge(*ei, G);
  }

  // Find random number in [0,1].
  double GetUniform() {
    thread_local std::default_random_engine re;
    std::uniform_real_distribution<double> Dist(0, 1);
    return Dist(re);
  }

  // Subsample points.
  void SampleWithoutReplacement(int populationSize, int sampleSize, std::vector<int>& samples) {
    int t = 0;
    int m = 0;
    double u;
    while (m < sampleSize) {
      u = GetUniform();
      if ((populationSize - t) * u >= sampleSize - m) {
        t++;
      } else {
        samples[m] = t;
        t++;
        m++;
      }
    }
  }

  // *******************************************************************************************************************
  // Utils.
  // *******************************************************************************************************************

 public:
  /** \brief Specifies whether the type of the output simplicial complex.
   *
   * @param[in] t std::string (either "GIC" or "Nerve").
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


  /** \brief Reads and stores the input point cloud from vector stored in memory.
   *
   * @param[in] point_cloud input vector representing the point cloud. Each row is a point and each coordinate is a vector.
   *
   */
  void set_point_cloud_from_range(const std::vector<std::vector<double> > & point_cloud) {
    this->num_points = point_cloud.size(); data_dimension = point_cloud[0].size();
    point_cloud_name = "cloud"; cover.resize(this->num_points);
    for(int i = 0; i < this->num_points; i++){
      boost::add_vertex(one_skeleton_OFF);
      vertices.push_back(boost::add_vertex(one_skeleton));
    }
    this->point_cloud = point_cloud;
  }

  /** \brief Reads and stores the input point cloud from .(n)OFF file.
   *
   * @param[in] off_file_name name of the input .OFF or .nOFF file.
   *
   */
  bool read_point_cloud(const std::string& off_file_name) {
    point_cloud_name = off_file_name;
    std::ifstream input(off_file_name);
    std::string line;

    char comment = '#';
    while (comment == '#') {
      std::getline(input, line);
      if (!line.empty() && !all_of(line.begin(), line.end(), (int (*)(int))isspace))
        comment = line[line.find_first_not_of(' ')];
    }
    if (strcmp((char*)line.c_str(), "nOFF") == 0) {
      comment = '#';
      while (comment == '#') {
        std::getline(input, line);
        if (!line.empty() && !all_of(line.begin(), line.end(), (int (*)(int))isspace))
          comment = line[line.find_first_not_of(' ')];
      }
      std::stringstream stream(line);
      stream >> data_dimension;
    } else {
      data_dimension = 3;
    }

    comment = '#';
    int numedges, numfaces, i, dim;
    while (comment == '#') {
      std::getline(input, line);
      if (!line.empty() && !all_of(line.begin(), line.end(), (int (*)(int))isspace))
        comment = line[line.find_first_not_of(' ')];
    }
    std::stringstream stream(line);
    stream >> this->num_points;
    stream >> numfaces;
    stream >> numedges;

    i = 0;
    while (i < this->num_points) {
      std::getline(input, line);
      if (!line.empty() && line[line.find_first_not_of(' ')] != '#' &&
          !all_of(line.begin(), line.end(), (int (*)(int))isspace)) {
        std::stringstream iss(line);
        std::vector<double> point;
        point.assign(std::istream_iterator<double>(iss), std::istream_iterator<double>());
        point_cloud.emplace_back(point.begin(), point.begin() + data_dimension);
        boost::add_vertex(one_skeleton_OFF);
        vertices.push_back(boost::add_vertex(one_skeleton));
        cover.emplace_back();
        i++;
      }
    }

    i = 0;
    while (i < numfaces) {
      std::getline(input, line);
      if (!line.empty() && line[line.find_first_not_of(' ')] != '#' &&
          !all_of(line.begin(), line.end(), (int (*)(int))isspace)) {
        std::vector<int> simplex;
        std::stringstream iss(line);
        simplex.assign(std::istream_iterator<int>(iss), std::istream_iterator<int>());
        dim = simplex[0];
        for (int j = 1; j <= dim; j++)
          for (int k = j + 1; k <= dim; k++)
            boost::add_edge(vertices[simplex[j]], vertices[simplex[k]], one_skeleton_OFF);
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
  void set_graph_from_file(const std::string& graph_file_name) {
    remove_edges(one_skeleton);
    int neighb;
    std::ifstream input(graph_file_name);
    std::string line;
    int source;
    while (std::getline(input, line)) {
      std::stringstream stream(line);
      stream >> source;
      while (stream >> neighb) boost::add_edge(vertices[source], vertices[neighb], one_skeleton);
    }
  }

 public:  // Set graph from OFF file.
          /** \brief Creates a graph G from the triangulation given by the input .OFF file.
           *
           */
  void set_graph_from_OFF() {
    remove_edges(one_skeleton);
    if (num_edges(one_skeleton_OFF))
      one_skeleton = one_skeleton_OFF;
    else
      std::cerr << "No triangulation read in OFF file!" << std::endl;
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
    remove_edges(one_skeleton);
    if (distances.size() == 0) compute_pairwise_distances(distance);
    for (int i = 0; i < this->num_points; i++) {
      for (int j = i + 1; j < this->num_points; j++) {
        if (distances[i][j] <= threshold) {
          boost::add_edge(vertices[i], vertices[j], one_skeleton);
          boost::put(boost::edge_weight, one_skeleton, boost::edge(vertices[i], vertices[j], one_skeleton).first,
                     distances[i][j]);
        }
      }
    }
  }

 public:
  void set_graph_weights() {
    Index_map index = boost::get(boost::vertex_index, one_skeleton);
    Weight_map weight = boost::get(boost::edge_weight, one_skeleton);
    boost::graph_traits<Graph>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(one_skeleton); ei != ei_end; ++ei)
      boost::put(weight, *ei,
                 distances[index[boost::source(*ei, one_skeleton)]][index[boost::target(*ei, one_skeleton)]]);
  }

 public:
  /** \brief Reads and stores the distance matrices from vector stored in memory.
   *
   * @param[in] distance_matrix input vector representing the distance matrix.
   *
   */
  void set_distances_from_range(const std::vector<std::vector<double> > & distance_matrix) {
    this->num_points = distance_matrix.size(); data_dimension = 0; point_cloud_name = "matrix";
    cover.resize(this->num_points); point_cloud.resize(this->num_points);
    for(int i = 0; i < this->num_points; i++){
      boost::add_vertex(one_skeleton_OFF);
      vertices.push_back(boost::add_vertex(one_skeleton));
    }
    distances = distance_matrix;
  }

 public:  // Pairwise distances.
          /** \private \brief Computes all pairwise distances.
           */
  template <typename Distance>
  void compute_pairwise_distances(Distance ref_distance) {
    std::vector<double> zeros(this->num_points);
    for (int i = 0; i < this->num_points; i++) distances.push_back(zeros);
    if (verbose) std::clog << "Computing distances..." << std::endl;
    for (int i = 0; i < this->num_points; i++) {
      int state = 100 * (i + 1) / this->num_points;
      if (verbose && state % 10 == 0) std::clog << "\r" << state << "%" << std::flush;
      for (int j = i; j < this->num_points; j++) {
        double dis = ref_distance(point_cloud[i], point_cloud[j]);
        distances[i][j] = dis;
        distances[j][i] = dis;
      }
    }
    if (verbose) std::clog << std::endl;
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
    int m = floor(this->num_points / std::exp((1 + rate_power) * std::log(std::log(this->num_points) / std::log(rate_constant))));
    m = (std::min)(m, this->num_points - 1);
    double delta = 0;

    if (verbose) std::clog << this->num_points << " points in R^" << data_dimension << std::endl;
    if (verbose) std::clog << "Subsampling " << m << " points" << std::endl;

    if (distances.size() == 0) compute_pairwise_distances(distance);

    #ifdef GUDHI_USE_TBB
    std::mutex deltamutex;
    tbb::parallel_for(0, N, [&](int i){
        std::vector<int> samples(m);
        SampleWithoutReplacement(this->num_points, m, samples);
        double hausdorff_dist = 0;
        for (int j = 0; j < this->num_points; j++) {
          double mj = distances[j][samples[0]];
          for (int k = 1; k < m; k++) mj = (std::min)(mj, distances[j][samples[k]]);
          hausdorff_dist = (std::max)(hausdorff_dist, mj);
        }
        deltamutex.lock();
        delta += hausdorff_dist / N;
        deltamutex.unlock();
      });
    #else
      for (int i = 0; i < N; i++) {
        std::vector<int> samples(m);
        SampleWithoutReplacement(this->num_points, m, samples);
        double hausdorff_dist = 0;
        for (int j = 0; j < this->num_points; j++) {
          double mj = distances[j][samples[0]];
          for (int k = 1; k < m; k++) mj = (std::min)(mj, distances[j][samples[k]]);
          hausdorff_dist = (std::max)(hausdorff_dist, mj);
        }
        delta += hausdorff_dist / N;
      }
    #endif

    if (verbose) std::clog << "delta = " << delta << std::endl;
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
    int i = 0;
    std::ifstream input(func_file_name);
    std::string line;
    double f;
    while (std::getline(input, line)) {
      std::stringstream stream(line);
      stream >> f;
      func.push_back(f);
      i++;
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
    if(point_cloud[0].size() > 0){
      for (int i = 0; i < this->num_points; i++) func.push_back(point_cloud[i][k]);
      functional_cover = true;
      cover_name = "coordinate " + std::to_string(k);
    }
    else{
      std::cerr << "Only pairwise distances provided---cannot access " << k << "th coordinate; returning null vector instead" << std::endl;
      for (int i = 0; i < this->num_points; i++) func.push_back(0.0);
      functional_cover = true;
      cover_name = "null";
    }
  }

 public:  // Set function from vector.
          /** \brief Creates the function f from a vector stored in memory.
           *
           * @param[in] function input vector of values.
           *
           */
  template <class InputRange>
  void set_function_from_range(InputRange const& function) {
    for (int i = 0; i < this->num_points; i++) func.push_back(function[i]);
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
      std::cerr << "Cover needs to come from the preimages of a function." << std::endl;
      return 0;
    }
    if (type != "Nerve" && type != "GIC") {
      std::cerr << "Type of complex needs to be specified." << std::endl;
      return 0;
    }

    double reso = 0;
    Index_map index = boost::get(boost::vertex_index, one_skeleton);

    if (type == "GIC") {
      boost::graph_traits<Graph>::edge_iterator ei, ei_end;
      for (boost::tie(ei, ei_end) = boost::edges(one_skeleton); ei != ei_end; ++ei)
        reso = (std::max)(reso, std::abs(func[index[boost::source(*ei, one_skeleton)]] -
                                       func[index[boost::target(*ei, one_skeleton)]]));
      if (verbose) std::clog << "resolution = " << reso << std::endl;
      resolution_double = reso;
    }

    if (type == "Nerve") {
      boost::graph_traits<Graph>::edge_iterator ei, ei_end;
      for (boost::tie(ei, ei_end) = boost::edges(one_skeleton); ei != ei_end; ++ei)
        reso = (std::max)(reso, std::abs(func[index[boost::source(*ei, one_skeleton)]] -
                                       func[index[boost::target(*ei, one_skeleton)]]) /
                                  gain);
      if (verbose) std::clog << "resolution = " << reso << std::endl;
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
      std::cerr << "Number and/or length of intervals not specified" << std::endl;
      return;
    }
    if (gain == -1) {
      std::cerr << "Gain not specified" << std::endl;
      return;
    }

    // Read function values and compute min and max
    double minf = (std::numeric_limits<float>::max)();
    double maxf = std::numeric_limits<float>::lowest();
    for (int i = 0; i < this->num_points; i++) {
      minf = (std::min)(minf, func[i]);
      maxf = (std::max)(maxf, func[i]);
    }
    if (verbose) std::clog << "Min function value = " << minf << " and Max function value = " << maxf << std::endl;

    // Compute cover of im(f)
    std::vector<std::pair<double, double> > intervals;
    int res;

    if (resolution_double == -1) {  // Case we use an integer for the number of intervals.
      double incr = (maxf - minf) / resolution_int;
      double x = minf;
      double alpha = (incr * gain) / (2 - 2 * gain);
      double y = minf + incr + alpha;
      std::pair<double, double> interm(x, y);
      intervals.push_back(interm);
      for (int i = 1; i < resolution_int - 1; i++) {
        x = minf + i * incr - alpha;
        y = minf + (i + 1) * incr + alpha;
        std::pair<double, double> inter(x, y);
        intervals.push_back(inter);
      }
      x = minf + (resolution_int - 1) * incr - alpha;
      y = maxf;
      std::pair<double, double> interM(x, y);
      intervals.push_back(interM);
      res = intervals.size();
      if (verbose) {
        for (int i = 0; i < res; i++)
          std::clog << "Interval " << i << " = [" << intervals[i].first << ", " << intervals[i].second << "]"
                    << std::endl;
      }
    } else {
      if (resolution_int == -1) {  // Case we use a double for the length of the intervals.
        double x = minf;
        double y = x + resolution_double;
        while (y <= maxf && maxf - (y - gain * resolution_double) >= resolution_double) {
          std::pair<double, double> inter(x, y);
          intervals.push_back(inter);
          x = y - gain * resolution_double;
          y = x + resolution_double;
        }
        std::pair<double, double> interM(x, maxf);
        intervals.push_back(interM);
        res = intervals.size();
        if (verbose) {
          for (int i = 0; i < res; i++)
            std::clog << "Interval " << i << " = [" << intervals[i].first << ", " << intervals[i].second << "]"
                      << std::endl;
        }
      } else {  // Case we use an integer and a double for the length of the intervals.
        double x = minf;
        double y = x + resolution_double;
        int count = 0;
        while (count < resolution_int && y <= maxf && maxf - (y - gain * resolution_double) >= resolution_double) {
          std::pair<double, double> inter(x, y);
          intervals.push_back(inter);
          count++;
          x = y - gain * resolution_double;
          y = x + resolution_double;
        }
        res = intervals.size();
        if (verbose) {
          for (int i = 0; i < res; i++)
            std::clog << "Interval " << i << " = [" << intervals[i].first << ", " << intervals[i].second << "]"
                      << std::endl;
        }
      }
    }

    // Sort points according to function values
    std::vector<int> points(this->num_points);
    for (int i = 0; i < this->num_points; i++) points[i] = i;
    std::sort(points.begin(), points.end(), [this](int p1, int p2){return (this->func[p1] < this->func[p2]);});

    int id = 0;
    int pos = 0;
    Index_map index = boost::get(boost::vertex_index, one_skeleton);  // int maxc = -1;
    std::map<int, std::vector<int> > preimages;
    std::map<int, double> funcstd;

    if (verbose) std::clog << "Computing preimages..." << std::endl;
    for (int i = 0; i < res; i++) {
      // Find points in the preimage
      std::pair<double, double> inter1 = intervals[i];
      int tmp = pos;
      double u, v;

      if (i != res - 1) {
        if (i != 0) {
          std::pair<double, double> inter3 = intervals[i - 1];
          while (func[points[tmp]] < inter3.second && tmp != this->num_points) {
            preimages[i].push_back(points[tmp]);
            tmp++;
          }
          u = inter3.second;
        } else {
          u = inter1.first;
        }

        std::pair<double, double> inter2 = intervals[i + 1];
        while (func[points[tmp]] < inter2.first && tmp != this->num_points) {
          preimages[i].push_back(points[tmp]);
          tmp++;
        }
        v = inter2.first;
        pos = tmp;
        while (func[points[tmp]] < inter1.second && tmp != this->num_points) {
          preimages[i].push_back(points[tmp]);
          tmp++;
        }

      } else {
        std::pair<double, double> inter3 = intervals[i - 1];
        while (func[points[tmp]] < inter3.second && tmp != this->num_points) {
          preimages[i].push_back(points[tmp]);
          tmp++;
        }
        while (tmp != this->num_points) {
          preimages[i].push_back(points[tmp]);
          tmp++;
        }
        u = inter3.second;
        v = inter1.second;
      }

      funcstd[i] = 0.5 * (u + v);
    }

    #ifdef GUDHI_USE_TBB
      if (verbose) std::clog << "Computing connected components (parallelized)..." << std::endl;
      std::mutex covermutex, idmutex;
      tbb::parallel_for(0, res, [&](int i){
        // Compute connected components
        Graph G = one_skeleton.create_subgraph();
        int num = preimages[i].size();
        std::vector<int> component(num);
        for (int j = 0; j < num; j++) boost::add_vertex(index[vertices[preimages[i][j]]], G);
        boost::connected_components(G, &component[0]);
        int max = 0;

        // For each point in preimage
        for (int j = 0; j < num; j++) {
          // Update number of components in preimage
          if (component[j] > max) max = component[j];

          // Identify component with Cantor polynomial N^2 -> N
          int identifier = ((i + component[j])*(i + component[j]) + 3 * i + component[j]) / 2;

          // Update covers
          covermutex.lock();
          cover[preimages[i][j]].push_back(identifier);
          cover_back[identifier].push_back(preimages[i][j]);
          cover_fct[identifier] = i;
          cover_std[identifier] = funcstd[i];
          cover_color[identifier].second += func_color[preimages[i][j]];
          cover_color[identifier].first += 1;
          covermutex.unlock();
        }

        // Maximal dimension is total number of connected components
        idmutex.lock();
        id += max + 1;
        idmutex.unlock();
      });
    #else
      if (verbose) std::clog << "Computing connected components..." << std::endl;
      for (int i = 0; i < res; i++) {
        // Compute connected components
        Graph G = one_skeleton.create_subgraph();
        int num = preimages[i].size();
        std::vector<int> component(num);
        for (int j = 0; j < num; j++) boost::add_vertex(index[vertices[preimages[i][j]]], G);
        boost::connected_components(G, &component[0]);
        int max = 0;

        // For each point in preimage
        for (int j = 0; j < num; j++) {
          // Update number of components in preimage
          if (component[j] > max) max = component[j];

          // Identify component with Cantor polynomial N^2 -> N
          int identifier = (std::pow(i + component[j], 2) + 3 * i + component[j]) / 2;

          // Update covers
          cover[preimages[i][j]].push_back(identifier);
          cover_back[identifier].push_back(preimages[i][j]);
          cover_fct[identifier] = i;
          cover_std[identifier] = funcstd[i];
          cover_color[identifier].second += func_color[preimages[i][j]];
          cover_color[identifier].first += 1;
        }

        // Maximal dimension is total number of connected components
        id += max + 1;
      }
    #endif

    maximal_dim = id - 1;
    for (std::map<int, std::pair<int, double> >::iterator iit = cover_color.begin(); iit != cover_color.end(); iit++)
      iit->second.second /= iit->second.first;
  }

 public:  // Set cover from file.
  /** \brief Creates the cover C from a file containing the cover elements of each point (the order has to be the same
  * as in the input file!).
  *
  * @param[in] cover_file_name name of the input cover file.
  *
  */
  void set_cover_from_file(const std::string& cover_file_name) {
    int i = 0;
    int cov;
    std::vector<int> cov_elts, cov_number;
    std::ifstream input(cover_file_name);
    std::string line;
    while (std::getline(input, line)) {
      cov_elts.clear();
      std::stringstream stream(line);
      while (stream >> cov) {
        cov_elts.push_back(cov);
        cov_number.push_back(cov);
        cover_fct[cov] = cov;
        cover_color[cov].second += func_color[i];
        cover_color[cov].first++;
        cover_back[cov].push_back(i);
      }
      cover[i] = cov_elts;
      i++;
    }

    std::sort(cov_number.begin(), cov_number.end());
    std::vector<int>::iterator it = std::unique(cov_number.begin(), cov_number.end());
    cov_number.resize(std::distance(cov_number.begin(), it));

    maximal_dim = cov_number.size() - 1;
    for (int i = 0; i <= maximal_dim; i++) cover_color[i].second /= cover_color[i].first;
    cover_name = cover_file_name;
  }

 public:  // Set cover from range.
  /** \brief Creates the cover C from a vector of assignments stored in memory. The assignments, or clusters IDs, are vectors of integers.
  *
  * @tparam AssignmentRange random access range of input ranges of integers.
  * @param[in] assignments the cover assignments, specifying for each point the list of clusters it belongs to.
  *
  */
  template <class AssignmentRange>
  void set_cover_from_range(AssignmentRange const& assignments) {    
    std::vector<int> cov_elts, cov_number;
    for(int i=0; i < static_cast<int>(assignments.size()); i++){
      cov_elts.clear();
      for (int cov : assignments[i]){
        cov_elts.push_back(cov);
        cov_number.push_back(cov);
        cover_fct[cov] = cov;
        auto& cc = cover_color[cov];
        cc.second += func_color[i];
        cc.first++;
        cover_back[cov].push_back(i);
      }
      cover[i] = cov_elts;
    }

    std::sort(cov_number.begin(), cov_number.end());
    std::vector<int>::iterator it = std::unique(cov_number.begin(), cov_number.end());
    cov_number.resize(std::distance(cov_number.begin(), it));

    maximal_dim = cov_number.size() - 1;
    for (int i = 0; i <= maximal_dim; i++) cover_color[i].second /= cover_color[i].first;
  }

 public:  // Set cover from Voronoi
          /** \brief Creates the cover C from the Voronoï cells of a subsampling of the point cloud.
          *
          * @param[in] distance distance between the points.
          * @param[in] m number of points in the subsample.
          *
          */
  template <typename Distance>
  void set_cover_from_Voronoi(Distance distance, int m = 100) {
    voronoi_subsamples.resize(m);
    SampleWithoutReplacement(this->num_points, m, voronoi_subsamples);
    if (distances.size() == 0) compute_pairwise_distances(distance);
    set_graph_weights();
    Weight_map weight = boost::get(boost::edge_weight, one_skeleton);
    Index_map index = boost::get(boost::vertex_index, one_skeleton);
    std::vector<double> mindist(this->num_points);
    for (int j = 0; j < this->num_points; j++) mindist[j] = (std::numeric_limits<double>::max)();

    // Compute the geodesic distances to subsamples with Dijkstra
    #ifdef GUDHI_USE_TBB
      if (verbose) std::clog << "Computing geodesic distances (parallelized)..." << std::endl;
      std::mutex coverMutex; std::mutex mindistMutex;
      tbb::parallel_for(0, m, [&](int i){
        int seed = voronoi_subsamples[i];
        std::vector<double> dmap(this->num_points);
        boost::dijkstra_shortest_paths(
            one_skeleton, vertices[seed],
            boost::weight_map(weight).distance_map(boost::make_iterator_property_map(dmap.begin(), index)));

        coverMutex.lock(); mindistMutex.lock();
        for (int j = 0; j < this->num_points; j++)
          if (mindist[j] > dmap[j]) {
            mindist[j] = dmap[j];
            if (cover[j].size() == 0)
              cover[j].push_back(i);
            else
              cover[j][0] = i;
          }
        coverMutex.unlock(); mindistMutex.unlock();
      });
    #else
      for (int i = 0; i < m; i++) {
        if (verbose) std::clog << "Computing geodesic distances to seed " << i << "..." << std::endl;
        int seed = voronoi_subsamples[i];
        std::vector<double> dmap(this->num_points);
        boost::dijkstra_shortest_paths(
            one_skeleton, vertices[seed],
            boost::weight_map(weight).distance_map(boost::make_iterator_property_map(dmap.begin(), index)));

        for (int j = 0; j < this->num_points; j++)
          if (mindist[j] > dmap[j]) {
            mindist[j] = dmap[j];
            if (cover[j].size() == 0)
              cover[j].push_back(i);
            else
              cover[j][0] = i;
          }
      }
    #endif

    for (int i = 0; i < this->num_points; i++) {
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
  const std::vector<int>& subpopulation(int c) { return cover_back[c]; }


 public:
          /** \brief Returns the mean color corresponding to a specific node of the created complex.
          *
          * @param[in] c ID of the node.
          * @result cover_color(c) mean color value.
          *
          */
  double subcolor(int c) { return cover_color[c].second; }

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
  void set_color_from_file(const std::string& color_file_name) {
    int i = 0;
    std::ifstream input(color_file_name);
    std::string line;
    double f;
    while (std::getline(input, line)) {
      std::stringstream stream(line);
      stream >> f;
      func_color.push_back(f);
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
    if(point_cloud[0].size() > 0){
      for (int i = 0; i < this->num_points; i++) func_color.push_back(point_cloud[i][k]);
      color_name = "coordinate ";
      color_name.append(std::to_string(k));
    }
    else{
      std::cerr << "Only pairwise distances provided---cannot access " << k << "th coordinate; returning null vector instead" << std::endl;
      for (int i = 0; i < this->num_points; i++) func.push_back(0.0);
      functional_cover = true;
      cover_name = "null";
    }
  }

 public:  // Set color from vector.
  /** \brief Computes the function used to color the nodes of the simplicial complex from a vector stored in memory.
   *
   * @param[in] color input vector of values.
   *
   */
  void set_color_from_range(std::vector<double> color) {
    for (unsigned int i = 0; i < color.size(); i++) func_color.push_back(color[i]);
  }

 public:  // Create a .dot file that can be compiled with neato to produce a .pdf file.
  /** \brief Creates a .dot file called SC.dot for neato (part of the graphviz package) once the simplicial complex is
   * computed to get a visualization
   * of its 1-skeleton in a .pdf file.
   */
  void plot_DOT() {
    std::string mapp = point_cloud_name + "_sc.dot";
    std::ofstream graphic(mapp);

    double maxv = std::numeric_limits<double>::lowest();
    double minv = (std::numeric_limits<double>::max)();
    for (std::map<int, std::pair<int, double> >::iterator iit = cover_color.begin(); iit != cover_color.end(); iit++) {
      maxv = (std::max)(maxv, iit->second.second);
      minv = (std::min)(minv, iit->second.second);
    }

    int k = 0;
    std::vector<int> nodes;
    nodes.clear();

    graphic << "graph GIC {" << std::endl;
    int id = 0;
    for (std::map<int, std::pair<int, double> >::iterator iit = cover_color.begin(); iit != cover_color.end(); iit++) {
      if (iit->second.first > mask) {
        nodes.push_back(iit->first);
        name2id[iit->first] = id;
        name2idinv[id] = iit->first;
        id++;
        graphic << name2id[iit->first] << "[shape=circle fontcolor=black color=black label=\"" << name2id[iit->first]
                << ":" << iit->second.first << "\" style=filled fillcolor=\""
                << (1 - (maxv - iit->second.second) / (maxv - minv)) * 0.6 << ", 1, 1\"]" << std::endl;
        k++;
      }
    }
    int ke = 0;
    int num_simplices = simplices.size();
    for (int i = 0; i < num_simplices; i++)
      if (simplices[i].size() == 2) {
        if (cover_color[simplices[i][0]].first > mask && cover_color[simplices[i][1]].first > mask) {
          graphic << "  " << name2id[simplices[i][0]] << " -- " << name2id[simplices[i][1]] << " [weight=15];"
                  << std::endl;
          ke++;
        }
      }
    graphic << "}";
    graphic.close();
    std::clog << mapp << " file generated. It can be visualized with e.g. neato." << std::endl;
  }

 public:  // Create a .txt file that can be compiled with KeplerMapper.
          /** \brief Creates a .txt file called SC.txt describing the 1-skeleton, which can then be plotted with e.g.
           * KeplerMapper.
           */
  void write_info() {
    int num_simplices = simplices.size();
    int num_edges = 0;
    std::string mapp = point_cloud_name + "_sc.txt";
    std::ofstream graphic(mapp);

    for (int i = 0; i < num_simplices; i++)
      if (simplices[i].size() == 2)
        if (cover_color[simplices[i][0]].first > mask && cover_color[simplices[i][1]].first > mask) num_edges++;

    graphic << point_cloud_name << std::endl;
    graphic << cover_name << std::endl;
    graphic << color_name << std::endl;
    graphic << resolution_double << " " << gain << std::endl;
    graphic << cover_color.size() << " " << num_edges << std::endl;

    int id = 0;
    for (std::map<int, std::pair<int, double> >::iterator iit = cover_color.begin(); iit != cover_color.end(); iit++) {
      graphic << id << " " << iit->second.second << " " << iit->second.first << std::endl;
      name2id[iit->first] = id;
      name2idinv[id] = iit->first;
      id++;
    }

    for (int i = 0; i < num_simplices; i++)
      if (simplices[i].size() == 2)
        if (cover_color[simplices[i][0]].first > mask && cover_color[simplices[i][1]].first > mask)
          graphic << name2id[simplices[i][0]] << " " << name2id[simplices[i][1]] << std::endl;
    graphic.close();
    std::clog << mapp
              << " generated. It can be visualized with e.g. python KeplerMapperVisuFromTxtFile.py and firefox."
              << std::endl;
  }

 public:  // Create a .off file that can be visualized (e.g. with Geomview).
          /** \brief Creates a .off file called SC.off for 3D visualization, which contains the 2-skeleton of the GIC.
           * This function assumes that the cover has been computed with Voronoi. If data points are in 1D or 2D,
           * the remaining coordinates of the points embedded in 3D are set to 0.
           */
  void plot_OFF() {
    assert(cover_name == "Voronoi");

    int m = voronoi_subsamples.size();
    int numedges = 0;
    int numfaces = 0;
    std::vector<std::vector<int> > edges, faces;
    int numsimplices = simplices.size();

    std::string mapp = point_cloud_name + "_sc.off";
    std::ofstream graphic(mapp);

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
    std::clog << mapp << " generated. It can be visualized with e.g. geomview." << std::endl;
  }

  // *******************************************************************************************************************
  // Extended Persistence Diagrams.
  // *******************************************************************************************************************

 public:
  /** \brief Computes the extended persistence diagram of the complex.
   *
   */
  Persistence_diagram compute_PD() {
    Simplex_tree st;

    // Compute max and min
    double maxf = std::numeric_limits<double>::lowest();
    double minf = (std::numeric_limits<double>::max)();
    for (std::map<int, double>::iterator it = cover_std.begin(); it != cover_std.end(); it++) {
      maxf = (std::max)(maxf, it->second);
      minf = (std::min)(minf, it->second);
    }

    // Build filtration
    for (auto const& simplex : simplices) {
      std::vector<int> splx = simplex;
      splx.push_back(-2);
      st.insert_simplex_and_subfaces(splx, -3);
    }

    for (std::map<int, double>::iterator it = cover_std.begin(); it != cover_std.end(); it++) {
      int vertex = it->first; float val = it->second;
      int vert[] = {vertex}; int edge[] = {vertex, -2};
      if(st.find(vert) != st.null_simplex()){
        st.assign_filtration(st.find(vert), -2 + (val - minf)/(maxf - minf));
        st.assign_filtration(st.find(edge),  2 - (val - minf)/(maxf - minf));
      }
    }
    st.make_filtration_non_decreasing();

    // Compute PD
    Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Gudhi::persistent_cohomology::Field_Zp> pcoh(st); pcoh.init_coefficients(2);
    pcoh.compute_persistent_cohomology();

    // Output PD
    int max_dim = st.dimension();
    for (int i = 0; i < max_dim; i++) {
      std::vector<std::pair<double, double> > bars = pcoh.intervals_in_dimension(i);
      int num_bars = bars.size(); if(i == 0)  num_bars -= 1;
      if(verbose)  std::clog << num_bars << " interval(s) in dimension " << i << ":" << std::endl;
      for (int j = 0; j < num_bars; j++) {
        double birth = bars[j].first;
        double death = bars[j].second;
        if (i == 0 && std::isinf(death)) continue;
        if (birth < 0)
          birth = minf + (birth + 2) * (maxf - minf);
        else
          birth = minf + (2 - birth) * (maxf - minf);
        if (death < 0)
          death = minf + (death + 2) * (maxf - minf);
        else
          death = minf + (2 - death) * (maxf - minf);
        PD.push_back(std::pair<double, double>(birth, death));
        if (verbose) std::clog << "  [" << birth << ", " << death << "]" << std::endl;
      }
    }
    return PD;
  }

 public:
  /** \brief Computes bootstrapped distances distribution.
   *
   * @param[in] N number of bootstrap iterations.
   *
   */
  void compute_distribution(unsigned int N = 100) {
    unsigned int sz = distribution.size();
    if (sz < N) {
      for (unsigned int i = 0; i < N - sz; i++) {
        if (verbose)  std::clog << "Computing " << i << "th bootstrap, bottleneck distance = ";

        Cover_complex Cboot; Cboot.num_points = this->num_points; Cboot.data_dimension = this->data_dimension; Cboot.type = this->type; Cboot.functional_cover = true;

        std::vector<int> boot(this->num_points);
        for (int j = 0; j < this->num_points; j++) {
          double u = GetUniform();
          int id = std::floor(u * (this->num_points)); boot[j] = id;
          Cboot.point_cloud.push_back(this->point_cloud[id]); Cboot.cover.emplace_back(); Cboot.func.push_back(this->func[id]);
          boost::add_vertex(Cboot.one_skeleton_OFF); Cboot.vertices.push_back(boost::add_vertex(Cboot.one_skeleton));
        }
        Cboot.set_color_from_range(Cboot.func);

        for (int j = 0; j < this->num_points; j++) {
          std::vector<double> dist(this->num_points);
          for (int k = 0; k < this->num_points; k++) dist[k] = distances[boot[j]][boot[k]];
          Cboot.distances.push_back(dist);
        }

        Cboot.set_graph_from_automatic_rips(Gudhi::Euclidean_distance());
        Cboot.set_gain();
        Cboot.set_automatic_resolution();
        Cboot.set_cover_from_function();
        Cboot.find_simplices();
        Cboot.compute_PD();
#ifdef GUDHI_GIC_USE_CGAL
        double db = Gudhi::persistence_diagram::bottleneck_distance(this->PD, Cboot.PD);
#elif defined GUDHI_GIC_USE_HERA
        double db = hera::bottleneckDistExact(this->PD, Cboot.PD);
#else
        double db;
        throw std::logic_error("This function requires CGAL or Hera for the bottleneck distance.");
#endif
        if (verbose)  std::clog << db << std::endl;
        distribution.push_back(db);
      }

      std::sort(distribution.begin(), distribution.end());
    }
  }

 public:
  /** \brief Computes the bottleneck distance threshold corresponding to a specific confidence level.
   *
   * @param[in] alpha Confidence level.
   *
   */
  double compute_distance_from_confidence_level(double alpha) {
    unsigned int N = distribution.size();
    double d = distribution[std::floor(alpha * N)];
    if (verbose)  std::clog << "Distance corresponding to confidence " << alpha << " is " << d << std::endl;
    return d;
  }

 public:
  /** \brief Computes the confidence level of a specific bottleneck distance threshold.
   *
   * @param[in] d Bottleneck distance.
   *
   */
  double compute_confidence_level_from_distance(double d) {
    unsigned int N = distribution.size();
    double level = 1;
    for (unsigned int i = 0; i < N; i++)
      if (distribution[i] >= d){ level = i * 1.0 / N; break; }
    if (verbose)  std::clog << "Confidence level of distance " << d << " is " << level << std::endl;
    return level;
  }

 public:
  /** \brief Computes the p-value, i.e. the opposite of the confidence level of the largest bottleneck
   * distance preserving the points in the persistence diagram of the output simplicial complex.
   *
   */
  double compute_p_value() {
    double distancemin = (std::numeric_limits<double>::max)(); int N = PD.size();
    for (int i = 0; i < N; i++) distancemin = (std::min)(distancemin, 0.5 * std::abs(PD[i].second - PD[i].first));
    double p_value = 1 - compute_confidence_level_from_distance(distancemin);
    if (verbose)  std::clog << "p value = " << p_value << std::endl;
    return p_value;
  }

  // *******************************************************************************************************************
  // Computation of simplices.
  // *******************************************************************************************************************

 public:
  /** \brief Creates the simplicial complex.
   *
   * @param[in] complex SimplicialComplex to be created.
   *
   */
  template <typename SimplicialComplex>
  void create_complex(SimplicialComplex& complex) {
    unsigned int dimension = 0;
    for (auto const& simplex : simplices) {
      int numvert = simplex.size();
      double filt = std::numeric_limits<double>::lowest();
      for (int i = 0; i < numvert; i++) filt = (std::max)(cover_color[simplex[i]].second, filt);
      complex.insert_simplex_and_subfaces(simplex, filt);
      if (dimension < simplex.size() - 1) dimension = simplex.size() - 1;
    }
  }

 public:
  /** \brief Computes the simplices of the simplicial complex.
   */
  void find_simplices() {
    if (type != "Nerve" && type != "GIC") {
      std::cerr << "Type of complex needs to be specified." << std::endl;
      return;
    }

    if (type == "Nerve") {
      simplices = cover;
      std::sort(simplices.begin(), simplices.end());
      std::vector<std::vector<int> >::iterator it = std::unique(simplices.begin(), simplices.end());
      simplices.resize(std::distance(simplices.begin(), it));
    }

    if (type == "GIC") {
      Index_map index = boost::get(boost::vertex_index, one_skeleton);

      if (functional_cover) {
        // Computes the simplices in the GIC by looking at all the edges of the graph and adding the
        // corresponding edges in the GIC if the images of the endpoints belong to consecutive intervals.

        if (gain >= 0.5)
          throw std::invalid_argument(
              "the output of this function is correct ONLY if the cover is minimal, i.e. the gain is less than 0.5.");

        // Loop on all edges.
        boost::graph_traits<Graph>::edge_iterator ei, ei_end;
        for (boost::tie(ei, ei_end) = boost::edges(one_skeleton); ei != ei_end; ++ei) {
          int nums = cover[index[boost::source(*ei, one_skeleton)]].size();
          for (int i = 0; i < nums; i++) {
            int vs = cover[index[boost::source(*ei, one_skeleton)]][i];
            int numt = cover[index[boost::target(*ei, one_skeleton)]].size();
            for (int j = 0; j < numt; j++) {
              int vt = cover[index[boost::target(*ei, one_skeleton)]][j];
              if (cover_fct[vs] == cover_fct[vt] + 1 || cover_fct[vt] == cover_fct[vs] + 1) {
                std::vector<int> edge(2);
                edge[0] = (std::min)(vs, vt);
                edge[1] = (std::max)(vs, vt);
                simplices.push_back(edge);
                goto afterLoop;
              }
            }
          }
        afterLoop:;
        }
        std::sort(simplices.begin(), simplices.end());
        std::vector<std::vector<int> >::iterator it = std::unique(simplices.begin(), simplices.end());
        simplices.resize(std::distance(simplices.begin(), it));

      } else {
        // Find edges to keep
        Simplex_tree st;
        boost::graph_traits<Graph>::edge_iterator ei, ei_end;
        for (boost::tie(ei, ei_end) = boost::edges(one_skeleton); ei != ei_end; ++ei)
          if (!(cover[index[boost::target(*ei, one_skeleton)]].size() == 1 &&
                cover[index[boost::target(*ei, one_skeleton)]] == cover[index[boost::source(*ei, one_skeleton)]])) {
            std::vector<int> edge(2);
            edge[0] = index[boost::source(*ei, one_skeleton)];
            edge[1] = index[boost::target(*ei, one_skeleton)];
            st.insert_simplex_and_subfaces(edge);
          }

        // st.insert_graph(one_skeleton);

        // Build the Simplex Tree corresponding to the graph
        st.expansion(maximal_dim);

        // Find simplices of GIC
        simplices.clear();
        for (auto simplex : st.complex_simplex_range()) {
          if (!st.has_children(simplex)) {
            std::vector<int> simplx;
            for (auto vertex : st.simplex_vertex_range(simplex)) {
              unsigned int sz = cover[vertex].size();
              for (unsigned int i = 0; i < sz; i++) {
                simplx.push_back(cover[vertex][i]);
              }
            }
            std::sort(simplx.begin(), simplx.end());
            std::vector<int>::iterator it = std::unique(simplx.begin(), simplx.end());
            simplx.resize(std::distance(simplx.begin(), it));
            simplices.push_back(simplx);
          }
        }
        std::sort(simplices.begin(), simplices.end());
        std::vector<std::vector<int> >::iterator it = std::unique(simplices.begin(), simplices.end());
        simplices.resize(std::distance(simplices.begin(), it));
      }
    }
  }
};

}  // namespace cover_complex

}  // namespace Gudhi

#endif  // GIC_H_

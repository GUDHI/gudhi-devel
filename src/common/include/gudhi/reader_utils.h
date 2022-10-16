/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clement Maria, Pawel Dlotko, Clement Jamin
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef READER_UTILS_H_
#define READER_UTILS_H_

#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/Debug_utils.h>

# include <boost/iterator/function_output_iterator.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <iostream>
#include <fstream>
#include <map>
#include <limits>  // for numeric_limits
#include <string>
#include <vector>
#include <utility>  // for pair
#include <tuple>  // for std::make_tuple

namespace Gudhi {

// Keep this file tag for Doxygen to parse the code, otherwise, functions are not documented.
// It is required for global functions and variables.

/** @file
 * @brief This file includes common file reader for GUDHI
 */

/**
 * @brief Read a set of points to turn it into a vector< vector<double> > by filling points.
 *
 * File format: 1 point per line<br>
 * X11 X12 ... X1d<br>
 * X21 X22 ... X2d<br>
 * etc<br>
 */
inline void read_points(std::string file_name, std::vector<std::vector<double>>& points) {
  std::ifstream in_file(file_name.c_str(), std::ios::in);
  if (!in_file.is_open()) {
    std::cerr << "Unable to open file " << file_name << std::endl;
    return;
  }

  std::string line;
  double x;
  while (getline(in_file, line)) {
    std::vector<double> point;
    std::istringstream iss(line);
    while (iss >> x) {
      point.push_back(x);
    }
    // Check for empty lines
    if (!point.empty()) points.push_back(point);
  }
  in_file.close();
}

/**
 * @brief Read a graph from a file.
 *
 * \tparam Graph_t Type for the return graph. Must be constructible from iterators on pairs of Vertex_handle
 * \tparam Filtration_value Type for the value of the read filtration
 * \tparam Vertex_handle Type for the value of the read vertices
 *
 * File format: 1 simplex per line<br>
 * Dim1 X11 X12 ... X1d Fil1<br>
 * Dim2 X21 X22 ... X2d Fil2<br>
 * etc<br>
 *
 * The vertices must be labeled from 0 to n-1.
 * Every simplex must appear exactly once.
 * Simplices of dimension more than 1 are ignored.
 */
template <typename Graph_t, typename Filtration_value, typename Vertex_handle>
Graph_t read_graph(std::string file_name) {
  std::ifstream in_(file_name.c_str(), std::ios::in);
  if (!in_.is_open()) {
    std::string error_str("read_graph - Unable to open file ");
    error_str.append(file_name);
    std::cerr << error_str << std::endl;
    throw std::invalid_argument(error_str);
  }

  typedef std::pair<Vertex_handle, Vertex_handle> Edge_t;
  std::vector<Edge_t> edges;
  std::vector<Filtration_value> edges_fil;
  std::map<Vertex_handle, Filtration_value> vertices;

  std::string line;
  int dim;
  Vertex_handle u, v, max_h = -1;
  Filtration_value fil;
  while (getline(in_, line)) {
    std::istringstream iss(line);
    while (iss >> dim) {
      switch (dim) {
        case 0: {
          iss >> u;
          iss >> fil;
          vertices[u] = fil;
          if (max_h < u) {
            max_h = u;
          }
          break;
        }
        case 1: {
          iss >> u;
          iss >> v;
          iss >> fil;
          edges.push_back(Edge_t(u, v));
          edges_fil.push_back(fil);
          break;
        }
        default: { break; }
      }
    }
  }
  in_.close();

  if ((size_t)(max_h + 1) != vertices.size()) {
    std::cerr << "Error: vertices must be labeled from 0 to n-1 \n";
  }

  Graph_t skel_graph(edges.begin(), edges.end(), edges_fil.begin(), vertices.size());
  auto vertex_prop = boost::get(vertex_filtration_t(), skel_graph);

  typename boost::graph_traits<Graph_t>::vertex_iterator vi, vi_end;
  auto v_it = vertices.begin();
  for (std::tie(vi, vi_end) = boost::vertices(skel_graph); vi != vi_end; ++vi, ++v_it) {
    boost::put(vertex_prop, *vi, v_it->second);
  }

  return skel_graph;
}

/**
 * @brief Read a face from a file.
 *
 * File format: 1 simplex per line<br>
 * Dim1 X11 X12 ... X1d Fil1<br>
 * Dim2 X21 X22 ... X2d Fil2<br>
 * etc<br>
 *
 * The vertices must be labeled from 0 to n-1.
 * Every simplex must appear exactly once.
 * Simplices of dimension more than 1 are ignored.
 */
template <typename Vertex_handle, typename Filtration_value>
bool read_simplex(std::istream& in_, std::vector<Vertex_handle>& simplex, Filtration_value& fil) {
  int dim = 0;
  if (!(in_ >> dim)) return false;
  Vertex_handle v;
  for (int i = 0; i < dim + 1; ++i) {
    if (!(in_ >> v)) return false;
    simplex.push_back(v);
  }
  if (!(in_ >> fil)) return false;
  in_.ignore((std::numeric_limits<std::streamsize>::max)(), '\n');  // ignore until the carriage return
  return true;
}

/**
 * @brief Read a hasse simplex from a file.
 *
 * File format: 1 simplex per line<br>
 * Dim1 k11 k12 ... k1Dim1 Fil1<br>
 * Dim2 k21 k22 ... k2Dim2 Fil2<br>
 * etc<br>
 *
 * The key of a simplex is its position in the filtration order and also the number of its row in the file.
 * Dimi ki1 ki2 ... kiDimi Fili means that the ith simplex in the filtration has dimension Dimi, filtration value
 * fil1 and simplices with key ki1 ... kiDimi in its boundary.*/
template <typename Simplex_key, typename Filtration_value>
bool read_hasse_simplex(std::istream& in_, std::vector<Simplex_key>& boundary, Filtration_value& fil) {
  int dim;
  if (!(in_ >> dim)) return false;
  if (dim == 0) {
    in_ >> fil;
    return true;
  }
  Simplex_key key;
  for (int i = 0; i < dim + 1; ++i) {
    in_ >> key;
    boundary.push_back(key);
  }
  in_ >> fil;
  return true;
}

/**
 * @brief Read a lower triangular distance matrix from a csv file. We assume that the .csv store the whole
 * (square) matrix.
 *
 * @author Pawel Dlotko
 *
 * Square matrix file format:<br>
 * 0;D12;...;D1j<br>
 * D21;0;...;D2j<br>
 * ...<br>
 * Dj1;Dj2;...;0<br>
 *
 * lower matrix file format:<br>
 * 0<br>
 * D21;<br>
 * D31;D32;<br>
 * ...<br>
 * Dj1;Dj2;...;Dj(j-1);<br>
 *
 **/
template <typename Filtration_value>
std::vector<std::vector<Filtration_value>> read_lower_triangular_matrix_from_csv_file(const std::string& filename,
                                                                                      const char separator = ';') {
#ifdef DEBUG_TRACES
  std::clog << "Using procedure read_lower_triangular_matrix_from_csv_file \n";
#endif  // DEBUG_TRACES
  std::vector<std::vector<Filtration_value>> result;
  std::ifstream in;
  in.open(filename.c_str());
  if (!in.is_open()) {
    return result;
  }

  std::string line;

  // the first line is empty, so we ignore it:
  std::getline(in, line);
  std::vector<Filtration_value> values_in_this_line;
  result.push_back(values_in_this_line);

  int number_of_line = 0;

  // first, read the file line by line to a string:
  while (std::getline(in, line)) {
    // if line is empty, break
    if (line.size() == 0) break;

    // if the last element of a string is comma:
    if (line[line.size() - 1] == separator) {
      // then shrink the string by one
      line.pop_back();
    }

    // replace all commas with spaces
    std::replace(line.begin(), line.end(), separator, ' ');

    // put the new line to a stream
    std::istringstream iss(line);
    // and now read the doubles.

    int number_of_entry = 0;
    std::vector<Filtration_value> values_in_this_line;
    while (iss.good()) {
      double entry;
      iss >> entry;
      if (number_of_entry <= number_of_line) {
        values_in_this_line.push_back(entry);
      }
      ++number_of_entry;
    }
    if (!values_in_this_line.empty()) result.push_back(values_in_this_line);
    ++number_of_line;
  }
  in.close();

#ifdef DEBUG_TRACES
  std::clog << "Here is the matrix we read : \n";
  for (size_t i = 0; i != result.size(); ++i) {
    for (size_t j = 0; j != result[i].size(); ++j) {
      std::clog << result[i][j] << " ";
    }
    std::clog << std::endl;
  }
#endif  // DEBUG_TRACES

  return result;
}  // read_lower_triangular_matrix_from_csv_file

/**
Reads a file containing persistence intervals.
Each line might contain 2, 3 or 4 values: [[field] dimension] birth death
The output iterator `out` is used this way: `*out++ = std::make_tuple(dim, birth, death);`
where `dim` is an `int`, `birth` a `double`, and `death` a `double`.
Note: the function does not check that birth <= death.
**/
template <typename OutputIterator>
void read_persistence_intervals_and_dimension(std::string const& filename, OutputIterator out) {
#ifdef DEBUG_TRACES
  std::clog << "read_persistence_intervals_and_dimension - " << filename << std::endl;
#endif  // DEBUG_TRACES
  std::ifstream in(filename);
  if (!in.is_open()) {
    std::string error_str("read_persistence_intervals_and_dimension - Unable to open file ");
    error_str.append(filename);
    std::cerr << error_str << std::endl;
    throw std::invalid_argument(error_str);
  }

  while (!in.eof()) {
    std::string line;
    getline(in, line);
    if (line.length() != 0 && line[0] != '#') {
      double numbers[4];
      int n = sscanf(line.c_str(), "%lf %lf %lf %lf", &numbers[0], &numbers[1], &numbers[2], &numbers[3]);
#ifdef DEBUG_TRACES
      std::clog << "[" << n << "] = ";
      for (int i = 0; i < n; i++) {
        std::clog << numbers[i] << ",";
      }
      std::clog << std::endl;
#endif  // DEBUG_TRACES
      if (n >= 2) {
        int dim = (n >= 3 ? static_cast<int>(numbers[n - 3]) : -1);
        *out++ = std::make_tuple(dim, numbers[n - 2], numbers[n - 1]);
      }
    }
  }
}

/**
Reads a file containing persistence intervals.
Each line might contain 2, 3 or 4 values: [[field] dimension] birth death
The return value is an `std::map<dim, std::vector<std::pair<birth, death>>>`
where `dim` is an `int`, `birth` a `double`, and `death` a `double`.
Note: the function does not check that birth <= death.
**/
inline std::map<int, std::vector<std::pair<double, double>>> read_persistence_intervals_grouped_by_dimension(
    std::string const& filename) {
  std::map<int, std::vector<std::pair<double, double>>> ret;
  read_persistence_intervals_and_dimension(
      filename, boost::make_function_output_iterator([&ret](std::tuple<int, double, double> t) {
        ret[get<0>(t)].push_back(std::make_pair(get<1>(t), get<2>(t)));
      }));
  return ret;
}

/**
Reads a file containing persistence intervals.
Each line might contain 2, 3 or 4 values: [[field] dimension] birth death
If `only_this_dim` = -1, dimension is ignored and all lines are returned.
If `only_this_dim` is >= 0, only the lines where dimension = `only_this_dim`
(or where dimension is not specified) are returned.
The return value is an `std::vector<std::pair<birth, death>>`
where `dim` is an `int`, `birth` a `double`, and `death` a `double`.
Note: the function does not check that birth <= death.
**/
inline std::vector<std::pair<double, double>> read_persistence_intervals_in_dimension(std::string const& filename,
                                                                                      int only_this_dim = -1) {
  std::vector<std::pair<double, double>> ret;
  read_persistence_intervals_and_dimension(
      filename, boost::make_function_output_iterator([only_this_dim, &ret](std::tuple<int, double, double> t) {
        if (only_this_dim == get<0>(t) || only_this_dim == -1) ret.emplace_back(get<1>(t), get<2>(t));
      }));
  return ret;
}

}  // namespace Gudhi

#endif  // READER_UTILS_H_

/*    This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014  INRIA Sophia Antipolis-Méditerranée (France)
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

#ifndef READER_UTILS_H_
#define READER_UTILS_H_

#include <gudhi/graph_simplicial_complex.h>

#include <boost/graph/adjacency_list.hpp>

#include <iostream>
#include <fstream>
#include <map>
#include <limits>  // for numeric_limits<>
#include <string>
#include <vector>

/**
 * \brief Read a set of points to turn it
 * into a vector< vector<double> > by filling points
 *
 * File format: 1 point per line
 * X11 X12 ... X1d 
 * X21 X22 ... X2d
 * etc
 */
inline void read_points(std::string file_name, std::vector< std::vector< double > > & points) {
  std::ifstream in_file(file_name.c_str(), std::ios::in);
  if (!in_file.is_open()) {
    std::cerr << "Unable to open file " << file_name << std::endl;
    return;
  }

  std::string line;
  double x;
  while (getline(in_file, line)) {
    std::vector< double > point;
    std::istringstream iss(line);
    while (iss >> x) {
      point.push_back(x);
    }
    // Check for empty lines
    if (!point.empty())
      points.push_back(point);
  }
  in_file.close();
}

/**
 * \brief Read a graph from a file.
 *
 * File format: 1 simplex per line
 * Dim1 X11 X12 ... X1d Fil1 
 * Dim2 X21 X22 ... X2d Fil2
 * etc
 *
 * The vertices must be labeled from 0 to n-1.
 * Every simplex must appear exactly once.
 * Simplices of dimension more than 1 are ignored.
 */
template< typename Graph_t, typename Edge_t, typename Filtration_value, typename Vertex_handle >
inline Graph_t read_graph(std::string file_name) {
  std::ifstream in_(file_name.c_str(), std::ios::in);
  if (!in_.is_open()) {
    std::cerr << "Unable to open file " << file_name << std::endl;
  }

  std::vector< Edge_t > edges;
  std::vector< Filtration_value > edges_fil;
  std::map< Vertex_handle, Filtration_value > vertices;

  std::string line;
  int dim;
  Vertex_handle u, v, max_h = -1;
  Filtration_value fil;
  while (getline(in_, line)) {
    std::istringstream iss(line);
    while (iss >> dim) {
      switch (dim) {
        case 0:
        {
          iss >> u;
          iss >> fil;
          vertices[u] = fil;
          if (max_h < u) {
            max_h = u;
          }
          break;
        }
        case 1:
        {
          iss >> u;
          iss >> v;
          iss >> fil;
          edges.push_back(Edge_t(u, v));
          edges_fil.push_back(fil);
          break;
        }
        default:
        {
          break;
        }
      }
    }
  }
  in_.close();

  if ((size_t) (max_h + 1) != vertices.size()) {
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
 * \brief Read a face from a file.
 *
 * File format: 1 simplex per line
 * Dim1 X11 X12 ... X1d Fil1 
 * Dim2 X21 X22 ... X2d Fil2
 * etc
 *
 * The vertices must be labeled from 0 to n-1.
 * Every simplex must appear exactly once.
 * Simplices of dimension more than 1 are ignored.
 */
template< typename Vertex_handle, typename Filtration_value >
bool read_simplex(std::istream & in_, std::vector< Vertex_handle > & simplex, Filtration_value & fil) {
  int dim = 0;
  if (!(in_ >> dim)) return false;
  Vertex_handle v;
  for (int i = 0; i < dim + 1; ++i) {
    in_ >> v;
    simplex.push_back(v);
  }
  in_ >> fil;
  in_.ignore((std::numeric_limits<std::streamsize>::max)(), '\n');  // ignore until the carriage return
  return true;
}

/**
 * \brief Read a hasse simplex from a file.
 *
 * File format: 1 simplex per line
 * Dim1 k11 k12 ... k1Dim1 Fil1 
 * Dim2 k21 k22 ... k2Dim2 Fil2
 * etc
 *
 * The key of a simplex is its position in the filtration order
 * and also the number of its row in the file.
 * Dimi ki1 ki2 ... kiDimi Fili means that the ith simplex in the 
 * filtration has dimension Dimi, filtration value fil1 and simplices with 
 * key ki1 ... kiDimi in its boundary.*/
template< typename Simplex_key, typename Filtration_value >
bool read_hasse_simplex(std::istream & in_, std::vector< Simplex_key > & boundary, Filtration_value & fil) {
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

#endif  // READER_UTILS_H_

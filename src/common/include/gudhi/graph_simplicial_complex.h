/*    This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014  INRIA
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

#ifndef GRAPH_SIMPLICIAL_COMPLEX_H_
#define GRAPH_SIMPLICIAL_COMPLEX_H_

#include <boost/graph/adjacency_list.hpp>

#include <utility>  // for pair<>
#include <vector>
#include <map>

/* Edge tag for Boost PropertyGraph. */
struct edge_filtration_t {
  typedef boost::edge_property_tag kind;
};

/* Vertex tag for Boost PropertyGraph. */
struct vertex_filtration_t {
  typedef boost::vertex_property_tag kind;
};

#endif  // GRAPH_SIMPLICIAL_COMPLEX_H_

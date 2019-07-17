/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef IO_BUILD_MESH_FROM_CELL_COMPLEX_H_
#define IO_BUILD_MESH_FROM_CELL_COMPLEX_H_

#include <gudhi/IO/Mesh_medit.h>

struct Configuration {
  bool toggle_edges = true,
    toggle_triangles = true,
    toggle_tetrahedra = true;
  std::size_t ref_edges = 1,
    ref_triangles = 1,
    ref_tetrahedra = 1;
};


template <class Cell_complex>
Mesh_medit build_primal_mesh(const Cell_complex& cell_complex,
			     Configuration configuration = Configuration()) {
  using Hasse_cell = typename Cell_complex::Hasse_cell;
  Mesh_medit output;
  std::map<Hasse_cell*, std::size_t> ci_map;
  std::size_t index = 1; // current size of output.vertex_points

  for (auto cp_pair: cell_complex.cell_point_map()) {
    ci_map.emplace(std::make_pair(cp_pair.first, index++));
    output.vertex_points.push_back(cp_pair.second);
  }

  return output;
}

#endif

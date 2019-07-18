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

namespace Gudhi {

namespace coxeter_triangulation {

struct Configuration {
  bool toggle_edges = true,
    toggle_triangles = true,
    toggle_tetrahedra = true;
  std::size_t ref_edges = 1,
    ref_triangles = 1,
    ref_tetrahedra = 1;
};


template <class Cell_complex>
Mesh_medit build_mesh_from_cell_complex(const Cell_complex& cell_complex,
					Configuration configuration = Configuration()) {
  using Hasse_cell = typename Cell_complex::Hasse_cell;
  using Mesh_element_vertices = Mesh_medit::Mesh_elements::value_type::first_type;
  Mesh_medit output;
  std::map<Hasse_cell*, std::size_t> vi_map, ci_map; // one for vertices, other for 2d-cells
  std::vector<Hasse_cell*> edge_cells, polygon_cells, polytope_cells;
  std::size_t index = 1; // current size of output.vertex_points

  if (cell_complex.cell_point_map().empty())
    return output;
  std::size_t amb_d = std::min((int) cell_complex.cell_point_map().begin()->second.size(), 3);
  
  for (const auto& cp_pair: cell_complex.cell_point_map()) {
    vi_map.emplace(std::make_pair(cp_pair.first, index++));
    output.vertex_points.push_back(cp_pair.second);
  }
  if (cell_complex.intrinsic_dimension() >= 2) 
    for (const auto& sc_pair: cell_complex.simplex_cell_map(2)) {
      Eigen::VectorXd barycenter = Eigen::VectorXd::Zero(amb_d);
      std::set<std::size_t> vertex_indices;
      Hasse_cell* cell = sc_pair.second;
      for (const auto& ei_pair: cell->get_boundary())
	for (const auto& vi_pair: ei_pair.first->get_boundary())
	  vertex_indices.emplace(vi_map[vi_pair.first]);
      for (const std::size_t& v: vertex_indices)
	barycenter += output.vertex_points[v-1];
      ci_map.emplace(std::make_pair(cell, index++));
      output.vertex_points.emplace_back((1./vertex_indices.size()) * barycenter);    
    }

  if (configuration.toggle_edges && cell_complex.intrinsic_dimension() >= 1)
    for (const auto& sc_map: cell_complex.simplex_cell_map(1)) {
      Hasse_cell* edge_cell = sc_map.second;
      Mesh_element_vertices edge;
      for (const auto& vi_pair: edge_cell->get_boundary())
	edge.push_back(vi_map[vi_pair.first]);
      output.edges.emplace_back(std::make_pair(edge, configuration.ref_edges));
    }
  
  if (configuration.toggle_triangles && cell_complex.intrinsic_dimension() >= 2)
    for (const auto& sc_pair: cell_complex.simplex_cell_map(2))
      for (const auto& ei_pair: sc_pair.second->get_boundary()) {
	Mesh_element_vertices triangle(1, ci_map[sc_pair.second]);
	for (const auto& vi_pair: ei_pair.first->get_boundary())
	  triangle.push_back(vi_map[vi_pair.first]);
	output.triangles.emplace_back(std::make_pair(triangle, configuration.ref_triangles));
      }
  
  if (configuration.toggle_tetrahedra && cell_complex.intrinsic_dimension() >= 3)
    for (const auto& sc_pair: cell_complex.simplex_cell_map(3)) {
      Eigen::VectorXd barycenter = Eigen::VectorXd::Zero(amb_d);
      std::set<std::size_t> vertex_indices;
      Hasse_cell* cell = sc_pair.second;
      for (const auto& ci_pair: cell->get_boundary())
	for (const auto& ei_pair: ci_pair.first->get_boundary())
	  for (const auto& vi_pair: ei_pair.first->get_boundary())
	    vertex_indices.emplace(vi_map[vi_pair.first]);
      for (const std::size_t& v: vertex_indices)
	barycenter += output.vertex_points[v-1];
      output.vertex_points.emplace_back((1./vertex_indices.size()) * barycenter);

      for (const auto& ci_pair: cell->get_boundary())
	for (const auto& ei_pair: ci_pair.first->get_boundary()) {
	  Mesh_element_vertices tetrahedron = {index, ci_map[sc_pair.second]};
	  for (const auto& vi_pair: ei_pair.first->get_boundary())
	    tetrahedron.push_back(vi_map[vi_pair.first]);
	  output.tetrahedra.emplace_back(std::make_pair(tetrahedron, configuration.ref_tetrahedra));
	}
      index++;
    }
      
  return output;
}

} // namespace coxeter_triangulation 

} // namespace Gudhi

#endif

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


template <class Hasse_cell,
	  class Simplex_cell_map>
void populate_mesh(Mesh_medit& output,
		   Simplex_cell_map& sc_map,
		   Configuration configuration,
		   std::size_t amb_d,
		   std::map<Hasse_cell*, std::size_t> vi_map) {
  using Mesh_element_vertices = Mesh_medit::Mesh_elements::value_type::first_type;
  std::map<Hasse_cell*, std::size_t> ci_map;
  std::size_t index = vi_map.size() + 1; // current size of output.vertex_points
  if (sc_map.size() >= 3) 
    for (const auto& sc_pair: sc_map[2]) {
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
#ifdef GUDHI_COX_OUTPUT_TO_HTML
      std::string vlist;
      for (const std::size_t& v: vertex_indices)
	vlist += " " + std::to_string(v);
      cell_vlist_map.emplace(std::make_pair(to_string(cell), vlist));
#endif
    }

  if (configuration.toggle_edges && sc_map.size() >= 2)
    for (const auto& sc_map: sc_map[1]) {
      Hasse_cell* edge_cell = sc_map.second;
      Mesh_element_vertices edge;
      for (const auto& vi_pair: edge_cell->get_boundary())
	edge.push_back(vi_map[vi_pair.first]);
      output.edges.emplace_back(std::make_pair(edge, configuration.ref_edges));
#ifdef GUDHI_COX_OUTPUT_TO_HTML
      std::string vlist;
      for (const std::size_t& v: edge)
	vlist += " " + std::to_string(v);
      cell_vlist_map.emplace(std::make_pair(to_string(edge_cell), vlist));
#endif
    }
  
  if (configuration.toggle_triangles && sc_map.size() >= 3)
    for (const auto& sc_pair: sc_map[2]) {
      for (const auto& ei_pair: sc_pair.second->get_boundary()) {
	Mesh_element_vertices triangle(1, ci_map[sc_pair.second]);
	for (const auto& vi_pair: ei_pair.first->get_boundary())
	  triangle.push_back(vi_map[vi_pair.first]);
	output.triangles.emplace_back(std::make_pair(triangle, configuration.ref_triangles));
      }
    }
  
  if (configuration.toggle_tetrahedra && sc_map.size() >= 4)
    for (const auto& sc_pair: sc_map[3]) {
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
#ifdef GUDHI_COX_OUTPUT_TO_HTML
      std::string vlist;
      for (const std::size_t& v: vertex_indices)
	vlist += " " + std::to_string(v);
      cell_vlist_map.emplace(std::make_pair(to_string(cell), vlist));
#endif

      for (const auto& ci_pair: cell->get_boundary())
	for (const auto& ei_pair: ci_pair.first->get_boundary()) {
	  Mesh_element_vertices tetrahedron = {index, ci_map[sc_pair.second]};
	  for (const auto& vi_pair: ei_pair.first->get_boundary())
	    tetrahedron.push_back(vi_map[vi_pair.first]);
	  output.tetrahedra.emplace_back(std::make_pair(tetrahedron, configuration.ref_tetrahedra));
	}
      index++;
    }  
}

template <class Cell_complex>
Mesh_medit build_mesh_from_cell_complex(const Cell_complex& cell_complex,
					Configuration i_configuration = Configuration(),
					Configuration b_configuration = Configuration()) {
  using Hasse_cell = typename Cell_complex::Hasse_cell;
  Mesh_medit output;
  std::map<Hasse_cell*, std::size_t> vi_map; // one for vertices, other for 2d-cells
  std::size_t index = 1; // current size of output.vertex_points

  if (cell_complex.cell_point_map().empty())
    return output;
  std::size_t amb_d = std::min((int) cell_complex.cell_point_map().begin()->second.size(), 3);
  
  for (const auto& cp_pair: cell_complex.cell_point_map()) {
#ifdef GUDHI_COX_OUTPUT_TO_HTML
    std::string vlist;
    vlist += " " + std::to_string(index);
    cell_vlist_map.emplace(std::make_pair(to_string(cp_pair.first), vlist));
#endif
    vi_map.emplace(std::make_pair(cp_pair.first, index++));
    output.vertex_points.push_back(cp_pair.second);
    output.vertex_points.back().conservativeResize(amb_d);
  }
  

  populate_mesh(output, cell_complex.interior_simplex_cell_maps(), i_configuration, amb_d, vi_map);
#ifdef GUDHI_COX_OUTPUT_TO_HTML
  for (const auto& sc_map: cell_complex.interior_simplex_cell_maps())    
    for (const auto& sc_pair: sc_map) {
      std::string simplex = "I" + to_string(sc_pair.first);
      std::string cell = to_string(sc_pair.second);
      std::string vlist = cell_vlist_map.at(cell).substr(1);
      simplex_vlist_map.emplace(std::make_pair(simplex, vlist));
    }
#endif  
  populate_mesh(output, cell_complex.boundary_simplex_cell_maps(), b_configuration, amb_d, vi_map);  
#ifdef GUDHI_COX_OUTPUT_TO_HTML
  for (const auto& sc_map: cell_complex.boundary_simplex_cell_maps())    
    for (const auto& sc_pair: sc_map) {
      std::string simplex = "B" + to_string(sc_pair.first);
      std::string cell = to_string(sc_pair.second);
      std::string vlist = cell_vlist_map.at(cell).substr(1);
      simplex_vlist_map.emplace(std::make_pair(simplex, vlist));
    }
#endif  
  return output;
}

} // namespace coxeter_triangulation 

} // namespace Gudhi

#endif

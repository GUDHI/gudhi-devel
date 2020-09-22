/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef IO_OUTPUT_MESHES_TO_MEDIT_H_
#define IO_OUTPUT_MESHES_TO_MEDIT_H_

#include <gudhi/IO/Mesh_medit.h>

#include <Eigen/Dense>

#include <cstdlib>  // for std::size_t
#include <fstream>  // for std::ofstream
#include <vector>
#include <type_traits>  // for std::enable_if
#include <tuple>  // for std::get
#include <utility>  // for std::make_pair

namespace Gudhi {

namespace coxeter_triangulation {

using Vertex_points = Mesh_medit::Vertex_points;
using Mesh_elements = Mesh_medit::Mesh_elements;
using Scalar_field_range = Mesh_medit::Scalar_field_range;

template <std::size_t I = 0,
	  typename... Meshes>
typename std::enable_if<I == sizeof... (Meshes), void>::type
fill_meshes (Vertex_points& vertex_points,
	     Mesh_elements& edges,
	     Mesh_elements& triangles,
	     Mesh_elements& tetrahedra,
	     Scalar_field_range& triangles_scalar_range,
	     Scalar_field_range& tetrahedra_scalar_range,
	     std::size_t index,
	     const Meshes&... meshes) {
}

template <std::size_t I = 0,
	  typename... Meshes>
typename std::enable_if<I != sizeof... (Meshes), void>::type
fill_meshes (Vertex_points& vertex_points,
	     Mesh_elements& edges,
	     Mesh_elements& triangles,
	     Mesh_elements& tetrahedra,
	     Scalar_field_range& triangles_scalar_range,
	     Scalar_field_range& tetrahedra_scalar_range,
	     std::size_t index,
	     const Meshes&... meshes) {
  auto mesh = std::get<I>(std::forward_as_tuple(meshes...));
  for (const auto& v: mesh.vertex_points)
    vertex_points.push_back(v);
  for (const auto& e: mesh.edges) {
    std::vector<std::size_t> edge;
    for (const auto& v_i: e.first)
      edge.push_back(v_i + index);
    edges.emplace_back(std::make_pair(edge, e.second));
  }
  for (const auto& t: mesh.triangles) {
    std::vector<std::size_t> triangle;
    for (const auto& v_i: t.first)
      triangle.push_back(v_i + index);
    triangles.emplace_back(std::make_pair(triangle, t.second));
  }
  for (const auto& t: mesh.tetrahedra) {
    std::vector<std::size_t> tetrahedron;
    for (const auto& v_i: t.first)
      tetrahedron.push_back(v_i + index);
    tetrahedra.emplace_back(std::make_pair(tetrahedron, t.second));
  }
  for (const auto& b: mesh.triangles_scalar_range)
    triangles_scalar_range.push_back(b);
  for (const auto& b: mesh.tetrahedra_scalar_range)
    tetrahedra_scalar_range.push_back(b);
  fill_meshes<I+1, Meshes...>(vertex_points,
			      edges,
			      triangles,
			      tetrahedra,
			      triangles_scalar_range,
			      tetrahedra_scalar_range,
			      index + mesh.vertex_points.size(),
			      meshes...);
}

/** \brief Outputs a text file with specified meshes that can be visualized in Medit.
 *  
 *  @param[in] amb_d Ambient dimension. Can be 2 or 3.
 *  @param[in] file_name The name of the output file.
 *  @param[in] meshes A pack of meshes to be specified separated by commas.
 */
template <typename... Meshes>
void output_meshes_to_medit(std::size_t amb_d, std::string file_name, const Meshes&... meshes) {
  Vertex_points vertex_points;
  Mesh_elements edges, triangles, tetrahedra;
  Scalar_field_range triangles_scalar_range, tetrahedra_scalar_range;
  fill_meshes(vertex_points,
	      edges,
	      triangles,
	      tetrahedra,
	      triangles_scalar_range,
	      tetrahedra_scalar_range,
	      0,
	      meshes...);
  
  std::ofstream ofs (file_name + ".mesh", std::ofstream::out);
  std::ofstream ofs_bb (file_name + ".bb", std::ofstream::out);
  
  if (amb_d == 2) {
    ofs << "MeshVersionFormatted 1\nDimension 2\n";
    ofs_bb << "2 1 ";
    ofs << "Vertices\n" << vertex_points.size() << "\n";
    for (auto p: vertex_points) {
      ofs << p[0] << " " << p[1] << " 2\n";
    }
    ofs << "Edges " << edges.size() << "\n";
    for (auto e: edges) {
      for (auto v: e.first)      
    	ofs << v << " ";
      ofs << e.second << std::endl;
    }
    ofs << "Triangles " << triangles.size() << "\n";
    for (auto s: triangles) {
      for (auto v: s.first) {
    	ofs << v << " ";
      }
      ofs << s.second << std::endl;
    }

    ofs_bb << triangles_scalar_range.size() << " 1\n";
    for (auto& b: triangles_scalar_range)
      ofs_bb << b << "\n";

  }
  else {
    ofs << "MeshVersionFormatted 1\nDimension 3\n";
    ofs_bb << "3 1 ";
    ofs << "Vertices\n" << vertex_points.size() << "\n";
    for (auto p: vertex_points) {
      ofs << p[0] << " " << p[1] << " " << p[2] << " 2\n";
    }
    ofs << "Edges " << edges.size() << "\n";
    for (auto e: edges) {
      for (auto v: e.first)      
    	ofs << v << " ";
      ofs << e.second << std::endl;
    }
    ofs << "Triangles " << triangles.size() << "\n";
    for (auto s: triangles) {
      for (auto v: s.first) {
    	ofs << v << " ";
      }
      ofs << s.second << std::endl;
    }
    ofs << "Tetrahedra " << tetrahedra.size() << "\n";
    for (auto s: tetrahedra) {
      for (auto v: s.first) {
    	ofs << v << " ";
      }
      ofs << s.second << std::endl;
    }

    ofs_bb << triangles_scalar_range.size() + tetrahedra_scalar_range.size() << " 1\n";
    for (auto& b: triangles_scalar_range)
      ofs_bb << b << "\n";
    for (auto& b: tetrahedra_scalar_range)
      ofs_bb << b << "\n";
  }

  
  ofs.close();
  ofs_bb.close();

}

} // namespace coxeter_triangulation 

} // namespace Gudhi

#endif

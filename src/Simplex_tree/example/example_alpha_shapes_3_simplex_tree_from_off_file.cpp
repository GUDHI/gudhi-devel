/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Simplex_tree.h>
#include <gudhi/Points_3D_off_io.h>

#include <boost/variant.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Alpha_shape_vertex_base_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#include <CGAL/iterator.h>

#include <fstream>
#include <cmath>
#include <string>
#include <tuple>  // for tuple<>
#include <map>
#include <utility>  // for pair<>
#include <list>
#include <vector>

// Alpha_shape_3 templates type definitions
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Alpha_shape_vertex_base_3<Kernel> Vb;
typedef CGAL::Alpha_shape_cell_base_3<Kernel> Fb;
typedef CGAL::Triangulation_data_structure_3<Vb, Fb> Tds;
typedef CGAL::Delaunay_triangulation_3<Kernel, Tds> Triangulation_3;
typedef CGAL::Alpha_shape_3<Triangulation_3> Alpha_shape_3;

// From file type definition
typedef Kernel::Point_3 Point;

// filtration with alpha values needed type definition
typedef Alpha_shape_3::FT Alpha_value_type;
typedef CGAL::Object Object;
typedef CGAL::Dispatch_output_iterator<
CGAL::cpp11::tuple<Object, Alpha_value_type>,
CGAL::cpp11::tuple<std::back_insert_iterator< std::vector<Object> >,
    std::back_insert_iterator< std::vector<Alpha_value_type> > > > Dispatch;
typedef Alpha_shape_3::Cell_handle Cell_handle;
typedef Alpha_shape_3::Facet Facet;
typedef Alpha_shape_3::Edge Edge;
typedef std::list<Alpha_shape_3::Vertex_handle> Vertex_list;

// gudhi type definition
typedef Gudhi::Simplex_tree<> Simplex_tree;
typedef Simplex_tree::Vertex_handle Simplex_tree_vertex;
typedef std::map<Alpha_shape_3::Vertex_handle, Simplex_tree_vertex > Alpha_shape_simplex_tree_map;
typedef std::pair<Alpha_shape_3::Vertex_handle, Simplex_tree_vertex> Alpha_shape_simplex_tree_pair;
typedef std::vector< Simplex_tree_vertex > Simplex_tree_vector_vertex;

Vertex_list from(const Cell_handle& ch) {
  Vertex_list the_list;
  for (auto i = 0; i < 4; i++) {
#ifdef DEBUG_TRACES
    std::clog << "from cell[" << i << "]=" << ch->vertex(i)->point() << std::endl;
#endif  // DEBUG_TRACES
    the_list.push_back(ch->vertex(i));
  }
  return the_list;
}

Vertex_list from(const Facet& fct) {
  Vertex_list the_list;
  for (auto i = 0; i < 4; i++) {
    if (fct.second != i) {
#ifdef DEBUG_TRACES
      std::clog << "from facet=[" << i << "]" << fct.first->vertex(i)->point() << std::endl;
#endif  // DEBUG_TRACES
      the_list.push_back(fct.first->vertex(i));
    }
  }
  return the_list;
}

Vertex_list from(const Edge& edg) {
  Vertex_list the_list;
  for (auto i = 0; i < 4; i++) {
    if ((edg.second == i) || (edg.third == i)) {
#ifdef DEBUG_TRACES
      std::clog << "from edge[" << i << "]=" << edg.first->vertex(i)->point() << std::endl;
#endif  // DEBUG_TRACES
      the_list.push_back(edg.first->vertex(i));
    }
  }
  return the_list;
}

Vertex_list from(const Alpha_shape_3::Vertex_handle& vh) {
  Vertex_list the_list;
#ifdef DEBUG_TRACES
  std::clog << "from vertex=" << vh->point() << std::endl;
#endif  // DEBUG_TRACES
  the_list.push_back(vh);
  return the_list;
}

int main(int argc, char * const argv[]) {
  // program args management
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0]
        << " path_to_off_file \n";
    return 0;
  }

  // Read points from file
  std::string offInputFile(argv[1]);
  // Read the OFF file (input file name given as parameter) and triangulate points
  Gudhi::Points_3D_off_reader<Point> off_reader(offInputFile);
  // Check the read operation was correct
  if (!off_reader.is_valid()) {
    std::cerr << "Unable to read file " << argv[1] << std::endl;
    return 0;
  }
  // Retrieve the triangulation
  std::vector<Point> lp = off_reader.get_point_cloud();

  // alpha shape construction from points. CGAL has a strange behavior in REGULARIZED mode.
  Alpha_shape_3 as(lp.begin(), lp.end(), 0, Alpha_shape_3::GENERAL);
#ifdef DEBUG_TRACES
  std::clog << "Alpha shape computed in GENERAL mode" << std::endl;
#endif  // DEBUG_TRACES

  // filtration with alpha values from alpha shape
  std::vector<Object> the_objects;
  std::vector<Alpha_value_type> the_alpha_values;

  Dispatch disp = CGAL::dispatch_output<Object, Alpha_value_type>(std::back_inserter(the_objects),
      std::back_inserter(the_alpha_values));

  as.filtration_with_alpha_values(disp);
#ifdef DEBUG_TRACES
  std::clog << "filtration_with_alpha_values returns : " << the_objects.size() << " objects" << std::endl;
#endif  // DEBUG_TRACES

  Alpha_shape_3::size_type count_vertices = 0;
  Alpha_shape_3::size_type count_edges = 0;
  Alpha_shape_3::size_type count_facets = 0;
  Alpha_shape_3::size_type count_cells = 0;

  // Loop on objects vector
  Vertex_list vertex_list;
  Simplex_tree simplex_tree;
  Alpha_shape_simplex_tree_map map_cgal_simplex_tree;
  std::vector<Alpha_value_type>::iterator the_alpha_value_iterator = the_alpha_values.begin();
  for (auto object_iterator : the_objects) {
    // Retrieve Alpha shape vertex list from object
    if (const Cell_handle * cell = CGAL::object_cast<Cell_handle>(&object_iterator)) {
      vertex_list = from(*cell);
      count_cells++;
    } else if (const Facet * facet = CGAL::object_cast<Facet>(&object_iterator)) {
      vertex_list = from(*facet);
      count_facets++;
    } else if (const Edge * edge = CGAL::object_cast<Edge>(&object_iterator)) {
      vertex_list = from(*edge);
      count_edges++;
    } else if (const Alpha_shape_3::Vertex_handle * vertex =
              CGAL::object_cast<Alpha_shape_3::Vertex_handle>(&object_iterator)) {
      count_vertices++;
      vertex_list = from(*vertex);
    }
    // Construction of the vector of simplex_tree vertex from list of alpha_shapes vertex
    Simplex_tree_vector_vertex the_simplex_tree;
    for (auto the_alpha_shape_vertex : vertex_list) {
      Alpha_shape_simplex_tree_map::iterator the_map_iterator = map_cgal_simplex_tree.find(the_alpha_shape_vertex);
      if (the_map_iterator == map_cgal_simplex_tree.end()) {
        // alpha shape not found
        Simplex_tree_vertex vertex = map_cgal_simplex_tree.size();
#ifdef DEBUG_TRACES
        std::clog << "vertex [" << the_alpha_shape_vertex->point() << "] not found - insert_simplex " << vertex << "\n";
#endif  // DEBUG_TRACES
        the_simplex_tree.push_back(vertex);
        map_cgal_simplex_tree.insert(Alpha_shape_simplex_tree_pair(the_alpha_shape_vertex, vertex));
      } else {
        // alpha shape found
        Simplex_tree_vertex vertex = the_map_iterator->second;
#ifdef DEBUG_TRACES
        std::clog << "vertex [" << the_alpha_shape_vertex->point() << "] found in " << vertex << std::endl;
#endif  // DEBUG_TRACES
        the_simplex_tree.push_back(vertex);
      }
    }
    // Construction of the simplex_tree
#ifdef DEBUG_TRACES
    std::clog << "filtration = " << *the_alpha_value_iterator << std::endl;
#endif  // DEBUG_TRACES
    simplex_tree.insert_simplex(the_simplex_tree, std::sqrt(*the_alpha_value_iterator));
    if (the_alpha_value_iterator != the_alpha_values.end())
      ++the_alpha_value_iterator;
    else
      std::cerr << "This shall not happen" << std::endl;
  }
#ifdef DEBUG_TRACES
  std::clog << "vertices \t\t" << count_vertices << std::endl;
  std::clog << "edges \t\t" << count_edges << std::endl;
  std::clog << "facets \t\t" << count_facets << std::endl;
  std::clog << "cells \t\t" << count_cells << std::endl;


  std::clog << "Information of the Simplex Tree:\n";
  std::clog << "  Number of vertices = " << simplex_tree.num_vertices() << " ";
  std::clog << "  Number of simplices = " << simplex_tree.num_simplices() << std::endl << std::endl;
#endif  // DEBUG_TRACES

#ifdef DEBUG_TRACES
  std::clog << "Iterator on vertices: \n";
  for (auto vertex : simplex_tree.complex_vertex_range()) {
    std::clog << vertex << " ";
  }
#endif  // DEBUG_TRACES

  std::clog << simplex_tree << std::endl;

#ifdef DEBUG_TRACES
  std::clog << std::endl << std::endl << "Iterator on simplices:\n";
  for (auto simplex : simplex_tree.complex_simplex_range()) {
    std::clog << "   ";
    for (auto vertex : simplex_tree.simplex_vertex_range(simplex)) {
      std::clog << vertex << " ";
    }
    std::clog << std::endl;
  }
#endif  // DEBUG_TRACES
#ifdef DEBUG_TRACES
  std::clog << std::endl << std::endl << "Iterator on Simplices in the filtration, with [filtration value]:\n";
  for (auto f_simplex : simplex_tree.filtration_simplex_range()) {
    std::clog << "   " << "[" << simplex_tree.filtration(f_simplex) << "] ";
    for (auto vertex : simplex_tree.simplex_vertex_range(f_simplex)) {
      std::clog << vertex << " ";
    }
    std::clog << std::endl;
  }
#endif  // DEBUG_TRACES
#ifdef DEBUG_TRACES
  std::clog << std::endl << std::endl << "Iterator on Simplices in the filtration, and their boundary simplices:\n";
  for (auto f_simplex : simplex_tree.filtration_simplex_range()) {
    std::clog << "   " << "[" << simplex_tree.filtration(f_simplex) << "] ";
    for (auto vertex : simplex_tree.simplex_vertex_range(f_simplex)) {
      std::clog << vertex << " ";
    }
    std::clog << std::endl;

    for (auto b_simplex : simplex_tree.boundary_simplex_range(f_simplex)) {
      std::clog << "      " << "[" << simplex_tree.filtration(b_simplex) << "] ";
      for (auto vertex : simplex_tree.simplex_vertex_range(b_simplex)) {
        std::clog << vertex << " ";
      }
      std::clog << std::endl;
    }
  }
#endif  // DEBUG_TRACES

  return 0;
}

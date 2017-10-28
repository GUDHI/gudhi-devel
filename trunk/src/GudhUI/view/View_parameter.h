/* This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014  INRIA Sophia Antipolis-Mediterranee (France)
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
 * 
 */

#ifndef VIEW_VIEW_PARAMETER_H_
#define VIEW_VIEW_PARAMETER_H_

#include <iostream>

/**
 * Different parameters for the view such as the camera angle,
 * the light, options for vertices/edges/triangles.
 */
class View_parameter {
 public:
  bool light;
  bool relative_light;

  double size_vertices;
  double size_edges;
  double light_edges;  // in 0-1
  double light_triangles;  // in 0-1

  /**
   * light angle
   */
  double theta;
  double phi;

  enum VERTEX_MODE {
    V_NONE, V_SIMPLE, V_COUNT
  };

  enum EDGE_MODE {
    E_NONE, E_SIMPLE, E_COUNT
  };

  enum TRIANGLE_MODE {
    T_NONE, T_SIMPLE, T_COUNT
  };

  VERTEX_MODE vertex_mode;
  EDGE_MODE edge_mode;
  TRIANGLE_MODE triangle_mode;

  void change_vertex_mode() {
    int current_value = vertex_mode;
    vertex_mode = static_cast<VERTEX_MODE> (++current_value % V_COUNT);
    std::cout << "Vertex mode : ";
    switch (vertex_mode) {
      case V_NONE:
        std::cout << "empty\n";
        break;
      case V_SIMPLE:
        std::cout << "simple\n";
        break;
      default:
        break;
    }
  }

  void change_vertex_mode(int new_mode) {
    vertex_mode = static_cast<VERTEX_MODE> (new_mode % V_COUNT);
  }

  void change_edge_mode() {
    int current_value = edge_mode;
    edge_mode = static_cast<EDGE_MODE> (++current_value % E_COUNT);
  }

  void change_edge_mode(int new_mode) {
    edge_mode = static_cast<EDGE_MODE> (new_mode % E_COUNT);
  }

  void change_triangle_mode() {
    int current_value = triangle_mode;
    triangle_mode = static_cast<TRIANGLE_MODE> (++current_value % T_COUNT);
  }

  View_parameter() {
    light = true;
    relative_light = true;
    vertex_mode = V_SIMPLE;
    edge_mode = E_SIMPLE;
    triangle_mode = T_NONE;

    size_vertices = 3;
    size_edges = 2;

    light_edges = 0.3;
    light_triangles = 0.85;
    theta = 0;
    phi = 0;
  }

  friend std::ostream& operator<<(std::ostream& stream, const View_parameter& param) {
    stream << param.light << " ";
    stream << param.relative_light << " ";
    stream << param.vertex_mode << " ";
    stream << param.edge_mode << " ";
    stream << param.triangle_mode << " ";
    stream << param.size_vertices << " ";
    stream << param.size_edges << " ";
    stream << param.light_edges << " ";
    stream << param.light_triangles << " ";
    stream << param.theta << " ";
    stream << param.phi << " ";
    return stream;
  }

  friend std::istream& operator>>(std::istream& stream, View_parameter& param) {
    stream >> param.light;
    stream >> param.relative_light;
    int a;
    stream >> a;
    param.vertex_mode = static_cast<VERTEX_MODE> (a % V_COUNT);
    stream >> a;
    param.edge_mode = static_cast<EDGE_MODE> (a % E_COUNT);
    stream >> a;
    param.triangle_mode = static_cast<TRIANGLE_MODE> (a % T_COUNT);
    stream >> a;
    stream >> param.size_vertices;
    stream >> param.size_edges;
    stream >> param.light_edges;
    stream >> param.light_triangles;
    stream >> param.theta;
    stream >> param.phi;
    return stream;
  }
};

#endif  // VIEW_VIEW_PARAMETER_H_

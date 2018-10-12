/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2014 Inria
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

#ifndef ALPHA_COMPLEX_3D_HELPER_H_
#define ALPHA_COMPLEX_3D_HELPER_H_

template <class Vertex_list, class Cell_handle>
Vertex_list from_cell(const Cell_handle& ch) {
  Vertex_list the_list;
  for (auto i = 0; i < 4; i++) {
#ifdef DEBUG_TRACES
    std::cout << "from cell[" << i << "]=" << ch->vertex(i)->point() << std::endl;
#endif  // DEBUG_TRACES
    the_list.push_back(ch->vertex(i));
  }
  return the_list;
}

template <class Vertex_list, class Facet>
Vertex_list from_facet(const Facet& fct) {
  Vertex_list the_list;
  for (auto i = 0; i < 4; i++) {
    if (fct.second != i) {
#ifdef DEBUG_TRACES
      std::cout << "from facet=[" << i << "]" << fct.first->vertex(i)->point() << std::endl;
#endif  // DEBUG_TRACES
      the_list.push_back(fct.first->vertex(i));
    }
  }
  return the_list;
}

template <class Vertex_list, class Edge_3>
Vertex_list from_edge(const Edge_3& edg) {
  Vertex_list the_list;
  for (auto i : {edg.second, edg.third}) {
#ifdef DEBUG_TRACES
    std::cout << "from edge[" << i << "]=" << edg.first->vertex(i)->point() << std::endl;
#endif  // DEBUG_TRACES
    the_list.push_back(edg.first->vertex(i));
  }
  return the_list;
}

template <class Vertex_list, class Vertex_handle>
Vertex_list from_vertex(const Vertex_handle& vh) {
  Vertex_list the_list;
#ifdef DEBUG_TRACES
  std::cout << "from vertex=" << vh->point() << std::endl;
#endif  // DEBUG_TRACES
  the_list.push_back(vh);
  return the_list;
}

#endif  // ALPHA_COMPLEX_3D_HELPER_H_

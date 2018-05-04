/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2015  INRIA Sophia Antipolis-Méditerranée (France)
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

#ifndef COLLAPSE_H_
#define COLLAPSE_H_

#include <boost/graph/adjacency_list.hpp>

/* NOTE: The elementary collapse is lazy. The coface is not removed physically, but all its adjacencies are
 * removed. Hence, the collapsed coface does not feature anywhere else in the algorithm.
 */
template <class Map,
          class Inv_map,
          class Graph,
          class CollapseTraits>
void elementary_collapse(typename Map::iterator facet_it,
                         Map* faces,
                         Inv_map* inv_faces,
                         Map* cofaces,
                         Inv_map* inv_cofaces,
                         Graph& face_coface_graph, 
                         const CollapseTraits& collapse_traits) {
  using Cell_comparison = typename CollapseTraits::Cell_comparison; 
  using Graph_v = typename Graph::vertex_descriptor;
  
  // check if there is one coface
  if (boost::in_degree(facet_it->second, face_coface_graph) == 1)
    return;

  // check if it is valid wrt traits
  Graph_v facet_v = facet_it->second;
  Graph_v cofacet_v = boost::source(*boost::in_edges(facet_v, face_coface_graph).first,
                                    face_coface_graph);
  typename Inv_map::iterator inv_cofacet_it = inv_cofaces->find(cofacet_v);
  typename Map::iterator cofacet_it = inv_cofacet_it->second; 
  if (!collapse_traits.is_valid_collapse(facet_it,
                                         cofacet_it))
    return;

  typename Graph::out_edge_iterator out_edge_it, out_edge_end;
  std::tie(out_edge_it, out_edge_end) = boost::out_edges(cofacet_it->second, face_coface_graph);
  while (out_edge_it != out_edge_end) {
    Graph_v another_facet_v = boost::target(*out_edge_it, face_coface_graph);
    boost::remove_edge(*out_edge_it++, face_coface_graph);
    auto another_facet_it = inv_faces.find(another_facet_v);
    if (another_facet_v != facet_it->second &&
        Cell_comparison()(another_facet_it->first, facet_it->first))
      elementary_collapse(another_facet_it, faces, inv_faces, cofaces, inv_cofaces, collapse_traits);
  }

  inv_cofaces->erase(inv_cofaces->find(facet_v));
  cofaces->erase(facet_it);
  boost::remove_vertex(facet_v, face_coface_graph);
  
  inv_cofaces->erase(inv_cofacet_it);
  cofaces->erase(cofacet_it);
  boost::remove_vertex(cofacet_v, face_coface_graph);
  
  // typedef typename Simplex_with_cofaces::List_of_facets List_of_facets;
  // Map_iterator coface_it = facet_it->second.first.coface().first;
  // // std::cout << "EColapse: facet = " << facet_it->first << ", coface = " << coface_it->first << "\n";
  // typename List_of_facets::iterator col_it = facet_it->second.first.coface().second;
  // List_of_facets& facet_list = coface_it->second.first.facets();
  // auto list_it = facet_list.begin();
  // while (list_it != facet_list.end())
  //   if (list_it != col_it) {
  //       (*list_it)->second.first.remove_coface(coface_it);
  //       if (Coface_compare<Simplex_with_cofaces>()(*list_it, *col_it))
  //         if ((*list_it)->second.first.number_of_cofaces() == 1 && (*list_it)->second.second == (*list_it)->second.first.coface().first->second.second)
  //           elementary_collapse(*list_it++);
  //         else 
  //           list_it++;
  //       else
  //         list_it++;
  //     }
  //     else
  //       list_it++;
  //   facets_->erase(facet_it);
  //   cofaces_->erase(coface_it);
}


template <class InputRange,
          class OutputIterator,
          class CollapseTraits>
void collapse(const InputRange& input_range,
              OutputIterator& output_it,
              const CollapseTraits& collapse_traits) {
  using Graph = boost::adjacency_list< boost::listS,
                                       boost::listS,
                                       boost::bidirectionalS >;
  using Graph_v = typename Graph::vertex_descriptor;
  using Cell_type = typename CollapseTraits::Cell_type;
  using Cell_comparison = typename CollapseTraits::Cell_comparison;
  using Map = std::map<Cell_type, Graph_v, Cell_comparison>;
  using Inv_map = std::map<Graph_v, typename Map::iterator>;  
  Graph face_coface_graph;
  Map* faces, cofaces;
  Inv_map* inv_faces, inv_cofaces;
  typename InputRange::iterator current_it = input_range.begin();
  if (current_it == input_range.end())
    return;
  int d = collapse_traits.dimension(*current_it);
  cofaces = new Map();
  inv_cofaces = new Inv_map();
  for (int curr_dim = d; curr_dim > 0; curr_dim--) {
    while (current_it != input_range.end() && collapse_traits.dimension(*current_it) == curr_dim) 
      if (cofaces->find(*current_it) == cofaces->end()) {
        Graph_v v = boost::add_vertex(face_coface_graph);
        auto cofacet_it = cofaces->emplace(std::make_pair(*current_it, v)).first;
        inv_cofaces->emplace(std::make_pair(v, cofacet_it));
      }
    faces = new Map();
    inv_faces = new Inv_map();
    for (auto cf_pair: *cofaces) {
      for (auto facet: collapse_traits.facets(cf_pair.first)) {
        auto facet_it = faces->find(facet);
        if (facet_it == faces->end()) {
          Graph_v v = boost::add_vertex(face_coface_graph);
          facet_it = faces->emplace(std::make_pair(facet, v)).first;
          inv_cofaces->emplace(std::make_pair(v, facet_it));
        }
        boost::add_edge(cf_pair.second, facet_it->second, face_coface_graph);
      }
    }
    auto facet_it = faces->begin();
    while (facet_it != faces->end())
      elementary_collapse(facet_it++, faces, inv_faces, cofaces, inv_cofaces, collapse_traits);
    for (auto cf_pair: *cofaces) {
      boost::clear_in_edges(cf_pair.second, face_coface_graph);
      boost::remove_vertex(cf_pair.second, face_coface_graph);
      output_it++ = cf_pair.first;
    }
    delete cofaces;
    delete inv_cofaces;
    cofaces = faces;
    inv_cofaces = inv_faces;
  }
  for (auto cf_pair: *cofaces)
    output_it++ = cf_pair.first;
  delete cofaces;
}

#endif

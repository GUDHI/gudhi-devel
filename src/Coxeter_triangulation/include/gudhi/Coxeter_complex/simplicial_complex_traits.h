/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2018  INRIA Sophia Antipolis-Méditerranée (France)
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

#ifndef SIMPLICIAL_COMPLEX_TRAITS_H_
#define SIMPLICIAL_COMPLEX_TRAITS_H_

#include <vector>

/* The input traits for pairs (simplex,filtration)
 */
template <class Simplex_filtration_pair>
struct Simplex_filtration_pair_input_traits {
  int dimension(const Simplex_filtration_pair& sfp) {
    return sfp.first.size()-1;
  }
  typename Simplex_filtration_pair::second_type filtration(const Simplex_filtration_pair& sfp) {
    return sfp.second;
  }
  typename Simplex_filtration_pair::first_type cell(const Simplex_filtration_pair& sfp) {
    return sfp.first;
  }
};

/* The traits for a filtered simplicial complex collapse.
 * A simplex is assumed to be given as a vector of unsigned integers.
 * The filtrations are assumed to be stored in the maps in the collapse function.
 */
struct Simplicial_complex_collapse_traits {
  typedef std::vector<size_t> Cell_type;
  typedef std::less<Cell_type> Cell_comparison;

  struct Boundary_element : public Cell_type {
    double f;
    Boundary_element(int size, double filtration)
      : Cell_type(size), f(filtration) {}
  };
  // An iterator that traverses through facets of a simplex given as a vector of vertices
  class Simplex_facet_iterator : public boost::iterator_facade<Simplex_facet_iterator,
                                                               Boundary_element const,
                                                               boost::forward_traversal_tag> {
    const Cell_type& simplex_;
    Boundary_element value_;
    unsigned i_;
    friend class boost::iterator_core_access;
    void update_value() {
      if (i_ != simplex_.size()) {
        for (unsigned j = 0; j < i_; ++j)
          value_[j] = simplex_[j];
        for (unsigned j = i_+1; j < simplex_.size(); ++j) 
          value_[j-1] = simplex_[j];
      }
    }
    bool equal(Simplex_facet_iterator const& other) const {
      return i_ == other.i_;
    }
    Boundary_element const& dereference() const {
      return value_;
    }
    void increment() {
      i_++;
      update_value();
    }
  public:
    Simplex_facet_iterator(const Cell_type& simplex, int i, double f)
      : simplex_(simplex), value_(simplex.size()-1, f), i_(i) {
      update_value();
    }
  };
  typedef boost::iterator_range<Simplex_facet_iterator> Boundary_range;

  template <class MapIterator>
  Boundary_range boundary(const MapIterator& cofacet_it) const {
    return Boundary_range(Simplex_facet_iterator(cofacet_it->first,
                                                 0,
                                                 cofacet_it->second.f),
                          Simplex_facet_iterator(cofacet_it->first,
                                                 cofacet_it->first.size(),
                                                 cofacet_it->second.f));
  }

  double facet_filtration(const Boundary_element& b) const {
    return b.f;
  }

  Cell_type facet_cell(const Boundary_element& b) const {
    return (Cell_type)b;
  }

  template <class MapIterator>
  bool is_valid_collapse(const MapIterator& facet_it,
                         const MapIterator& cofacet_it) const {
    return facet_it->second.f == cofacet_it->second.f;
  }
};

/* A structure to be used as an argument boost::make_function_output_iterator
 */
template <class SimplexTree>
struct simplex_tree_inserter {
  void operator() (const std::pair<std::vector<std::size_t>, double>& simplex_and_filtration) {
    st_.insert_simplex(simplex_and_filtration.first,
                       simplex_and_filtration.second);
  }

  simplex_tree_inserter(SimplexTree& st)
    : st_(st) {}

  SimplexTree& st_;
};


/* A function that relabels the vertices in a Simplex tree
 */
template <class SimplexTree>
void relabel_simplex_tree(SimplexTree& stree, SimplexTree& new_stree) {
  using Vertex_handle = typename SimplexTree::Vertex_handle;
  std::map<Vertex_handle, Vertex_handle> label_map;
  Vertex_handle max_index = 0;
  for (auto v: stree.complex_vertex_range())
    label_map.emplace(std::make_pair(v, max_index++));
  // stree.rec_copy(stree.root(),
  //                rec_relabel_nodes(stree.root(), (Siblings*)nullptr, stree, label_map));
  // rec_relabel_nodes(stree.root(), stree, label_map);
  std::vector<int> vertices;
  for (auto sh: stree.complex_simplex_range()) {
    vertices.clear();
    for (auto v: stree.simplex_vertex_range(sh))
      vertices.push_back(label_map[v]);
    new_stree.insert_simplex(vertices);
  }
}

#endif

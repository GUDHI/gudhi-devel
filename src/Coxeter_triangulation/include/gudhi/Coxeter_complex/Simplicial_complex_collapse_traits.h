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

#ifndef SIMPLICIAL_COMPLEX_COLLAPSE_TRAITS_H_
#define SIMPLICIAL_COMPLEX_COLLAPSE_TRAITS_H_

#include <vector>

  /* The traits for a (non-filtered) simplicial complex collapse.
   */
  struct Simplicial_complex_collapse_traits {
    typedef std::vector<size_t> Cell_type;
    typedef std::less<Cell_type> Cell_comparison;

    int dimension(const Cell_type& simplex) const {
      return simplex.size()-1;
    }

    class Simplex_facet_iterator : public boost::iterator_facade<Simplex_facet_iterator,
                                                                 Cell_type const,
                                                                 boost::forward_traversal_tag> {
      const Cell_type& simplex_;
      Cell_type value_;
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
      Cell_type const& dereference() const {
        return value_;
      }
      void increment() {
        i_++;
        update_value();
      }
    public:
      Simplex_facet_iterator(const Cell_type& simplex, int i)
        : simplex_(simplex), value_(simplex.size()-1), i_(i) {
        update_value();
      }
    };
    using Facet_range = boost::iterator_range<Simplex_facet_iterator>;
    
    Facet_range facets(const Cell_type& simplex) const {
      return Facet_range(Simplex_facet_iterator(simplex, 0),
                         Simplex_facet_iterator(simplex, simplex.size()));
    }

    template <class MapIterator>
    bool is_valid_collapse(const MapIterator& facet_it,
                           const MapIterator& cofacet_it) const {
      return true;
    }
  };

/* A structure to be used as an argument boost::make_function_output_iterator
 */
template <class SimplexTree>
struct simplex_tree_non_filtered_inserter {
  void operator() (const std::vector<std::size_t>& simplex) {
    // std::vector<std::size_t> relabeled_simplex;
    // for (std::size_t v: simplex) {
    //   auto m_it = label_map_.find(v);
    //   if (m_it == label_map_.end()) {
    //     label_map_.emplace(std::make_pair(v, max_index_));
    //     relabeled_simplex.push_back(max_index_++);
    //   }
    //   else
    //     relabeled_simplex.push_back(m_it->second);
    // }
    // st_.insert_simplex_and_subfaces(relabeled_simplex);
    st_.insert_simplex(simplex);
  }

  simplex_tree_non_filtered_inserter(SimplexTree& st)
    : st_(st), max_index_(0) {}

  SimplexTree& st_;
  // std::map<std::size_t, std::size_t> label_map_; // Used to output contiguous labels for vertices
  std::size_t max_index_ = 0;
};

#endif

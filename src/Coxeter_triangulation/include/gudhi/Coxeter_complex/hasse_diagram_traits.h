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

#ifndef HASSE_DIAGRAM_TRAITS_H_
#define HASSE_DIAGRAM_TRAITS_H_

#include <functional>   // for std::less
#include <utility>      // for std::pair
#include <vector>       // for std::vector

/* The input traits for pairs simplex-filtration
 */
template <class Hasse_cell>
struct Hasse_cell_input_traits {
  int dimension(const Hasse_cell* cell) const {
    return cell->get_dimension();
  }
  using Filtration_type = typename Hasse_cell::Filtration_type;
  Filtration_type filtration(const Hasse_cell* cell) const {
    return cell->get_filtration();
  }
  Hasse_cell* cell(Hasse_cell* cell) const {
    return cell;
  }
};

/* The traits for a Hasse diagram collapse.
 * Hasse_diagram template parameter should have remove_cell function, which
 * marks a cell for deletion instead of deleting it.
 * The traits mark the collapsed cells for deletion.
 */
template <class Hasse_cell,
          class Hasse_diagram>
struct Hasse_diagram_collapse_traits {
  typedef Hasse_cell* Cell_type;
  typedef std::less<Cell_type> Cell_comparison;

  using Incidence_type = typename Hasse_cell::Incidence_type;
  typedef std::pair<Hasse_cell*,Incidence_type> Boundary_element;
  typedef std::vector<Boundary_element> Boundary_range;
  
  template <class MapIterator>
  Boundary_range boundary(const MapIterator& cofacet_it) const {
    return cofacet_it->first->get_boundary();
  }

  using Filtration_type = typename Hasse_cell::Filtration_type;
  Filtration_type facet_filtration(const Boundary_element& b) const {
    return b.first->get_filtration();
  }

  Cell_type facet_cell(const Boundary_element& b) const {
    return b.first;
  }
  
  template <class MapIterator>
  bool is_valid_collapse(const MapIterator& facet_it,
                         const MapIterator& cofacet_it) {
    bool is_valid = ignore_filtration_ || facet_it->second.f == cofacet_it->second.f;
    if (is_valid) {
      hasse_diagram_.remove_cell(facet_it->first);
      hasse_diagram_.remove_cell(cofacet_it->first);
    }
    return is_valid;
  }

  Hasse_diagram_collapse_traits(Hasse_diagram& hasse_diagram, bool ignore_filtration = false) 
    : hasse_diagram_(hasse_diagram), ignore_filtration_(ignore_filtration) {}

protected:
  Hasse_diagram& hasse_diagram_;
  bool ignore_filtration_;
};


#endif

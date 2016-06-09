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

#ifndef DIM_LISTS_H_
#define DIM_LISTS_H_

#include <boost/container/flat_map.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <algorithm>
#include <utility>
#include "gudhi/reader_utils.h"
#include "gudhi/distance_functions.h"
#include <gudhi/Dim_list_iterator.h>
#include <vector>
#include <list>
#include <set>
#include <queue>
#include <limits>
#include <math.h>
#include <ctime>
#include <iostream>

#include <boost/tuple/tuple.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/iterator_range.hpp>

// Needed for the adjacency graph in bad link search
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

namespace Gudhi {

namespace witness_complex {

  /** \addtogroup simplex_tree
   *  Witness complex is a simplicial complex defined on two sets of points in \f$\mathbf{R}^D\f$:
   *  \f$W\f$ set of witnesses and \f$L \subseteq W\f$ set of landmarks. The simplices are based on points in \f$L\f$
   *  and a simplex belongs to the witness complex if and only if it is witnessed (there exists a point \f$w \in W\f$ such that
   *  w is closer to the vertices of this simplex than others) and all of its faces are witnessed as well. 
   */
template< class Simplicial_complex >
class Dim_lists {

private:

  typedef typename Simplicial_complex::Simplex_handle Simplex_handle;
  typedef typename Gudhi::witness_complex::Dim_lists_iterator<Simplex_handle> Iterator;
  
  std::vector<std::list<Simplex_handle>> table_;
  Simplicial_complex& sc_;

public:

  Dim_lists(Simplicial_complex & sc, int dim, double alpha_max = 100)
    : sc_(sc)
  {
    table_ = std::vector<std::list<Simplex_handle>>(dim+1);
    for (auto sh: sc.filtration_simplex_range()) {
      if (sc_.filtration(sh) < alpha_max)
        table_[sc.dimension(sh)].push_front(sh);
    }
    auto t_it = table_.rbegin();
    while (t_it->empty()) {
      t_it++;
      table_.pop_back();
    }
  }

  Iterator begin()
  {
    return Iterator(table_.rbegin()->begin(), table_.rbegin(), table_);
  }

  Iterator end()
  {
    return Iterator(table_[0].end(), table_.rend(), table_);
  }
  
  unsigned size()
  {
    unsigned curr_size = 0;
    for (auto l: table_)
      curr_size += l.size();
    return curr_size;
  }
  
  bool is_face(Simplex_handle face, Simplex_handle coface)
  {
    // vertex range is sorted in decreasing order
    auto fvr = sc_.simplex_vertex_range(face);
    auto cfvr = sc_.simplex_vertex_range(coface);
    auto fv_it = fvr.begin();
    auto cfv_it = cfvr.begin();
    while (fv_it != fvr.end() && cfv_it != cfvr.end()) {
      if (*fv_it < *cfv_it)
        ++cfv_it;
      else if (*fv_it == *cfv_it) {
        ++cfv_it;
        ++fv_it;
      }
      else
        return false;
      
    }
    return (fv_it == fvr.end());
  }

  void output_simplices() {
    std::cout << "Size of vector: " << size() << std::endl;
    for (auto line_it = table_.rbegin(); line_it != table_.rend(); ++line_it)
      for (auto sh: *line_it) {
        std::cout << sc_.dimension(sh) << " ";
        for (auto v : sc_.simplex_vertex_range(sh))
          std::cout << v << " ";
        std::cout << sc_.filtration(sh) << "\n";
      }
  }
  
  void collapse()
  {
    auto coface_list_it = table_.rbegin();
    auto face_list_it = table_.rbegin()+1;
    for ( ;
          face_list_it != table_.rend();
          ++face_list_it) {
      auto face_it = face_list_it->begin();
      while (face_it != face_list_it->end() && sc_.filtration(*face_it) != 0) {
        int coface_count = 0;
        auto reduced_coface = coface_list_it->begin();
        for (auto coface_it = coface_list_it->begin(); coface_it != coface_list_it->end() && sc_.filtration(*coface_it) != 0; ++coface_it)
          if (is_face(*face_it, *coface_it)) {
            coface_count++;
            if (coface_count == 1)
              reduced_coface = coface_it;
            else
              break;
          }
        if (coface_count == 1) {
          /*
          std::cout << "Erase ( ";
          for (auto v: sc_.simplex_vertex_range(*reduced_coface))
            std::cout << v << " ";        
          */
          coface_list_it->erase(reduced_coface);
          /*
          std::cout << ") and then ( ";
          for (auto v: sc_.simplex_vertex_range(*face_it))
            std::cout << v << " ";
          std::cout << ")\n";
          */
          face_list_it->erase(face_it);
          face_it = face_list_it->begin();
        } 
        else
          face_it++;
      }
      if ((coface_list_it++)->empty())
        table_.pop_back();
    }
  }
  
  bool is_pseudomanifold()
  {

    return true;
  }
  
}; //class Dim_lists

} // namespace witness_complex
  
} // namespace Gudhi

#endif

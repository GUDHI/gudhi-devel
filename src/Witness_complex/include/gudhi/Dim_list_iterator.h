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

#ifndef DIM_LISTS_ITERATOR_H_
#define DIM_LISTS_ITERATOR_H_

#include <boost/container/flat_map.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <algorithm>
#include <utility>
#include "gudhi/reader_utils.h"
#include "gudhi/distance_functions.h"
#include "gudhi/Simplex_tree.h"
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
template< class Simplex_handle >
class Dim_lists_iterator {

private:

  typedef typename std::list<Simplex_handle>::iterator List_iterator;
  typedef typename std::vector<std::list<Simplex_handle>>::reverse_iterator Vector_iterator;
  typedef typename Gudhi::witness_complex::Dim_lists_iterator<Simplex_handle> Iterator;
  
  
  List_iterator sh_;
  Vector_iterator curr_list_;
  typename std::vector<std::list<Simplex_handle>>& table_;

public:

  Dim_lists_iterator(List_iterator sh,
                     Vector_iterator curr_list,
                     typename std::vector<std::list<Simplex_handle>>& table)
    : sh_(sh), curr_list_(curr_list), table_(table)
  {
  }

  Simplex_handle operator*()
  {
    return *sh_;
  }
  
  Iterator operator++()
  {
    increment();
    return (*this);
  }

  Iterator operator++(int)
  {
    Iterator prev_it(sh_, curr_list_, table_);
    increment();
    return prev_it;
  }

  Iterator dim_begin()
  {
    return Iterator(curr_list_->begin(), curr_list_, table_);
  }  

  Iterator dim_end()
  {
    return Iterator(curr_list_->end(), curr_list_, table_);
  }
  
  Iterator dimp1_begin()
  {
    return Iterator((curr_list_-1)->begin(), curr_list_-1, table_); 
  }

  Iterator dimp1_end()
  {
    return Iterator((curr_list_-1)->end(), curr_list_-1, table_); 
  }
  
  bool operator==(const Iterator& it2) const
  {
    return (sh_ == it2.sh_);
  }

  bool operator!=(const Iterator& it2) const
  {
    return (sh_ != it2.sh_);
  }

  void remove_incr()
  {
    
  }
  
private:
  
  void increment()
  {
    if (++sh_ == curr_list_->end())
      if (++curr_list_ != table_.rend())
        sh_ = curr_list_->begin();
    // The iterator of the end of the table is the end of the last list
  }

  
}; //class Dim_lists_iterator

} // namespace witness_complex
  
} // namespace Gudhi

#endif

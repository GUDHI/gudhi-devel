/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2016  INRIA Sophia Antipolis-Méditerranée (France)
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

#ifndef ACTIVE_WITNESS_H_
#define ACTIVE_WITNESS_H_

#include "Active_witness_iterator.h"
#include <vector>
#include <utility>

namespace Gudhi {

namespace witness_complex {

  /** \brief Class representing a list of nearest neighbors to a given witness.
   *  \detail Every element is a pair of a landmark identifier and the squared distance to it.
   */
template< typename Id_distance_pair,
          typename INS_range >
class Active_witness {
public:  
  typedef Active_witness<Id_distance_pair, INS_range> ActiveWitness;
  typedef typename INS_range::iterator INS_iterator;
  typedef Active_witness_iterator< ActiveWitness, Id_distance_pair, INS_iterator > iterator;
  typedef typename std::list<Id_distance_pair> Table;

  Table end_element_table_ = {Id_distance_pair(-1,0)};
  typename Table::iterator end_pointer = end_element_table_.begin();

  Table nearest_landmark_table_;
  INS_range    search_range_;
  INS_iterator iterator_last_;
  INS_iterator iterator_end_;

  Active_witness(INS_range search_range)
    : search_range_(search_range), iterator_last_(search_range.begin()), iterator_end_(search_range.end())
  {
    nearest_landmark_table_.push_back(*iterator_last_);
  }
 
  iterator begin()
  {
    return iterator(this, nearest_landmark_table_.begin());
  }

  iterator end()
  {
    return iterator(this);
  }
};

}
}
  
#endif

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

#ifndef ACTIVE_WITNESS_ITERATOR_H_
#define ACTIVE_WITNESS_ITERATOR_H_

//#include "Active_witness.h"
#include <boost/iterator/iterator_facade.hpp>
#include <vector>

namespace Gudhi {

namespace witness_complex {

  //  /** \brief Iterator in the nearest landmark list.
  //   *  \details After the iterator reaches the end of the list,
  //   *          the list is augmented by a (nearest landmark, distance) pair if possible.
  //   *          If all the landmarks are present in the list, iterator returns the specific end value
  //   *          of the corresponding 'Active_witness' object.
  //   */
  
template< typename Active_witness,
          typename Id_distance_pair,
          typename INS_iterator >
class Active_witness_iterator
  : public boost::iterator_facade< Active_witness_iterator <Active_witness, Id_distance_pair, INS_iterator>
, Id_distance_pair const
, boost::forward_traversal_tag
, Id_distance_pair const> {
  friend class boost::iterator_core_access;
  
  //typedef Active_witness<Id_distance_pair, INS_iterator> Active_witness;
  typedef typename std::list<Id_distance_pair>::iterator Pair_iterator;
  typedef typename Gudhi::witness_complex::Active_witness_iterator<Active_witness, Id_distance_pair, INS_iterator> Iterator;
  
  
  Active_witness *aw_;
  Pair_iterator lh_; // landmark handle
  //INS_iterator iterator_last;
  //INS_iterator iterator_end;

public:
  Active_witness_iterator(Active_witness* aw)
    : aw_(aw), lh_(aw_->end_pointer)
  {
  }

  Active_witness_iterator(Active_witness* aw, Pair_iterator lh)
    : aw_(aw), lh_(lh)
  {
  }
  
private :

  Id_distance_pair& dereference() const
  {
    return *lh_;
  }

  bool equal(const Iterator& other) const
  {
    return (lh_ == other.lh_);
  }
  
  void increment()
  {
    // if neighbor search is at its end, check if lh_++ is end
    if (aw_->iterator_last_ == aw_->iterator_end_) {
      if (lh_++ == aw_->nearest_landmark_table_.end()) {
        lh_ = aw_->end_pointer;
        return;
      }
      return;
    }
    // if the id of the current landmark is the same as the last one
    if (lh_->first == aw_->iterator_last_->first) {
      // if the next iterator is end, lh_it = end pointer
      if (++(aw_->iterator_last_) == aw_->iterator_end_) {
        lh_ = aw_->end_pointer;
        return;
      }
      else 
        aw_->nearest_landmark_table_.push_back(*(aw_->iterator_last_));
    }
    lh_++;
  }

};

}
}
  
#endif

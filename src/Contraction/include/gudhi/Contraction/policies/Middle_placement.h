/*    This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author(s):       David Salinas
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

#ifndef CONTRACTION_POLICIES_MIDDLE_PLACEMENT_H_
#define CONTRACTION_POLICIES_MIDDLE_PLACEMENT_H_

#include <gudhi/Contraction/policies/Placement_policy.h>

namespace Gudhi {

namespace contraction {

template< typename EdgeProfile>
class Middle_placement : public Placement_policy<EdgeProfile> {
 public:
  typedef typename EdgeProfile::Point Point;
  typedef typename EdgeProfile::Edge_handle Edge_handle;
  typedef typename EdgeProfile::Graph_vertex Graph_vertex;

  typedef typename Placement_policy<EdgeProfile>::Placement_type Placement_type;

  Placement_type operator()(const EdgeProfile& profile) const override {
    // todo compute the middle
    return Placement_type(profile.p0());
  }
};

}  // namespace contraction

}  // namespace Gudhi

#endif  // CONTRACTION_POLICIES_MIDDLE_PLACEMENT_H_

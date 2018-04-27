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

#ifndef CONTRACTION_POLICIES_PLACEMENT_POLICY_H_
#define CONTRACTION_POLICIES_PLACEMENT_POLICY_H_

#include <boost/optional.hpp>

namespace Gudhi {

namespace contraction {

/**
 *@brief Policy to specify where the merged point had to be placed after an edge contraction. 
 *@ingroup contr
 */
template< typename EdgeProfile>
class Placement_policy {
 public:
  typedef typename EdgeProfile::Point Point;
  typedef boost::optional<Point> Placement_type;

  virtual Placement_type operator()(const EdgeProfile& profile) const = 0;

  virtual ~Placement_policy() { }
};

}  // namespace contraction

}  // namespace Gudhi

#endif  // CONTRACTION_POLICIES_PLACEMENT_POLICY_H_

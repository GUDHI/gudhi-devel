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

#ifndef CONTRACTION_POLICIES_VALID_CONTRACTION_POLICY_H_
#define CONTRACTION_POLICIES_VALID_CONTRACTION_POLICY_H_

namespace Gudhi {

namespace contraction {

/**
 *@brief Policy to specify if an edge contraction is valid or not. 
 *@ingroup contr
 */
template< typename EdgeProfile>
class Valid_contraction_policy {
 public:
  typedef typename EdgeProfile::Point Point;
  typedef typename EdgeProfile::Edge_handle Edge_handle;
  typedef typename EdgeProfile::Graph_vertex Graph_vertex;

  virtual bool operator()(const EdgeProfile& profile, const boost::optional<Point>& placement) const = 0;

  virtual ~Valid_contraction_policy() { }
};

}  // namespace contraction

}  // namespace Gudhi


#endif  // CONTRACTION_POLICIES_VALID_CONTRACTION_POLICY_H_

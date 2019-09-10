/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
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

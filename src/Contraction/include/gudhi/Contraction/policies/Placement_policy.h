/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
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

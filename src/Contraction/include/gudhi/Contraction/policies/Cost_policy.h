/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CONTRACTION_POLICIES_COST_POLICY_H_
#define CONTRACTION_POLICIES_COST_POLICY_H_

#include <boost/optional.hpp>

namespace Gudhi {

namespace contraction {

/**
 *@brief Policy to specify the cost of contracting an edge.
 *@ingroup contr
 */
template< typename EdgeProfile>
class Cost_policy {
 public:
  typedef typename EdgeProfile::Point Point;
  typedef typename EdgeProfile::Graph_vertex Graph_vertex;

  typedef boost::optional<double> Cost_type;

  virtual Cost_type operator()(const EdgeProfile& profile, const boost::optional<Point>& placement) const = 0;

  virtual ~Cost_policy() { }
};

}  // namespace contraction

}  // namespace Gudhi

#endif  // CONTRACTION_POLICIES_COST_POLICY_H_

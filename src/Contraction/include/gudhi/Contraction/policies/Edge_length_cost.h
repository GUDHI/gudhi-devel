/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CONTRACTION_POLICIES_EDGE_LENGTH_COST_H_
#define CONTRACTION_POLICIES_EDGE_LENGTH_COST_H_

#include <gudhi/Contraction/policies/Cost_policy.h>

namespace Gudhi {

namespace contraction {

/**
 * @brief return a cost corresponding to the squared length of the edge
 */
template< typename EdgeProfile>
class Edge_length_cost : public Cost_policy<EdgeProfile> {
 public:
  typedef typename Cost_policy<EdgeProfile>::Cost_type Cost_type;
  typedef typename EdgeProfile::Point Point;

  Cost_type operator()(const EdgeProfile& profile, const boost::optional<Point>& placement) const override {
    double res = 0;
    auto p0_coord = profile.p0().begin();
    auto p1_coord = profile.p1().begin();
    for (; p0_coord != profile.p0().end(); p0_coord++, p1_coord++) {
      res += (*p0_coord - *p1_coord) * (*p0_coord - *p1_coord);
    }
    return res;
  }
};

}  // namespace contraction

}  // namespace Gudhi

#endif  // CONTRACTION_POLICIES_EDGE_LENGTH_COST_H_

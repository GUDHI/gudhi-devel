/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CONTRACTION_POLICIES_FIRST_VERTEX_PLACEMENT_H_
#define CONTRACTION_POLICIES_FIRST_VERTEX_PLACEMENT_H_

#include <gudhi/Contraction/policies/Placement_policy.h>

namespace Gudhi {

namespace contraction {

/**
 * @brief Places the contracted point onto the first point of the edge
 */
template< typename EdgeProfile>
class First_vertex_placement : public Placement_policy<EdgeProfile> {
 public:
  typedef typename EdgeProfile::Point Point;
  typedef typename EdgeProfile::Edge_handle Edge_handle;

  typedef typename Placement_policy<EdgeProfile>::Placement_type Placement_type;

  Placement_type operator()(const EdgeProfile& profile) const override {
    return Placement_type(profile.p0());
  }
};

}  // namespace contraction

}  // namespace Gudhi

#endif  // CONTRACTION_POLICIES_FIRST_VERTEX_PLACEMENT_H_

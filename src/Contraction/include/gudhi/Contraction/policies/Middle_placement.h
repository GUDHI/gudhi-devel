/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
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

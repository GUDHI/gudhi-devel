/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CONTRACTION_POLICIES_LINK_CONDITION_VALID_CONTRACTION_H_
#define CONTRACTION_POLICIES_LINK_CONDITION_VALID_CONTRACTION_H_

#include <gudhi/Contraction/policies/Valid_contraction_policy.h>
#include <gudhi/Debug_utils.h>


namespace Gudhi {

namespace contraction {

/**
 *@brief Policy that only accept edges verifying the link condition (and therefore whose contraction preserving homotopy type). 
 *@ingroup contr
 */
template< typename EdgeProfile>
class Link_condition_valid_contraction : public Valid_contraction_policy<EdgeProfile> {
 public:
  typedef typename EdgeProfile::Edge_handle Edge_handle;
  typedef typename EdgeProfile::Point Point;
  // typedef typename EdgeProfile::Edge_handle Edge_handle;

  bool operator()(const EdgeProfile& profile, const boost::optional<Point>& placement) const override {
    Edge_handle edge(profile.edge_handle());
    DBGMSG("Link_condition_valid_contraction:", profile.complex().link_condition(edge));
    return profile.complex().link_condition(edge);
  }
};

}  // namespace contraction

}  // namespace Gudhi

#endif  // CONTRACTION_POLICIES_LINK_CONDITION_VALID_CONTRACTION_H_

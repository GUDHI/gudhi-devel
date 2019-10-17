/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CONTRACTION_POLICIES_DUMMY_VALID_CONTRACTION_H_
#define CONTRACTION_POLICIES_DUMMY_VALID_CONTRACTION_H_

#include <gudhi/Contraction/policies/Valid_contraction_policy.h>

namespace Gudhi {

namespace contraction {

/**
 *@brief Policy that accept all edge contraction.
 */
template< typename EdgeProfile>
class Dummy_valid_contraction : public Valid_contraction_policy<EdgeProfile> {
 public:
  typedef typename EdgeProfile::Point Point;

  bool operator()(const EdgeProfile& profile, const boost::optional<Point>& placement) {
    return true;
  }
};

}  // namespace contraction

}  // namespace Gudhi

#endif  // CONTRACTION_POLICIES_DUMMY_VALID_CONTRACTION_H_

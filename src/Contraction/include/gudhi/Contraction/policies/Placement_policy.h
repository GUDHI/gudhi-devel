/*
 * Placement_policy.h
 *
 *  Created on: Feb 13, 2014
 *      Author: dsalinas
 */

#ifndef GUDHI_PLACEMENT_POLICY_H_
#define GUDHI_PLACEMENT_POLICY_H_

#include <boost/optional.hpp>

namespace Gudhi {
namespace contraction {

template< typename EdgeProfile> class Placement_policy{
public:
	typedef typename EdgeProfile::Point Point;
	typedef boost::optional<Point> Placement_type;

	virtual Placement_type operator()(const EdgeProfile& profile) const=0;
	virtual ~Placement_policy(){};
};


}  // namespace contraction
}  // namespace GUDHI

#endif /* GUDHI_PLACEMENT_POLICY_H_ */

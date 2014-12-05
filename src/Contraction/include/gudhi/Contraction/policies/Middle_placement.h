/*
 * Middle_placement.h
 *
 *  Created on: Feb 13, 2014
 *      Author: dsalinas
 */

#ifndef GUDHI_MIDDLE_PLACEMENT_H_
#define GUDHI_MIDDLE_PLACEMENT_H_

#include "Placement_policy.h"


namespace Gudhi{

namespace contraction {

template< typename EdgeProfile> class Middle_placement : public Placement_policy<EdgeProfile>{

public:
	typedef typename EdgeProfile::Point Point;
	typedef typename EdgeProfile::Edge_handle Edge_handle;
	typedef typename EdgeProfile::Graph_vertex Graph_vertex;

	typedef typename Placement_policy<EdgeProfile>::Placement_type Placement_type;

	Placement_type operator()(const EdgeProfile& profile) const override{
		//todo compute the middle
		return Placement_type(profile.p0());
	}
};
}  // namespace contraction
}  // namespace GUDHI

#endif /* GUDHI_MIDDLE_PLACEMENT_H_ */

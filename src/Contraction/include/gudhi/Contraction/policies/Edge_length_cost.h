/*
 * Edge_length_cost.h
 *
 *  Created on: Feb 13, 2014
 *      Author: dsalinas
 */

#ifndef GUDHI_EDGE_LENGTH_COST_H_
#define GUDHI_EDGE_LENGTH_COST_H_

#include "Cost_policy.h"

namespace Gudhi{

namespace contraction {


/**
 * @brief return a cost corresponding to the squared length of the edge
 */
template< typename EdgeProfile> class Edge_length_cost : public Cost_policy<EdgeProfile>{
public:
	typedef typename Cost_policy<EdgeProfile>::Cost_type Cost_type;
	typedef typename EdgeProfile::Point Point;
	Cost_type operator()(const EdgeProfile& profile, const boost::optional<Point>& placement) const override{
		double res = 0;
		auto p0_coord = profile.p0().begin();
		auto p1_coord = profile.p1().begin();
		for(; p0_coord != profile.p0().end(); p0_coord++, p1_coord++){
			res += (*p0_coord - *p1_coord) * (*p0_coord - *p1_coord);
		}
		return res;
	}

};

}  // namespace contraction

}  // namespace GUDHI

#endif /* GUDHI_EDGE_LENGTH_COST_H_ */

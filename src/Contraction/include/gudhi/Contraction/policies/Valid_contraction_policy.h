/*
 * Valid_contraction_policy.h
 *
 *  Created on: Feb 13, 2014
 *      Author: dsalinas
 */

#ifndef GUDHI_VALID_CONTRACTION_POLICY_H_
#define GUDHI_VALID_CONTRACTION_POLICY_H_

namespace Gudhi {
namespace contraction {
template< typename EdgeProfile> class Valid_contraction_policy{
public:
	typedef typename EdgeProfile::Point Point;
	typedef typename EdgeProfile::Edge_handle Edge_handle;
	typedef typename EdgeProfile::Graph_vertex Graph_vertex;

	virtual bool operator()(const EdgeProfile& profile,const boost::optional<Point>& placement) const =0;
	virtual ~Valid_contraction_policy(){};
};

}  // namespace contraction
}  // namespace GUDHI


#endif /* GUDHI_VALID_CONTRACTION_POLICY_H_ */

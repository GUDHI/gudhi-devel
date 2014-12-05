/*
 * Dummy_valid_contraction.h
 *
 *  Created on: Feb 13, 2014
 *      Author: dsalinas
 */

#ifndef GUDHI_DUMMY_VALID_CONTRACTION_H_
#define GUDHI_DUMMY_VALID_CONTRACTION_H_

#include "Valid_contraction_policy.h"

namespace Gudhi{

namespace contraction {




template< typename EdgeProfile> class Dummy_valid_contraction : public Valid_contraction_policy<EdgeProfile>{
public:
	typedef typename EdgeProfile::Point Point;
	bool operator()(const EdgeProfile& profile,const boost::optional<Point>& placement){
		return true;
	}
};

}  // namespace contraction

}  // namespace GUDHI



#endif /* GUDHI_DUMMY_VALID_CONTRACTION_H_ */

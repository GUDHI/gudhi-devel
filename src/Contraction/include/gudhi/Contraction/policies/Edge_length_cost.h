 /*    This file is part of the Gudhi Library. The Gudhi library 
  *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
  *    library for computational topology.
  *
  *    Author(s):       David Salinas
  *
  *    Copyright (C) 2014  INRIA Sophia Antipolis-Mediterranee (France)
  *
  *    This program is free software: you can redistribute it and/or modify
  *    it under the terms of the GNU General Public License as published by
  *    the Free Software Foundation, either version 3 of the License, or
  *    (at your option) any later version.
  *
  *    This program is distributed in the hope that it will be useful,
  *    but WITHOUT ANY WARRANTY; without even the implied warranty of
  *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *    GNU General Public License for more details.
  *
  *    You should have received a copy of the GNU General Public License
  *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
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

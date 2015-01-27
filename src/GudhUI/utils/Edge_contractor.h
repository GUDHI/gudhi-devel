/*
 * Contractor.h
 *
 *  Created on: Sep 25, 2014
 *      Author: dsalinas
 */

#ifndef EDGE_CONTRACTOR_H_
#define EDGE_CONTRACTOR_H_


#include "gudhi/Skeleton_blocker_contractor.h"

#include "gudhi/Contraction/Edge_profile.h"
#include "gudhi/Contraction/policies/Cost_policy.h"


/**
 * Iteratively puts every vertex at the center of its neighbors
 */
template<typename SkBlComplex> class Edge_contractor{
private:
	SkBlComplex& complex_;
	unsigned num_contractions_;

	/**
	 * @brief return a cost corresponding to the squared length of the edge
	 */
	template< typename EdgeProfile> class Length_cost : public contraction::Cost_policy<EdgeProfile>{
	public:
		typedef typename contraction::Cost_policy<EdgeProfile>::Cost_type Cost_type;
		typedef typename EdgeProfile::Point Point;
		Cost_type operator()(const EdgeProfile& profile, const boost::optional<Point>& placement) const override{
			Cost_type res;
			if(!placement) return res;
			return Geometry_trait::Squared_distance_d()(profile.p0(),profile.p1()); //not working??
		}
	};

	public:
	typedef typename SkBlComplex::Vertex_handle Vertex_handle;
	typedef typename SkBlComplex::Edge_handle Edge_handle;

	/**
	 * @brief Modify complex to be the expansion of the k-nearest neighbor
	 * symetric graph.
	 */
	Edge_contractor(SkBlComplex& complex,unsigned num_contractions):
		complex_(complex),num_contractions_(num_contractions)
	{
		typedef typename contraction::Edge_profile<Complex> Profile;
		num_contractions = (std::min)((int)num_contractions,(int)(complex_.num_vertices()-1));
		contraction::Skeleton_blocker_contractor<Complex> contractor(
				complex_,
				new Length_cost<contraction::Edge_profile<Complex>>(),
				contraction::make_first_vertex_placement<Profile>(),
				contraction::make_link_valid_contraction<Profile>(),
				contraction::make_remove_popable_blockers_visitor<Profile>()
		);
		contractor.contract_edges(num_contractions);
	}


};



#endif /* EDGE_CONTRACTOR_H_ */

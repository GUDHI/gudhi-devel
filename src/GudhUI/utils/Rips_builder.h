/*
 * Rips_builder.h
 *
 *  Created on: Sep 10, 2014
 *      Author: dsalinas
 */

#ifndef RIPS_BUILDER_H_
#define RIPS_BUILDER_H_

#include <boost/iterator/iterator_facade.hpp>

#include "utils/UI_utils.h"
#include "model/Complex_typedefs.h"
#include <CGAL/Euclidean_distance.h>

#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_d.h>

template<typename SkBlComplex> class Rips_builder{
private:
	SkBlComplex& complex_;
public:

	/**
	 * @brief Modify complex to be the Rips complex
	 * of its points with offset alpha.
	 */
	Rips_builder(SkBlComplex& complex,double alpha):complex_(complex){
		complex.keep_only_vertices();
		if (alpha<=0) return;
		compute_edges(alpha);
	}

private:


	double squared_eucl_distance(const Point& p1,const Point& p2) const{
		return Geometry_trait::Squared_distance_d()(p1,p2);
	}

	void compute_edges(double alpha){
		auto vertices = complex_.vertex_range();
		for(auto p = vertices.begin(); p!= vertices.end(); ++p){
			std::cout << *p << " "; std::cout.flush();
			for (auto q = p; ++q != vertices.end(); /**/)
				if (squared_eucl_distance(complex_.point(*p),complex_.point(*q)) < 4*alpha*alpha)
					complex_.add_edge(*p,*q);
		}
		std::cout << std::endl;
	}

};


#endif /* RIPS_BUILDER_H_ */

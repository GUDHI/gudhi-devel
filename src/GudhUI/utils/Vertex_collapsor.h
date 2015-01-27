/*
 * Vertex_collapsor.h
 *
 *  Created on: Sep 25, 2014
 *      Author: dsalinas
 */

#ifndef VERTEX_COLLAPSOR_H_
#define VERTEX_COLLAPSOR_H_

#include "utils/Edge_contractor.h"
#include "utils/Furthest_point_epsilon_net.h"
#include "utils/UI_utils.h"
/**
 * Iteratively puts every vertex at the center of its neighbors
 */
template<typename SkBlComplex> class Vertex_collapsor{
private:
	SkBlComplex& complex_;
	size_t num_collapses_;
public:
	typedef typename SkBlComplex::Vertex_handle Vertex_handle;
	typedef typename SkBlComplex::Edge_handle Edge_handle;

	/**
	 * @brief Modify complex to be the expansion of the k-nearest neighbor
	 * symetric graph.
	 */
	Vertex_collapsor(SkBlComplex& complex, size_t num_collapses) :
		complex_(complex),num_collapses_(num_collapses)
	{
//		std::list<Vertex_handle> vertices;
//		vertices.insert(vertices.begin(),complex_.vertex_range().begin(),complex_.vertex_range().end());
//		UIDBG("Collapse vertices");
//		collapse_vertices(vertices);

		std::list<Vertex_handle> vertices;

		UIDBG("Compute eps net");
		Furthest_point_epsilon_net<Complex> eps_net(complex_);

		for(auto vh : eps_net.net_filtration_)
			vertices.push_back(vh.vertex_handle);

		UIDBG("Collapse vertices");
		collapse_vertices(vertices);



	}

private:


	void collapse_vertices(std::list<Vertex_handle>& vertices){
		while(!vertices.empty() && num_collapses_--){
			Vertex_handle current_vertex = vertices.front();
			vertices.pop_front();
			if(is_link_reducible(current_vertex))
				complex_.remove_vertex(current_vertex);
		}
	}

	bool is_link_reducible(Vertex_handle v){
		auto link = complex_.link(v);
		if(link.empty()) return false;
		if(link.is_cone()) return true;
		if(link.num_connected_components()>1) return false;
		Edge_contractor<Complex> contractor(link,link.num_vertices()-1);
		return (link.num_vertices()==1);
	}

};


#endif /* VERTEX_COLLAPSOR_H_ */

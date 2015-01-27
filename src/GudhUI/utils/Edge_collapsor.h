/*
 * Collapsor.h
 *
 *  Created on: Sep 25, 2014
 *      Author: dsalinas
 */

#ifndef COLLAPSOR_H_
#define COLLAPSOR_H_

#include <list>
#include "utils/Edge_contractor.h"

/**
 * Iteratively puts every vertex at the center of its neighbors
 */
template<typename SkBlComplex> class Edge_collapsor{
private:
	SkBlComplex& complex_;
	unsigned num_collapses_;
public:
	typedef typename SkBlComplex::Vertex_handle Vertex_handle;
	typedef typename SkBlComplex::Edge_handle Edge_handle;

	/**
	 * @brief Modify complex to be the expansion of the k-nearest neighbor
	 * symetric graph.
	 */
	Edge_collapsor(SkBlComplex& complex,unsigned num_collapses):
		complex_(complex),num_collapses_(num_collapses)
	{
		std::list<Edge_handle> edges;
		edges.insert(edges.begin(),complex_.edge_range().begin(),complex_.edge_range().end());

		edges.sort(
				[&](Edge_handle e1,Edge_handle e2){
			return squared_edge_length(e1) < squared_edge_length(e2);
		});

		collapse_edges(edges);

	}

private:


	void collapse_edges(std::list<Edge_handle>& edges){
		while(!edges.empty() && num_collapses_--){
			Edge_handle current_edge = edges.front();
			edges.pop_front();
			if(is_link_reducible(current_edge))
				complex_.remove_edge(current_edge);
		}

	}

	bool is_link_reducible(Edge_handle e){
		auto link = complex_.link(e);

		if(link.empty()) return false;

		if(link.is_cone()) return true;

		if(link.num_connected_components()>1) return false;

		Edge_contractor<Complex> contractor(link,link.num_vertices()-1);

		return (link.num_vertices()==1);
	}


	double squared_edge_length(Edge_handle e) const{
		return squared_eucl_distance(complex_.point(complex_.first_vertex(e)),complex_.point(complex_.second_vertex(e)));
	}

	double squared_eucl_distance(const Point& p1,const Point& p2) const{
		return Geometry_trait::Squared_distance_d()(p1,p2);
	}

};



#endif /* COLLAPSOR_H_ */

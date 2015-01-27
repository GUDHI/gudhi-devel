/*
 * Furthest_point_epsilon_net.h
 *
 *  Created on: Sep 26, 2014
 *      Author: dsalinas
 */

#ifndef FURTHEST_POINT_EPSILON_NET_H_
#define FURTHEST_POINT_EPSILON_NET_H_

#include "utils/UI_utils.h"
#include <vector>

/**
 * Computes an epsilon net with furthest point strategy.
 */
template<typename SkBlComplex> class Furthest_point_epsilon_net{
private:
	SkBlComplex& complex_;
	typedef typename SkBlComplex::Vertex_handle Vertex_handle;
	typedef typename SkBlComplex::Edge_handle Edge_handle;

	/**
	 * Let V be the set of vertices.
	 * Initially v0 is one arbitrarly vertex and the set V0 is {v0}.
	 * Then Vk is computed as follows.
	 * First we compute the vertex pk that is the furthest from Vk
	 * then Vk = Vk \cup pk.
	 * The radius of pk is its distance to Vk and its meeting vertex
	 * is the vertex of Vk for which this distance is achieved.
	 */
	struct Net_filtration_vertex{
		Vertex_handle vertex_handle;
		Vertex_handle meeting_vertex;
		double radius;


		Net_filtration_vertex(
				Vertex_handle vertex_handle_,
				Vertex_handle meeting_vertex_,
				double radius_):
					vertex_handle(vertex_handle_),meeting_vertex(meeting_vertex_),radius(radius_)
		{}

		bool operator<(const Net_filtration_vertex& other ) const{
			return radius < other.radius;
		}

	};

public:


	std::vector<Net_filtration_vertex> net_filtration_;

	/**
	 * @brief Modify complex to be the expansion of the k-nearest neighbor
	 * symetric graph.
	 */
	Furthest_point_epsilon_net(SkBlComplex& complex):
		complex_(complex)
	{
		if(!complex.empty()){
			init_filtration();
			for(int k = 2; k < net_filtration_.size(); ++k){
				update_radius_value(k);
			}
		}
	}

	//xxx does not work if complex not full
	double radius(Vertex_handle v){
		return net_filtration_[v.vertex].radius;
	}




private:

	void init_filtration(){
		Vertex_handle v0 = *(complex_.vertex_range().begin());
		net_filtration_.reserve(complex_.num_vertices());
		for(auto v : complex_.vertex_range()){
			if(v != v0)
				net_filtration_.push_back(
						Net_filtration_vertex(v,
								Vertex_handle(-1),
								squared_eucl_distance(v,v0))
				);
		}
		net_filtration_.push_back(Net_filtration_vertex(v0,Vertex_handle(-1),1e10));
		auto n = net_filtration_.size();
		sort_filtration(n-1);
	}

	void update_radius_value(int k){
		int n = net_filtration_.size();
		int index_to_update = n-k;
		for(int i = 0; i< index_to_update; ++i){
			net_filtration_[i].radius = (std::min)(net_filtration_[i].radius ,
					squared_eucl_distance(
							net_filtration_[i].vertex_handle,
							net_filtration_[index_to_update].vertex_handle
					)
			);
		}
		sort_filtration(n-k);
	}

	/**
	 * sort all i first elements.
	 */
	void sort_filtration(int i){
		std::sort(net_filtration_.begin(),net_filtration_.begin()+ i);
	}

	double squared_eucl_distance(Vertex_handle v1,Vertex_handle v2) const{
		return std::sqrt(Geometry_trait::Squared_distance_d()(
				complex_.point(v1),complex_.point(v2))
		);
	}

	void print_filtration() const{
		for(auto v : net_filtration_){
			std::cerr <<"v="<<v.vertex_handle<<"-> d="<<v.radius<<std::endl;
		}
	}

};



#endif /* FURTHEST_POINT_EPSILON_NET_H_ */

/*
 * ComplexVisitor.h
 *
 *  Created on: Dec 11, 2013
 *      Author: dsalinas
 */

#ifndef GUDHI_COMPLEXVISITOR_H_
#define GUDHI_COMPLEXVISITOR_H_


#include "gudhi/Skeleton_blocker/Skeleton_blocker_simplex.h"

namespace Gudhi{

namespace skbl {
// todo rajouter les const

/**
 *@class Skeleton_blocker_complex_visitor
 *@brief Interface for a visitor of a simplicial complex.
 */
template <typename Vertex_handle>
class Skeleton_blocker_complex_visitor {
public:
	virtual ~Skeleton_blocker_complex_visitor(){};

	virtual void on_add_vertex(Vertex_handle) = 0;
	virtual void on_remove_vertex(Vertex_handle) = 0;

	virtual void on_add_edge(Vertex_handle a,Vertex_handle b) = 0;
	virtual void on_remove_edge(Vertex_handle a,Vertex_handle b) = 0;

	/**
	 * @brief Called when an edge changes, during the contraction of
	 * an edge
	 */
	virtual void on_changed_edge(Vertex_handle a,Vertex_handle b) = 0;

	/**
	 * @brief Called when performing an edge contraction when
	 * an edge bx is replaced by an edge ax (not already present).
	 * Precisely, this methods is called this way in contract_edge :
	 * add_edge(a,x)
	 * on_swaped_edge(a,b,x)
	 * remove_edge(b,x)
	 */
	virtual void on_swaped_edge(Vertex_handle a,Vertex_handle b,Vertex_handle x)=0;
	virtual void on_add_blocker(const Skeleton_blocker_simplex<Vertex_handle>&) = 0;
	virtual void on_delete_blocker(const Skeleton_blocker_simplex<Vertex_handle>*) = 0;
};


/**
 *@class Dummy_complex_visitor
 *@brief A dummy visitor of a simplicial complex that does nothing
 *
 */
template <typename Vertex_handle>
class Dummy_complex_visitor  : public Skeleton_blocker_complex_visitor<Vertex_handle>{
public:
	void on_add_vertex(Vertex_handle) {}
	void on_remove_vertex(Vertex_handle){}
	void on_add_edge(Vertex_handle a,Vertex_handle b){}
	void on_remove_edge(Vertex_handle a,Vertex_handle b){}
	void on_changed_edge(Vertex_handle a,Vertex_handle b){}
	void on_swaped_edge(Vertex_handle a,Vertex_handle b,Vertex_handle x){}
	void on_add_blocker(const Skeleton_blocker_simplex<Vertex_handle>&){}
	void on_delete_blocker(const Skeleton_blocker_simplex<Vertex_handle>*){}

};



/**
 *@class Print_complex_visitor
 *@brief A visitor that prints the visit information.
 *
 */
template <typename Vertex_handle>
class Print_complex_visitor  : public  Skeleton_blocker_complex_visitor<Vertex_handle>{
public:
	void on_add_vertex(Vertex_handle v) {
		std::cerr << "on_add_vertex:"<<v<<std::endl;
	}
	void on_remove_vertex(Vertex_handle v){
		std::cerr << "on_remove_vertex:"<<v<<std::endl;
	}
	void on_add_edge(Vertex_handle a,Vertex_handle b){
		std::cerr << "on_add_edge:"<<a<<","<<b<<std::endl;
	}
	void on_remove_edge(Vertex_handle a,Vertex_handle b){
		std::cerr << "on_remove_edge:"<<a<<","<<b<<std::endl;
	}
	void on_changed_edge(Vertex_handle a,Vertex_handle b){
		std::cerr << "on_changed_edge:"<<a<<","<<b<<std::endl;
	}
	void on_swaped_edge(Vertex_handle a,Vertex_handle b,Vertex_handle x){
			std::cerr << "on_swaped_edge:"<<a<<","<<b<<","<<x<<std::endl;
	}
	void on_add_blocker(const Skeleton_blocker_simplex<Vertex_handle>& b){
		std::cerr << "on_add_blocker:"<<b<<std::endl;
	}
	void on_delete_blocker(const Skeleton_blocker_simplex<Vertex_handle>* b){
		std::cerr << "on_delete_blocker:"<<b<<std::endl;
	}
};

}

}  // namespace GUDHI

#endif /* GUDHI_COMPLEXVISITOR_H_ */

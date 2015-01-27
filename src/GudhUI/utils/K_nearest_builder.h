/*
 * K_nearest_builder.h
 *
 *  Created on: Sep 10, 2014
 *      Author: dsalinas
 */

#ifndef K_NEAREST_BUILDER_H_
#define K_NEAREST_BUILDER_H_

#include <unordered_map>
#include <boost/iterator/iterator_facade.hpp>
#include <CGAL/Euclidean_distance.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_d.h>

#include "utils/UI_utils.h"
#include "model/Complex_typedefs.h"

template<typename SkBlComplex> class K_nearest_builder{
private:

	//todo uggh no virtual delete operator in Point
	// so do a composition asap
	class Point_d_with_id : public Point{
		typedef Point Base;
	public:
		Complex::Vertex_handle vertex_handle;
		Point_d_with_id(int d=0) : Point(d) {}
//		Point_d_with_id(int d, const Origin &o) : Base(d,o) {}

		Point_d_with_id(int a, int b, int c = 1) :
			Base(RT(a),RT(b),RT(c)) {}
		Point_d_with_id(int a, int b, int c, int d) :
			Base(RT(a),RT(b),RT(c),RT(d)) {}

		template <class InputIterator>
		Point_d_with_id (int d, InputIterator first, InputIterator last)
		: Base (d, first, last) {}
		template <class InputIterator>

		Point_d_with_id(const Point_d_with_id &p) : Base(p) {}
		Point_d_with_id(const Base& p) : Base(p) {}
		Point_d_with_id(const Base& p,Complex::Vertex_handle v) : Base(p),vertex_handle(v)	{}

	};

	struct Kernel_with_id : public Geometry_trait{
		typedef Point_d_with_id Point_d;
	};


	/**
	 * TODO wrap vertex handle in class passed to tree
	 */
	typedef Kernel_with_id K;
	typedef typename K::Point_d Point_d;
	typedef typename CGAL::Search_traits_d<K> TreeTraits;
	typedef typename CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
	typedef typename Neighbor_search::Tree Tree;

	SkBlComplex& complex_;
public:

	/**
	 * @brief Modify complex to be the expansion of the k-nearest neighbor
	 * symetric graph.
	 */
	K_nearest_builder(SkBlComplex& complex,unsigned k):complex_(complex){
		complex.keep_only_vertices();
		compute_edges(k);
	}

private:



	double squared_eucl_distance(const Point& p1,const Point& p2) const{
		return Geometry_trait::Squared_distance_d()(p1,p2);
	}

	void compute_edges(unsigned k){

		std::list<Point_d_with_id> points_with_id;
		for(auto v: complex_.vertex_range()){
			Point_d_with_id point_v(complex_.point(v),v);
			points_with_id.push_back(point_v);
		}

		Tree tree(points_with_id.begin(),points_with_id.end());


		for (auto p : complex_.vertex_range()){
			Neighbor_search search(tree, complex_.point(p),k+1);
			for(auto it = ++search.begin(); it != search.end(); ++it){
				auto q = it->first.vertex_handle;
				if (p != q && complex_.contains_vertex(p) && complex_.contains_vertex(q))
					complex_.add_edge(p,q);
			}
		}
	}


};

#endif /* K_NEAREST_BUILDER_H_ */

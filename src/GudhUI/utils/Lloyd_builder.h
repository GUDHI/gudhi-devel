/*
 * Lloyd.h
 *
 *  Created on: Sep 25, 2014
 *      Author: dsalinas
 */

#ifndef LLOYD_H_
#define LLOYD_H_


#include <vector>

/**
 * Iteratively puts every vertex at the center of its neighbors
 */
template<typename SkBlComplex> class Lloyd_builder{
private:
	SkBlComplex& complex_;
	int dim;
public:
	typedef typename SkBlComplex::Vertex_handle Vertex_handle;
	/**
	 * @brief Modify complex to be the expansion of the k-nearest neighbor
	 * symetric graph.
	 */
	Lloyd_builder(SkBlComplex& complex,unsigned num_iterations):complex_(complex),dim(-1){
		if(!complex_.empty()){
			dim = get_dimension();
			while(num_iterations--){
				std::list<Point> new_points;
				for(auto v : complex.vertex_range())
					new_points.push_back(barycenter_neighbors(v));

				auto new_points_it = new_points.begin();
				for(auto v : complex.vertex_range())
					complex_.point(v) = *(new_points_it++);
			}
		}
	}

private:

	int get_dimension(){
		assert(!complex_.empty());
		for(auto v : complex_.vertex_range())
			return complex_.point(v).dimension();
		return -1;
	}

	Point barycenter_neighbors(Vertex_handle v) const{
		if(complex_.degree(v)==0)
			return complex_.point(v);

		std::vector<double> res(dim,0);
		unsigned num_points = 0;
		for(auto nv : complex_.vertex_range(v)){
			++num_points;
			const Point& point = complex_.point(nv);
			assert(point.dimension() == dim);
			for(int i = 0; i < point.dimension() ; ++i)
				res[i] += point[i];
		}
		for(auto& x : res)
			x/= num_points;
		return Point(dim,res.begin(),res.end());
	}



	double squared_eucl_distance(const Point& p1,const Point& p2) const{
		return Geometry_trait::Squared_distance_d()(p1,p2);
	}

};


#endif /* LLOYD_H_ */

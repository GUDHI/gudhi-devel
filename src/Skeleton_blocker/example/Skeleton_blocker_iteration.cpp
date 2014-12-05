#include <boost/timer/timer.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <sstream>


#include "gudhi/Skeleton_blocker.h"
//#include "gudhi/Skeleton_blocker_complex.h"
//#include "gudhi/Skeleton_blocker/Skeleton_blocker_simple_traits.h"

using namespace std;
using namespace Gudhi;
using namespace skbl;

typedef Skeleton_blocker_complex<Skeleton_blocker_simple_traits> Complex;
typedef Complex::Vertex_handle Vertex_handle;
typedef Complex::Simplex_handle Simplex;


Complex build_complete_complex(int n){
	// build a full complex with 10 vertices and 2^n-1 simplices
		Complex complex;
		for(int i=0;i<n;i++)
			complex.add_vertex();
		for(int i=0;i<n;i++)
			for(int j=0;j<i;j++)
				//note that add_edge, add the edge and all its cofaces
				complex.add_edge(Vertex_handle(i),Vertex_handle(j));
		return complex;
}

int main (int argc, char *argv[]){
	boost::timer::auto_cpu_timer t;

	const int n = 15;

	// build a full complex with 10 vertices and 2^n-1 simplices
	Complex complex(build_complete_complex(n));

	// this is just to illustrate iterators, to count number of vertices
	// or edges, complex.num_vertices() and complex.num_edges() are
	// more appropriated!
	unsigned num_vertices = 0;
	for(auto v : complex.vertex_range()){
		if(v==v); //do something with v as removing a ennoying warning.
		++num_vertices;
	}

	unsigned num_edges = 0;
	for(auto e : complex.edge_range()){
		if(e==e);
		++num_edges;
	}

	unsigned euler = 0;
	unsigned num_simplices = 0;
	// we use a reference to a simplex instead of a copy
	// value here because a simplex is a set of integers
	// and copying it cost time
	for(const Simplex & s : complex.simplex_range()){
		++num_simplices;
		if(s.dimension()%2 == 0) 
			euler += 1;
		else 
			euler -= 1;
	}
	std::cout << "Saw "<<num_vertices<<" vertices, "<<num_edges<<" edges and "<<num_simplices<<" simplices"<<std::endl;
	std::cout << "The Euler Characteristic is "<<euler<<std::endl;
	return EXIT_SUCCESS;
}


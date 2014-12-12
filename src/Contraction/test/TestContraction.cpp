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
#include <ctime>
#include <list>

#include "combinatorics/Skeleton_blocker/Skeleton_blocker_simple_traits.h"
#include "geometry/Skeleton_blocker_simple_geometric_traits.h"
//#include "Skeleton_blocker/Simplex.h"
#include "contraction/Skeleton_blocker_contractor.h"
#include "Utils.h"
#include "iofile.h"
#include "Test.h"
#include "Skeleton_blocker_geometric_complex.h"
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>


#include "contraction/Edge_profile.h"

#include "contraction/policies/Cost_policy.h"
#include "contraction/policies/Edge_length_cost.h"
#include "contraction/policies/Placement_policy.h"
#include "contraction/policies/Middle_placement.h"

#include "contraction/policies/Valid_contraction_policy.h"
#include "contraction/policies/Dummy_valid_contraction.h"
#include "contraction/policies/Link_condition_valid_contraction.h"


using namespace std;

using namespace Gudhi;

using namespace skbl;

struct Geometry_trait{
	typedef std::vector<double> Point;
};

typedef Geometry_trait::Point Point;


typedef Skeleton_blocker_complex<Skeleton_blocker_simple_traits> AbstractComplex;
typedef Skeleton_blocker_simple_geometric_traits<Geometry_trait> Complex_geometric_traits;


typedef Skeleton_blocker_geometric_complex< Complex_geometric_traits > Complex;

typedef Complex::Vertex_handle Vertex_handle;
typedef Complex::Simplex_handle Simplex_handle;

typedef Complex::Root_vertex_handle Root_vertex_handle;

using namespace contraction;

typedef Skeleton_blocker_contractor<Complex> Complex_contractor;

typedef Edge_profile<Complex> Profile;

// compute the distance todo utiliser Euclidean_geometry a la place
template<typename Point>
double eucl_distance(const Point& a,const Point& b){
	double res = 0;
	auto a_coord = a.begin();
	auto b_coord = b.begin();
	for(; a_coord != a.end(); a_coord++, b_coord++){
		res += (*a_coord - *b_coord) * (*a_coord - *b_coord);
	}
	return sqrt(res);
}

// build the Rips complex todo utiliser Euclidean_geometry a la place
template<typename ComplexType>
void build_rips(ComplexType& complex, double offset){
	if (offset<=0) return;
	auto vertices = complex.vertex_range();
	for (auto p = vertices.begin(); p != vertices.end(); ++p)
		for (auto q = p; ++q != vertices.end(); /**/)
		if (eucl_distance(complex.point(*p), complex.point(*q)) < 2 * offset){
				complex.add_edge(*p,*q);
			}
}




void test_contraction_rips(string name_file, double offset){
	Complex complex;
	// load the points
	if (!read_off_file<Complex>(name_file,complex,true)){
		std::cerr << "Unable to read file:"<<name_file<<std::endl;
		std::cerr << "current path : ";
		system("pwd");
		std::cerr<<endl;
		return;
	}

	clock_t time = clock();

	TEST("build the Rips complex");

	build_rips(complex,offset);

	std::cerr << "Rips contruction took "<< ( (float)(clock()-time))/CLOCKS_PER_SEC << " seconds\n";

	TESTMSG("Initial number of vertices :",complex.num_vertices());
	TESTMSG("Initial number of edges :",complex.num_edges());
	TESTMSG("Initial number of blockers:",complex.num_blockers());

	time = clock();

	Complex_contractor contractor(complex,
			new Edge_length_cost<Profile>,
			contraction::make_first_vertex_placement<Profile>(),
			contraction::make_link_valid_contraction<Profile>(),
			contraction::make_remove_popable_blockers_visitor<Profile>());
	contractor.contract_edges();

	TESTVALUE(complex.to_string());

	TESTVALUE(complex.num_vertices());
	TESTVALUE(complex.num_edges());
	TESTVALUE(complex.num_blockers());

	std::cerr << "Edge contractions took "<< ( (float)(clock()-time))/CLOCKS_PER_SEC << " seconds\n";

}


void test_geometric_link(){

	Complex complex;
	std::vector<double> p0(2,0);
	std::vector<double> p1(2,0); p1[0] = 1.;
	std::vector<double> p2(2,1);
	complex.add_vertex(p0);
	complex.add_vertex(p1);
	complex.add_vertex(p2);

	complex.add_edge(Vertex_handle(0),Vertex_handle(1));
	complex.add_edge(Vertex_handle(1),Vertex_handle(2));
	complex.add_edge(Vertex_handle(2),Vertex_handle(0));



	cerr << "complex points:" <<endl;
	for(auto v : complex.vertex_range()){
		cerr <<v <<" -> ";
		DBGCONT(complex.point(v));
	}

	cerr << "complex : "<<complex.to_string()<<endl;




	auto link = complex.link(Vertex_handle(0));


	cerr << "link of 0 points:" <<endl;
	for(auto v : link.vertex_range()){
		cerr <<v <<" -> ";
		DBGCONT(link.point(v));
	}

	cerr << "link : "<<link.to_string()<<endl;
}




int main (int argc, char *argv[])
{
	if (argc!=3){
		std::cerr << "Usage "<<argv[0]<<" GUDHIPATH/src/data/sphere3D.off 0.1 to load the file GUDHIPATH/src/data/sphere3D.off and contract the Rips complex built with paremeter 0.2.\n";
		return -1;
	}

	std::string name_file(argv[1]);
	test_contraction_rips(name_file,atof(argv[2]));
}



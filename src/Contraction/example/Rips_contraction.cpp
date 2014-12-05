/*
 * Rips_contraction.cpp
 *
 *  Created on: Nov 28, 2014
 *      Author: dsalinas
 */
#include <list>
#include <boost/timer/timer.hpp>
#include <iostream>
#include "gudhi/Edge_contraction.h"
#include "gudhi/Skeleton_blocker.h"
#include "gudhi/Off_reader.h"


using namespace std;
using namespace Gudhi;
using namespace skbl;
using namespace contraction;

struct Geometry_trait{
	typedef std::vector<double> Point;
};

typedef Geometry_trait::Point Point;
typedef Skeleton_blocker_simple_geometric_traits<Geometry_trait> Complex_geometric_traits;
typedef Skeleton_blocker_geometric_complex< Complex_geometric_traits > Complex;
typedef Edge_profile<Complex> Profile;
typedef Skeleton_blocker_contractor<Complex> Complex_contractor;

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

template<typename ComplexType>
void build_rips(ComplexType& complex, double offset){
	if (offset<=0) return;
	auto vertices = complex.vertex_range();
	for (auto p = vertices.begin(); p != vertices.end(); ++p)
		for (auto q = p; ++q != vertices.end(); /**/)
			if (eucl_distance(complex.point(*p), complex.point(*q)) < 2 * offset)
				complex.add_edge(*p,*q);
}




void test_contraction_rips(string name_file, double offset){
	boost::timer::auto_cpu_timer t;
	Complex complex;

	// load the points
	Skeleton_blocker_off_reader<Complex> off_reader(name_file,complex,true);
	if(!off_reader.is_valid()){
		std::cerr << "Unable to read file:"<<name_file<<std::endl;
		return;
	}
	std::cerr << "build the Rips complex"<<std::endl;

	build_rips(complex,offset);


	std::cerr << "Initial complex has "<<
			complex.num_vertices()<<" vertices, and "<<
			complex.num_edges()<<" edges."<<std::endl;

	Complex_contractor contractor(complex,
			new Edge_length_cost<Profile>, //todo make_edge_length_cost
			contraction::make_first_vertex_placement<Profile>(),
			contraction::make_link_valid_contraction<Profile>(),
			contraction::make_remove_popable_blockers_visitor<Profile>());
	contractor.contract_edges();

	std::cerr << "Resulting complex has "<<
			complex.num_vertices()<<" vertices, "<<
			complex.num_edges()<<"edges and "<<
			complex.num_blockers()<<" blockers"<<std::endl;


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



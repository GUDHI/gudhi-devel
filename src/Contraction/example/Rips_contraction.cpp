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

int main (int argc, char *argv[])
{
	if (argc!=3){
		std::cerr << "Usage "<<argv[0]<<" ../../../data/meshes/SO3_10000.off 0.3 to load the file ../../data/SO3_10000.off and contract the Rips complex built with paremeter 0.3.\n";
		return -1;
	}

	Complex complex;

	// load the points
	Skeleton_blocker_off_reader<Complex> off_reader(argv[1],complex,true);
	if(!off_reader.is_valid()){
		std::cerr << "Unable to read file:"<<argv[1]<<std::endl;
		return EXIT_FAILURE;
	}
	std::cout << "Build the Rips complex"<<std::endl;

	build_rips(complex,atof(argv[2]));

	boost::timer::auto_cpu_timer t;

	std::cout << "Initial complex has "<<
			complex.num_vertices()<<" vertices and "<<
			complex.num_edges()<<" edges"<<std::endl;

	Complex_contractor contractor(complex,
			new Edge_length_cost<Profile>,
			contraction::make_first_vertex_placement<Profile>(),
			contraction::make_link_valid_contraction<Profile>(),
			contraction::make_remove_popable_blockers_visitor<Profile>());
	contractor.contract_edges();

	std::cout << "Counting final number of simplices \n";
	unsigned num_simplices = std::distance(complex.simplex_range().begin(),complex.simplex_range().end());

	std::cout << "Final complex has "<<
			complex.num_vertices()<<" vertices, "<<
			complex.num_edges()<<" edges, "<<
			complex.num_blockers()<<" blockers and "<<
			num_simplices<<" simplices"<<std::endl;


	std::cout << "Time to simplify and enumerate simplices:\n";

	return EXIT_SUCCESS;
}



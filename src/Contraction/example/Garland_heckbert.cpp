/*
 * Garland_heckbert.h
 *  Created on: Feb 10, 2015
 * This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014  INRIA Sophia Antipolis-Méditerranée (France)
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
 * 
 */


#ifndef GARLAND_HECKBERT_H_
#define GARLAND_HECKBERT_H_


#include <iostream>
#include "Garland_heckbert/Point.h"
#include "Garland_heckbert/Error_quadric.h"
#include "Garland_heckbert/Garland_heckbert_policies.h"

#include "gudhi/Edge_contraction.h"
#include "gudhi/Skeleton_blocker.h"
#include "gudhi/Off_reader.h"

using namespace std;
using namespace Gudhi;
using namespace skbl;
using namespace contraction;


struct Geometry_trait{
	typedef Point_d Point;
};

struct Garland_heckbert_traits : public Skeleton_blocker_simple_geometric_traits<Geometry_trait> {

public:
	// the vertex stored in the complex contains a quadric
	struct Garland_heckbert_vertex : public Simple_geometric_vertex{
		Error_quadric<Geometry_trait::Point> quadric;
	};
	typedef Garland_heckbert_vertex Graph_vertex;

};

typedef Skeleton_blocker_geometric_complex< Garland_heckbert_traits > Complex;
typedef Edge_profile<Complex> Profile;
typedef Skeleton_blocker_contractor<Complex> Complex_contractor;




int main(int argc, char *argv[]){
	if (argc!=3){
		std::cerr << "Usage "<<argv[0]<<" ../../../data/meshes/test.off N to load the file ../../data/test.off and contract N edges.\n";
		return -1;
	}

	Complex complex;

	// load the points
	Skeleton_blocker_off_reader<Complex> off_reader(argv[1],complex);
	if(!off_reader.is_valid()){
		std::cerr << "Unable to read file:"<<argv[1]<<std::endl;
		return EXIT_FAILURE;
	}

	int num_contractions = atoi(argv[2]);

	std::cout << "Initial complex has "<<
			complex.num_vertices()<<" vertices, "<<
			complex.num_blockers()<<" blockers, "<<
			complex.num_edges()<<" edges and" <<
			complex.num_triangles()<<" triangles.";

	Complex_contractor contractor(complex,
			new GH_cost<Profile>(complex),
			new GH_placement<Profile>(complex),
			contraction::make_link_valid_contraction<Profile>(),
			new GH_visitor<Profile>(complex)
	);

	std::cout<<"Contract "<<num_contractions<<" edges\n";
	contractor.contract_edges(num_contractions);

	std::cout << "Final complex has "<<
			complex.num_vertices()<<" vertices, "<<
			complex.num_edges()<<" edges and" <<
			complex.num_triangles()<<" triangles.";

	Skeleton_blocker_off_writer<Complex> off_writer("reduced.off",complex);


	return EXIT_SUCCESS;
}



#endif /* GARLAND_HECKBERT_H_ */

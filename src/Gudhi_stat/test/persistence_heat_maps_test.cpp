/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2015  INRIA Sophia-Saclay (France)
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



#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "gudhi_stat"
#include <boost/test/unit_test.hpp>
#include <gudhi/reader_utils.h>
#include <gudhi/abstract_classes/Abs_Topological_data.h>
#include <gudhi/concretizations/Persistence_heat_maps.h>

#include <iostream>



using namespace Gudhi;
using namespace Gudhi::Gudhi_stat;

using namespace std;


	

BOOST_AUTO_TEST_CASE(check_construction_of_landscape) 
{	
	std::vector< std::vector<double> > filter = create_Gaussian_filter(100,1);		
	Persistence_heat_maps p( "data/file_with_diagram" , filter ,  constant_function, false , 1000 , 0 , 1 );
	p.write_to_file( "data/persistence_heat_map_from_file_with_diagram" );
	
	Persistence_heat_maps q;
	q.load_from_file( "data/persistence_heat_map_from_file_with_diagram" );	
	
	
	
	BOOST_CHECK( p == q );
}

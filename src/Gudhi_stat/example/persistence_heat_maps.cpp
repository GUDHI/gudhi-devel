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



#include <gudhi/reader_utils.h>
#include <gudhi/abstract_classes/Abs_Topological_data.h>
#include <gudhi/concretizations/Persistence_heat_maps.h>

#include <iostream>



using namespace Gudhi;
using namespace Gudhi::Gudhi_stat;


double epsilon = 0.0000005;

using namespace std;

	
int main( int argc , char** argv )
{
	//if ( argc != 2 )
	//{
	//	cout << "To run this program, please provide the name of a file with persistence diagram \n";
	//	return 1;
	//}
	
/*
	std::vector< std::pair< double,double > > intervals;
	intervals.push_back( std::make_pair(0.5,0.5) );	
	std::vector< std::vector<double> > filter = create_Gaussian_filter(5,1);		
	Persistence_heat_maps p( intervals , filter ,  constant_function, false , 100 , 0 , 1 );
	p.plot( "heat_map_1" );

	
	std::vector< std::pair< double,double > > intervals2;
	intervals2.push_back( std::make_pair(7,12) );		
	Persistence_heat_maps q( intervals2 , filter ,  constant_function, false , 100 , 0 , 10 );
	q.plot( "heat_map_2" );
*/
/*
	std::vector< std::pair< double,double > > intervals;
	intervals.push_back( std::make_pair(0.5,0.5) );	
	std::vector< std::vector<double> > filter = create_Gaussian_filter(5,1);		
	Persistence_heat_maps p( intervals , filter ,  constant_function, false , 10 , 0 , 1 );
	p.write_to_file( "aaa" );
	
	Persistence_heat_maps q;
	q.load_from_file( "aaa" );
	
	cerr << ( p == q ) << endl;
	*/

/*
	std::vector< std::vector<double> > filter = create_Gaussian_filter(30,1);		
	Persistence_heat_maps p( "file_with_diagram" , filter ,  constant_function, false , 100 , 0 , 1 );
	p.plot( "heat_map_1" );
*/ 
	 
/*   
   //test to construct persistence heat map:
    std::vector< std::vector<double> > filter = create_Gaussian_filter(100,1);		
	Persistence_heat_maps p( "file_with_diagram" , filter ,  constant_function, false , 1000 , 0 , 1 );
	p.write_to_file( "persistence_heat_map_from_file_with_diagram" );
	
	Persistence_heat_maps q;
	q.load_from_file( "persistence_heat_map_from_file_with_diagram" );	
	
	cerr << (p == q) << endl;
*/	
	
	

	
	

	
	return 0;
}








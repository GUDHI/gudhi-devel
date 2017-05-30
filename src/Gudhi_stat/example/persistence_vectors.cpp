/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2015  INRIA (France)
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



#include <gudhi/Persistence_vectors.h>
#include <iostream>


#include <gudhi/reader_utils.h>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace Gudhi;
using namespace Gudhi::Persistence_representations;




int main( int argc , char** argv )
{
	//create two simple vectors with birth--death pairs:
	
	std::vector< std::pair< double , double > > persistence1;
	std::vector< std::pair< double , double > > persistence2;
	
	persistence1.push_back( std::make_pair(1,2) );
	persistence1.push_back( std::make_pair(6,8) );
	persistence1.push_back( std::make_pair(0,4) );
	persistence1.push_back( std::make_pair(3,8) );
	
	persistence2.push_back( std::make_pair(2,9) );
	persistence2.push_back( std::make_pair(1,6) );
	persistence2.push_back( std::make_pair(3,5) );
	persistence2.push_back( std::make_pair(6,10) );
	
	//create two persistence vectors based on persistence1 and persistence2:
	Vector_distances_in_diagram<Euclidean_distance > v1 = Vector_distances_in_diagram<Euclidean_distance >( persistence1 , std::numeric_limits<size_t>::max() );
	Vector_distances_in_diagram<Euclidean_distance > v2 = Vector_distances_in_diagram<Euclidean_distance >( persistence2 , std::numeric_limits<size_t>::max() );
	
	//writing to a stream:
	std::cout << "v1 : " << v1 << std::endl;
	std::cout << "v2 : " << v2 << std::endl;
	
	//averages:
	Vector_distances_in_diagram<Euclidean_distance > average;
	average.compute_average( {&v1,&v2} );
	std::cout << "Average : " << average << std::endl;
	
	//computations of distances:
	std::cout << "l^1 distance : " << v1.distance( v2 ) << std::endl;
    
   //computations of scalar product:
    std::cout << "Scalar product of l1 and l2 : " << v1.compute_scalar_product( v2 ) << std::endl;
    
    //create a file with a gnuplot script:
    v1.plot( "plot_of_vector_representation" );

	return 0;
}

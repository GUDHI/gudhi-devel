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



#include <gudhi/concretizations/Vector_distances_in_diagram.h>
#include <iostream>


#include <gudhi/reader_utils.h>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace Gudhi;
using namespace Gudhi::Gudhi_stat;

using namespace std;


double epsilon = 0.000005;



int main( int argc , char** argv )
{
	if ( argc < 2 )
	{
		cout << "To run this program, please provide the name of a file with persistence diagram. If you provide two files, we will do distance, scalar produc and average computations \n";
		return 1;
	}
	
	Vector_distances_in_diagram< euclidean_distance<double> > p( argv[1] , 100 );
	cout << "This is a vector corresponding to the input persistence diagram : \n";	
	cout << p << endl;



	if ( argc == 3 )
	{
		Vector_distances_in_diagram< euclidean_distance<double> > p_prime( argv[2] , 100);	
		
		cout << "p_prime : " <<p_prime << endl;
			
		cout << "Distance between input persistence diagrams : " << p.distance( (Abs_Topological_data_with_distances*)(&p_prime) ) << endl;	
		std::vector< Abs_Topological_data_with_averages* > to_average;
		to_average.push_back( (Abs_Topological_data_with_averages*)(&p) );
		to_average.push_back( (Abs_Topological_data_with_averages*)(&p_prime) );
	
		Vector_distances_in_diagram< euclidean_distance<double> > average; 
		average.compute_average( to_average );
	
		cout << "Here is an average : " << average << endl;
	}
	return 0;
}



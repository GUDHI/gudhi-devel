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


#include <gudhi/permutation_test.h>
#include <gudhi/persistence_representations/Persistence_landscape.h>
#include <iostream>
#include <cstring>

using namespace Gudhi;
using namespace Gudhi::Gudhi_stat;

int main( int argc , char** argv )
{	
	
	std::cout << "This program require the following parameters: \n";
	std::cout << "(1-2) Names of files each of them contains the names of files with persistence diagrams. The diagrams from a single file are assumed to be in the same group \n";
	std::cout << "Third parameter is an integer being the number of permutations to be made \n";
	std::cout << "The last parameter is a double, the power of a distance \n";
	if ( argc != 5 )
	{
		std::cout << "Wrong number of parameters, the program will now terminat \n";
		return 1;
	}
	std::cout << "We will now read the data from files : " << argv[1] << " and " << argv[2] << std::endl;
	size_t number_of_permutations = (size_t)(atoi(argv[3]));
	size_t exponent = (double)atof( argv[4] );

	std::vector< std::string > first_group = readFileNames( argv[1] );
	std::vector< std::string > second_group =readFileNames( argv[2] );
	
	std::cout << "Here are the filenames in the first group :\n";
	for ( size_t i = 0 ; i != first_group.size() ; ++i )
	{
		std::cout << first_group[i] << std::endl;
	}
	std::cout << "Here are the filenames in the second group :\n";
	for ( size_t i = 0 ; i != second_group.size() ; ++i )
	{
		std::cout << second_group[i] << std::endl;
	}
		
	std::vector< Persistence_landscape* > first_collection( first_group.size() );
	for ( size_t i = 0 ; i != first_group.size() ; ++i )
	{
		std::vector< std::pair< double , double > > diag = read_persistence_intervals_in_one_dimension_from_file( first_group[i].c_str() );
		Persistence_landscape* l = new Persistence_landscape( diag );			
		first_collection[i] = l;		
	}
	
	std::vector< Persistence_landscape* > second_collection( second_group.size() );
	for ( size_t i = 0 ; i != second_group.size() ; ++i )
	{
		std::vector< std::pair< double , double > > diag = read_standard_persistence_file( second_group[i].c_str() );
		Persistence_landscape* l = new Persistence_landscape( diag );		
		second_collection[i] = l;
	}	
		
	std::cout << "The p-value form a permutation test is : " << permutation_test<Persistence_landscape>( first_collection , second_collection , number_of_permutations , exponent  ) << std::endl;
	

	
	return 0;
}

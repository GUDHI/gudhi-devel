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

#ifndef PERMUTATION_TEST_H
#define PERMUTATION_TEST_H


//concretizations
#include <vector>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <gudhi/concretizations/read_persitence_from_file.h>

using namespace std;

namespace Gudhi 
{
namespace Gudhi_stat 
{

template <typename Representation_of_persistence>
double permutation_test( const std::vector<Representation_of_persistence*>& data_1 , const std::vector<Representation_of_persistence*>& data_2 , size_t number_of_permutations , double exponent = 1 )
{	
	try
	{
		Representation_of_persistence* av_1 = new Representation_of_persistence;		
		av_1->compute_average( data_1 );						
		Representation_of_persistence* av_2 = new Representation_of_persistence;			
	    av_2->compute_average( data_2 );	    	
		double initial_distance = av_1->distance( *av_2 , exponent );
		double counter = 0;
		
		for  ( size_t i = 0 ; i != number_of_permutations ; ++i )
		{					
			std::vector<Representation_of_persistence*> all_data;
			all_data.insert(all_data.end() , data_1.begin() , data_1.end() );
			all_data.insert(all_data.end() , data_2.begin() , data_2.end() );
			std::random_shuffle( all_data.begin() , all_data.end() );
			
			std::vector<Representation_of_persistence*> first_part(&all_data[0],&all_data[data_1.size()]);
			std::vector<Representation_of_persistence*> second_part(&all_data[data_1.size()],&all_data[all_data.size()]);
		
			av_1->compute_average( first_part );
			av_2->compute_average( second_part );			
			double distance_after_permutations = av_1->distance( *av_2 , exponent );
			if ( distance_after_permutations > initial_distance )++counter;
		}
		return counter / (double)number_of_permutations;
	}
	catch (...)
	{
		std::cout << "The data structure do not support the operations that are neccessay for a permutation test (averaging, distances) \n";
	}
	return 0;
}//permutation_test

}//Gudhi_stat
}//Gudhi

#endif

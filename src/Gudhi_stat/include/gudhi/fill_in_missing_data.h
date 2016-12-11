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

#ifndef FILL_IN_MISSING_DATA_H
#define FILL_IN_MISSING_DATA_H

namespace Gudhi
{
namespace Gudhi_stat
{

/**
 * Quite often in biological sciences we are facing a problem of missing data. We may have for instance a number of sequences of observations made in between times A and B in a discrete 
 * collection of times A = t1, t2,...,tn = B. But quite typically some of the observations may be missing. Then quite often it is hard to estimate the values in the missing times. 
 * The procedure below assumes that we compute some topological descriptor of the observations we are given. Most typically this will be some type of persistence homology representation.
 * The the values in the missing points are filled in by linear interpolation and extrapolation. The procedure below have minimal requirements: it assumes that we have two elements that 
 * filled in. The rest will be filled in by using them.
 * The fill-in process is done based on the idea of linear approximation. Let us assume that we have two positions A and B which are filled in with proper objects of a type Representation_of_topology.
 * Any intermediate time step can be interpolated by taking a linear combination t*A + (1-t)*B, where t \in [0,1]. 
 * If the missing data are located at the beginning of at the end of vector of Representation_of_topology, then we use the following extrapolation scheme:
 * First we make sure that there are no missing data except at the very beginning and at the very end of a vector of data. Then, we pick two constitutive closest filled-in data A and B and we extrapolate
 * by using a formula: t*A-(t-1)*B, where t is a natural number > 1. 
 * Note that both vector of data and the vector is_the_position_filled_in are modified by this procedure. Upon successful termination of the procedure, the vector is_the_position_filled_in has only 'true' entries,
 * and the vector data do not have missing data. 
**/  

template < typename Representation_of_topology >
void fill_in_missing_data( std::vector< Representation_of_topology* >& data , std::vector< bool >& is_the_position_filled_in )
{
	bool dbg = false;
	
	//first check if at least two positions are filled in:
	size_t number_of_positions_that_are_filled_in = 0;
	for ( size_t i = 0 ; i != is_the_position_filled_in.size() ; ++i )
	{
		if ( is_the_position_filled_in[i] )++number_of_positions_that_are_filled_in;
	}
	
	if ( number_of_positions_that_are_filled_in < 2 )
	{
		std::cerr << "There are too few positions filled in to do extrapolation / interpolation. The program will now terminate.\n";
		throw "There are too few positions filled in to do extrapolation / interpolation. The program will now terminate.\n";
	}
	
	for ( size_t i = 0 ; i != data.size() ; ++i )
	{
		if ( !is_the_position_filled_in[i] )
		{				
			//This position is not filled in. Find the next position which is nonzero.
			size_t j = 1;
			while ( (is_the_position_filled_in[i+j] == false) && ( i+j != data.size() ) )++j;
			if ( dbg )
			{
				std::cout << "The position number : " << i << " is not filled in. The next filled-in position is : " << i+j << std::endl;
			}
			if ( i != 0 )
			{
				//this is not the first position of the data:
				if ( i + j != data.size() )
				{
					//this is not the last position of the data either
					for ( size_t k = 0 ; k != j ; ++k )
					{
						double weight1 = double(j-k)/(double)(j+1);
						double weight2 = double(k+1)/(double)(j+1);							
						data[i+k] = new Representation_of_topology(weight1*(*data[i-1]) + weight2*(*data[i+j]));
						is_the_position_filled_in[i+k] = true;
						if ( dbg )
						{
							std::cerr << "We fill in a position : " << i+k << " with: position " << i-1 << " with weight " << weight1 << " and position " << i+j << " with weight " << weight2 << std::endl;
						}
					}
				}					
			}
			else
			{
				//this is the first position of the data, i.e. i == 0.
				while ( is_the_position_filled_in[i] == 0 )++i;
			}
			
		}
	}					
		
	
	if ( is_the_position_filled_in[0] == false )
	{
		//find the first nonzero (then, we know that the second one will be nonzero too, since we filled it in above):
		size_t i = 0;
		while ( is_the_position_filled_in[i] == false )++i;
		//the data at position i is declared. Since, we made sure that all other positions, except maybe a few first and a few last, are declared too. Therefore, if we find a 
		//first declared, then the next one will be declared too. So, we can do a telescopic declaration backward and forward.
				
		for ( size_t j = i ; j != 0 ; --j )
		{
			data[j-1] = new Representation_of_topology(2*(*data[j]) + (-1)*(*data[j+1]));
			is_the_position_filled_in[j-1] = true;
			if ( dbg )
			{
				std::cerr << "Filling in a position : " << j-1 << " by using: " << j << " and " << j+1 <<std::endl;
			}
		}
	}
	
	if ( is_the_position_filled_in[ is_the_position_filled_in.size()-1 ] == false )
	{
		//first the first nonzero (then, we know that the second one will be nonzero too):
		size_t i = is_the_position_filled_in.size()-1;
		while ( is_the_position_filled_in[i] == false )--i;
		
		for ( size_t j = i+1 ; j != data.size() ; ++j )
		{
			data[j] = new Representation_of_topology(2*(*data[j-1]) + (-1)*(*data[j-2]));
			is_the_position_filled_in[j] = true;
			if ( dbg )
			{
				std::cerr << "Filling in a position : " << j << " by using: " << j-1 << " and " << j-2 <<std::endl;
			}
		}
	}
}//fill_in_missing_data

}//namespace Gudhi_stat
}//namespace Gudhi

#endif

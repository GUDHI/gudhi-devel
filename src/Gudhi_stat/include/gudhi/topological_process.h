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

#ifndef TOPOLOGICAL_PROCESS_H
#define TOPOLOGICAL_PROCESS_H


//concretizations
#include <gudhi/concretizations/Vector_distances_in_diagram.h>
#include <gudhi/concretizations/Persistence_landscape.h>
#include <gudhi/concretizations/Persistence_landscape_on_grid.h>
#include <vector>

//extras
#include <gudhi/common_gudhi_stat.h>

namespace Gudhi 
{
namespace Gudhi_stat 
{
	

template <typename Representation>
std::vector< Representation* > construct_representation_from_file( const char* filename )	
{
	bool dbg = false;
	std::vector< std::string > files = readFileNames( filename );	
	
	std::cout << "Here are the filenames in the file : " << filename << std::endl;
	for ( size_t i = 0 ; i != files.size() ; ++i )
	{
		std::cout << files[i] << std::endl;
	}
	
	std::vector< Representation* > result( files.size() );
	for ( size_t i = 0 ; i != files.size() ; ++i )
	{
		std::vector< std::pair< double , double > > diag = read_standard_file( files[i].c_str() );	
		
		if ( dbg )
		{
			std::cerr << "Here is a diagram from a file : " << files[i].c_str() << std::endl;
			for ( size_t aa = 0  ; aa != diag.size() ; ++aa )
			{
				std::cout << diag[aa].first << " " << diag[aa].second << std::endl;
			}
			getchar();
		}
						
		Representation* l = new Representation( diag );			
		result[i] = l;		
	}
	return result;
}
	

template <typename Representation>
class Topological_process
{
public:
	Topological_process();
	Topological_process( const std::vector< Representation* >& data_ ):data(data_){}
	double distance( const Topological_process& second )
	{
		if ( this->data.size() != second.data.size() )
		{
			throw "Incompatible lengths of topological processes, we cannot compute the distance in this case \n";
		}
		double result = 0;
		for ( size_t i = 0 ; i != this->data.size() ; ++i )
		{
			result += this->data[i]->distance( *second.data[i] );
		}
		return result;
	}	
	
	void compute_average( const std::vector< Representation* >& to_average )
	{
		//since we will substitute whatever data we have in this object with an average, we clear the data in this object first:
		this->data.clear();
		//then we need to check if all the topological processes in the vector to_average have the same length.
		if ( to_average.size() == 0 )return;
		for ( size_t i = 1 ; i != to_average.size() ; ++i )
		{
			if ( to_average[0].data.size() != to_average[i].data.size )
			{
				throw "Incompatible lengths of topological processes, the averages cannot be computed \n";
			}
		}
	}
	//scalar products?
	//confidence bounds?
private:
	std::vector< Representation* > data;
};



}//Gudhi_stat
}//Gudhi

#endif

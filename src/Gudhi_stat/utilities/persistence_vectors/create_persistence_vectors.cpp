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

#include <gudhi/persistence_vectors.h>



using namespace Gudhi;
using namespace Gudhi::Persistence_representations;

#include <iostream>
#include <sstream>


int main( int argc , char** argv )
{
	std::cout << "This program creates persistence vectors of diagrams provided as an input. The first parameter of this program is a dimension of persistence ";
	std::cout << " that will be used in creation of the persistence vectors. If our input files contain persistence pairs of various dimension, as a second parameter of the ";
	std::cout << " procedure please provide the dimension of persistence you want to use. If in your file there are only birth-death pairs of the same dimension, set the first parameter to -1." << std::endl;
	std::cout << "The remaining parameters are the names of files with persistence diagrams. \n";
	int dim = atoi( argv[1] );
	unsigned dimension = std::numeric_limits<unsigned>::max();
	if ( dim >= 0 )
	{
		dimension = (unsigned)dim;
	}
	
	std::vector< const char* > filenames;
	for ( int i = 2 ; i < argc ; ++i )
	{
		filenames.push_back( argv[i] );
	}
	
	for ( size_t i = 0 ; i != filenames.size() ; ++i )
	{
		std::cerr << "Creatign persistence vectors based on a file : " << filenames[i] << std::endl;
		//std::vector< std::pair< double , double > > persistence_pairs = read_gudhi_persistence_file_in_one_dimension( filenames[i] , size_t dimension = 0 )
		Vector_distances_in_diagram< Euclidean_distance > l( filenames[i] , dimension );				
		std::stringstream ss;
		ss << filenames[i] << ".vect";
		l.print_to_file( ss.str().c_str() );
	}
	std::cout << "Done \n";
	return 0;
}


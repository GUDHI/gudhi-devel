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


#include <gudhi/Hausdorff_distances.h>
#include <gudhi/multiplicative_bootstrap.h>
#include <gudhi/read_persistence_from_file.h>
#include <gudhi/persistence_vectors.h>

using namespace Gudhi;
using namespace Gudhi::Gudhi_stat;
using namespace Gudhi::Persistence_representations;



int main( int argc , char** argv )
{
	std::cout << "The parameters of this program are : " << std::endl;
	std::cout << "(a) a name of a file with names of files with persistence diagrams," << std:: endl;
	std::cout << "(b) a number of repetitions of bootstrap (integer)," << std::endl;
	std::cout << "(c) a quantile (real number between 0 and 1. If you do not know what to set, set it to 0.95." << std::endl;
	if ( argc != 4 )
	{
		std::cerr << "Wrong number of parameters, the program will now terminate.\n";
		return 1;
	}
	
	const char* file_with_filenames = argv[1];
	size_t number_of_repetitions_of_bootstrap = (size_t)atoi( argv[2] );
	double quantile = atof( argv[3] );
	
	std::vector< std::string > filenames = readFileNames( file_with_filenames );	
	std::vector< Persistence_landscape* > collection_of_landscapes( filenames.size() );
	for ( size_t i = 0 ; i != filenames.size() ; ++i )
	{
		std::vector< std::pair< double , double > > diag = read_persistence_intervals_in_one_dimension_from_file( filenames[i].c_str() );
		collection_of_landscapes[i] = new Persistence_landscape( diag );		
	}

	//now we can run the bootstrap:
	difference_of_objects<Persistence_landscape> diff;
	norm_of_objects<Persistence_landscape> norm;
		
	double result = 
	multiplicative_bootstrap< Persistence_landscape , difference_of_objects<Persistence_landscape> , norm_of_objects<Persistence_landscape> >
	( collection_of_landscapes , number_of_repetitions_of_bootstrap , diff , norm , quantile , 1 );
	
	std::cout << "result of bootstrap : " << result << std::endl;
	return 0;	
}


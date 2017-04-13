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

#include <gudhi/persistence_representations/Persistence_landscape.h>


using namespace Gudhi;
using namespace Gudhi::Gudhi_stat;

#include <iostream>
#include <sstream>


int main( int argc , char** argv )
{
	std::cout << "This program creates persistence landscapes of diagrams provided as an input. \n";
	std::cout << "The first parameter of the program is the dimension of persistence to be used to construct persistence landscapes. If your file contains ";
	std::cout << "the information about dimension of persistence pairs, please provide here the dimension of persistence pairs you want to use. If your input files consist only ";
	std::cout << "of birth-death pairs, please set this first parameter to -1 \n";	
	std::cout << "The remaining parameters of the program are the names of files with persistence diagrams. \n";
	std::vector< const char* > filenames;
	int dim = atoi(argv[1]);
	unsigned dimension = std::numeric_limits<unsigned>::max();
	if ( dim >= 0 )
	{
		dimension = (unsigned)dim;
	}
	for ( int i = 2 ; i < argc ; ++i )
	{
		filenames.push_back( argv[i] );
	}
	
	std::cout << "Creating persistence landscapes...\n";	
	for ( size_t i = 0 ; i != filenames.size() ; ++i )
	{
		Persistence_landscape l( filenames[i] , dimension );
		std::stringstream ss;
		ss << filenames[i] << ".land";
		l.print_to_file( ss.str().c_str() );
	}
	std::cout << "Done \n";
	return 0;
}

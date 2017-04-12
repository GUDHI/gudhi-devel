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


#include <gudhi/persistence_representations/Persistence_heat_maps.h>



using namespace Gudhi;
using namespace Gudhi::Gudhi_stat;

#include <iostream>
#include <sstream>



int main( int argc , char** argv )
{
	std::cout << "This program creates persistence heat map of diagrams provided as an input. The Gaussian kernels are weighted by the square of distance of a center from the diagonal.\n";
	std::cout << "The first parameter of a program is an integer, a size of a grid.\n";
	std::cout << "The second and third parameters are min and max of the grid. If you want those numbers to be computed based on the data, set them both to -1 \n";
	std::cerr << "The fourth parameter is an integer, the standard deviation of a gaussian kernel expressed in a number of pixels \n";
	std::cout << "The fifth parameter of this program is a dimension of persistence that will be used in creation of the persistence heat maps.";
	std::cout << "If our input files contain persistence pairs of various dimension, as a fifth parameter of the procedure please provide the dimension of persistence you want to use.";
	std::cout << "If in your file there are only birth-death pairs of the same dimension, set the first parameter to -1." << std::endl;				
	std::cout << "The remaining parameters are the names of files with persistence diagrams. \n";

	if ( argc < 5 )
	{
		std::cout << "Wrong parameter list, the program will now terminate \n";
		return 1;
	}
	
	size_t size_of_grid = (size_t)atoi( argv[1] );
	double min_ = atof( argv[2] );
	double max_ = atof( argv[3] );
	size_t stdiv = atof( argv[4] );

	unsigned dimension = std::numeric_limits<unsigned>::max();
	int dim = atoi( argv[5] );	 
	if ( (dim != -1) && (dim >= 0) )
	{
		dimension = (unsigned)dim;
	}
	
	std::vector< const char* > filenames;
	for ( int i = 6 ; i < argc ; ++i )
	{
		filenames.push_back( argv[i] );
	}
	
	std::cout << "Creating persistence heat maps...\n";	
	std::vector< std::vector<double> > filter = create_Gaussian_filter(stdiv,1);
	
	for ( size_t i = 0 ; i != filenames.size() ; ++i )
	{		
		std::cout << "Creating a heat map based on a file : " << filenames[i] << std::endl;
		Persistence_heat_maps<squared_distance_from_diagonal_scaling> l( filenames[i] , filter, false , size_of_grid , min_ , max_ , dimension );		
		
		std::stringstream ss;
		ss << filenames[i] << ".mps";
		l.print_to_file( ss.str().c_str() );
	}
	std::cout << "Done \n";
	return 0;
}

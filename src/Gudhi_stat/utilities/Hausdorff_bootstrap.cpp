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
#include <gudhi/bootstrap.h>
#include <gudhi/read_persitence_from_file.h>
#include <gudhi/concretizations/Vector_distances_in_diagram.h>

using namespace std;
using namespace Gudhi;
using namespace Gudhi::Gudhi_stat;



int main( int argc , char** argv )
{
	std::cout << "The parameters of this program are : " << std::endl;
	std::cout << "(a) a name of a file with points," << std:: endl;
	std::cout << "(b) a number of repetitions of bootstrap (integer)," << std::endl;
	std::cout << "(c) a size of subsample (integer, smaller than the number of points," << std::endl;
	std::cout << "(d) a quantile (real number between 0 and 1. If you do not know what to set, set it to 0.95." << std::endl;
	if ( argc != 5 )
	{
		std::cerr << "Wrong number of parameters, the program will now terminate.\n";
		return 1;
	}
	
	const char* filename = argv[1];
	size_t number_of_repetitions_of_bootstrap = (size_t)atoi( argv[2] );
	size_t size_of_subsample = (size_t)atoi( argv[3] );
	double quantile = atof( argv[4] );
	
	std::cout << "Now we will read points from the file : " << filename << " and then perform " << number_of_repetitions_of_bootstrap << " times the bootstrap on it by choosing subsample of a size " << size_of_subsample << std::endl;
	
	std::vector< std::vector< double > > points = read_numbers_from_file_line_by_line( filename );
	
	std::cout << "Read : " << points.size() << " points.\n";
	
	//comute all-to-all distance matrix:
	std::vector< std::vector<double> > all_to_all_distance_matrix_between_points = compute_all_to_all_distance_matrix_between_points< std::vector<double> , Euclidean_distance<double> >( points );
	Hausdorff_distance_between_subspace_and_the_whole_metric_space distance( all_to_all_distance_matrix_between_points );
	identity< std::vector<size_t> > identity_char;
	
	//and now we can run the real bootstrap.
	//template < typename PointCloudCharacteristics , typename CharacteristicFunction , typename DistanceBetweenPointsCharacteristics >
	//In this case, the PointCloudCharacteristics is just a vector of numbers of points (in a order fixed on points vector). 
	//CharacteristicFunction is just identity, transforming std::vector< size_t > to itself.
	//DistanceBetweenPointsCharacteristics is the place were all happens. This class hace the information about the coordinates of the points, and allows to compute a Hausdorff distance between 
	//the collection of all points, and the subsample. 
	double result = bootstrap< 
							   std::vector< size_t > , //PointCloudCharacteristics
							   identity< std::vector<size_t> > , //CharacteristicFunction
							   Hausdorff_distance_between_subspace_and_the_whole_metric_space //DistanceBetweenPointsCharacteristics. This function have the information about point's coordinates. 
							   >
	( points.size() ,  identity_char , distance , number_of_repetitions_of_bootstrap , size_of_subsample , quantile );
	
	std::cout << "result of bootstrap : " << result << std::endl;
	
	
	return 0;	
}

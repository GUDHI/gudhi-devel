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

//stat part:
#include <gudhi/Hausdorff_distances.h>
#include <gudhi/bootstrap.h>
#include <gudhi/concretizations/Persistence_landscape.h>
#include <gudhi/read_persitence_from_file.h>
#include <gudhi/concretizations/Vector_distances_in_diagram.h>
//persistence part:
#include <gudhi/reader_utils.h>
#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>



using namespace Gudhi;
using namespace Gudhi::Gudhi_stat;
using namespace Gudhi::persistent_cohomology;

typedef int Vertex_handle;
typedef double Filtration_value;


//if this variable is -1, then the infinite interals are ignored. If not, they infinite values are replaced with what_to_replace_infinite_intervals_with:
double what_to_replace_infinite_intervals_with = -1;



class compute_persistence_landscape_of_a_point_cloud_in_certain_dimension
{
public:
	compute_persistence_landscape_of_a_point_cloud_in_certain_dimension( std::vector< std::vector< double > >& points_ , int dimension , double threshold_ , int coeficient_field_ = 11 , double min_persistence_ = 0 ):dim( dimension ),points(points_),threshold(threshold_),coeficient_field(coeficient_field_),min_persistence(min_persistence_){}
	//This function takes a vector of indices (numbers_to_sample). It will select the points from this->points having those indices, construct Rips complex and persistence intervals based on this.
	//Then it will filter the intervals to find only those in the dimension this->dim, and construct a persistence landascape based on this. Thie will be the result of the procedure.
	Persistence_landscape operator()( std::vector< size_t > numbers_to_sample )
	{
		bool dbg = false;
		//take the subsampled points:
		std::vector< std::vector< double > > points_in_subsample;
		points_in_subsample.reserve( numbers_to_sample.size() );
		for ( size_t i = 0 ; i != numbers_to_sample.size() ; ++i )
		{
			points_in_subsample.push_back( this->points[ numbers_to_sample[i] ] );
		}
		//construct a Rips complex based on it and compute its persistence:
		Graph_t prox_graph = compute_proximity_graph(points_in_subsample, this->threshold , euclidean_distance< std::vector< double > >);
		// Construct the Rips complex in a Simplex Tree		
		Simplex_tree<Simplex_tree_options_fast_persistence> st;
		// insert the proximity graph in the simplex tree
		st.insert_graph(prox_graph);
		// expand the graph until dimension dim_max
		st.expansion(this->dim + 1);		
		// Sort the simplices in the order of the filtration
		st.initialize_filtration();
		// Compute the persistence diagram of the complex
		persistent_cohomology::Persistent_cohomology<Simplex_tree<Simplex_tree_options_fast_persistence>, Field_Zp > pcoh(st);
		// initializes the coefficient field for homology
		pcoh.init_coefficients( this->coeficient_field );
		pcoh.compute_persistent_cohomology(this->min_persistence);
		auto persistence_pairs = pcoh.get_persistent_pairs();				
		//From the persistence take only this in the dimension this->dim:
		
		if ( dbg )std::cerr << "Here are the persistence pairs :\n";
		std::vector< std::pair< double,double > > persistence_in_fixed_dimension;
		for ( size_t i = 0 ; i != persistence_pairs.size() ; ++i )
		{
			if ( st.dimension( std::get<0>(persistence_pairs[i]) ) == this->dim )
			{
				double birth = st.filtration( std::get<0>(persistence_pairs[i]) );
				double death = st.filtration( std::get<1>(persistence_pairs[i]) );				
								
				if ( std::get<1>(persistence_pairs[i]) != st.null_simplex() )
				{
					//finite interval
					persistence_in_fixed_dimension.push_back( std::pair<double,double>( birth , death ) );
					if (dbg){std::cout << "birth : " << birth << " , death : " << death << std::endl;}
				}
				else
				{
					//infinite interval
					if ( what_to_replace_infinite_intervals_with != -1 )
					{
						persistence_in_fixed_dimension.push_back( std::pair<double,double>( birth , what_to_replace_infinite_intervals_with ) );
						if (dbg){std::cout << "birth : " << birth << " , death : " << what_to_replace_infinite_intervals_with << std::endl;}
					}
				}				
			}
		}		
		if ( dbg )std::cerr << "Persistence pairs computed \n";
		//Construct and return the persistence landscape:		
		return Persistence_landscape( persistence_in_fixed_dimension );
	}
private:
	int dim;
	std::vector< std::vector< double > >& points;
	double threshold;
	int coeficient_field;
	double min_persistence;
};

class distance_between_landscapes
{
public:
	distance_between_landscapes( double exponent_ ):exponent(exponent_){}
	double operator()( const Persistence_landscape& first , const Persistence_landscape& second )
	{
		return first.distance( second, this->exponent );
	}
private:
	double exponent;
};


int main( int argc , char** argv )
{
	std::cout << "The parameters of this program are : " << std::endl;
	std::cout << "(1) a name of a file with points," << std:: endl;
	std::cout << "(2) a number of repetitions of bootstrap (integer)," << std::endl;
	std::cout << "(3) a size of subsample (integer, smaller than the number of points. " << std::endl;
	std::cout << "(4) An real value p such that L^p distance is going to be computed. \n";
	std::cout << "(5) A dimension of persistence that is to be taken into account (positive integer) \n";
	std::cout << "(6) A maximal diameter to which complex is to be grown (positive integer) \n";
	std::cout << "(d) a quantile (real number between 0 and 1. If you do not know what to set, set it to 0.95." << std::endl;
	if ( argc != 8 )
	{
		std::cerr << "Wrong number of parameters, the program will now terminate.\n";
		return 1;
	}
	
	const char* filename = argv[1];
	size_t number_of_repetitions_of_bootstrap = (size_t)atoi( argv[2] );
	size_t size_of_subsample = (size_t)atoi( argv[3] );
	double p = atoi( argv[4] );
	int dimension = atoi( argv[5] );
	double threshold = atof( argv[6] );
	double quantile = atof( argv[7] );
	
	std::cout << "Now we will read points from the file : " << filename << " and then perform " << number_of_repetitions_of_bootstrap << " times the bootstrap on it by choosing subsample of a size " << size_of_subsample << std::endl;
	
	std::vector< std::vector< double > > points = read_numbers_from_file_line_by_line( filename );
	
	std::cout << "Read : " << points.size() << " points.\n";
	
	distance_between_landscapes distance( p );//L^p distance.
	compute_persistence_landscape_of_a_point_cloud_in_certain_dimension characteristic_fun( points , dimension , threshold );
		
		
	//and now we can run the real bootstrap.
	//template < typename PointCloudCharacteristics , typename CharacteristicFunction , typename DistanceBetweenPointsCharacteristics >
	//In this case, the PointCloudCharacteristics is just a vector of numbers of points (in a order fixed on points vector). 
	//CharacteristicFunction is just identity, transforming std::vector< size_t > to itself.
	//DistanceBetweenPointsCharacteristics is the place were all happens. This class hace the information about the coordinates of the points, and allows to compute a Hausdorff distance between 
	//the collection of all points, and the subsample. 
	double result = bootstrap< 
							   Persistence_landscape , //PointCloudCharacteristics, persistence landascapes constructed based on vector of 
																			 //pairs of birth--death values in a cartain dimension.
							   compute_persistence_landscape_of_a_point_cloud_in_certain_dimension , //CharacteristicFunction, in this case, we will need to compute persistence in a certain dimension.
							   distance_between_landscapes //DistanceBetweenPointsCharacteristics. In this case
							   >
	( points.size() , characteristic_fun , distance , number_of_repetitions_of_bootstrap , size_of_subsample , quantile );
	
	std::cout << "result of bootstrap : " << result << std::endl;
	
	
	return 0;	
}

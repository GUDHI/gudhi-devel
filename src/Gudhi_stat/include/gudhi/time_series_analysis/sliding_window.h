/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2017  Swansea University (UK)
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

#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>

//for persistent homology computations.
#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>

using namespace std;

namespace Gudhi 
{
namespace Gudhi_stat 
{
	
using namespace Gudhi;
using namespace Gudhi::persistent_cohomology;

/**
 * This is an implementation of a sliding window embedding class. Sliding window embedding can be constructed either based on a time series (given as a vector of values), or a real-valued
 * function which is sampled on a nunmber of points (which also give a vector of real values). 
 * Given a vector of real values V=[v0,v1,v2,...,v(n-1),vn] and a positive integer M (assumed in this example to be 2) the sliding window embedding of V into R^M is the following point cloud:
 * (v0,v1)
 * (v1,v2)
 * (v2,v3)
 * ...
 * (v(n-1),vn)
 * The construction of sliding window embedding allows to build a point-cloud representation of a time series and later analyze this representation by using topological methods. 
 * The class presented below build a representation, and compute persistence of it. 
**/ 

//construction from time series
//construction from real value function
//higher dimensional analaogous?? Shall we design periodicity testing in general, multi-d setting??

typedef Simplex_tree<Simplex_tree_options_fast_persistence> ST;//for persistence computations

class Sliding_window_embedding
{
public:
	
	/**
	 * This is a basic constructor of a Sliding_window_embedding class. It take as an input a vector of time series values and the dimension of embedding and 
	 * given those construct an object of a Sliding_window_embedding class. 
	**/ 
	Sliding_window_embedding( const std::vector< double >& ts , unsigned _dimension_of_embedding ):dimension_of_embedding(_dimension_of_embedding)
	{
		this->time_series = std::vector<double>(ts);		
	}
	
	/**
	 * This is a constructor of a Sliding_window_embedding class. It take as an input a c-style pointer to a function, range of parameters for which the sliding window will be constructed, 
	 * number of steps in the construction (i.e. number of equally distributed points in the desired domain of the function on which the function will be sampled) and the dimension of the 
	 * embedding. Given those information it construct an object of a Sliding_window_embedding class. 
	**/ 
	Sliding_window_embedding( double (*function)(double) , double x_min , double x_max , unsigned number_of_steps , unsigned _dimension_of_embedding ):dimension_of_embedding(_dimension_of_embedding)
	{
		assert( x_min < x_max );
		this->time_series.reserve( number_of_steps+1 );
		double x = x_min;
		double dx = ( x_max-x_min )/(double)number_of_steps;
		for ( size_t i = 0 ; i <= number_of_steps ; ++i )
		{
			this->time_series.push_back( function(x) );
			x += dx;	
		}		
	}
	
	/**
	 * This is a constructor of a Sliding_window_embedding class. It take as an input a filename of file with the values and the dimension of embedding and 
	 * given those construct an object of a Sliding_window_embedding class. 
	 * We assume that the input file is a text file containing numbers separated with a white space. 
	**/ 
	Sliding_window_embedding( const char* filename , unsigned _dimension_of_embedding ):dimension_of_embedding(_dimension_of_embedding)
	{
		ifstream inputFile( filename );
		//read the file and populate this->time_series
		if (inputFile) 
		{        
			double value;
			while ( inputFile >> value ) 
			{
				this->time_series.push_back(value);
			}
		}
		inputFile.close();
	}
	
	/**
	 * This procedure create a point cloud from a sliding window embedding.
	**/ 
	std::vector< std::vector<double> > create_point_cloud()
	{		
		bool dbg = false;
		std::vector< std::vector<double> > result;
		if ( this->time_series.size() < this->dimension_of_embedding )return result;
		result.reserve( this->time_series.size() - this->dimension_of_embedding );
		for ( size_t i =  0 ; i <= this->time_series.size() - this->dimension_of_embedding ; ++i )
		{						
			std::vector<double> point( this->dimension_of_embedding );
			for ( size_t j = 0 ; j != this->dimension_of_embedding ; ++j )
			{
				point[j] = this->time_series[i+j];
			}
			result.push_back( point );
			if ( dbg )
			{
				std::cerr << "Adding point : ";
				for ( size_t k = 0 ; k != point.size() ; ++k )
				{
					std::cerr << point[k] << " ";
				}
				std::cerr << std::endl;
			}
		}
		return result;
	}
	
	/**
	 * This procedure store sliding window embedding point cloud in a file.
	**/ 
	void create_point_cloud( const char* filename )
	{
		ofstream out( filename );
		std::vector< std::vector<double> > point_cloud = this->create_point_cloud();
		for ( size_t i = 0 ; i != point_cloud.size() ; ++i )
		{
			for ( size_t j = 0 ; j != point_cloud[i].size() ; ++j )
			{
				out << point_cloud[i][j] << " ";
			}
			out << endl;
		}
		out.close();
	}
	
	/**
	 * This procedure compute persistence of the sliding window embedding point cloud by using Vietoris-Rips filtration.
	**/
	persistent_cohomology::Persistent_cohomology<ST, Field_Zp > compute_persistence_of_Vietoris_Rips_complex( double threshold , unsigned dim_max , unsigned field_coef = 2 , double min_persistence = 0 ) 
	{	
		//compute points of the sliding window embedding.
		std::vector< std::vector< double > > points = this->create_point_cloud();
		
		// Compute the proximity graph of the points.
		Graph_t prox_graph = compute_proximity_graph(points, threshold , euclidean_distance< std::vector< double > >);

		// Construct the Rips complex in a Simplex Tree.		
		ST st;
		// insert the proximity graph in the simplex tree.
		st.insert_graph(prox_graph);
		// expand the graph until dimension dim_max.
		st.expansion(dim_max);

		std::cout << "The complex contains " << st.num_simplices() << " simplices \n";
		std::cout << "   and has dimension " << st.dimension() << " \n";

		// Sort the simplices in the order of the filtration.
		st.initialize_filtration();

		// Compute the persistence diagram of the complex.
		persistent_cohomology::Persistent_cohomology<ST, Field_Zp > pcoh(st);
		// initializes the coefficient field for homology.
		pcoh.init_coefficients( field_coef );

		//compute persistent cohomology.
		pcoh.compute_persistent_cohomology(min_persistence);
		return pcoh;
		pcoh.output_diagram();
	}//compute_persistence_of_Vietoris_Rips_complex
	
	/**
	 * This procedure compute persistence of the sliding window embedding point cloud by using Alpha complex filtration.
	**/ 
	
	/*persistent_cohomology::Persistent_cohomology<ST, Field_Zp >*/
	/*void compute_persistence_of_Alpha_complex( double threshold , unsigned dim_max , unsigned field_coef = 2 , double min_persistence = 0 ) 
	{	
		std::string off_file_points;
		std::string output_file_diag;
		Filtration_value alpha_square_max_value;
		int coeff_field_characteristic;
		Filtration_value min_persistence;

		program_options(argc, argv, off_file_points, output_file_diag, alpha_square_max_value,
			  coeff_field_characteristic, min_persistence);

		// ----------------------------------------------------------------------------
		// Init of an alpha complex from an OFF file
		// ----------------------------------------------------------------------------
		using Kernel = CGAL::Epick_d< CGAL::Dynamic_dimension_tag >;
		Gudhi::alpha_complex::Alpha_complex<Kernel> alpha_complex_from_file(off_file_points, alpha_square_max_value);

		// ----------------------------------------------------------------------------
		// Display information about the alpha complex
		// ----------------------------------------------------------------------------
		std::cout << "Alpha complex is of dimension " << alpha_complex_from_file.dimension() <<
		" - " << alpha_complex_from_file.num_simplices() << " simplices - " <<
		alpha_complex_from_file.num_vertices() << " vertices." << std::endl;

		// Sort the simplices in the order of the filtration
		alpha_complex_from_file.initialize_filtration();

		std::cout << "Simplex_tree dim: " << alpha_complex_from_file.dimension() << std::endl;
		// Compute the persistence diagram of the complex
		Gudhi::persistent_cohomology::Persistent_cohomology< Gudhi::alpha_complex::Alpha_complex<Kernel>,
		Gudhi::persistent_cohomology::Field_Zp > pcoh(alpha_complex_from_file);
		// initializes the coefficient field for homology
		pcoh.init_coefficients(coeff_field_characteristic);

		pcoh.compute_persistent_cohomology(min_persistence);

		// Output the diagram in filediag
		if (output_file_diag.empty()) {
		pcoh.output_diagram();
		} else {
		std::cout << "Result in file: " << output_file_diag << std::endl;
		std::ofstream out(output_file_diag);
		pcoh.output_diagram(out);
		out.close();
		}
	}//compute_persistence_of_Alpha_complex
		*/

	
private:
	std::vector< double > time_series;
	unsigned dimension_of_embedding;
};

}//namespace Gudhi_stat 
}//namespace Gudhi 

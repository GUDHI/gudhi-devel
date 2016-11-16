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
#include <gudhi/concretizations/Persistence_heat_maps.h>
#include <vector>

//extras
#include <gudhi/common_gudhi_stat.h>

namespace Gudhi 
{
namespace Gudhi_stat 
{
	
//over here we will need a few version of construct_representation_from_file procedure, since different representations may require different parameters. This is a procedure that in my 
//oppinion cannot be standarize, since construction of representation cannot. But, the remaining part of the code in my opinion is free from any details of representation.


//the reason I am separating the process of getting the intervals from the process of consstructing the represnetation is the following: for soem representations, we may need to ahve the same
//scale. To determine the scale, we need to know all the intervals before. So, we read the intervals first, then if needed, process them, to get the parametersof the representation,
//and only then construct the representations. 
std::vector< std::vector< std::pair< double , double > > > read_persistence_pairs_from_file_with_names_of_files( const char* filename )
{
	bool dbg = false;
	std::vector< std::vector< std::pair< double , double > > > result;
	
	std::vector< std::string > files = readFileNames( filename );	
	
	std::cout << "Here are the filenames in the file : " << filename << std::endl;
	for ( size_t i = 0 ; i != files.size() ; ++i )
	{
		std::cout << files[i] << std::endl;
	}
	
	for ( size_t i = 0 ; i != files.size() ; ++i )
	{
		std::vector< std::pair< double , double > > diag = read_standard_file( files[i].c_str() );	
		result.push_back( diag );		
		if ( dbg )
		{
			std::cerr << "Here is a diagram from a file : " << files[i].c_str() << std::endl;
			for ( size_t aa = 0  ; aa != diag.size() ; ++aa )
			{
				std::cout << diag[aa].first << " " << diag[aa].second << std::endl;
			}
			getchar();
		}
	}
	return result;
}

//When workign with time varying data, we tpically have a collection of files with intervals to produce one process. The intervals from all the files that constitute a single process can be read with 
//the read_persistence_pairs_from_file_with_names_of_files procedure. 
//But then, we also may need to read the persistence intervals of collection of proesses we may want to read the intervals first. This is what the procedre 
std::vector< std::vector< std::vector< std::pair< double , double > > > > read_persistence_pairs_from_files_with_names_of_files( const std::vector<const char*>& filenames )
{
	std::vector< std::vector< std::vector< std::pair< double , double > > > > result;
	result.reserve( filenames.size() );
	for ( size_t file_no = 0 ; file_no != filenames.size() ; ++file_no )
	{
		result.push_back( read_persistence_pairs_from_file_with_names_of_files( filenames[file_no] ) );
	}
	return result;
}//read_persistence_pairs_from_files_with_names_of_files


std::pair< std::pair< double,double > , std::pair< double,double > > find_x_and_y_ranges_of_intervals( const std::vector< std::vector< std::pair< double , double > > >& intervals_from_file )
{
	double min_x = std::numeric_limits< double >::max();
	double max_x = -std::numeric_limits< double >::max();
	double min_y = std::numeric_limits< double >::max();
	double max_y = -std::numeric_limits< double >::max();
	for ( size_t i = 0 ; i != intervals_from_file.size() ; ++i )
	{	
		for ( size_t j = 0 ; j != intervals_from_file[i].size() ; ++j )
		{
			if ( min_x > intervals_from_file[i][j].first )min_x = intervals_from_file[i][j].first;
			if ( max_x < intervals_from_file[i][j].first )max_x = intervals_from_file[i][j].first;
			
			if ( min_y > intervals_from_file[i][j].second )min_y = intervals_from_file[i][j].second;
			if ( max_y < intervals_from_file[i][j].second )max_y = intervals_from_file[i][j].second;
		}		
	}
	return std::make_pair( std::make_pair( min_x,max_x ) , std::make_pair( min_y,max_y ) );
}//find_x_and_y_ranges_of_intervals


std::pair< std::pair< double,double > , std::pair< double,double > > find_x_and_y_ranges_of_intervals( const std::vector< std::vector< std::vector< std::pair< double , double > > > >& intervals_from_files )
{
	double min_x = std::numeric_limits< double >::max();
	double max_x = -std::numeric_limits< double >::max();
	double min_y = std::numeric_limits< double >::max();
	double max_y = -std::numeric_limits< double >::max();
	for ( size_t i = 0 ; i != intervals_from_files.size() ; ++i )
	{
		std::pair< std::pair< double,double > , std::pair< double,double > > ranges = find_x_and_y_ranges_of_intervals( intervals_from_files[i] );
		if ( min_x > ranges.first.first ) min_x = ranges.first.first;
		if ( max_x < ranges.first.second ) max_x = ranges.first.second;
		if ( min_y > ranges.second.first ) min_y = ranges.second.first;
		if ( max_y < ranges.second.second ) max_y = ranges.second.second;
	}
	return std::make_pair( std::make_pair( min_x,max_x ) , std::make_pair( min_y,max_y ) );
}


template <typename Representation>
std::vector< Representation* > construct_representation_from_file( const char* filename )	
{
	std::vector< std::vector< std::pair< double , double > > > intervals_from_file = read_persistence_pairs_from_file_with_names_of_files( filename );
	std::vector< Representation* > result( intervals_from_file.size() );
	for ( size_t i = 0 ; i != intervals_from_file.size() ; ++i )
	{										
		Representation* l = new Representation( intervals_from_file[i] );			
		result[i] = l;		
	}
	return result;
}

//this one can be use for Persistence_intervals.h and Persistence_landscape.h
template <typename Representation>
std::vector< Representation* > construct_representation_from_file( const std::vector< std::vector< std::pair< double , double > > >& intervals_from_file)	
{
	
	std::vector< Representation* > result( intervals_from_file.size() );
	for ( size_t i = 0 ; i != intervals_from_file.size() ; ++i )
	{										
		Representation* l = new Representation( intervals_from_file[i] );			
		result[i] = l;		
	}
	return result;
}

//this one can be use for Persistence_heat_maps.h
template <typename Representation>
std::vector< Representation* > construct_representation_from_file( const std::vector< std::vector< std::pair< double , double > > >& intervals_from_file , 
																   std::vector< std::vector<double> > filter = create_Gaussian_filter(5,1), 
																   bool erase_below_diagonal = false , 
																   size_t number_of_pixels = 1000 , 
																   double min_ = -1 , double max_ = -1  )	
{	
	std::vector< Representation* > result( intervals_from_file.size() );
	for ( size_t i = 0 ; i != intervals_from_file.size() ; ++i )
	{										
		Representation* l = new Representation( intervals_from_file[i] , filter , erase_below_diagonal , number_of_pixels , min_ , max_ );			
		result[i] = l;		
	}
	return result;
}

//this one can be use for Persistence_landscape_on_grid.h
template <typename Representation>
std::vector< Representation* > construct_representation_from_file( const std::vector< std::vector< std::pair< double , double > > >& intervals_from_file, 
																   double grid_min_ , double grid_max_ , size_t number_of_points_ 
)	
{	
	std::vector< Representation* > result( intervals_from_file.size() );
	for ( size_t i = 0 ; i != intervals_from_file.size() ; ++i )
	{										
		Representation* l = new Representation( intervals_from_file[i] , grid_min_ , grid_max_ , number_of_points_ );			
		result[i] = l;		
	}
	return result;
}


//this one can be use for Vector_distances_in_diagram.h
template <typename Representation>
std::vector< Representation* > construct_representation_from_file( const std::vector< std::vector< std::pair< double , double > > >& intervals_from_file, 
																   size_t where_to_cut
)	
{
	std::vector< Representation* > result( intervals_from_file.size() );
	for ( size_t i = 0 ; i != intervals_from_file.size() ; ++i )
	{										
		Representation* l = new Representation( intervals_from_file[i] , where_to_cut );			
		result[i] = l;		
	}
	return result;
}
	

template <typename Representation>
class Topological_process
{
public:
	Topological_process(){};
	~Topological_process()
	{		
		for ( size_t i = 0  ; i != this->data.size() ; ++i )
		{
			delete this->data[i];
		}
	}
	Topological_process( const Topological_process& org )
	{
		this->data = std::vector< Representation* >( org.data.size() );
		for ( size_t i = 0 ; i != org.data.size() ; ++i )
		{
			this->data[i] = new Representation( *org.data[i] );
		}		
	}
	Topological_process& operator = ( const Topological_process& rhs )
	{
		for ( size_t i = 0 ; i != rhs.data.size() ; ++i )
		{
			this->data[i] = new Representation( *rhs.data[i] );
		}	
		return *this;	
	}
	
	Topological_process( const std::vector< Representation* >& data_ ):data(data_){}
	double distance( const Topological_process& second , double exponent = 1 )
	{
		if ( this->data.size() != second.data.size() )
		{
			throw "Incompatible lengths of topological processes, we cannot compute the distance in this case \n";
		}
		double result = 0;
		for ( size_t i = 0 ; i != this->data.size() ; ++i )
		{
			result += this->data[i]->distance( *second.data[i] , exponent );
		}
		return result;
	}	
	
	void compute_average( const std::vector< Topological_process* >& to_average )
	{
		//since we will substitute whatever data we have in this object with an average, we clear the data in this object first:
		this->data.clear();
		//then we need to check if all the topological processes in the vector to_average have the same length.
		if ( to_average.size() == 0 )return;
		for ( size_t i = 1 ; i != to_average.size() ; ++i )
		{
			if ( to_average[0]->data.size() != to_average[i]->data.size() )
			{
				throw "Incompatible lengths of topological processes, the averages cannot be computed \n";
			}
		}
		
		this->data.reserve( to_average[0]->data.size() );
		
		for ( size_t level = 0 ; level != to_average[0]->data.size() ; ++level )
		{
			//now we will be averaging the level level:
			std::vector< Representation* > to_average_at_this_level;
			to_average_at_this_level.reserve( to_average.size() );
			for ( size_t i = 0 ; i != to_average.size() ; ++i )
			{
				to_average_at_this_level.push_back( to_average[i]->data[level] );
			}
			Representation* average_representation_on_this_level = new Representation;
			average_representation_on_this_level->compute_average( to_average_at_this_level );
			this->data.push_back( average_representation_on_this_level );
		}				
	}
	
	std::pair< double , double > gimme_x_range()const
	{
		double min_x = std::numeric_limits< double >::max();		
		double max_x = -std::numeric_limits< double >::max();
		for ( size_t i = 0 ; i != this->data.size() ; ++i )
		{
			std::pair< double , double > xrange = this->data[i]->gimme_x_range();			
			if ( min_x > xrange.first )min_x = xrange.first;
			if ( max_x < xrange.second )max_x = xrange.second;			
		}
		return std::make_pair( min_x , max_x );
	}
	
	std::pair< double , double > gimme_y_range()const
	{
		double min_y = std::numeric_limits< double >::max();		
		double max_y = -std::numeric_limits< double >::max();
		for ( size_t i = 0 ; i != this->data.size() ; ++i )
		{
			std::pair< double , double > yrange = this->data[i]->gimme_y_range();			
			if ( min_y > yrange.first )min_y = yrange.first;
			if ( max_y < yrange.second )max_y = yrange.second;			
		}
		return std::make_pair( min_y , max_y );
	}
	
	/**
	 * The procedure checks if the ranges of data are the same for all of them.
	**/ 
	bool are_the_data_aligned()const
	{
		if ( this->data.size() == 0 )return true;//empty collection is aligned 
		std::pair< double , double > x_range = this->data[0]->gimme_x_range();
		std::pair< double , double > y_range = this->data[0]->gimme_y_range();
		for ( size_t i = 1 ; i != this->data.size() ; ++i )
		{
			if ( (x_range != this->data[i]->gimme_x_range()) || (y_range != this->data[i]->gimme_y_range()) )
			{
				return false;
			}
		}
		return true;
	}
	
	
	//scalar products?
	//confidence bounds?
	
	void plot( const char* filename , size_t delay = 30 , double min_x = -1 , double max_x = -1 , double min_y = -1 , double max_y = -1 )
	{				
		std::vector< std::string > filenames;		
		//over here we need to			
		for ( size_t i = 0 ; i != this->data.size() ; ++i )
		{
			std::stringstream ss;
			ss << filename << "_" << i;
			if ( ( min_x != max_x ) && ( min_y != max_y ) )
			{
				//in this case, we set up the uniform min and max values for pciture
				//this->data[i]->plot( ss.str().c_str() , min_x , max_x , min_y , max_y );
			}
			else
			{
				//in this case, we DO NOT set up the uniform min and max values for pciture
				this->data[i]->plot( ss.str().c_str() );
			}
			ss << "_GnuplotScript";
			filenames.push_back( ss.str() );									
		}
		//and now we have to call it somehow for all the files. We will create a script that will call all the other scripts and ceate a sequence of jpg/png files. 

	
		std::stringstream gif_file_name;
		gif_file_name << filename << ".gif";
		std::stringstream gif_gnuplot_script_file_name;
		gif_gnuplot_script_file_name << filename << "_gif_gnuplot_script";
		
		ofstream out;
		out.open( gif_gnuplot_script_file_name.str().c_str() );
		out << "set terminal gif animate delay " << delay << std::endl;
		out << "set output '" << gif_file_name.str() << "'" << std::endl;
		if ( min_x != max_x )
		{
			out << "set xrange [" << min_x << ":" << max_x << "]" << std::endl;
			out << "set yrange [" << min_y << ":" << max_y << "]" << std::endl;		
		}
		
		for ( size_t i = 0 ; i != filenames.size() ; ++i )
		{
			out << " load '" << filenames[i] << "'" << std::endl;
		}
		out.close();		
		
		std::cout << std::endl << std::endl << std::endl << "Open gnuplot terminal and type load '" << gif_gnuplot_script_file_name.str() << "' to create an animated gif \n";
	}//plot
	
	
	std::vector< Representation* > gimme_data(){return this->data;}
private:
	std::vector< Representation* > data;
};



}//Gudhi_stat
}//Gudhi

#endif

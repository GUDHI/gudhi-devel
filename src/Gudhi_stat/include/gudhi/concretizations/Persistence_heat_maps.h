/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2015  INRIA Sophia-Saclay (France)
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
 


//standard include
#include <vector>
#include <sstream>
#include <iostream>

//gudhi include
#include "read_persitence_from_file.h"
#include "Distances_of_points_in_diagram.h"


using namespace std;

std::vector< std::vector<double> > create_Gaussian_filter( size_t pixel_radius , double sigma )
{
	//we are computing the kernel mask to 2 standard deviations away from the center. We discretize it in a grid of a size 2*pixel_radius times 2*pixel_radius.
	
    double r = 0;
    double sigma_sqr = sigma * sigma;
 
    // sum is for normalization
    double sum = 0;
    
    //initialization of a kernel:
    std::vector< std::vector<double> > kernel( 2*pixel_radius  );
    for ( size_t i = 0 ; i != 2*pixel_radius ; ++i )
    {
		std::vector<double> v( 2*pixel_radius , 0 );
		kernel[i] = v;
	}
 
    for (int x = -pixel_radius; x <= pixel_radius; x++)
    {
        for(int y = -pixel_radius; y <= pixel_radius; y++)
        {
			double real_x = 2*sigma*x/pixel_radius;
			double real_y = 2*sigma*y/pixel_radius;						
            r = sqrt(real_x*real_x + real_y*real_y);
            kernel[x + pixel_radius][y + pixel_radius] = (exp(-(r*r)/sigma_sqr))/(3.141592 * sigma_sqr);
            sum += gKernel[x + pixel_radius][y + pixel_radius];
        }
    }
 
    // normalize the kernel
    for(int i = 0; i != kernel.size() ; ++i)
    {
        for(int j = 0; j != kernel[i].size() ; ++j)
        {
            gKernel[i][j] /= sum;
		}
            
    }
    return kernel;
}

double constant_function( std::pair< double , double >& point_in_diagram )
{
	return 1;
}



class Persistence_heat_maps
{
public:
	Persistence_heat_maps()
	{
		this->scalling_function_with_respect_to_distance_from_diagonal = constant_function;
		this->erase_below_diagonal = false;
		this->filter = create_Gaussian_filter(5,1);
		this->min_ = this->max_ = 0;
	};
    Persistence_heat_maps( const std::vector< std::pair< double,double > >  & interval , std::vector< std::vector<double> > filter = create_Gaussian_filter(5,1) ,  double (*scalling_function_with_respect_to_distance_from_diagonal)( std::pair< double , double >& point_in_diagram ) = constant_function, bool erase_below_diagonal = false , int number_of_pixels = 1000 , double min_ = -1 , double max_ = -1  );
    Persistence_heat_maps( const char* name_of_file_with_names_of_files_with_interval , std::vector< std::vector<double> > filter = create_Gaussian_filter(5,1) , double (*scalling_function_with_respect_to_distance_from_diagonal)( std::pair< double , double >& point_in_diagram ) = constant_function, bool erase_below_diagonal = false , int number_of_pixels = 1000 , double min_ = -1 , double max_ = -1  );

	
	void load_mean( const std::vector<Persistence_heat_maps*>& maps );
	void load_median( const std::vector<Persistence_heat_maps*>& maps );	
	void load_percentage_of_active( const std::vector<Persistence_heat_maps*>& maps );
    
    //put to file subroutine
    void write_to_file( const char* filename );
    
    //read from file subroutine
    Persistence_heat_maps( const char* filename );
    
    void plot( const char* filename );
private:
	//private methods
	std::vector< std::vector<double> > check_and_initialize_maps( const std::vector<Persistence_heat_maps*>& maps );
    void construct( const std::vector< std::pair<double,double> >& intervals_ , double min_ = -1 , double max_ = -1);
    
    //data
    std::vector< std::vector<double> > filter;      
    double (*scalling_function_with_respect_to_distance_from_diagonal)( std::pair< double , double >& point_in_diagram );    
    bool erase_below_diagonal;
    double min_;
    double max_;
    std::vector< std::vector< double > > heat_map;    
};


//if min_ == max_, then the program is requested to set up the values itself based on persistence intervals
void Persistence_heat_maps::construct( const std::vector< std::pair<double,doubl> >& intervals_ , double min_ , double max_ )
{
    bool dbg = false;

    if ( min_ == max_ )
    {
        //in this case, we want the program to set up the min_ and max_ values by itself.
        min_ = INT_MAX;
        max_ = -INT_MAX;
        
        
        for ( size_t i = 0 ; i != intervals_.size() ; ++i )
        {
			if ( intervals_[i].first < min_ )min_ = min_max.first;
			if ( intervals_[i].second > max_ )max_ = min_max.second;
		}
        //now we have the structure filled in, and moreover we know min_ and max_ values of the interval, so we know the range.

        //add some more space:
        min_ -= fabs(max_ - min_)/100;
        max_ += fabs(max_ - min_)/100;
    }

	if ( dng )
	{
		cerr << "min_ : " << min_ << endl;
		cerr << "max_ : " << max_ << endl;
	}

    this->min_ = min_;
    this->max_ = max_;


    //initialization of the structure heat_map
    std::vector< std::vector<double> > heat_map_;
    for ( size_t i = 0 ; i != this->number_of_pixels ; ++i )
    {       
		std::vector<double> v( this->number_of_pixels , 0 );
        heat_map_.push_back( v );
    }
    this->heat_map = heat_map_;

    if (dbg)cerr << "Done creating of the heat map, now we will fill in the structure \n";

	//we will use this extra data structure to store the information about points in the diagram and their weights with respect to the weighting function.
	std::vector< std::vector< std::vector< double > > > weights_of_points;
	for ( aaa )
	

	for ( size_t i = 0 ; i != this->number_of_pixels ; ++i )
	{
		for ( size_t j = 0 ; j != this->number_of_pixels ; ++j )
		{
			//compute the value of a heat map at a point [i,j].
			aaa, todo
		}
	}
}//construct

todo
Persistence_heat_maps::Persistence_heat_maps( const std::vector< std::pair< double,double > >  & interval , std::vector< std::vector<double> > filter = create_Gaussian_filter(5,1) ,  double (*scalling_function_with_respect_to_distance_from_diagonal)( std::pair< double , double >& point_in_diagram ) = constant_function, bool erase_below_diagonal = false , int number_of_pixels = 1000 , double min_ = -1 , double max_ = -1  )
{
    this->construct( intervals_ , min_ , max_ );
}

todo
Persistence_heat_maps::Persistence_heat_maps( const char* name_of_file_with_names_of_files_with_interval , std::vector< std::vector<double> > filter = create_Gaussian_filter(5,1) , double (*scalling_function_with_respect_to_distance_from_diagonal)( std::pair< double , double >& point_in_diagram ) = constant_function, bool erase_below_diagonal = false , int number_of_pixels = 1000 , double min_ = -1 , double max_ = -1  )
{    
    std::vector< std::pair< double , double > > interval = read_standard_file( name_of_file_with_names_of_files_with_interval );   
    this->construct( interval , min_ , max_);
}


double vector_median(std::vector<double> vec)
{
    if(vec.empty()) return 0;
    else
    {
        std::sort(vec.begin(), vec.end());
        if(vec.size() % 2 == 0)
        {
            return (vec[vec.size()/2 - 1] + vec[vec.size()/2]) / 2;
        }
        else
        {
            return vec[vec.size()/2];
        }
    }
}


std::vector< std::vector<double> > Persistence_heat_maps::check_and_initialize_maps( const std::vector<Persistence_heat_maps*>& maps )
{
	//checking if all the heat maps are of the same size:
	for ( size_t i = 0 ; i != maps.heat_maps.size() ; ++i )
    {
		if ( maps.heat_maps[i].size() != maps.heat_maps[0].size() )
		{
			std::cerr << "Sizes of Persistence_heat_maps are not compatible. The program will terminate now \n";
			throw "Sizes of Persistence_heat_maps are not compatible. The program will terminate now \n";
		}
		if ( maps.heat_maps[i][0].size() != maps.heat_maps[0][0].size() )
		{
			std::cerr << "Sizes of Persistence_heat_maps are not compatible. The program will terminate now \n";
			throw "Sizes of Persistence_heat_maps are not compatible. The program will terminate now \n";			
		}
	} 	
	std::vector< std::vector<double> > heat_maps( maps.heat_maps.size() );
	for ( size_t i = 0 ; i != maps.heat_maps.size() ; ++i )
    {
		 std::vector<double> v( maps.heat_maps[0].size() , 0 );
		 heat_maps[i] = v;
	}
	return heat_maps;
}


void Persistence_heat_maps::load_median( const std::vector<Persistence_heat_maps*>& maps )
{
	std::vector< std::vector<double> > heat_maps = this->check_and_initialize_maps( maps );
	
	std::vector<double> to_compute_median( maps.size() );
    for ( size_t i = 0 ; i != heat_map.size() ; ++i )
    {
        for ( size_t j = 0 ; j != heat_map[i].size() ; ++j )
        {
			for ( size_t aa = 0 ; aa != maps.size() ; ++aa )
			{
				to_compute_median[aa] = maps[i][j];
			}
            heat_maps[i][j] =  std::nth_element(to_compute_median.begin(), to_compute_median.begin() + to_compute_median.size()/2, to_compute_median.end());
        }       
    }
	this->heat_map = heat_maps;
}



void Persistence_heat_maps::load_mean( const std::vector<Persistence_heat_maps*>& maps )
{
	std::vector< std::vector<double> > heat_maps = this->check_and_initialize_maps( maps );
	
	std::vector<double> to_compute_median( maps.size() );
    for ( size_t i = 0 ; i != heat_map.size() ; ++i )
    {
        for ( size_t j = 0 ; j != heat_map[i].size() ; ++j )
        {
			double mean = 0;
			for ( size_t aa = 0 ; aa != maps.size() ; ++aa )
			{
				mean += maps[i][j];
			}
            heat_maps[i][j] =  mean/(double)maps.size();
        }       
    }
	this->heat_map = heat_maps;
}





void Persistence_heat_maps::load_percentage_of_active( const std::vector<Persistence_heat_maps*> maps , size_t cutoff )
{
	std::vector< std::vector<double> > heat_maps = this->check_and_initialize_maps( maps );

    for ( size_t i = 0 ; i != this->heat_map.size() ; ++i )
    {
        for ( size_t j = 0 ; j != this->heat_map[i].size() ; ++j )
        {
            int number_of_active_levels = 0;
            for ( size_t k = 0 ; k != this->heat_map[j][i].size() ; ++k )
            {
                if ( this->heat_map[j][i][k] ) number_of_active_levels++;
            }
            if ( number_of_active_levels > cutoff )
            {
                heat_maps[i][j] =number_of_active_levels;
            }
            else
            {
                heat_maps[i][j] = 0;
            }
        }       
    }
    this->heat_map = heat_maps;
}



void Persistence_heat_maps::plot( const char* filename )
{
	 ofstream out;
    out.open( filename );
	out << "plot      '-' matrix with image" << std::endl;
    for ( size_t i = 0 ; i != this->heat_map.size() ; ++i )
    {
        for ( size_t j = 0 ; j != this->heat_map[i].size() ; ++j )
        {
            out << this->heat_map[i][j] << " ";
        }
        out << endl;
    }

    out.close();
}



void Persistence_heat_maps::write_to_file( const char* filename )
{
	ofstream out;
	out.open( filename );
	for ( size_t i = 0 ; i != this->heat_map.size() ; ++i )
    {
        for ( size_t j = 0 ; j != this->heat_map[i].size() ; ++j )
        {
            out << this->heat_map[i][j] << " ";
        }
        out << endl;
    }
	out.close();
}


Persistence_heat_maps::Persistence_heat_maps( const char* filename )
{
	bool dbg = true;
	
	ifstream in;
	in.open( filename );
	
	//checking if the file exist / if it was open. 
	if ( !( access( filename, F_OK ) != -1 ) )
	{
		cerr << "The file : " << filename << " do not exist. The program will now terminate \n";
		throw "The file from which you are trying to read the persistence landscape do not exist. The program will now terminate \n";
	}
	
	//now we read the file one by one. 
	std::string line;
	 
	while (!in.eof())
	{
        getline(in,line);
        std::stringstream lineSS;
        lineSS << line;
        
        std::vector< std::vector<double> > line_of_heat_map;    
        while ( lineSS.good() )
        {
			double point;
			lineSS >> point;
			line_of_heat_map.push_back( point );
			if ( dbg )
			{
				std::cout << point << " ";
			}
		}
		if ( dbg )
		{
			std::cout << std::endl;
		}

		this->heat_map.push_back( line_of_heat_map );
	}	
	in.close();
	if ( dbg )std::cout << "Done \n";
}

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

#pragma once
#ifndef PSSK_H
#define PSSK_H

//gudhi include
#include <gudhi/concretizations/Persistence_heat_maps.h>


namespace Gudhi
{
namespace Gudhi_stat
{

/**
* This is a version of a representation presented in https://arxiv.org/abs/1412.6821
* In that paper the authors are using the representation just to compute kernel. Over here, we extend the usability by far.
* Note that the version presented here is not exact, since we are discretizing the kernel.
* The only difference with respect to the original class is the method of creation. We have full (square) image, and for every point (p,q), we add a kernel at (p,q) and the negative kernel
* at (q,p)
**/

class PSSK : public Persistence_heat_maps<constant_scaling_function>
{
public:
	PSSK():Persistence_heat_maps(){}
	
	PSSK(const std::vector< std::pair< double,double > >  & interval , std::vector< std::vector<double> > filter = create_Gaussian_filter(5,1) , size_t number_of_pixels = 1000 , double min_ = -1 , double max_ = -1 )
	:Persistence_heat_maps()
	{
		this->construct( interval  , filter , number_of_pixels , min_ , max_ );
	}
	
	 
	 PSSK( const char* filename , std::vector< std::vector<double> > filter = create_Gaussian_filter(5,1) , size_t number_of_pixels = 1000 , double min_ = -1 , double max_ = -1  ):
	    Persistence_heat_maps()
	    {
			std::vector< std::pair< double , double > > intervals_ = read_standard_file( filename );         
			this->construct( intervals_ , filter , number_of_pixels , min_ , max_ );
		}

protected:
	void construct( const std::vector< std::pair<double,double> >& intervals_  , 
					std::vector< std::vector<double> > filter = create_Gaussian_filter(5,1),
                    size_t number_of_pixels = 1000 , double min_ = -1 , double max_ = -1 );
};

//if min_ == max_, then the program is requested to set up the values itself based on persistence intervals
void PSSK::construct( const std::vector< std::pair<double,double> >& intervals_  ,  
									   std::vector< std::vector<double> > filter,
									   size_t number_of_pixels , double min_ , double max_ )
{	
    bool dbg = false;       
    if ( dbg ){std::cerr << "Entering construct procedure \n";getchar();}

    if ( min_ == max_ )
    {
        //in this case, we want the program to set up the min_ and max_ values by itself.
        min_ = std::numeric_limits<int>::max();
        max_ = -std::numeric_limits<int>::max();
        
        
        for ( size_t i = 0 ; i != intervals_.size() ; ++i )
        {
			if ( intervals_[i].first < min_ )min_ = intervals_[i].first;
			if ( intervals_[i].second > max_ )max_ = intervals_[i].second;
		}
        //now we have the structure filled in, and moreover we know min_ and max_ values of the interval, so we know the range.

        //add some more space:
        min_ -= fabs(max_ - min_)/100;
        max_ += fabs(max_ - min_)/100;
    }

	if ( dbg )
	{
		std::cerr << "min_ : " << min_ << std::endl;
		std::cerr << "max_ : " << max_ << std::endl;
		std::cerr << "number_of_pixels : " << number_of_pixels << std::endl;
		getchar();
	}

    this->min_ = min_;
    this->max_ = max_;    



    //initialization of the structure heat_map
    std::vector< std::vector<double> > heat_map_;
    for ( size_t i = 0 ; i != number_of_pixels ; ++i )
    {       
		std::vector<double> v( number_of_pixels , 0 );
        heat_map_.push_back( v );
    }
    this->heat_map = heat_map_;

    if (dbg)std::cerr << "Done creating of the heat map, now we will fill in the structure \n";

	for ( size_t pt_nr = 0 ; pt_nr != intervals_.size() ; ++pt_nr )
	{
		//compute the value of intervals_[pt_nr] in the grid:
		int x_grid = (int)((intervals_[pt_nr].first - this->min_)/( this->max_-this->min_ )*number_of_pixels);
		int y_grid = (int)((intervals_[pt_nr].second - this->min_)/( this->max_-this->min_ )*number_of_pixels);
		
		if ( dbg )
		{
			std::cerr << "point : " << intervals_[pt_nr].first << " , " << intervals_[pt_nr].second << std::endl;
			std::cerr << "x_grid : " << x_grid << std::endl;
			std::cerr << "y_grid : " << y_grid << std::endl;
		}
		
		//x_grid and y_grid gives a center of the kernel. We want to have its lower left cordner. To get this, we need to shift x_grid and y_grid by a grid diameter.		
		x_grid -= filter.size()/2;
		y_grid -= filter.size()/2;
		//note that the numbers x_grid and y_grid may be negative. 
		
		if ( dbg )
		{
			std::cerr << "After shift : \n";;
			std::cerr << "x_grid : " << x_grid << std::endl;
			std::cerr << "y_grid : " << y_grid << std::endl;
			std::cerr << "filter.size() : " << filter.size() << std::endl;
			getchar();
		}					

		
		for ( size_t i = 0 ; i != filter.size() ; ++i )
		{
			for ( size_t j = 0 ; j != filter.size() ; ++j )
			{			
				//if the point (x_grid+i,y_grid+j) is the correct point in the grid.						
				if ( 
					  ((x_grid+i)>=0) && (x_grid+i<this->heat_map.size()) 
					  &&
					  ((y_grid+j)>=0) && (y_grid+j<this->heat_map.size()) 
				   )
				{
					if ( dbg ){std::cerr << y_grid+j << " " <<  x_grid+i << std::endl;}
					this->heat_map[ y_grid+j ][ x_grid+i ] += filter[i][j];
					this->heat_map[ x_grid+i ][ y_grid+j ] += -filter[i][j];
				}
			}
		}
		
	}
}//construct


#endif

}//namespace Gudhi_stat
}//namespace Gudhi

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

#ifndef Persistence_heat_maps_H
#define Persistence_heat_maps_H
 
//standard include
#include <vector>
#include <sstream>
#include <iostream>
#include <cmath>
#include <limits>
#include <vector>
#include <algorithm>

//gudhi include
#include <gudhi/read_persistence_from_file.h>
#include <gudhi/common_gudhi_stat.h>




namespace Gudhi
{
namespace Gudhi_stat
{


/**
 * This is a simple procedure to create n by n (or 2*pixel_radius times 2*pixel_radius cubical approximation of a Gaussian kernel.
**/ 
std::vector< std::vector<double> > create_Gaussian_filter( size_t pixel_radius , double sigma )
{
	bool dbg = false;
	//we are computing the kernel mask to 2 standard deviations away from the center. We discretize it in a grid of a size 2*pixel_radius times 2*pixel_radius.
	
    double r = 0;
    double sigma_sqr = sigma * sigma;
 
    // sum is for normalization
    double sum = 0;
    
    //initialization of a kernel:
    std::vector< std::vector<double> > kernel( 2*pixel_radius +1 );
    for ( size_t i = 0 ; i != kernel.size() ; ++i )
    {
		std::vector<double> v( 2*pixel_radius +1 , 0 );
		kernel[i] = v;
	}
	
	if ( dbg )
	{
		std::cerr << "Kernel initalize \n";	
		std::cerr << "pixel_radius : " << pixel_radius << std::endl; 
		std::cerr << "kernel.size() : " << kernel.size() << std::endl;
		getchar();
	}
 
    for (int x = -pixel_radius; x <= (int)pixel_radius; x++)
    {
        for(int y = -pixel_radius; y <= (int)pixel_radius; y++)
        {			
			double real_x = 2*sigma*x/pixel_radius;
			double real_y = 2*sigma*y/pixel_radius;						
            r = sqrt(real_x*real_x + real_y*real_y);						
            kernel[x + pixel_radius][y + pixel_radius] = (exp(-(r*r)/sigma_sqr))/(3.141592 * sigma_sqr);
            sum += kernel[x + pixel_radius][y + pixel_radius];            
        }
    }      
 
    // normalize the kernel
    for( size_t i = 0; i != kernel.size() ; ++i)
    {
        for( size_t j = 0; j != kernel[i].size() ; ++j)
        {
            kernel[i][j] /= sum;
		}
            
    }
    
    if ( dbg )
    {
		std::cerr << "Here is the kernel : \n";
		for( size_t i = 0; i != kernel.size() ; ++i)
		{
			for( size_t j = 0; j != kernel[i].size() ; ++j)
			{
				std::cerr << kernel[i][j] << " ";
			}
			std::cerr << std::endl;
		}
	}
    return kernel;
}


/*
* There are various options to scale the poits depending on their location. One can for instance:
* (1) do nothing (scale all of them with the weight 1), as in the function constant_function
* (2) Scale them by the distance to the diagonal. This is implemented in function
* (3) Scale them with the square of their distance to diagonal. This is implemented in function
* (4) Scale them with 
*/ 


/**
 * This is one of a scaling functions used to weight poits depending on their persistence and/or location in the diagram. 
 * This particular functiona is a finction which always assign value 1 to a point in the diagram.
**/ 
class constant_scaling_function
{
public:
	double operator()( const std::pair< double , double >& point_in_diagram )
	{	
		return 1;
	}
};


/**
 * This is one of a scaling functions used to weight poits depending on their persistence and/or location in the diagram. 
 * The scaling given by this function to a point (b,d) is Euclidean distance of (b,d) from diagonal.
**/ 
class distance_from_diagonal_scaling
{
public:
	double operator()( const std::pair< double , double >& point_in_diagram )
	{
		//(point_in_diagram.first+point_in_diagram.second)/2.0
		return sqrt( pow((point_in_diagram.first-(point_in_diagram.first+point_in_diagram.second)/2.0),2) + pow((point_in_diagram.second-(point_in_diagram.first+point_in_diagram.second)/2.0),2)  );
	}
};

/**
 * This is one of a scaling functions used to weight poits depending on their persistence and/or location in the diagram. 
 * The scaling given by this function to a point (b,d) is a square of Euclidean distance of (b,d) from diagonal.
**/ 
class squared_distance_from_diagonal_scaling
{
public:
	double operator()( const std::pair< double , double >& point_in_diagram )
	{
		return pow((point_in_diagram.first-(point_in_diagram.first+point_in_diagram.second)/2.0),2) + pow((point_in_diagram.second-(point_in_diagram.first+point_in_diagram.second)/2.0),2);
	}
};

/**
 * This is one of a scaling functions used to weight poits depending on their persistence and/or location in the diagram. 
 * The scaling given by this function to a point (b,d) is an arctan of a persistence of a point (i.e. arctan( b-d ).
**/ 
class arc_tan_of_persistence_of_point
{
public:
	double operator()( const std::pair< double , double >& point_in_diagram )
	{
		return atan( point_in_diagram.second - point_in_diagram.first );
	}
};

/**
 * This is one of a scaling functions used to weight poits depending on their persistence and/or location in the diagram. 
 * This scaling function do not only depend on a point (p,d) in the diagram, but it depends on the whole diagram. 
 * The longest persistence pair get a scaling 1. Any other pair get a scaling belong to [0,1], which is proportional 
 * to the persistence of that pair.
**/ 
class weight_by_setting_maximal_interval_to_have_length_one
{
public:
	weight_by_setting_maximal_interval_to_have_length_one( double len ):letngth_of_maximal_interval(len){}
	double operator()( const std::pair< double , double >& point_in_diagram )
	{
		return (point_in_diagram.second-point_in_diagram.first)/this->letngth_of_maximal_interval;
	}
private:
	double letngth_of_maximal_interval;
};


/**
 * This class implements the following concepts: Vectorized_topological_data, Topological_data_with_distances, Real_valued_topological_data, Topological_data_with_averages, Topological_data_with_scalar_product
**/ 
template <typename Scalling_of_kernels = constant_scaling_function>
class Persistence_heat_maps
{
public:
    /**
     * The default constructor. A scaling function from the diagonal is set up to a constant function. The image is not erased below the diagonal. The gaussian have diameter 5. 	 
    **/ 
	Persistence_heat_maps()
	{
		Scalling_of_kernels f;
		this->f = f;
		this->erase_below_diagonal = false;		
		this->min_ = this->max_ = 0;
		this->set_up_parameters_for_basic_classes();
	};
	
	/**
	 * Construction that takes at the input the following parameters:
	 * (1) A vector of pairs of doubles (representing persistence intervals). All other parameters are optional. They are:
	 * (2) a Gausian filter generated by create_Gaussian_filter filter (the default value of this vaiable is a Gaussian filter of a radius 5), 
	 * (3) a boolean value which determines if the area of image below diagonal should, or should not be erased (it will be erased by default). 
	 * (4) a number of pixels in each direction (set to 1000 by default). 
	 * (5) a min x and y value of points that are to be taken into account. By default it is set to std::numeric_limits<double>::max(), in which case the program compute the values based on the data,
	 * (6) a max x and y value of points that are to be taken into account. By default it is set to std::numeric_limits<double>::max(), in which case the program compute the values based on the data.
	**/  
    Persistence_heat_maps( const std::vector< std::pair< double,double > >  & interval , std::vector< std::vector<double> > filter = create_Gaussian_filter(5,1) , bool erase_below_diagonal = false , size_t number_of_pixels = 1000 , double min_ = std::numeric_limits<double>::max() , double max_ = std::numeric_limits<double>::max()  );
    
    /**
	 * Construction that takes at the input a name of a file with persistence intervals, a filter (radius 5 by default), a scaling function (constant by default), a boolean value which determines if the area of image below diagonal should, or should not be erased (should by default). The next parameter is the number of pixels in each direction (set to 1000 by default). and min and max values of images (both set to std::numeric_limits<double>::max() by defaulet. If this is the case, the program will pick the right values based on the data).
	**/  
	/**
	 * Construction that takes at the input the following parameters:
	 * (1) A a name of a file with persistence intervals. The file shold be readable by the function read_standard_persistence_file. All other parameters are optional. They are:
	 * (2) a Gausian filter generated by create_Gaussian_filter filter (the default value of this vaiable is a Gaussian filter of a radius 5), 	
	 * (3) a boolean value which determines if the area of image below diagonal should, or should not be erased (it will be erased by default). 
	 * (4) a number of pixels in each direction (set to 1000 by default). 
	 * (5) a min x and y value of points that are to be taken into account. By default it is set to std::numeric_limits<double>::max(), in which case the program compute the values based on the data,
	 * (6) a max x and y value of points that are to be taken into account. By default it is set to std::numeric_limits<double>::max(), in which case the program compute the values based on the data.
	**/ 	
    Persistence_heat_maps( const char* filename , std::vector< std::vector<double> > filter = create_Gaussian_filter(5,1), bool erase_below_diagonal = false , size_t number_of_pixels = 1000 , double min_ = std::numeric_limits<double>::max() , double max_ = std::numeric_limits<double>::max()  );

	
	/**
	 * Compute a mean value of a collection of heat maps and store it in the current object. Note that all the persistence maps send in a vector to this procedure need to have the same parameters. 
	 * If this is not the case, the program will throw an exception. 
	**/ 
	void compute_mean( const std::vector<Persistence_heat_maps*>& maps );
	
	/**
	 * Compute a median value of a collection of heat maps and store it in the current object. Note that all the persistence maps send in a vector to this procedure need to have the same parameters. 
	 * If this is not the case, the program will throw an exception. 
	**/ 
	void compute_median( const std::vector<Persistence_heat_maps*>& maps );	
	
	/**
	 * Compute a percentage of active (i.e) values above the cutoff of a collection of heat maps. 
	**/ 
	void compute_percentage_of_active( const std::vector<Persistence_heat_maps*>& maps , size_t cutoff = 1 );
    
    //put to file subroutine
    /**
     * The function outputs the perssitence image to a text file. The format as follow:
     * In the first line, the values min and max of the image are stored
     * In the next lines, we have the persistence images in a form of a bitmap image. 
    **/
    void print_to_file( const char* filename )const;
    
    /**
     * A function that load a heat map from file to the current object (and arase qhatever was stored in the current object before).
    **/
    void load_from_file( const char* filename );
    
    
    /**
     * The procedure checks if min_, max_ and this->heat_maps sizes are the same. 
    **/ 
    inline bool check_if_the_same( const Persistence_heat_maps& second  )const
    {
		bool dbg = false;
		if ( this->heat_map.size() != second.heat_map.size() )
		{
			if ( dbg )std::cerr << "this->heat_map.size() : " << this->heat_map.size()  << " \n second.heat_map.size() : " << second.heat_map.size() << std::endl;
			return false;
		}
		if ( this->min_ != second.min_ )
		{
			if ( dbg )std::cerr << "this->min_ : " << this->min_ << ", second.min_ : " << second.min_ << std::endl;
			return false;
		}
		if ( this->max_ != second.max_ )
		{
			if ( dbg )std::cerr << "this->max_ : " << this->max_ << ", second.max_ : " << second.max_ << std::endl;
			return false;
		}
		//in the other case we may assume that the persistence images are defined on the same domain.
		return true;
	}
	
	
	/**
	 * Return minimal range value of persistent image.
	**/ 
	inline double get_min()const{return this->min_;}

	/**
	 * Return maximal range value of persistent image.
	**/ 
	inline double get_max()const{return this->max_;}
	
	
	/**
	 * Operator == to check if to persistence heat maps are the same.
	**/ 
	bool operator == ( const Persistence_heat_maps& rhs )const
	{
		bool dbg = false;
		if ( !this->check_if_the_same(rhs) )
		{
			if ( dbg )std::cerr << "The domains are not the same \n";
			return false;//in this case, the domains are not the same, so the maps cannot be the same.
		}
		for ( size_t i = 0 ; i != this->heat_map.size() ; ++i )
		{
			for ( size_t j = 0 ; j != this->heat_map[i].size() ; ++j )
			{
				if ( !almost_equal(this->heat_map[i][j] , rhs.heat_map[i][j]) )
				{
					if ( dbg )
					{
						std::cerr << "this->heat_map[" << i << "][" << j << "] = " << this->heat_map[i][j] << std::endl;
						std::cerr << "rhs.heat_map[" << i << "][" << j << "] = " << rhs.heat_map[i][j] << std::endl;						
					}
					return false;
				}
			}
		}
		return true;
	}
	
	/**
	 * Operator != to check if to persistence heat maps are different.
	**/ 
	bool operator != ( const Persistence_heat_maps& rhs )const
	{
		return !( (*this) == rhs );
	}
    
    
    /**
     * A function to generate a gnuplot script to vizualize the persistent image.  
    **/
    void plot( const char* filename )const;
    
       
    template<typename Operation_type>
    friend Persistence_heat_maps operation_on_pair_of_heat_maps( const Persistence_heat_maps& first ,  const Persistence_heat_maps& second , Operation_type operation )
    {
		//first check if the heat maps are compatible 
		if ( !first.check_if_the_same( second ) )
		{
			std::cerr << "Sizes of the heat maps are not compatible. The program will now terminate \n";
			throw "Sizes of the heat maps are not compatible. The program will now terminate \n";
		}
		Persistence_heat_maps result;
		result.min_ = first.min_;
		result.max_ = first.max_;
		result.heat_map.reserve( first.heat_map.size() );
		for ( size_t i = 0 ; i != first.heat_map.size() ; ++i )
		{ 
			std::vector< double > v;
			v.reserve( first.heat_map[i].size() );
			for ( size_t j = 0 ; j != first.heat_map[i].size() ; ++j )
			{
				v.push_back( operation( first.heat_map[i][j] , second.heat_map[i][j] ) );
			}
			result.heat_map.push_back( v );
		}
		return result;
	}//operation_on_pair_of_heat_maps
	
	
	/**
	 * Multiplication of Persistence_heat_maps by scalar (so that all values of the heat map gets multiplied by that scalar). 
	**/ 
	Persistence_heat_maps multiply_by_scalar( double scalar )const 
    {
		Persistence_heat_maps result;
		result.min_ = this->min_;
		result.max_ = this->max_;
		result.heat_map.reserve( this->heat_map.size() );
		for ( size_t i = 0 ; i != this->heat_map.size() ; ++i )
		{ 
			std::vector< double > v;
			v.reserve( this->heat_map[i].size() );
			for ( size_t j = 0 ; j != this->heat_map[i].size() ; ++j )
			{
				v.push_back( this->heat_map[i][j] * scalar );
			}
			result.heat_map.push_back( v );
		}
		return result;
	}
	
    /**
     * This function computes a sum of two objects of a type Persistence_heat_maps.
    **/    
    friend Persistence_heat_maps operator+( const Persistence_heat_maps& first , const Persistence_heat_maps& second )
    {
		return operation_on_pair_of_heat_maps( first , second , std::plus<double>() );
	}	
	/**
     * This function computes a difference of two objects of a type Persistence_heat_maps.
    **/
    friend Persistence_heat_maps operator-( const Persistence_heat_maps& first , const Persistence_heat_maps& second )
    {
		return operation_on_pair_of_heat_maps( first , second , std::minus<double>() );
	}
	/**
     * This function computes a product of an object of a type Persistence_heat_maps with real number.
    **/
	friend Persistence_heat_maps operator*( double scalar , const Persistence_heat_maps& A )
	{
		return A.multiply_by_scalar( scalar );
	}
	/**
     * This function computes a product of an object of a type Persistence_heat_maps with real number.
    **/
	friend Persistence_heat_maps operator*( const Persistence_heat_maps& A , double scalar )
	{
		return A.multiply_by_scalar( scalar );
	}
	/**
     * This function computes a product of an object of a type Persistence_heat_maps with real number.
    **/
	Persistence_heat_maps operator*( double scalar )
	{
		return this->multiply_by_scalar( scalar );
	}
	/**
	 * += operator for Persistence_heat_maps.
	**/ 
    Persistence_heat_maps operator += ( const Persistence_heat_maps& rhs )
    {
        *this = *this + rhs;
        return *this;
    }    
    /**
	 * -= operator for Persistence_heat_maps.
	**/ 
    Persistence_heat_maps operator -= ( const Persistence_heat_maps& rhs )
    {
        *this = *this - rhs;
        return *this;
    }
    /**
	 * *= operator for Persistence_heat_maps.
	**/ 
    Persistence_heat_maps operator *= ( double x )
    {
        *this = *this*x;
        return *this;
    }
    /**
	 * /= operator for Persistence_heat_maps.
	**/ 
    Persistence_heat_maps operator /= ( double x )
    {
        if ( x == 0 )throw( "In operator /=, division by 0. Program terminated." );
        *this = *this * (1/x);
        return *this;
    }
    
    
    //Implementations of functions for various concepts.  
    
    /**
	 * This function produce a vector of doubles based on a persisence heat map. It is required in a concept Vectorized_topological_data
	*/
    std::vector<double> vectorize( int number_of_function )const;    
    /**
	 * This function return the number of functions that allows vectorization of persistence heat map. It is required in a concept Vectorized_topological_data.
	 **/ 
	size_t number_of_vectorize_functions()const
	{
		return this->number_of_functions_for_vectorization;	
	}
	
	/**	 
	 * This function is required by the Real_valued_topological_data concept. It returns various projections od the persistence heat map to a real line.
	**/ 	
	double project_to_R( int number_of_function )const;	
	/**
	 * The function gives the number of possible projections to R. This function is required by the Real_valued_topological_data concept.
	**/ 
	size_t number_of_projections_to_R()const
	{
		return this->number_of_functions_for_projections_to_reals;
	}
	
	/**
	* A function to compute distance between persistence heat maps.
	* The parameter of this function is a const reference to an object of a class Persistence_heat_maps.
	* This function is required in Topological_data_with_distances concept.
    * For max norm distance, set power to std::numeric_limits<double>::max()
	**/
	double distance( const Persistence_heat_maps& second_ , double power = 1)const;
	
	/**
	 * A function to compute averaged persistence heat map, based on vector of persistence heat maps.
	 * This function is required by Topological_data_with_averages concept.
	**/ 	
	void compute_average( const std::vector< Persistence_heat_maps* >& to_average );	
	
	/**
	* A function to compute scalar product of persistence heat maps.
	* The parameter of this functionis a const reference to an object of a class Persistence_heat_maps.
	* This function is required in Topological_data_with_scalar_product concept.
	**/
	double compute_scalar_product( const Persistence_heat_maps& second_ )const;
	
	//end of implementation of functions needed for concepts.
	
	
	/**
	 * The x-range of the persistence heat map. 
	**/ 
	std::pair< double , double > get_x_range()const
	{
		return std::make_pair( this->min_ , this->max_ );
	}
	
	/**
	 * The y-range of the persistence heat map. 
	**/ 
	std::pair< double , double > get_y_range()const
	{
		return this->get_x_range();
	}
	
	
	
	
protected:
	//private methods
	std::vector< std::vector<double> > check_and_initialize_maps( const std::vector<Persistence_heat_maps*>& maps );
	size_t number_of_functions_for_vectorization;
	size_t number_of_functions_for_projections_to_reals;
    void construct( const std::vector< std::pair<double,double> >& intervals_  , 
					std::vector< std::vector<double> > filter = create_Gaussian_filter(5,1),
					                    
                    bool erase_below_diagonal = false , size_t number_of_pixels = 1000 , double min_ = std::numeric_limits<double>::max() , double max_ = std::numeric_limits<double>::max() );
                    
	void set_up_parameters_for_basic_classes()
	{
		this->number_of_functions_for_vectorization = 1;
		this->number_of_functions_for_projections_to_reals = 1;
	}
    
    //data    
    //double (*scalling_function_with_respect_to_distance_from_diagonal)( const std::pair< double , double >& point_in_diagram );    
    Scalling_of_kernels f;
    bool erase_below_diagonal;
    double min_;
    double max_;
    std::vector< std::vector< double > > heat_map;    
};


//if min_ == max_, then the program is requested to set up the values itself based on persistence intervals
template <typename Scalling_of_kernels>
void Persistence_heat_maps<Scalling_of_kernels>::construct( const std::vector< std::pair<double,double> >& intervals_  ,  
									   std::vector< std::vector<double> > filter,									   
									   bool erase_below_diagonal , size_t number_of_pixels , double min_ , double max_ )
{
    bool dbg = false;       
    if ( dbg )std::cerr << "Entering construct procedure \n";
    Scalling_of_kernels f;
    this->f = f;

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
		}
				
		double scaling_value = this->f(intervals_[pt_nr]);

		
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
					this->heat_map[ y_grid+j ][ x_grid+i ] += scaling_value * filter[i][j];
					if ( dbg )
					{
						std::cerr << "Position : (" << x_grid+i << "," << y_grid+j << ") got increased by the value : " << filter[i][j] << std::endl;
					}					
				}
			}
		}
		
	}
	
	//now it remains to cut everything below diagonal if the user wants us to. 
	if ( erase_below_diagonal )
	{
		for ( size_t i = 0 ; i != this->heat_map.size() ; ++i )
		{
			for ( size_t j = i ; j != this->heat_map.size() ; ++j )
			{
				this->heat_map[i][j] = 0;
			}
		}
	}
}//construct

template <typename Scalling_of_kernels>
Persistence_heat_maps<Scalling_of_kernels>::Persistence_heat_maps( const std::vector< std::pair< double,double > >  & interval , 
											  std::vector< std::vector<double> > filter,  
											  bool erase_below_diagonal , size_t number_of_pixels , double min_ , double max_ )
{
    this->construct( interval ,  filter , erase_below_diagonal , number_of_pixels , min_ , max_ );
    this->set_up_parameters_for_basic_classes();
}


template <typename Scalling_of_kernels>
Persistence_heat_maps<Scalling_of_kernels>::Persistence_heat_maps( const char* filename , 
											  std::vector< std::vector<double> > filter, 											  
											  bool erase_below_diagonal , size_t number_of_pixels  , double min_ , double max_ )
{    
    std::vector< std::pair< double , double > > intervals_ = read_standard_persistence_file( filename );           
    this->construct( intervals_ ,  filter, erase_below_diagonal , number_of_pixels , min_ , max_ );
    this->set_up_parameters_for_basic_classes();
}


template <typename Scalling_of_kernels>
std::vector< std::vector<double> > Persistence_heat_maps<Scalling_of_kernels>::check_and_initialize_maps( const std::vector<Persistence_heat_maps*>& maps )
{
	//checking if all the heat maps are of the same size:
	for ( size_t i = 0 ; i != maps.size() ; ++i )
    {
		if ( maps[i]->heat_map.size() != maps[0]->heat_map.size() )
		{
			std::cerr << "Sizes of Persistence_heat_maps are not compatible. The program will terminate now \n";
			throw "Sizes of Persistence_heat_maps are not compatible. The program will terminate now \n";
		}
		if ( maps[i]->heat_map[0].size() != maps[0]->heat_map[0].size() )
		{
			std::cerr << "Sizes of Persistence_heat_maps are not compatible. The program will terminate now \n";
			throw "Sizes of Persistence_heat_maps are not compatible. The program will terminate now \n";			
		}
	} 	
	std::vector< std::vector<double> > heat_maps( maps[0]->heat_map.size() );
	for ( size_t i = 0 ; i != maps[0]->heat_map.size() ; ++i )
    {
		 std::vector<double> v( maps[0]->heat_map[0].size() , 0 );
		 heat_maps[i] = v;
	}
	return heat_maps;
}

template <typename Scalling_of_kernels>
void Persistence_heat_maps<Scalling_of_kernels>::compute_median( const std::vector<Persistence_heat_maps*>& maps )
{
	std::vector< std::vector<double> > heat_maps = this->check_and_initialize_maps( maps );
	
	std::vector<double> to_compute_median( maps.size() );
    for ( size_t i = 0 ; i != heat_maps.size() ; ++i )
    {
        for ( size_t j = 0 ; j != heat_maps[i].size() ; ++j )
        {
			for ( size_t map_no = 0 ; map_no != maps.size() ; ++map_no )
			{
				to_compute_median[map_no] = maps[map_no]->heat_map[i][j];
			}
			std::nth_element(to_compute_median.begin(), to_compute_median.begin() + to_compute_median.size()/2, to_compute_median.end());
            heat_maps[i][j] = to_compute_median[ to_compute_median.size()/2 ];
        }       
    }
	this->heat_map = heat_maps;
	this->min_= maps[0]->min_;
	this->max_= maps[0]->max_;
}


template <typename Scalling_of_kernels>
void Persistence_heat_maps<Scalling_of_kernels>::compute_mean( const std::vector<Persistence_heat_maps*>& maps )
{
	std::vector< std::vector<double> > heat_maps = this->check_and_initialize_maps( maps );
		    for ( size_t i = 0 ; i != heat_maps.size() ; ++i )
    {
        for ( size_t j = 0 ; j != heat_maps[i].size() ; ++j )
        {
			double mean = 0;
			for ( size_t map_no = 0 ; map_no != maps.size() ; ++map_no )
			{	
				mean += maps[map_no]->heat_map[i][j];
			}		
            heat_maps[i][j] =  mean/(double)maps.size();
        }       
    }
	this->heat_map = heat_maps;
	this->min_ = maps[0]->min_;
	this->max_ = maps[0]->max_;
}



template <typename Scalling_of_kernels>
void Persistence_heat_maps<Scalling_of_kernels>::compute_percentage_of_active( const std::vector<Persistence_heat_maps*>& maps , size_t cutoff  )
{
	std::vector< std::vector<double> > heat_maps = this->check_and_initialize_maps( maps );

    for ( size_t i = 0 ; i != heat_maps.size() ; ++i )
    {
        for ( size_t j = 0 ; j != heat_maps[i].size() ; ++j )
        {
            size_t number_of_active_levels = 0;
            for ( size_t map_no = 0 ; map_no != maps.size() ; ++map_no )
            {
                if ( maps[map_no]->heat_map[i][j] ) number_of_active_levels++;
            }
            if ( number_of_active_levels > cutoff )
            {
                heat_maps[i][j] = number_of_active_levels;
            }
            else
            {
                heat_maps[i][j] = 0;
            }
        }       
    }
    this->heat_map = heat_maps;
    this->min_ = maps[0]->min_;
	this->max_ = maps[0]->max_;
}


template <typename Scalling_of_kernels>
void Persistence_heat_maps<Scalling_of_kernels>::plot( const char* filename )const
{
	std::ofstream out;
	std::stringstream ss;
	ss << filename << "_GnuplotScript";
	
    out.open( ss.str().c_str() );
	out << "plot      '-' matrix with image" << std::endl;
    for ( size_t i = 0 ; i != this->heat_map.size() ; ++i )
    {
        for ( size_t j = 0 ; j != this->heat_map[i].size() ; ++j )
        {
            out << this->heat_map[i][j] << " ";
        }
        out << std::endl;
    }
    out.close();
    std::cout << "Gnuplot script have been created. Open gnuplot and type load \'" << ss.str().c_str() << "\' to see the picture." << std::endl; 
}


template <typename Scalling_of_kernels>
void Persistence_heat_maps<Scalling_of_kernels>::print_to_file( const char* filename )const
{

	std::ofstream out;
	out.open( filename );
	
	//First we store this->min_ and this->max_ values:
	out << this->min_ << " " << this->max_ << std::endl;	
	for ( size_t i = 0 ; i != this->heat_map.size() ; ++i )
    {
        for ( size_t j = 0 ; j != this->heat_map[i].size() ; ++j )
        {
            out << this->heat_map[i][j] << " ";
        }
        out << std::endl;
    }
	out.close();
}

template <typename Scalling_of_kernels>
void Persistence_heat_maps<Scalling_of_kernels>::load_from_file( const char* filename )
{
	bool dbg = false;	
	
	std::ifstream in;
	in.open( filename );
	
	//checking if the file exist / if it was open. 
	if ( !( access( filename, F_OK ) != -1 ) )
	{
		std::cerr << "The file : " << filename << " do not exist. The program will now terminate \n";
		throw "The file from which you are trying to read the persistence landscape do not exist. The program will now terminate \n";
	}
	
	//now we read the file one by one. 
	
	 
	 
	in >> this->min_ >> this->max_;	    
	if ( dbg )
	{
		std::cerr << "Reading the following values of min and max : " << this->min_ << " , " << this->max_ << std::endl;
	}

	std::string temp;
	std::getline(in, temp);
         
	while (!in.eof())
	{
		std::getline(in, temp);
        std::stringstream lineSS;
        lineSS << temp;
        
        std::vector<double> line_of_heat_map;    
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
			getchar();
		}

		if ( in.good() )this->heat_map.push_back( line_of_heat_map );
	}	
	in.close();
	if ( dbg )std::cout << "Done \n";
}


//Concretizations of virtual methods:   
template <typename Scalling_of_kernels> 
std::vector<double> Persistence_heat_maps<Scalling_of_kernels>::vectorize( int number_of_function )const
{
	//convert this->heat_map into one large vector:
	size_t size_of_result = 0;
	for ( size_t i = 0 ; i != this->heat_map.size() ; ++i )
	{
		size_of_result += this->heat_map[i].size();
	}
	
	std::vector< double > result;
	result.reserve( size_of_result );
	
	for ( size_t i = 0 ; i != this->heat_map.size() ; ++i )
	{
		for ( size_t j = 0 ; j != this->heat_map[i].size() ; ++j )
		{
			result.push_back( this->heat_map[i][j] );
		}
	}
	
	return result;
}

template <typename Scalling_of_kernels>
double Persistence_heat_maps<Scalling_of_kernels>::distance( const Persistence_heat_maps& second , double power )const
{			
	//first we need to check if (*this) and second are defined on the same domain and have the same dimensions:
	if ( !this->check_if_the_same(second) )
	{
		std::cerr << "The persistence images are of noncompatible sizes. We cannot therefore compute distance between them. The program will now terminate";
		throw "The persistence images are of noncompatible sizes. We cannot therefore compute distance between them. The program will now terminate";
	}
	
	//if we are here, we know that the two persistence iomages are defined on the same domain, so we can start computing their distances:
	
	double distance = 0;
	if ( power < std::numeric_limits<double>::max() )
	{
		for ( size_t i = 0 ; i != this->heat_map.size() ; ++i )
		{
			for ( size_t j = 0 ; j != this->heat_map[i].size() ; ++j )
			{
				distance += pow( fabs(this->heat_map[i][j] - second.heat_map[i][j]) , power );
			}
		}
	}
	else
	{		
		//in this case, we compute max norm distance
		for ( size_t i = 0 ; i != this->heat_map.size() ; ++i )
		{
			for ( size_t j = 0 ; j != this->heat_map[i].size() ; ++j )
			{
				if ( distance < fabs(this->heat_map[i][j] - second.heat_map[i][j]) )
				{
					distance = fabs(this->heat_map[i][j] - second.heat_map[i][j]);
				}
			}
		}
	}
	return distance;
}

template <typename Scalling_of_kernels>
double Persistence_heat_maps<Scalling_of_kernels>::project_to_R( int number_of_function )const
{
	double result = 0;
	for ( size_t i = 0 ; i != this->heat_map.size() ; ++i )
	{
		for ( size_t j = 0 ; j != this->heat_map[i].size() ; ++j )
		{
			result += this->heat_map[i][j];
		}
	}
	return result;
}

template <typename Scalling_of_kernels>
void Persistence_heat_maps<Scalling_of_kernels>::compute_average( const std::vector< Persistence_heat_maps* >& to_average )
{				
	this->compute_mean( to_average );
}

template <typename Scalling_of_kernels>
double Persistence_heat_maps<Scalling_of_kernels>::compute_scalar_product( const Persistence_heat_maps& second )const
{		
	//first we need to check if (*this) and second are defined on the same domain and have the same dimensions:
	if ( !this->check_if_the_same(second) )
	{
		std::cerr << "The persistence images are of noncompatible sizes. We cannot therefore compute distance between them. The program will now terminate";
		throw "The persistence images are of noncompatible sizes. We cannot therefore compute distance between them. The program will now terminate";
	}
	
	//if we are here, we know that the two persistence iomages are defined on the same domain, so we can start computing their scalar product:
	double scalar_prod = 0;
	for ( size_t i = 0 ; i != this->heat_map.size() ; ++i )
	{
		for ( size_t j = 0 ; j != this->heat_map[i].size() ; ++j )
		{
			scalar_prod += this->heat_map[i][j]*second.heat_map[i][j];
		}
	}
	return scalar_prod;
}




}//namespace Gudhi_stat
}//namespace Gudhi


#endif

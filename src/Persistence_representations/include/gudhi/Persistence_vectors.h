/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2017  INRIA (France)
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


#ifndef PERSISTENCE_VECTORS_H_
#define PERSISTENCE_VECTORS_H_

#include <fstream>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <limits>
#include <functional>

//gudhi include
#include <gudhi/read_persistence_from_file.h>
#include <gudhi/common_persistence_representations.h>
#include <gudhi/distance_functions.h>


namespace Gudhi 
{
namespace Persistence_representations 
{

template <typename T>
struct Maximum_distance
{
    double operator() ( const std::pair< T,T >& f , const std::pair<T,T>& s )
    {
        return  std::max( fabs( f.first - s.first ) , fabs( f.second - s.second ) );
    }
};



/**
 * \class Vector_distances_in_diagram Persistence_vectors.h gudhi/Persistence_vectors.h
 * \brief A class implementing persistence vectors.
 *
 * \ingroup Persistence_representations
 *
 * \details
 * This is an implementation of idea presented in the paper <i>Stable Topological Signatures for Points on 3D
 * Shapes</i> \cite Carriere_Oudot_Ovsjanikov_top_signatures_3d .<br>
 * The parameter of the class is the class that computes distance used to construct the vectors. The typical function
 * is either Euclidean of maximum (Manhattan) distance.
 *
 * This class implements the following concepts: Vectorized_topological_data, Topological_data_with_distances,
 * Real_valued_topological_data, Topological_data_with_averages, Topological_data_with_scalar_product
 **/
template <typename F>
class Vector_distances_in_diagram 
{
public:
	/**
	* The default constructor.  
	**/
    Vector_distances_in_diagram(){};
        
    /**
	* The constructor that takes as an input a multiset of persistence intervals (given as vector of birth-death pairs). The second parameter is the desired length of the output vectors. 
	**/
    Vector_distances_in_diagram( const std::vector< std::pair< double , double > >& intervals , size_t where_to_cut );
    
    /**
	* The constructor taking as an input a file with birth-death pairs. The second parameter is the desired length of the output vectors. 
	**/
    Vector_distances_in_diagram( const char* filename , size_t where_to_cut , unsigned dimension = std::numeric_limits<unsigned>::max() );

	
	/**
	 *  Writing to a stream. 
	**/
	template <typename K>
    friend std::ostream& operator << ( std::ostream& out , const Vector_distances_in_diagram<K>& d )
    {
        for ( size_t i = 0 ; i != std::min( d.sorted_vector_of_distances.size() , d.where_to_cut) ; ++i )
        {
            out << d.sorted_vector_of_distances[i] << " ";
        }
        return out;
    }

	/**
	* This procedure gives the value of a vector on a given position.  
	**/ 
    inline double vector_in_position( size_t position )const
    {
        if ( position >= this->sorted_vector_of_distances.size() )throw("Wrong position in accessing Vector_distances_in_diagram::sorted_vector_of_distances\n");
        return this->sorted_vector_of_distances[position];
    }

	/**
	 * Return a size of a vector. 
	**/
    inline size_t size()const{return this->sorted_vector_of_distances.size();}
    
    /**
	 * Write a vector to a file.
	**/
    void write_to_file( const char* filename )const;
    
     /**
	 * Write a vector to a file.
	**/
    void print_to_file( const char* filename )const
    {
		this->write_to_file( filename );
	}  
    
    /**
	 * Loading a vector to a file.
	**/
    void load_from_file( const char* filename );
    
    /**
     * Comparison operators:
    **/ 
    bool operator == ( const Vector_distances_in_diagram& second )const
    {
		if ( this->sorted_vector_of_distances.size() != second.sorted_vector_of_distances.size() )return false;
		for ( size_t i = 0 ; i != this->sorted_vector_of_distances.size() ; ++i )
		{
			if ( !almost_equal(this->sorted_vector_of_distances[i] , second.sorted_vector_of_distances[i]) )return false;
		}
		return true;
	}
	
	bool operator != ( const Vector_distances_in_diagram& second )const
    {
		return !( *this == second );
	}
    
   //Implementations of functions for various concepts.
     /**
     * Compute projection to real numbers of persistence vector. This function is required by the Real_valued_topological_data concept
     * At the moment this function is not tested, since it is quite likely to be changed in the future. Given this, when using it, keep in mind that it
     * will be most likely changed in the next versions.
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
     * Compute a vectorization of a persistent vectors. It is required in a concept Vectorized_topological_data.
    **/
    std::vector<double> vectorize( int number_of_function )const;
     /**
	 * This function return the number of functions that allows vectorization of a persistence vector. It is required in a concept Vectorized_topological_data.
	 **/ 
	size_t number_of_vectorize_functions()const
	{
		return this->number_of_functions_for_vectorization;	
	}
    
    /**
     * Compute a average of two persistent vectors. This function is required by Topological_data_with_averages concept.
    **/
    void compute_average( const std::vector< Vector_distances_in_diagram* >& to_average );
    
    /**
     * Compute a distance of two persistent vectors. This function is required in Topological_data_with_distances concept.
     * For max norm distance, set power to std::numeric_limits<double>::max()
    **/
    double distance( const Vector_distances_in_diagram& second , double power = 1)const;
    
    /**
     * Compute a scalar product of two persistent vectors. This function is required in Topological_data_with_scalar_product concept.
    **/ 
    double compute_scalar_product( const Vector_distances_in_diagram& second )const;        
    //end of implementation of functions needed for concepts.
    
    
    /**
    * For visualization use output from vectorize and build histograms. 
    **/
    std::vector< double > output_for_visualization()const
    {
		return this->sorted_vector_of_distances;
	}
	
	
	/**
	 * Create a gnuplot script to visualize the data structure. 
	 **/ 
	void plot( const char* filename )const
	{
		std::stringstream gnuplot_script;
		gnuplot_script << filename << "_GnuplotScript";
		std::ofstream out;
		out.open( gnuplot_script.str().c_str() );
		out << "set style data histogram" << std::endl;
		out << "set style histogram cluster gap 1" << std::endl;
		out << "set style fill solid border -1" << std::endl;
		out << "plot '-' notitle" << std::endl;    
		for ( size_t i = 0 ; i != this->sorted_vector_of_distances.size() ; ++i )
		{
			out << this->sorted_vector_of_distances[i]  << std::endl;
		}		
		out <<std::endl;
		out.close();
		std::cout << "To visualize, open gnuplot and type: load \'" << gnuplot_script.str().c_str() << "\'" <<  std::endl;			
	}
	
	/**
	 * The x-range of the persistence vector. 
	**/ 
	std::pair< double , double > get_x_range()const
	{
		return std::make_pair( 0 , this->sorted_vector_of_distances.size() );
	}
	
	/**
	 * The y-range of the persistence vector. 
	**/ 
	std::pair< double , double > get_y_range()const
	{
		if ( this->sorted_vector_of_distances.size() == 0 )return std::make_pair(0,0);
		return std::make_pair( this->sorted_vector_of_distances[0] , 0);
	}
	
	//arithmetic operations:		
	template < typename Operation_type > 
	friend Vector_distances_in_diagram operation_on_pair_of_vectors( const Vector_distances_in_diagram& first ,  const Vector_distances_in_diagram& second , Operation_type opertion )
    {
		Vector_distances_in_diagram result;		
		//Operation_type operation;		
		result.sorted_vector_of_distances.reserve(std::max( first.sorted_vector_of_distances.size() , second.sorted_vector_of_distances.size() ) );
        for ( size_t i = 0 ; i != std::min( first.sorted_vector_of_distances.size() , second.sorted_vector_of_distances.size() ) ; ++i )
        {
			result.sorted_vector_of_distances.push_back( opertion( first.sorted_vector_of_distances[i] , second.sorted_vector_of_distances[i]) );
		}
		if ( first.sorted_vector_of_distances.size() == std::min( first.sorted_vector_of_distances.size() , second.sorted_vector_of_distances.size() ) )
		{
			for ( size_t i = std::min( first.sorted_vector_of_distances.size() , second.sorted_vector_of_distances.size() ) ; 
			i != std::max( first.sorted_vector_of_distances.size() , second.sorted_vector_of_distances.size() ) ; ++i )
			{
				result.sorted_vector_of_distances.push_back( opertion(0,second.sorted_vector_of_distances[i]) );
			}
		}
		else
		{
			for ( size_t i = std::min( first.sorted_vector_of_distances.size() , second.sorted_vector_of_distances.size() ) ; 
			i != std::max( first.sorted_vector_of_distances.size() , second.sorted_vector_of_distances.size() ) ; ++i )
			{
				result.sorted_vector_of_distances.push_back( opertion(first.sorted_vector_of_distances[i],0) );
			}
		}		
		return result;
    }//operation_on_pair_of_vectors
    
    /**
     * This function implements an operation of multiplying Vector_distances_in_diagram by a scalar. 
    **/ 
    Vector_distances_in_diagram multiply_by_scalar( double scalar )const 
    {
		Vector_distances_in_diagram result;
		result.sorted_vector_of_distances.reserve( this->sorted_vector_of_distances.size() );
        for ( size_t i = 0 ; i != this->sorted_vector_of_distances.size() ; ++i )
        {
			result.sorted_vector_of_distances.push_back( scalar * this->sorted_vector_of_distances[i] );
		}		
		return result;
    }//multiply_by_scalar
    
    
    
    /**
     * This function computes a sum of two objects of a type Vector_distances_in_diagram.
    **/    
    friend Vector_distances_in_diagram operator+( const Vector_distances_in_diagram& first , const Vector_distances_in_diagram& second )
    {
		return operation_on_pair_of_vectors( first , second , std::plus<double>() );
	}	
	/**
     * This function computes a difference of two objects of a type Vector_distances_in_diagram.
    **/
    friend Vector_distances_in_diagram operator-( const Vector_distances_in_diagram& first , const Vector_distances_in_diagram& second )
    {
		return operation_on_pair_of_vectors( first , second , std::minus<double>() );
	}
	/**
     * This function computes a product of an object of a type Vector_distances_in_diagram with real number.
    **/
	friend Vector_distances_in_diagram operator*( double scalar , const Vector_distances_in_diagram& A )
	{
		return A.multiply_by_scalar( scalar );
	}
	/**
     * This function computes a product of an object of a type Vector_distances_in_diagram with real number.
    **/
	friend Vector_distances_in_diagram operator*( const Vector_distances_in_diagram& A , double scalar )
	{
		return A.multiply_by_scalar( scalar );
	}
	/**
     * This function computes a product of an object of a type Vector_distances_in_diagram with real number.
    **/
	Vector_distances_in_diagram operator*( double scalar )
	{
		return this->multiply_by_scalar( scalar );
	}
	/**
	 * += operator for Vector_distances_in_diagram.
	**/ 
    Vector_distances_in_diagram operator += ( const Vector_distances_in_diagram& rhs )
    {
        *this = *this + rhs;
        return *this;
    }    
    /**
	 * -= operator for Vector_distances_in_diagram.
	**/ 
    Vector_distances_in_diagram operator -= ( const Vector_distances_in_diagram& rhs )
    {
        *this = *this - rhs;
        return *this;
    }
    /**
	 * *= operator for Vector_distances_in_diagram.
	**/ 
    Vector_distances_in_diagram operator *= ( double x )
    {
        *this = *this*x;
        return *this;
    }
    /**
	 * /= operator for Vector_distances_in_diagram.
	**/ 
    Vector_distances_in_diagram operator /= ( double x )
    {
        if ( x == 0 )throw( "In operator /=, division by 0. Program terminated." );
        *this = *this * (1/x);
        return *this;
    }
    

private:
    std::vector< std::pair< double , double > > intervals;
    std::vector< double > sorted_vector_of_distances;
    size_t number_of_functions_for_vectorization;
	size_t number_of_functions_for_projections_to_reals;
	size_t where_to_cut;

    void compute_sorted_vector_of_distances_via_heap( size_t where_to_cut );
    void compute_sorted_vector_of_distances_via_vector_sorting( size_t where_to_cut );
    
    Vector_distances_in_diagram( const std::vector< double >& sorted_vector_of_distances_ ):sorted_vector_of_distances(sorted_vector_of_distances_)
    {
		this->set_up_numbers_of_functions_for_vectorization_and_projections_to_reals();
	}
    
    void set_up_numbers_of_functions_for_vectorization_and_projections_to_reals()
	{		
		//warning, this function can be only called after filling in the intervals vector.
		this->number_of_functions_for_vectorization = this->sorted_vector_of_distances.size();
		this->number_of_functions_for_projections_to_reals = this->sorted_vector_of_distances.size();
	}
};


template <typename F>
Vector_distances_in_diagram<F>::Vector_distances_in_diagram( const std::vector< std::pair< double,double > >& intervals_ , size_t where_to_cut_ ):where_to_cut(where_to_cut_)
{    
    std::vector< std::pair< double,double > > i( intervals_ );
    this->intervals = i;
    //this->compute_sorted_vector_of_distances_via_heap( where_to_cut );
    this->compute_sorted_vector_of_distances_via_vector_sorting(where_to_cut);
    this->set_up_numbers_of_functions_for_vectorization_and_projections_to_reals();
}

template <typename F>
Vector_distances_in_diagram<F>::Vector_distances_in_diagram( const char* filename , size_t where_to_cut , unsigned dimension  ):where_to_cut(where_to_cut)
{
    std::vector< std::pair< double , double > > intervals;
    if ( dimension == std::numeric_limits<unsigned>::max() )
    {
        intervals = read_persistence_intervals_in_one_dimension_from_file( filename );
    }
    else
    {
        intervals = read_persistence_intervals_in_one_dimension_from_file( filename , dimension );
    }
    this->intervals = intervals;
    this->compute_sorted_vector_of_distances_via_heap( where_to_cut );
    //this->compute_sorted_vector_of_distances_via_vector_sorting( where_to_cut );
    set_up_numbers_of_functions_for_vectorization_and_projections_to_reals();  
}

template < typename F>
void Vector_distances_in_diagram<F>::compute_sorted_vector_of_distances_via_heap( size_t where_to_cut )
{

    bool dbg = false;
    if ( dbg )
    {
        std::cerr << "Here are the intervals : \n";
        for ( size_t i = 0 ; i != this->intervals.size() ; ++i )
        {
            std::cerr << this->intervals[i].first << " , " << this->intervals[i].second <<std::endl;
        }
    }       
	where_to_cut = std::min(where_to_cut , (size_t)(0.5 * this->intervals.size() * ( this->intervals.size() - 1 ) + this->intervals.size()));

    std::vector< double > heap( where_to_cut , std::numeric_limits<int>::max() );
    std::make_heap (heap.begin(),heap.end());
    F f;

	//for every pair of points in the diagram, compute the minimum of their distance, and distance of those points from diagonal	
    for ( size_t i = 0 ; i < this->intervals.size() ; ++i )
    {
        for ( size_t j = i+1 ; j < this->intervals.size() ; ++j )
        {			
            double value = std::min(
                                f( this->intervals[i] , this->intervals[j] ),                               
                                std::min(
                                        f( this->intervals[i] , std::make_pair( 0.5*(this->intervals[i].first+this->intervals[i].second) , 0.5*(this->intervals[i].first+this->intervals[i].second) ) ),
                                        f( this->intervals[j] , std::make_pair( 0.5*(this->intervals[j].first+this->intervals[j].second) , 0.5*(this->intervals[j].first+this->intervals[j].second) ) )                                
                                        )    
                                );

							if ( dbg )
							{
								std::cerr << "Value : " << value <<std::endl;
								std::cerr << "heap.front() : " << heap.front() <<std::endl;
								getchar();
							}

            if ( -value < heap.front() )
            {
                if ( dbg ){std::cerr << "Replacing : " << heap.front() << " with : " << -value <<std::endl;getchar();}
                //remove the first element from the heap
                std::pop_heap (heap.begin(),heap.end());
                //heap.pop_back();
                //and put value there instead:
                //heap.push_back(-value);
                heap[ where_to_cut-1 ] = -value;
                std::push_heap (heap.begin(),heap.end());
            }
        }
    }
    
    //now add distances of all points from diagonal   
    for ( size_t i = 0 ; i < this->intervals.size() ; ++i )
    {
		double value = f( this->intervals[i] , std::make_pair( 0.5*(this->intervals[i].first+this->intervals[i].second) , 0.5*(this->intervals[i].first+this->intervals[i].second) ) );
		if ( -value < heap.front() )
            {
                //remove the first element from the heap
                std::pop_heap (heap.begin(),heap.end());
                //heap.pop_back();
                //and put value there instead:
                //heap.push_back(-value);
                heap[ where_to_cut-1 ] = -value;
                std::push_heap (heap.begin(),heap.end());
            }
	}
	 
    
    std::sort_heap (heap.begin(),heap.end());
    for ( size_t i = 0 ; i != heap.size() ; ++i )
    {
        if ( heap[i] == std::numeric_limits<int>::max() )
        {
            heap[i] = 0;
        }
        else
        {
            heap[i] *= -1;
        }
    }
    
    if ( dbg )
    {
		std::cerr << "This is the heap after all the operations :\n";
		for ( size_t i = 0 ; i != heap.size() ; ++i )
		{
			std::cout << heap[i] << " ";
		} 
		std::cout <<std::endl;
	} 

    this->sorted_vector_of_distances = heap;
}




template < typename F>
void Vector_distances_in_diagram<F>::compute_sorted_vector_of_distances_via_vector_sorting( size_t where_to_cut )
{		
	std::vector< double > distances;	
	distances.reserve( (size_t)(0.5 * this->intervals.size() * ( this->intervals.size() - 1 ) + this->intervals.size()) );
    F f;

	//for every pair of points in the diagram, compute the minimum of their distance, and distance of those points from diagonal	
    for ( size_t i = 0 ; i < this->intervals.size() ; ++i )
    {
		//add distance of i-th point in the diagram from the diagonal to the distances vector
		distances.push_back( f( this->intervals[i] , std::make_pair( 0.5*(this->intervals[i].first+this->intervals[i].second) , 0.5*(this->intervals[i].first+this->intervals[i].second) ) )  );
        for ( size_t j = i+1 ; j < this->intervals.size() ; ++j )
        {
            double value = std::min(
                                f( this->intervals[i] , this->intervals[j] ),                               
                                std::min(
                                        f( this->intervals[i] , std::make_pair( 0.5*(this->intervals[i].first+this->intervals[i].second) , 0.5*(this->intervals[i].first+this->intervals[i].second) ) ),
                                        f( this->intervals[j] , std::make_pair( 0.5*(this->intervals[j].first+this->intervals[j].second) , 0.5*(this->intervals[j].first+this->intervals[j].second) ) )                                
                                        )    
                                );
            distances.push_back( value );
            
        }
    }
    std::sort( distances.begin() , distances.end() , std::greater<double>() );    
    if ( distances.size() > where_to_cut )distances.resize( where_to_cut );       
    
    this->sorted_vector_of_distances =  distances;
}



//Implementations of functions for various concepts.  
template <typename F>
double Vector_distances_in_diagram<F>::project_to_R( int number_of_function )const
{
	if ( (size_t)number_of_function > this->number_of_functions_for_projections_to_reals )throw "Wrong index of a function in a method Vector_distances_in_diagram<F>::project_to_R";
	if ( number_of_function < 0 )throw "Wrong index of a function in a method Vector_distances_in_diagram<F>::project_to_R";
	
	double result = 0;
	for ( size_t i = 0 ; i != (size_t)number_of_function ; ++i )
	{
		result += sorted_vector_of_distances[i];
	}
	return result;
}

template <typename F>
void Vector_distances_in_diagram<F>::compute_average( const std::vector< Vector_distances_in_diagram* >& to_average )
{
	
	if ( to_average.size() == 0 )
	{
		(*this) = Vector_distances_in_diagram<F>();
		return;
	}
	
	size_t maximal_length_of_vector = 0;
	for ( size_t i = 0 ; i != to_average.size() ; ++i )
	{		
		if ( to_average[i]->sorted_vector_of_distances.size() > maximal_length_of_vector )
		{
			maximal_length_of_vector = to_average[i]->sorted_vector_of_distances.size();
		}
	}	
	
	std::vector< double > av(  maximal_length_of_vector , 0 );
	for ( size_t i = 0 ; i != to_average.size() ; ++i )
	{		
		for ( size_t j = 0  ; j != to_average[i]->sorted_vector_of_distances.size() ; ++j )
		{
			av[j] += to_average[i]->sorted_vector_of_distances[j];
		}
	}
	
	for ( size_t i = 0 ; i != maximal_length_of_vector ; ++i )
	{
		av[i] /= (double)to_average.size();
	}
	this->sorted_vector_of_distances = av;
	this->where_to_cut = av.size();
}

template <typename F>
double Vector_distances_in_diagram<F>::distance( const Vector_distances_in_diagram& second_ , double power )const
{
	bool dbg = false;
	
	if ( dbg )
	{
		std::cerr << "Entering double Vector_distances_in_diagram<F>::distance( const Abs_Topological_data_with_distances* second , double power ) procedure \n";
		std::cerr << "Power : " << power << std::endl;
		std::cerr << "This : " << *this << std::endl;
		std::cerr << "second : " << second_ << std::endl;
	}		
	
	
	double result = 0;
	for ( size_t i = 0 ; i != std::min(this->sorted_vector_of_distances.size(), second_.sorted_vector_of_distances.size())  ; ++i )
	{
		if ( power == 1 )
		{
			if ( dbg )
			{
				std::cerr << "|" << this->sorted_vector_of_distances[i] << " -  " << second_.sorted_vector_of_distances[i] << " |  : " << fabs( this->sorted_vector_of_distances[i] - second_.sorted_vector_of_distances[i] ) <<std::endl;
			}
			result += fabs( this->sorted_vector_of_distances[i] - second_.sorted_vector_of_distances[i] );
		}
		else
		{
			if ( power < std::numeric_limits<double>::max() )
			{
				result += std::pow( fabs( this->sorted_vector_of_distances[i] - second_.sorted_vector_of_distances[i] ) , power );
			}
			else
			{
				// max norm
				if ( result < fabs( this->sorted_vector_of_distances[i] - second_.sorted_vector_of_distances[i] ) )result = fabs( this->sorted_vector_of_distances[i] - second_.sorted_vector_of_distances[i] );
			}
			if ( dbg )
			{
				std::cerr << "| " << this->sorted_vector_of_distances[i] << " - " << second_.sorted_vector_of_distances[i]  << " : " << fabs( this->sorted_vector_of_distances[i] - second_.sorted_vector_of_distances[i] ) <<std::endl;
			}
		}
	}
	if ( this->sorted_vector_of_distances.size() != second_.sorted_vector_of_distances.size() )
	{
		if ( this->sorted_vector_of_distances.size() > second_.sorted_vector_of_distances.size() )
		{
			for ( size_t i = second_.sorted_vector_of_distances.size() ; i != this->sorted_vector_of_distances.size() ; ++i )			
			{
				result += fabs( this->sorted_vector_of_distances[i] );
			}
		}
		else
		{
			//this->sorted_vector_of_distances.size() < second_.sorted_vector_of_distances.size()
			for ( size_t i = this->sorted_vector_of_distances.size() ; i != second_.sorted_vector_of_distances.size() ; ++i )			
			{
				result += fabs( second_.sorted_vector_of_distances[i] );
			}			
		}
	}
	
	
	if ( power != 1 )
	{
		result = std::pow( result , (1.0/power) );
	}
	return result;
}
    
template < typename F>
std::vector<double> Vector_distances_in_diagram<F>::vectorize( int number_of_function )const
{
	if ( (size_t)number_of_function > this->number_of_functions_for_vectorization )throw "Wrong index of a function in a method Vector_distances_in_diagram<F>::vectorize";
	if ( number_of_function < 0 )throw "Wrong index of a function in a method Vector_distances_in_diagram<F>::vectorize";
	
    std::vector< double > result( std::min( (size_t)number_of_function , this->sorted_vector_of_distances.size() ) );
    for ( size_t i = 0 ; i != std::min( (size_t)number_of_function , this->sorted_vector_of_distances.size() ) ; ++i )
    {
        result[i] = this->sorted_vector_of_distances[i];
    }
    return result;
}


template < typename F>
void Vector_distances_in_diagram<F>::write_to_file( const char* filename )const
{
	std::ofstream out;
	out.open( filename );
	
	for ( size_t i = 0 ; i != this->sorted_vector_of_distances.size() ; ++i )
	{
		out << this->sorted_vector_of_distances[i] << " ";
	}
	
	out.close();
}

template < typename F>
void Vector_distances_in_diagram<F>::load_from_file( const char* filename )
{
	std::ifstream in;
	in.open( filename );
	//check if the file exist.
	if ( !in.good() )
	{
		std::cerr << "The file : " << filename << " do not exist. The program will now terminate \n";
		throw "The file from which you are trying to read the persistence landscape do not exist. The program will now terminate \n";
	}	

	double number;
	while ( true )
	{		
		in >> number;
		if ( in.eof() )break;
		this->sorted_vector_of_distances.push_back(number);
	}	
	in.close();
}

template < typename F>
double Vector_distances_in_diagram<F>::compute_scalar_product( const Vector_distances_in_diagram& second_vector )const
{
	double result = 0;
	for ( size_t i = 0 ; i != std::min(this->sorted_vector_of_distances.size(),second_vector.sorted_vector_of_distances.size()) ; ++i )
	{
		result += this->sorted_vector_of_distances[i] * second_vector.sorted_vector_of_distances[i];
	}
	return result;
}

}//namespace Persistence_representations
}//namespace Gudhi


#endif  // PERSISTENCE_VECTORS_H_

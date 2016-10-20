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


#ifndef Vector_distances_in_diagram_H
#define Vector_distances_in_diagram_H

#include <fstream>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <limits>


#include <gudhi/abstract_classes/Abs_Vectorized_topological_data.h>
#include <gudhi/abstract_classes/Abs_Topological_data_with_averages.h>
#include <gudhi/abstract_classes/Abs_Topological_data_with_distances.h>
#include <gudhi/abstract_classes/Abs_Real_valued_topological_data.h>
#include <gudhi/abstract_classes/Abs_Topological_data_with_scalar_product.h>
#include <gudhi/concretizations/read_persitence_from_file.h>
#include <gudhi/common_gudhi_stat.h>
using namespace std;

namespace Gudhi 
{
namespace Gudhi_stat 
{


template <typename T>
struct euclidean_distance
{
    double operator() ( const std::pair< T,T >& f , const std::pair<T,T>& s )
    {
        return  sqrt( (f.first-s.first)*(f.first-s.first) + (f.second-s.second)*(f.second-s.second) );
    }
};

template <typename T>
struct maximum_distance
{
    double operator() ( const std::pair< T,T >& f , const std::pair<T,T>& s )
    {
        return  std::min( fabs( f.first - s.first ) , fabs( f.second - s.second ) );
    }
};


/**
* This is an implementation of idea presented in the paper by Steve, Matthew and Max. The parameter of the class is the class that computes distance used to construct the vectors. The typical function is
* either Eucludean of maximum (Manhattan) distance. 
**/

template <typename F>
class Vector_distances_in_diagram : 
									public Abs_Vectorized_topological_data, 
									public Abs_Topological_data_with_distances, 
									public Abs_Real_valued_topological_data, 
									public Abs_Topological_data_with_averages, 
									public Abs_Topological_data_with_scalar_product
{
public:
	/**
	* The default constructor.  
	**/
    Vector_distances_in_diagram(){};
        
    /**
	* The constructor that takes as an input a multiset of persistence intervals (given as vector of birth-death pairs). The second parameter is the desiered length of the output vectors. 
	**/
    Vector_distances_in_diagram( const std::vector< std::pair< double , double > >& intervals , size_t where_to_cut );
    
    /**
	* The constructor taking as an input a file with birth-death pairs. The second parameter is the desiered length of the output vectors. 
	**/
    Vector_distances_in_diagram( const char* filename , size_t where_to_cut );

	/**
	 *  Assignement operator.
	**/ 
    Vector_distances_in_diagram<F>& operator =( const Vector_distances_in_diagram& org );


	/**
	 *  Copy constructor. 
	**/
    Vector_distances_in_diagram( const Vector_distances_in_diagram& org );
	
	/**
	 *  Writing to a stream. 
	**/
	template <typename K>
    friend ostream& operator << ( ostream& out , const Vector_distances_in_diagram<K>& d )
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
    inline double vector_in_position( size_t position )
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
    void write_to_file( const char* filename );
    
     /**
	 * Write a vector to a file.
	**/
    void print_to_file( const char* filename )
    {
		this->write_to_file( filename );
	}  
    
    /**
	 * Loading a vector to a file.
	**/
    void load_from_file( const char* filename );
    
    /**
     * Comparision operators:
    **/ 
    bool operator == ( const Vector_distances_in_diagram& second )
    {
		if ( this->sorted_vector_of_distances.size() != second.sorted_vector_of_distances.size() )return false;
		for ( size_t i = 0 ; i != this->sorted_vector_of_distances.size() ; ++i )
		{
			if ( !almost_equal(this->sorted_vector_of_distances[i] , second.sorted_vector_of_distances[i]) )return false;
		}
		return true;
	}
	
	bool operator != ( const Vector_distances_in_diagram& second )
    {
		return !( *this == second );
	}
    
    //concretization of abstract methods:
     /**
     * Compute projection to real numbers of persistence vector. 
    **/
    double project_to_R( int number_of_function );
    
    /**
     * Compute a vectorization of a persistent vectors.
    **/
    std::vector<double> vectorize( int number_of_function );
    
    /**
     * Compute a average of two persistent vectors.
    **/
    void compute_average( const std::vector< Abs_Topological_data_with_averages* >& to_average );
    
    /**
     * Compute a distance of two persistent vectors.
    **/
    double distance( const Abs_Topological_data_with_distances* second , double power = 1);
    
    /**
     * Compute a scalar product of two persistent vectors.
    **/ 
    double compute_scalar_product( const Abs_Topological_data_with_scalar_product* second );
    
    /**
    * For visualization use output from vectorize and build histograms. 
    **/
    std::vector< double > output_for_visualization()
    {
		return this->sorted_vector_of_distances;
	}
	
	/**
	 * Create a gnuplot script to vizualize the data structure. 
	 **/ 
	void plot( const char* filename )
	{
		std::stringstream gnuplot_script;
		gnuplot_script << filename << "_Gnuplot_script";
		ofstream out;
		out.open( gnuplot_script.str().c_str() );
		out << "set style data histogram" << std::endl;
		out << "set style histogram cluster gap 1" << std::endl;
		out << "set style fill solid border -1" << std::endl;
		out << "plot '-' notitle" << std::endl;    
		for ( size_t i = 0 ; i != this->sorted_vector_of_distances.size() ; ++i )
		{
			out << this->sorted_vector_of_distances[i]  << std::endl;
		}		
		out << endl;
		out.close();
		std::cout << "To vizualize, open gnuplot and type: load \'" << gnuplot_script.str().c_str() << "\'" <<  std::endl;			
	}
    

private:
    std::vector< std::pair< double , double > > intervals;
    std::vector< double > sorted_vector_of_distances;

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
Vector_distances_in_diagram<F>::Vector_distances_in_diagram( const std::vector< std::pair< double,double > >& intervals_ , size_t where_to_cut ):Abs_Vectorized_topological_data(where_to_cut)
{    
    std::vector< std::pair< double,double > > i( intervals_ );
    this->intervals = i;
    //this->compute_sorted_vector_of_distances_via_heap( where_to_cut );
    this->compute_sorted_vector_of_distances_via_vector_sorting(where_to_cut);
    this->set_up_numbers_of_functions_for_vectorization_and_projections_to_reals();
}

template <typename F>
Vector_distances_in_diagram<F>::Vector_distances_in_diagram( const Vector_distances_in_diagram<F>& org )
{
    std::vector< std::pair< double,double > > inter( org.intervals );
    this->intervals = inter;
    std::vector< double > sorted_vector_of_distances( org.sorted_vector_of_distances );
    this->sorted_vector_of_distances = sorted_vector_of_distances;
    set_up_numbers_of_functions_for_vectorization_and_projections_to_reals();
}


template <typename F>
Vector_distances_in_diagram<F>& Vector_distances_in_diagram<F>::operator =( const Vector_distances_in_diagram& org )
{
    std::vector< std::pair< double , double > > inter( org.intervals );
    this->intervals = inter;
    std::vector< double > sorted_vector_of_distances( org.sorted_vector_of_distances );
    this->sorted_vector_of_distances = sorted_vector_of_distances;
    return *this;
}


template <typename F>
Vector_distances_in_diagram<F>::Vector_distances_in_diagram( const char* filename , size_t where_to_cut  ):Abs_Vectorized_topological_data(where_to_cut)
{   	
    //standard file with barcode
    std::vector< std::pair< double , double > > intervals = read_standard_file( filename );    
    //gudhi file with barcode
    //std::vector< std::pair< double , double > > intervals = read_gudhi_file( filename , dimension );    
     
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
        cerr << "Here are the intervals : \n";
        for ( size_t i = 0 ; i != this->intervals.size() ; ++i )
        {
            cerr << this->intervals[i].first << " , " << this->intervals[i].second << endl;
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
								cerr << "Value : " << value << endl;
								cerr << "heap.front() : " << heap.front() << endl;
								getchar();
							}

            if ( -value < heap.front() )
            {
                if ( dbg ){cerr << "Replacing : " << heap.front() << " with : " << -value << endl;getchar();}
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
                //cerr << "Replacing : " << heap.front() << " with : " << -value << endl;getchar();
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
		std::cout << endl;
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



//concretization of abstract methods:
template <typename F>
double Vector_distances_in_diagram<F>::project_to_R( int number_of_function )
{
	if ( number_of_function > this->number_of_functions_for_projections_to_reals )throw "Wrong index of a function in a method Vector_distances_in_diagram<F>::project_to_R";
	if ( number_of_function < 0 )throw "Wrong index of a function in a method Vector_distances_in_diagram<F>::project_to_R";
	
	double result = 0;
	for ( size_t i = 0 ; i != number_of_function ; ++i )
	{
		result += sorted_vector_of_distances[i];
	}
	return result;
}

template <typename F>
void Vector_distances_in_diagram<F>::compute_average( const std::vector< Abs_Topological_data_with_averages* >& to_average )
{
	
	if ( to_average.size() == 0 )
	{
		(*this) = Vector_distances_in_diagram<F>();
		return;
	}
	
	size_t maximal_length_of_vector = 0;
	for ( size_t i = 0 ; i != to_average.size() ; ++i )
	{
		Vector_distances_in_diagram<F>* current = (Vector_distances_in_diagram<F>*)to_average[i];
		if ( current->sorted_vector_of_distances.size() > maximal_length_of_vector )
		{
			maximal_length_of_vector = current->sorted_vector_of_distances.size();
		}
	}
	
	//cerr << "maximal_length_of_vector  : " << maximal_length_of_vector << endl;
	
	std::vector< double > av(  maximal_length_of_vector , 0 );
	for ( size_t i = 0 ; i != to_average.size() ; ++i )
	{
		Vector_distances_in_diagram<F>* current = (Vector_distances_in_diagram<F>*)to_average[i];
		for ( size_t j = 0  ; j != current->sorted_vector_of_distances.size() ; ++j )
		{
			av[j] += current->sorted_vector_of_distances[j];
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
double Vector_distances_in_diagram<F>::distance( const Abs_Topological_data_with_distances* second , double power )
{
	bool dbg = false;
	
	if ( dbg )
	{
		std::cerr << "Entering double Vector_distances_in_diagram<F>::distance( const Abs_Topological_data_with_distances* second , double power ) procedure \n";
		std::cerr << "Power : " << power << std::endl;
		std::cerr << "This : " << *this << std::endl;
		std::cerr << "second : " << *((Vector_distances_in_diagram<F>*)second) << std::endl;
	}
	
	Vector_distances_in_diagram<F>* second_ = ((Vector_distances_in_diagram<F>*)second);
	
	
	double result = 0;
	for ( size_t i = 0 ; i != std::min(this->sorted_vector_of_distances.size(), second_->sorted_vector_of_distances.size())  ; ++i )
	{
		if ( power == 1 )
		{
			if ( dbg )
			{
				cerr << "|" << this->sorted_vector_of_distances[i] << " -  " << second_->sorted_vector_of_distances[i] << " |  : " << fabs( this->sorted_vector_of_distances[i] - second_->sorted_vector_of_distances[i] ) << endl;
			}
			result += fabs( this->sorted_vector_of_distances[i] - second_->sorted_vector_of_distances[i] );
		}
		else
		{
			result += std::pow( fabs( this->sorted_vector_of_distances[i] - second_->sorted_vector_of_distances[i] ) , power );
			if ( dbg )
			{
				cerr << "| " << this->sorted_vector_of_distances[i] << " - " << second_->sorted_vector_of_distances[i]  << " : " << fabs( this->sorted_vector_of_distances[i] - second_->sorted_vector_of_distances[i] ) << endl;
			}
		}
	}
	if ( this->sorted_vector_of_distances.size() != second_->sorted_vector_of_distances.size() )
	{
		if ( this->sorted_vector_of_distances.size() > second_->sorted_vector_of_distances.size() )
		{
			for ( size_t i = second_->sorted_vector_of_distances.size() ; i != this->sorted_vector_of_distances.size() ; ++i )			
			{
				result += fabs( this->sorted_vector_of_distances[i] );
			}
		}
		else
		{
			//this->sorted_vector_of_distances.size() < second_->sorted_vector_of_distances.size()
			for ( size_t i = this->sorted_vector_of_distances.size() ; i != second_->sorted_vector_of_distances.size() ; ++i )			
			{
				result += fabs( second_->sorted_vector_of_distances[i] );
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
std::vector<double> Vector_distances_in_diagram<F>::vectorize( int number_of_function )
{
	if ( number_of_function > this->number_of_functions_for_vectorization )throw "Wrong index of a function in a method Vector_distances_in_diagram<F>::vectorize";
	if ( number_of_function < 0 )throw "Wrong index of a function in a method Vector_distances_in_diagram<F>::vectorize";
	
    std::vector< double > result( std::min( (size_t)number_of_function , this->sorted_vector_of_distances.size() ) );
    for ( size_t i = 0 ; i != std::min( (size_t)number_of_function , this->sorted_vector_of_distances.size() ) ; ++i )
    {
        result[i] = this->sorted_vector_of_distances[i];
    }
    return result;
}


template < typename F>
void Vector_distances_in_diagram<F>::write_to_file( const char* filename )
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
	//check if the file exist.
	if ( !( access( filename, F_OK ) != -1 ) )
	{
		cerr << "The file : " << filename << " do not exist. The program will now terminate \n";
		throw "The file from which you are trying to read the persistence landscape do not exist. The program will now terminate \n";
	}	
	std::ifstream in;
	in.open( filename );
	
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
double Vector_distances_in_diagram<F>::compute_scalar_product( const Abs_Topological_data_with_scalar_product* second )
{
	Vector_distances_in_diagram<F>* second_vector = (Vector_distances_in_diagram<F>*)second;
		
	double result = 0;
	for ( size_t i = 0 ; i != std::min(this->sorted_vector_of_distances.size(),second_vector->sorted_vector_of_distances.size()) ; ++i )
	{
		result += this->sorted_vector_of_distances[i] * second_vector->sorted_vector_of_distances[i];
	}
	return result;
}



}//namespace Gudhi_stat 
}//namespace Gudhi


#endif // Vector_distances_in_diagram_H

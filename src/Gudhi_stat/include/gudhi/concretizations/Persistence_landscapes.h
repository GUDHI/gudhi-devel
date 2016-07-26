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


#ifndef Persistence_landscapes_H
#define Persistence_landscapes_H

//standard include
#include <cmath>
#include <iostream>
#include <vector>
#include <limits>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <unistd.h>


//gudhi include
#include <gudhi/abstract_classes/Abs_Vectorized_topological_data.h>
#include <gudhi/abstract_classes/Abs_Topological_data_with_averages.h>
#include <gudhi/abstract_classes/Abs_Topological_data_with_distances.h>
#include <gudhi/abstract_classes/Abs_Real_valued_topological_data.h>
#include <gudhi/abstract_classes/Abs_Topological_data_with_scalar_product.h>
#include <gudhi/concretizations/read_persitence_from_file.h>
using namespace std;




namespace Gudhi
{
namespace Gudhi_stat
{




//double epsi = std::numeric_limits<double>::epsilon();
double epsi = 0.000005;


/**
 *  A procedure used to compare doubles. Typically gien two doubles A and B, comparing A == B is not good idea. In this case, we use the procedure almostEqual with the epsi defined at
 *  the top of the file. Setting up the epsi give the user a tolerance on what should be consider equal.
**/
inline bool almost_equal( double a , double b )
{
    if ( fabs(a-b) < epsi )
        return true;
    return false;
}


/**
 * Extra functions needed in construction of barcodes.
**/
double minus_length( std::pair<double,double> a )
{
    return a.first-a.second;
}
double birth_plus_deaths( std::pair<double,double> a )
{
    return a.first+a.second;
}



/**
 * Given two points in R^2, the procedure compute the parameters A and B of the line y = Ax + B that crosses those two points.
**/
std::pair<double,double> compute_parameters_of_a_line( std::pair<double,double> p1 , std::pair<double,double> p2 )
{
    double a = (p2.second-p1.second)/( p2.first - p1.first );
    double b = p1.second - a*p1.first;
    return std::make_pair(a,b);
}

/**
 * This procedure given two points which lies on the opposide sides of x axis, compute x for which the line connecting those two points crosses x axis.
**/
double find_zero_of_a_line_segment_between_those_two_points ( std::pair<double,double> p1, std::pair<double,double> p2 )
{
    if ( p1.first == p2.first )return p1.first;
    if ( p1.second*p2.second > 0 )
    {
        std::ostringstream errMessage;
        errMessage <<"In function find_zero_of_a_line_segment_between_those_two_points the agguments are: (" << p1.first << "," << p1.second << ") and (" << p2.first << "," << p2.second << "). There is no zero in line between those two points. Program terminated.";
        std::string errMessageStr = errMessage.str();
        const char* err = errMessageStr.c_str();
        throw(err);
    }
    //we assume here, that x \in [ p1.first, p2.first ] and p1 and p2 are points between which we will put the line segment
    double a = (p2.second - p1.second)/(p2.first - p1.first);
    double b = p1.second - a*p1.first;
    //cerr << "Line crossing points : (" << p1.first << "," << p1.second << ") oraz (" << p2.first << "," << p2.second << ") : \n";
    //cerr << "a : " << a << " , b : " << b << " , x : " << x << endl;
    return -b/a;
}




/**
 * Lexicographical ordering of points	.
**/
bool compare_points_sorting( std::pair<double,double> f, std::pair<double,double> s )
{
    if ( f.first < s.first )
    {
        return true;
    }
    else
    {//f.first >= s.first
        if ( f.first > s.first )
        {
            return false;
        }
        else
        {//f.first == s.first
            if ( f.second > s.second )
            {
                return true;
            }
            else
            {
                return false;
            }
        }
    }
}



/**
 * This procedure takes two points in R^2 and a double value x. It conputes the line pasing through those two points and return the value of that linear function at x.
**/
double function_value ( std::pair<double,double> p1, std::pair<double,double> p2 , double x )
{
    //we assume here, that x \in [ p1.first, p2.first ] and p1 and p2 are points between which we will put the line segment
    double a = (p2.second - p1.second)/(p2.first - p1.first);
    double b = p1.second - a*p1.first;
    return (a*x+b);
}





/**
 * A clas implementing persistence landascpes data structures. For theroretical desciritpion, please consult a paper ''Statistical topological data analysis using persistence landscapes'' by Peter Bubenik.
 * For details of algorithms, please consult ''A persistence landscapes toolbox for topological statistics'' by Peter Bubenik and Pawel Dlotko.
 * Persistence landscapes allow vertorization, computations of distances, computations of projections to Real, computations of averages and scalar products. Therefore they implement suitable interfaces.
**/
class Persistence_landscape :
								public Abs_Vectorized_topological_data ,
								public Abs_Topological_data_with_distances,
								public Abs_Real_valued_topological_data,
								public Abs_Topological_data_with_averages,
								public Abs_Topological_data_with_scalar_product
{
public:
	/**
	 * Default constructor.
	**/
    Persistence_landscape()
    {
		this->set_up_numbers_of_functions_for_vectorization_and_projections_to_reals();
	}

    /**
	 * Constructor that takes as an input a vector of birth-death pairs.
	**/
    Persistence_landscape( const std::vector< std::pair< double , double > >& p );

    /**
	 * Assignement operator.
	**/
    Persistence_landscape& operator=( const Persistence_landscape& org );

    /**
	 * Copy constructor.
	**/
    Persistence_landscape(const Persistence_landscape&);

    /**
	 * Constructor that reads persistence intervals from file and creates persistence landscape. The format of the input file is the following: in each line we put birth-death pair. Last line is assumed
	 * to be empty. Even if the points within a line are not ordered, they will be ordered while the input is read.
	**/
    Persistence_landscape(const char* filename , size_t dimension = 0);



    /**
     * This procedure loads a landscape from file. It erase all the data that was previously stored in this landscape.
    **/
    void load_landscape_from_file( const char* filename );


    /**
     * The procedure stores a landscape to a file. The file can be later used by a procedure load_landscape_from_file.
    **/
    void print_to_file( const char* filename )const;



	/**
	 * This function compute integral of the landscape (defined formally as sum of integrals on R of all landscape functions)
	**/
    double compute_integral_of_landscape()const;


    /**
	 * This function compute integral of the 'level'-level of a landscape.
	**/
    double compute_integral_of_landscape( size_t level )const;


    /**
	 * This function compute integral of the landscape p-th power of a landscape (defined formally as sum of integrals on R of p-th powers of all landscape functions)
	**/
	double compute_integral_of_landscape( double p )const;//this function compute integral of p-th power of landscape.


    /**
     * A function that computes the value of a landscape at a given point. The parameters of the function are: unsigned level and double x.
     * The procedure will compute the value of the level-landscape at the point x.
    **/
    double compute_value_at_a_given_point( unsigned level , double x )const;

    /**
     * Writing landscape into a stream. A i-th level landscape starts with a string "lambda_i". Then the discontinuity points of the landscapes follows.
     * Shall those points be joined with lines, we will obtain the i-th landscape function.
    **/
    friend std::ostream& operator<<(std::ostream& out, Persistence_landscape& land );





	/**
	 * A function that compute sum of two landscapes.
	**/
    friend Persistence_landscape add_two_landscapes ( const Persistence_landscape& land1 ,  const Persistence_landscape& land2 )
    {
        return operation_on_pair_of_landscapes< std::plus<double> >(land1,land2);
    }

    /**
	 * A function that compute difference of two landscapes.
	**/
    friend Persistence_landscape subtract_two_landscapes ( const Persistence_landscape& land1 ,  const Persistence_landscape& land2 )
    {
        return operation_on_pair_of_landscapes< std::minus<double> >(land1,land2);
    }

	/**
	 * An operator +, that compute sum of two landscapes.
	**/
    friend Persistence_landscape operator+( const Persistence_landscape& first , const Persistence_landscape& second )
    {
        return add_two_landscapes( first,second );
    }

	/**
	 * An operator -, that compute difference of two landscapes.
	**/
    friend Persistence_landscape operator-( const Persistence_landscape& first , const Persistence_landscape& second )
    {
        return subtract_two_landscapes( first,second );
    }

	/**
	 * An operator * that allows multipilication of a landscape by a real number.
	**/
    friend Persistence_landscape operator*( const Persistence_landscape& first , double con )
    {
        return first.multiply_lanscape_by_real_number_not_overwrite(con);
    }

	/**
	 * An operator * that allows multipilication of a landscape by a real number (order of parameters swapped).
	**/
    friend Persistence_landscape operator*( double con , const Persistence_landscape& first  )
    {
        return first.multiply_lanscape_by_real_number_not_overwrite(con);
    }

	/**
	 * Operator +=. The second parameter is persistnece landwscape.
	**/
    Persistence_landscape operator += ( const Persistence_landscape& rhs )
    {
        *this = *this + rhs;
        return *this;
    }

	/**
	 * Operator -=. The second parameter is persistnece landwscape.
	**/
    Persistence_landscape operator -= ( const Persistence_landscape& rhs )
    {
        *this = *this - rhs;
        return *this;
    }


	/**
	 * Operator *=. The second parameter is a real number by which the y values of all landscape functions are multiplied. The x-values remain unchanged. 
	**/
    Persistence_landscape operator *= ( double x )
    {
        *this = *this*x;
        return *this;
    }

	/**
	 * Operator /=. The second parameter is a real number.
	**/
    Persistence_landscape operator /= ( double x )
    {
        if ( x == 0 )throw( "In operator /=, division by 0. Program terminated." );
        *this = *this * (1/x);
        return *this;
    }

	/**
	 * An operator to compare two persistence landscapes.
	**/
    bool operator == ( const Persistence_landscape& rhs  )const;


    /**
	 * An operator to compare two persistence landscapes.
	**/
    bool operator != ( const Persistence_landscape& rhs  )const
    {
		return !((*this) == rhs);
	}


	/**
	 * Computations of maximum (y) value of landscape.
	**/
    double compute_maximum()const
    {
        double maxValue = 0;
        if ( this->land.size() )
        {
            maxValue = -std::numeric_limits<int>::max();
            for ( size_t i = 0 ; i != this->land[0].size() ; ++i )
            {
                if ( this->land[0][i].second > maxValue )maxValue = this->land[0][i].second;
            }
        }
        return maxValue;
    }

	/**
	 * Computations of a L^i norm of landscape, where i is the input parameter.
	**/
    double compute_norm_of_landscape( double i )
    {
        Persistence_landscape l;
        if ( i != -1 )
        {
            return compute_discance_of_landscapes(*this,l,i);
        }
        else
        {
            return compute_max_norm_discance_of_landscapes(*this,l);
        }
    }

	/**
 	 * An operator to compute the value of a landscape in the level 'level' at the argument 'x'.
	**/
    double operator()(unsigned level,double x)const{return this->compute_value_at_a_given_point(level,x);}

	/**
	 * Computations of L^{\infty} distance between two landscapes.
	**/
    friend double compute_max_norm_discance_of_landscapes( const Persistence_landscape& first, const Persistence_landscape& second );
    //friend double compute_max_norm_discance_of_landscapes( const Persistence_landscape& first, const Persistence_landscape& second , unsigned& nrOfLand , double&x , double& y1, double& y2 );


	/**
	 * Computations of L^{p} distance between two landscapes. p is the parameter of the procedure.
	**/
    friend double compute_discance_of_landscapes( const Persistence_landscape& first, const Persistence_landscape& second , int p );



	/**
	 * Function to compute absolute value of a PL function. The representation of persistence landscapes allow to store general PL-function. When computing distance betwen two landscapes, we compute difference between
	 * them. In this case, a general PL-function with negative value can appear as a result. Then in order to compute distance, we need to take its absolute value. This is the purpose of this procedure.
	**/
    Persistence_landscape abs();

	/**
	 * Computes the number of landscape functions.
	**/
    size_t size()const{return this->land.size(); }

	/**
	 *  Computate maximal value of lambda-level landscape.
	**/
    double find_max( unsigned lambda )const;

	/**
	 * Function to compute inner (scalar) product of two landscapes.
	**/
    friend double compute_inner_product( const Persistence_landscape& l1 , const Persistence_landscape& l2 );




	//concretization of abstract functions:

	/**
	 * The number of projections to R is defined to the number of nonzero landscape functions. I-th projection is an integral of i-th landscape function over whole R.
	**/
    double project_to_R( int number_of_function )
    {
		return this->compute_integral_of_landscape( (size_t)number_of_function );
	}

    std::vector<double> vectorize( int number_of_function )
    {
		//TODO, think of something smarter over here
		std::vector<double> v;
		if ( (size_t)number_of_function > this->land.size() )
		{
			return v;
		}
		v.reserve( this->land[number_of_function].size() );
		for ( size_t i = 0 ; i != this->land[number_of_function].size() ; ++i )
		{
			v.push_back( this->land[number_of_function][i].second );
		}
		return v;
	}
    void compute_average( std::vector< Abs_Topological_data_with_averages* > to_average )
    {
		bool dbg = false;

		std::vector< Persistence_landscape* > nextLevelMerge( to_average.size() );
        for ( size_t i = 0 ; i != to_average.size() ; ++i )
        {
            nextLevelMerge[i] = (Persistence_landscape*)to_average[i];
        }
        bool is_this_first_level = true;//in the loop, we will create dynamically a unmber of intermediate complexes. We have to clean that up, but we cannot erase the initial andscapes we have
										//to average. In this case, we simply check if the nextLevelMerge are the input landscapes or the ones created in that loop by usig this extra variable.

        while ( nextLevelMerge.size() != 1 )
        {
            if ( dbg ){cerr << "nextLevelMerge.size() : " << nextLevelMerge.size() << endl;}
            std::vector< Persistence_landscape* > nextNextLevelMerge;
            nextNextLevelMerge.reserve( to_average.size() );
            for ( size_t i = 0 ; i < nextLevelMerge.size() ; i=i+2 )
            {
                if ( dbg ){cerr << "i : " << i << endl;}
                Persistence_landscape* l = new Persistence_landscape;
                if ( i+1 != nextLevelMerge.size() )
                {
                    (*l) = (*nextLevelMerge[i])+(*nextLevelMerge[i+1]);
                }
                else
                {
                    (*l) = *nextLevelMerge[i];
                }
                nextNextLevelMerge.push_back( l );
            }
            if ( dbg ){cerr << "After this iteration \n";}

            if ( !is_this_first_level )
            {
				//deallocate the memory if the vector nextLevelMerge do not consist of the initial landscapes
				for ( size_t i = 0 ; i != nextLevelMerge.size() ; ++i )
				{
					delete nextLevelMerge[i];
				}
			}
			is_this_first_level = false;
            nextLevelMerge.swap(nextNextLevelMerge);
        }
        (*this) = (*nextLevelMerge[0]);
        (*this) *= 1/( (double)to_average.size() );
	}


    double distance( const Abs_Topological_data_with_distances* second , double power = 1 )
    {
		if ( power != -1 )
		{
			return compute_discance_of_landscapes( *this , *((Persistence_landscape*)second) , power );
		}
		else
		{
			return compute_max_norm_discance_of_landscapes( *this , *((Persistence_landscape*)second) );
		}
	}


	double compute_scalar_product( const Abs_Topological_data_with_scalar_product* second )
	{
		return compute_inner_product( (*this) , *((Persistence_landscape*)second) );
	}


	std::vector< std::vector< std::pair<double,double> > > output_for_visualization()
	{
		return this->land;
	}
	
	
	//a function used to create a gnuplot script for visualization of landscapes
	void plot( const char* filename ,int from = -1, int to = -1 ,  double xRangeBegin = -1 , double xRangeEnd = -1 , double yRangeBegin = -1 , double yRangeEnd = -1 );


private:
    std::vector< std::vector< std::pair<double,double> > > land;

    void construct_persistence_landscape_from_barcode( const std::vector< std::pair< double , double > > & p );
    Persistence_landscape multiply_lanscape_by_real_number_not_overwrite( double x )const;
    void multiply_lanscape_by_real_number_overwrite( double x );
    template < typename oper > friend Persistence_landscape operation_on_pair_of_landscapes ( const Persistence_landscape& land1 ,  const Persistence_landscape& land2 );
    friend double compute_maximal_distance_non_symmetric( const Persistence_landscape& pl1, const Persistence_landscape& pl2 );

    void set_up_numbers_of_functions_for_vectorization_and_projections_to_reals()
	{
		//warning, this function can be only called after filling in the intervals vector.
		this->number_of_functions_for_vectorization = this->land.size();
		this->number_of_functions_for_projections_to_reals = this->land.size();
	}
};


Persistence_landscape::Persistence_landscape(const Persistence_landscape& oryginal)
{
    //std::cerr << "Running copy constructor \n";
    std::vector< std::vector< std::pair<double,double> > > land( oryginal.land.size() );
    for ( size_t i = 0 ; i != oryginal.land.size() ; ++i )
    {
        land[i].insert( land[i].end() , oryginal.land[i].begin() , oryginal.land[i].end() );
    }
    this->land = land;
    this->set_up_numbers_of_functions_for_vectorization_and_projections_to_reals();
}





Persistence_landscape::Persistence_landscape(const char* filename , size_t dimension)
{
    bool dbg = false;

    if ( dbg )
    {
        std::cerr << "Using constructor : Persistence_landscape(char* filename)" << std::endl;
    }   
    //standard file with barcode
    //std::vector< std::pair< double , double > > barcode = read_standard_file( filename );    
    //gudhi file with barcode
    std::vector< std::pair< double , double > > barcode = read_gudhi_file( filename , dimension );        
	this->construct_persistence_landscape_from_barcode( barcode );
	this->set_up_numbers_of_functions_for_vectorization_and_projections_to_reals();
}


bool operatorEqualDbg = false;
bool Persistence_landscape::operator == ( const Persistence_landscape& rhs  )const
{
    if ( this->land.size() != rhs.land.size() )
    {
        if (operatorEqualDbg)std::cerr << "1\n";
        return false;
    }
    for ( size_t level = 0 ; level != this->land.size() ; ++level )
    {
        if ( this->land[level].size() != rhs.land[level].size() )
        {
            if (operatorEqualDbg)std::cerr << "this->land[level].size() : " << this->land[level].size() <<  "\n";
            if (operatorEqualDbg)std::cerr << "rhs.land[level].size() : " << rhs.land[level].size() <<  "\n";
            if (operatorEqualDbg)std::cerr << "2\n";
            return false;
        }
        for ( size_t i = 0 ; i != this->land[level].size() ; ++i )
        {
            if ( !( almost_equal(this->land[level][i].first , rhs.land[level][i].first) && almost_equal(this->land[level][i].second , rhs.land[level][i].second) ) )
            {
				//cerr<< this->land[level][i].first << " , " <<  rhs.land[level][i].first << " and " << this->land[level][i].second << " , " << rhs.land[level][i].second << endl;
                if (operatorEqualDbg)std::cerr << "this->land[level][i] : " << this->land[level][i].first << " " << this->land[level][i].second << "\n";
                if (operatorEqualDbg)std::cerr << "rhs.land[level][i] : " << rhs.land[level][i].first << " " << rhs.land[level][i].second	 << "\n";
                if (operatorEqualDbg)std::cerr << "3\n";
                return false;
            }
        }
    }
    return true;
}


Persistence_landscape& Persistence_landscape::operator=( const Persistence_landscape& oryginal )
{
    std::vector< std::vector< std::pair<double,double> > > land( oryginal.land.size() );
    for ( size_t i = 0 ; i != oryginal.land.size() ; ++i )
    {
        land[i].insert( land[i].end() , oryginal.land[i].begin() , oryginal.land[i].end() );
    }
    this->land = land;
    return *this;
}



Persistence_landscape::Persistence_landscape( const std::vector< std::pair< double , double > > & p )
{
    this->construct_persistence_landscape_from_barcode( p );
    this->set_up_numbers_of_functions_for_vectorization_and_projections_to_reals();
}


void Persistence_landscape::construct_persistence_landscape_from_barcode( const std::vector< std::pair< double , double > > & p )
{
	bool dbg = false;
    if ( dbg ){cerr << "Persistence_landscape::Persistence_landscape( const std::vector< std::pair< double , double > >& p )" << endl;}

	//this is a general algorithm to construct persistence landscapes.
	std::vector< std::pair<double,double> > bars;
	bars.insert( bars.begin() , p.begin() , p.end() );
	std::sort( bars.begin() , bars.end() , compare_points_sorting );

	if (dbg)
	{
		std::cerr << "Bars : \n";
		for ( size_t i = 0 ; i != bars.size() ; ++i )
		{
			std::cerr << bars[i].first << " " << bars[i].second << "\n";
		}
		getchar();
	}

	std::vector< std::pair<double,double> > characteristicPoints(p.size());
	for ( size_t i = 0 ; i != bars.size() ; ++i )
	{
		characteristicPoints[i] = std::make_pair((bars[i].first+bars[i].second)/2.0 , (bars[i].second - bars[i].first)/2.0);
	}
	std::vector< std::vector< std::pair<double,double> > > Persistence_landscape;
	while ( !characteristicPoints.empty() )
	{
		if(dbg)
		{
			for ( size_t i = 0 ; i != characteristicPoints.size() ; ++i )
			{
				std::cout << "("  << characteristicPoints[i].first << " " << characteristicPoints[i].second << ")\n";
			}
			std::cin.ignore();
		}

		std::vector< std::pair<double,double> > lambda_n;
		lambda_n.push_back( std::make_pair( -std::numeric_limits<int>::max() , 0 ) );
		lambda_n.push_back( std::make_pair(minus_length(characteristicPoints[0]),0) );
		lambda_n.push_back( characteristicPoints[0] );

		if (dbg)
		{
			std::cerr << "1 Adding to lambda_n : (" << -std::numeric_limits<int>::max() << " " << 0  << ") , (" << minus_length(characteristicPoints[0]) << " " << 0 << ") , (" << characteristicPoints[0].first << " " << characteristicPoints[0].second << ") \n";
		}

		size_t i = 1;
		std::vector< std::pair<double,double> >  newCharacteristicPoints;
		while ( i < characteristicPoints.size() )
		{
			size_t p = 1;
			if ( (minus_length(characteristicPoints[i]) >= minus_length(lambda_n[lambda_n.size()-1])) && (birth_plus_deaths(characteristicPoints[i]) > birth_plus_deaths(lambda_n[lambda_n.size()-1])) )
			{
				if ( minus_length(characteristicPoints[i]) < birth_plus_deaths(lambda_n[lambda_n.size()-1]) )
				{
					std::pair<double,double> point = std::make_pair( (minus_length(characteristicPoints[i])+birth_plus_deaths(lambda_n[lambda_n.size()-1]))/2 , (birth_plus_deaths(lambda_n[lambda_n.size()-1])-minus_length(characteristicPoints[i]))/2 );
					lambda_n.push_back( point );
					if (dbg)
					{
						std::cerr << "2 Adding to lambda_n : (" << point.first << " " << point.second << ")\n";
					}


					if ( dbg )
					{
						std::cerr << "characteristicPoints[i+p] : " << characteristicPoints[i+p].first << " " << characteristicPoints[i+p].second << "\n";
						std::cerr << "point : " << point.first << " " << point.second << "\n";
						getchar();
					}

					while ( (i+p < characteristicPoints.size() ) && ( almost_equal(minus_length(point),minus_length(characteristicPoints[i+p])) ) && ( birth_plus_deaths(point) <= birth_plus_deaths(characteristicPoints[i+p]) ) )
					{
						newCharacteristicPoints.push_back( characteristicPoints[i+p] );
						if (dbg)
						{
							std::cerr << "3.5 Adding to newCharacteristicPoints : (" << characteristicPoints[i+p].first << " " << characteristicPoints[i+p].second << ")\n";
							getchar();
						}
						++p;
					}


					newCharacteristicPoints.push_back( point );
					if (dbg)
					{
						std::cerr << "4 Adding to newCharacteristicPoints : (" << point.first << " " << point.second << ")\n";
					}


					while ( (i+p < characteristicPoints.size() ) && ( minus_length(point) <= minus_length(characteristicPoints[i+p]) ) && (birth_plus_deaths(point)>=birth_plus_deaths(characteristicPoints[i+p])) )
					{
						newCharacteristicPoints.push_back( characteristicPoints[i+p] );
						if (dbg)
						{
							std::cerr << "characteristicPoints[i+p] : " << characteristicPoints[i+p].first << " " << characteristicPoints[i+p].second << "\n";
							std::cerr << "point : " << point.first << " " << point.second << "\n";
							std::cerr << "characteristicPoints[i+p] birth and death : " << minus_length(characteristicPoints[i+p]) << " , " << birth_plus_deaths(characteristicPoints[i+p]) << "\n";
							std::cerr << "point birth and death : " << minus_length(point) << " , " << birth_plus_deaths(point) << "\n";

							std::cerr << "3 Adding to newCharacteristicPoints : (" << characteristicPoints[i+p].first << " " << characteristicPoints[i+p].second << ")\n";
							getchar();
						}
						++p;
					}

				}
				else
				{
					lambda_n.push_back( std::make_pair( birth_plus_deaths(lambda_n[lambda_n.size()-1]) , 0 ) );
					lambda_n.push_back( std::make_pair( minus_length(characteristicPoints[i]) , 0 ) );
					if (dbg)
					{
						std::cerr << "5 Adding to lambda_n : (" << birth_plus_deaths(lambda_n[lambda_n.size()-1]) << " " <<  0 << ")\n";
						std::cerr << "5 Adding to lambda_n : (" << minus_length(characteristicPoints[i]) << " " << 0  << ")\n";
					}
				}
				lambda_n.push_back( characteristicPoints[i] );
				if (dbg)
				{
					std::cerr << "6 Adding to lambda_n : (" << characteristicPoints[i].first << " " << characteristicPoints[i].second << ")\n";
				}
			}
			else
			{
				newCharacteristicPoints.push_back( characteristicPoints[i] );
				if (dbg)
				{
					std::cerr << "7 Adding to newCharacteristicPoints : (" << characteristicPoints[i].first << " " << characteristicPoints[i].second << ")\n";
				}
			}
			i = i+p;
		}
		lambda_n.push_back( std::make_pair(birth_plus_deaths(lambda_n[lambda_n.size()-1]),0) );
		lambda_n.push_back( std::make_pair( std::numeric_limits<int>::max() , 0 ) );

		characteristicPoints = newCharacteristicPoints;

		lambda_n.erase(std::unique(lambda_n.begin(), lambda_n.end()), lambda_n.end());
		this->land.push_back( lambda_n );
	}
}



//this function find maximum of lambda_n
double Persistence_landscape::find_max( unsigned lambda )const
{
    if ( this->land.size() < lambda )return 0;
    double maximum = -std::numeric_limits<int>::max();
    for ( size_t i = 0 ; i != this->land[lambda].size() ; ++i )
    {
        if ( this->land[lambda][i].second > maximum )maximum = this->land[lambda][i].second;
    }
    return maximum;
}


double Persistence_landscape::compute_integral_of_landscape()const
{
    double result = 0;
    for ( size_t i = 0 ; i != this->land.size() ; ++i )
    {
        for ( size_t nr = 2 ; nr != this->land[i].size()-1 ; ++nr )
        {
            //it suffices to compute every planar integral and then sum them ap for each lambda_n
            result += 0.5*( this->land[i][nr].first - this->land[i][nr-1].first )*(this->land[i][nr].second + this->land[i][nr-1].second);
        }
    }
    return result;
}

double Persistence_landscape::compute_integral_of_landscape( size_t  level )const
{
    double result = 0;
    if ( level >= this->land.size() )
    {
		//this landscape function is constantly equal 0, so is the intergral.
		return result;
	}
	//also negative landscapes are assumed to be zero.
	if ( level < 0 )return 0;

	for ( size_t nr = 2 ; nr != this->land[ level ].size()-1 ; ++nr )
	{
		//it suffices to compute every planar integral and then sum them ap for each lambda_n
		result += 0.5*( this->land[ level ][nr].first - this->land[ level ][nr-1].first )*(this->land[ level ][nr].second + this->land[ level ][nr-1].second);
	}

    return result;
}


bool compute_integral_of_landscapeDbg = false;
double Persistence_landscape::compute_integral_of_landscape( double p )const
{
    double result = 0;
    for ( size_t i = 0 ; i != this->land.size() ; ++i )
    {
        for ( size_t nr = 2 ; nr != this->land[i].size()-1 ; ++nr )
        {
            if (compute_integral_of_landscapeDbg)std::cout << "nr : " << nr << "\n";
            //In this interval, the landscape has a form f(x) = ax+b. We want to compute integral of (ax+b)^p = 1/a * (ax+b)^{p+1}/(p+1)
            std::pair<double,double> coef = compute_parameters_of_a_line( this->land[i][nr] , this->land[i][nr-1] );
            double a = coef.first;
            double b = coef.second;

            if (compute_integral_of_landscapeDbg)std::cout << "(" << this->land[i][nr].first << "," << this->land[i][nr].second << ") , " << this->land[i][nr-1].first << "," << this->land[i][nr].second << ")" << std::endl;
            if ( this->land[i][nr].first == this->land[i][nr-1].first )continue;
            if ( a != 0 )
            {
                result += 1/(a*(p+1)) * ( pow((a*this->land[i][nr].first+b),p+1) - pow((a*this->land[i][nr-1].first+b),p+1));
            }
            else
            {
                result += ( this->land[i][nr].first - this->land[i][nr-1].first )*( pow(this->land[i][nr].second,p) );
            }
            if ( compute_integral_of_landscapeDbg )
            {
                std::cout << "a : " <<a << " , b : " << b << std::endl;
                std::cout << "result : " << result << std::endl;
            }
        }
        //if (compute_integral_of_landscapeDbg) std::cin.ignore();
    }
    return result;
}


//this is O(log(n)) algorithm, where n is number of points in this->land.
double Persistence_landscape::compute_value_at_a_given_point( unsigned level , double x )const
{
	bool compute_value_at_a_given_pointDbg = false;
    //in such a case lambda_level = 0.
    if ( level > this->land.size() ) return 0;

    //we know that the points in this->land[level] are ordered according to x coordinate. Therefore, we can find the point by using bisection:
    unsigned coordBegin = 1;
    unsigned coordEnd = this->land[level].size()-2;

    if ( compute_value_at_a_given_pointDbg )
    {
        std::cerr << "Tutaj \n";
        std::cerr << "x : " << x << "\n";
        std::cerr << "this->land[level][coordBegin].first : " << this->land[level][coordBegin].first << "\n";
        std::cerr << "this->land[level][coordEnd].first : " << this->land[level][coordEnd].first << "\n";
    }

    //in this case x is outside the support of the landscape, therefore the value of the landscape is 0.
    if ( x <= this->land[level][coordBegin].first )return 0;
    if ( x >= this->land[level][coordEnd].first )return 0;

    if (compute_value_at_a_given_pointDbg)std::cerr << "Entering to the while loop \n";

    while ( coordBegin+1 != coordEnd )
    {
        if (compute_value_at_a_given_pointDbg)
        {
            std::cerr << "coordBegin : " << coordBegin << "\n";
            std::cerr << "coordEnd : " << coordEnd << "\n";
            std::cerr << "this->land[level][coordBegin].first : " << this->land[level][coordBegin].first << "\n";
            std::cerr << "this->land[level][coordEnd].first : " << this->land[level][coordEnd].first << "\n";
        }


        unsigned newCord = (unsigned)floor((coordEnd+coordBegin)/2.0);

        if (compute_value_at_a_given_pointDbg)
        {
            std::cerr << "newCord : " << newCord << "\n";
            std::cerr << "this->land[level][newCord].first : " << this->land[level][newCord].first << "\n";
            std::cin.ignore();
        }

        if ( this->land[level][newCord].first <= x )
        {
            coordBegin = newCord;
            if ( this->land[level][newCord].first == x )return this->land[level][newCord].second;
        }
        else
        {
            coordEnd = newCord;
        }
    }

    if (compute_value_at_a_given_pointDbg)
    {
        std::cout << "x : " << x << " is between : " << this->land[level][coordBegin].first << " a  " << this->land[level][coordEnd].first << "\n";
        std::cout << "the y coords are : " << this->land[level][coordBegin].second << " a  " << this->land[level][coordEnd].second << "\n";
        std::cerr << "coordBegin : " << coordBegin << "\n";
        std::cerr << "coordEnd : " << coordEnd << "\n";
        std::cin.ignore();
    }
    return function_value( this->land[level][coordBegin] , this->land[level][coordEnd] , x );
}

std::ostream& operator<<(std::ostream& out, Persistence_landscape& land )
{
    for ( size_t level = 0 ; level != land.land.size() ; ++level )
    {
        out << "Lambda_" << level << ":" << std::endl;
        for ( size_t i = 0 ; i != land.land[level].size() ; ++i )
        {
            if ( land.land[level][i].first == -std::numeric_limits<int>::max() )
            {
                out << "-inf";
            }
            else
            {
                if ( land.land[level][i].first == std::numeric_limits<int>::max() )
                {
                    out << "+inf";
                }
                else
                {
                    out << land.land[level][i].first;
                }
            }
            out << " , " << land.land[level][i].second << std::endl;
        }
    }
    return out;
}




void Persistence_landscape::multiply_lanscape_by_real_number_overwrite( double x )
{
    for ( size_t dim = 0 ; dim != this->land.size() ; ++dim )
    {
        for ( size_t i = 0 ; i != this->land[dim].size() ; ++i )
        {
             this->land[dim][i].second *= x;
        }
    }
}

bool AbsDbg = false;
Persistence_landscape Persistence_landscape::abs()
{
    Persistence_landscape result;
    for ( size_t level = 0 ; level != this->land.size() ; ++level )
    {
        if ( AbsDbg ){ std::cout << "level: " << level << std::endl; }
        std::vector< std::pair<double,double> > lambda_n;
        lambda_n.push_back( std::make_pair( -std::numeric_limits<int>::max() , 0 ) );
        for ( size_t i = 1 ; i != this->land[level].size() ; ++i )
        {
            if ( AbsDbg ){std::cout << "this->land[" << level << "][" << i << "] : " << this->land[level][i].first << " " << this->land[level][i].second << std::endl;}
            //if a line segment between this->land[level][i-1] and this->land[level][i] crosses the x-axis, then we have to add one landscape point t oresult
            if ( (this->land[level][i-1].second)*(this->land[level][i].second)  < 0 )
            {
                double zero = find_zero_of_a_line_segment_between_those_two_points( this->land[level][i-1] , this->land[level][i] );

                lambda_n.push_back( std::make_pair(zero , 0) );
                lambda_n.push_back( std::make_pair(this->land[level][i].first , fabs(this->land[level][i].second)) );
                if ( AbsDbg )
                {
                    std::cout << "Adding pair : (" << zero << ",0)" << std::endl;
                    std::cout << "In the same step adding pair : (" << this->land[level][i].first << "," << fabs(this->land[level][i].second) << ") " << std::endl;
                    std::cin.ignore();
                }
            }
            else
            {
                lambda_n.push_back( std::make_pair(this->land[level][i].first , fabs(this->land[level][i].second)) );
                if ( AbsDbg )
                {
                    std::cout << "Adding pair : (" << this->land[level][i].first << "," << fabs(this->land[level][i].second) << ") " << std::endl;
                    std::cin.ignore();
                }
            }
        }
        result.land.push_back( lambda_n );
    }
    return result;
}


Persistence_landscape Persistence_landscape::multiply_lanscape_by_real_number_not_overwrite( double x )const
{
    std::vector< std::vector< std::pair<double,double> > > result(this->land.size());
    for ( size_t dim = 0 ; dim != this->land.size() ; ++dim )
    {
        std::vector< std::pair<double,double> > lambda_dim( this->land[dim].size() );
        for ( size_t i = 0 ; i != this->land[dim].size() ; ++i )
        {
            lambda_dim[i] = std::make_pair( this->land[dim][i].first , x*this->land[dim][i].second );
        }
        result[dim] = lambda_dim;
    }
    Persistence_landscape res;
    //CHANGE
    //res.land = result;
    res.land.swap(result);
    return res;
}//multiply_lanscape_by_real_number_overwrite


void Persistence_landscape::print_to_file( const char* filename )const
{
    std::ofstream write;
    write.open(filename);
    for ( size_t dim = 0 ; dim != this->land.size() ; ++dim )
    {
        write << "#lambda_" << dim << std::endl;
        for ( size_t i = 1 ; i != this->land[dim].size()-1 ; ++i )
        {
            write << this->land[dim][i].first << "  " << this->land[dim][i].second << std::endl;
        }
    }
    write.close();
}

void Persistence_landscape::load_landscape_from_file( const char* filename )
{
	bool dbg = false;
	//removing the current content of the persistence landscape.
	this->land.clear();


	//this constructor reads persistence landscape form a file. This file have to be created by this software beforehead
    std::ifstream in;
    in.open( filename );
	if ( !( access( filename, F_OK ) != -1 ) )
	{
		cerr << "The file : " << filename << " do not exist. The program will now terminate \n";
		throw "The file from which you are trying to read the persistence landscape do not exist. The program will now terminate \n";
	}

    std::string line;
    std::vector< std::pair<double,double> > landscapeAtThisLevel;

    bool isThisAFirsLine = true;
    while ( !in.eof() )
    {
        getline(in,line);
        if ( !(line.length() == 0 || line[0] == '#') )
        {
            std::stringstream lineSS;
            lineSS << line;
            double beginn, endd;
            lineSS >> beginn;
            lineSS >> endd;
            landscapeAtThisLevel.push_back( std::make_pair( beginn , endd ) );
            if (dbg){std::cerr << "Reading a pont : " << beginn << " , " << endd << std::endl;}
        }
        else
        {
            if (dbg)
            {
                std::cout << "IGNORE LINE\n";
                getchar();
            }
            if ( !isThisAFirsLine )
            {
                landscapeAtThisLevel.push_back( std::make_pair( std::numeric_limits<int>::max() , 0 ) );
                this->land.push_back(landscapeAtThisLevel);
                std::vector< std::pair<double,double> > newLevelOdLandscape;
                landscapeAtThisLevel.swap(newLevelOdLandscape);
            }
            landscapeAtThisLevel.push_back( std::make_pair( -std::numeric_limits<int>::max() , 0 ) );
            isThisAFirsLine = false;
        }
	}
	if ( landscapeAtThisLevel.size() > 1 )
    {
        //seems that the last line of the file is not finished with the newline sign. We need to put what we have in landscapeAtThisLevel to the constructed landscape.
        landscapeAtThisLevel.push_back( std::make_pair( std::numeric_limits<int>::max() , 0 ) );
        this->land.push_back(landscapeAtThisLevel);
    }

    in.close();
}


template < typename T >
Persistence_landscape operation_on_pair_of_landscapes ( const Persistence_landscape& land1 ,  const Persistence_landscape& land2 )
{
	bool operation_on_pair_of_landscapesDBG = false;
    if ( operation_on_pair_of_landscapesDBG ){std::cout << "operation_on_pair_of_landscapes\n";std::cin.ignore();}
    Persistence_landscape result;
    std::vector< std::vector< std::pair<double,double> > > land( std::max( land1.land.size() , land2.land.size() ) );
    result.land = land;
    T oper;

    for ( size_t i = 0 ; i != std::min( land1.land.size() , land2.land.size() ) ; ++i )
    {
        std::vector< std::pair<double,double> > lambda_n;
        size_t p = 0;
        size_t q = 0;
        while ( (p+1 < land1.land[i].size()) && (q+1 < land2.land[i].size()) )
        {
            if ( operation_on_pair_of_landscapesDBG )
            {
                std::cerr << "p : " << p << "\n";
                std::cerr << "q : " << q << "\n";
                std::cout << "land1.land[i][p].first : " << land1.land[i][p].first << "\n";
                std::cout << "land2.land[i][q].first : " << land2.land[i][q].first << "\n";
            }

            if ( land1.land[i][p].first < land2.land[i][q].first )
            {
                if ( operation_on_pair_of_landscapesDBG )
                {
                    std::cout << "first \n";
                    std::cout << " function_value(land2.land[i][q-1],land2.land[i][q],land1.land[i][p].first) : "<<  function_value(land2.land[i][q-1],land2.land[i][q],land1.land[i][p].first) << "\n";
                    //std::cout << "oper( " << land1.land[i][p].second <<"," << function_value(land2.land[i][q-1],land2.land[i][q],land1.land[i][p].first) << " : " << oper( land1.land[i][p].second , function_value(land2.land[i][q-1],land2.land[i][q],land1.land[i][p].first) ) << "\n";
                }
                lambda_n.push_back( 
									std::make_pair( 
									              land1.land[i][p].first , 
									              oper( (double)land1.land[i][p].second , function_value(land2.land[i][q-1],land2.land[i][q],land1.land[i][p].first) ) 
												  ) 
								  );
                ++p;
                continue;
            }
            if ( land1.land[i][p].first > land2.land[i][q].first )
            {
                if ( operation_on_pair_of_landscapesDBG )
                {
                    std::cout << "Second \n";
                    std::cout << "function_value("<< land1.land[i][p-1].first << " " << land1.land[i][p-1].second <<" ,"<< land1.land[i][p].first << " " << land1.land[i][p].second <<", " << land2.land[i][q].first<<" ) : " << function_value( land1.land[i][p-1] , land1.land[i][p-1] ,land2.land[i][q].first ) << "\n";
                    std::cout << "oper( " << function_value( land1.land[i][p] , land1.land[i][p-1] ,land2.land[i][q].first ) <<"," << land2.land[i][q].second <<" : " << oper( land2.land[i][q].second , function_value( land1.land[i][p] , land1.land[i][p-1] ,land2.land[i][q].first ) ) << "\n";
                }
                lambda_n.push_back( std::make_pair( land2.land[i][q].first , oper( function_value( land1.land[i][p] , land1.land[i][p-1] ,land2.land[i][q].first ) , land2.land[i][q].second )  )  );
                ++q;
                continue;
            }
            if ( land1.land[i][p].first == land2.land[i][q].first )
            {
                if (operation_on_pair_of_landscapesDBG)std::cout << "Third \n";
                lambda_n.push_back( std::make_pair( land2.land[i][q].first , oper( land1.land[i][p].second , land2.land[i][q].second ) ) );
                ++p;++q;
            }
            if (operation_on_pair_of_landscapesDBG){std::cout << "Next iteration \n";getchar();}
        }
        while ( (p+1 < land1.land[i].size())&&(q+1 >= land2.land[i].size()) )
        {
            if (operation_on_pair_of_landscapesDBG)
            {
                std::cout << "New point : " << land1.land[i][p].first << "  oper(land1.land[i][p].second,0) : " <<  oper(land1.land[i][p].second,0) << std::endl;
            }
            lambda_n.push_back( std::make_pair(land1.land[i][p].first , oper(land1.land[i][p].second,0) ) );
            ++p;
        }
        while ( (p+1 >= land1.land[i].size())&&(q+1 < land2.land[i].size()) )
        {
            if (operation_on_pair_of_landscapesDBG)
            {
                std::cout << "New point : " << land2.land[i][q].first << " oper(0,land2.land[i][q].second) : " <<  oper(0,land2.land[i][q].second) << std::endl;
            }
            lambda_n.push_back( std::make_pair(land2.land[i][q].first , oper(0,land2.land[i][q].second) ) );
            ++q;
        }
        lambda_n.push_back( std::make_pair( std::numeric_limits<int>::max() , 0 ) );
        //CHANGE
        //result.land[i] = lambda_n;
        result.land[i].swap(lambda_n);
    }
    if ( land1.land.size() > std::min( land1.land.size() , land2.land.size() ) )
    {
        if (operation_on_pair_of_landscapesDBG){std::cout << "land1.land.size() > std::min( land1.land.size() , land2.land.size() )" << std::endl;}
        for ( size_t i = std::min( land1.land.size() , land2.land.size() ) ; i != std::max( land1.land.size() , land2.land.size() ) ; ++i )
        {
            std::vector< std::pair<double,double> > lambda_n( land1.land[i] );
            for ( size_t nr = 0 ; nr != land1.land[i].size() ; ++nr )
            {
                lambda_n[nr] = std::make_pair( land1.land[i][nr].first , oper( land1.land[i][nr].second , 0 ) );
            }
            //CHANGE
            //result.land[i] = lambda_n;
            result.land[i].swap(lambda_n);
        }
    }
    if ( land2.land.size() > std::min( land1.land.size() , land2.land.size() ) )
    {
        if (operation_on_pair_of_landscapesDBG){std::cout << "( land2.land.size() > std::min( land1.land.size() , land2.land.size() ) ) " << std::endl;}
        for ( size_t i = std::min( land1.land.size() , land2.land.size() ) ; i != std::max( land1.land.size() , land2.land.size() ) ; ++i )
        {
            std::vector< std::pair<double,double> > lambda_n( land2.land[i] );
            for ( size_t nr = 0 ; nr != land2.land[i].size() ; ++nr )
            {
                lambda_n[nr] = std::make_pair( land2.land[i][nr].first , oper( 0 , land2.land[i][nr].second ) );
            }
            //CHANGE
            //result.land[i] = lambda_n;
            result.land[i].swap(lambda_n);
        }
    }
    if ( operation_on_pair_of_landscapesDBG ){std::cout << "operation_on_pair_of_landscapes\n";std::cin.ignore();}
    return result;
}//operation_on_pair_of_landscapes



double compute_maximal_distance_non_symmetric( const Persistence_landscape& pl1, const Persistence_landscape& pl2 )
{
    bool dbg = false;
    if (dbg)std::cerr << " compute_maximal_distance_non_symmetric \n";
    //this distance is not symmetric. It compute ONLY distance between inflection points of pl1 and pl2.
    double maxDist = 0;
    size_t minimalNumberOfLevels = std::min( pl1.land.size() , pl2.land.size() );
    for ( size_t level = 0 ; level != minimalNumberOfLevels ; ++ level )
    {
        if (dbg)
        {
            std::cerr << "Level : " << level << std::endl;
            std::cerr << "PL1 : \n";
            for ( size_t i = 0 ; i  != pl1.land[level].size() ; ++i )
            {
                std::cerr << "(" <<pl1.land[level][i].first << "," << pl1.land[level][i].second << ") \n";
            }
            std::cerr << "PL2 : \n";
            for ( size_t i = 0 ; i  != pl2.land[level].size() ; ++i )
            {
                std::cerr << "(" <<pl2.land[level][i].first << "," << pl2.land[level][i].second << ") \n";
            }
            std::cin.ignore();
        }

        int p2Count = 0;
        for ( size_t i = 1 ; i != pl1.land[level].size()-1 ; ++i ) //w tym przypadku nie rozwarzam punktow w nieskocznosci
        {
            while ( true )
            {
                if (  (pl1.land[level][i].first>=pl2.land[level][p2Count].first) && (pl1.land[level][i].first<=pl2.land[level][p2Count+1].first)  )break;
                p2Count++;
            }
            double val = fabs( function_value( pl2.land[level][p2Count] , pl2.land[level][p2Count+1] , pl1.land[level][i].first ) - pl1.land[level][i].second);
            if ( maxDist <= val )maxDist = val;

            if (dbg)
            {
                std::cerr << pl1.land[level][i].first <<"in [" << pl2.land[level][p2Count].first << "," <<  pl2.land[level][p2Count+1].first <<"] \n";
                std::cerr << "pl1[level][i].second : " << pl1.land[level][i].second << std::endl;
                std::cerr << "function_value( pl2[level][p2Count] , pl2[level][p2Count+1] , pl1[level][i].first ) : " << function_value( pl2.land[level][p2Count] , pl2.land[level][p2Count+1] , pl1.land[level][i].first ) << std::endl;
                std::cerr << "val : "  << val << std::endl;
                std::cin.ignore();
            }
        }
    }

    if (dbg)std::cerr << "minimalNumberOfLevels : " << minimalNumberOfLevels << std::endl;

    if ( minimalNumberOfLevels < pl1.land.size() )
    {
        for ( size_t level = minimalNumberOfLevels ; level != pl1.land.size() ; ++ level )
        {
            for ( size_t i = 0 ; i != pl1.land[level].size() ; ++i )
            {
                if (dbg)std::cerr << "pl1[level][i].second  : " << pl1.land[level][i].second << std::endl;
                if ( maxDist < pl1.land[level][i].second )maxDist = pl1.land[level][i].second;
            }
        }
    }
    return maxDist;
}




double compute_discance_of_landscapes( const Persistence_landscape& first, const Persistence_landscape& second , int p )
{
    //This is what we want to compute: (\int_{- \infty}^{+\infty}| first-second |^p)^(1/p). We will do it one step at a time:

    //first-second :
    Persistence_landscape lan = first-second;

    //| first-second |:
    lan = lan.abs();
	if ( p != -1 )
	{
		//\int_{- \infty}^{+\infty}| first-second |^p
		double result;
		if ( p != 1 )
		{
			result = lan.compute_integral_of_landscape( (double)p );
		}
		else
		{
			result = lan.compute_integral_of_landscape();
		}

		//(\int_{- \infty}^{+\infty}| first-second |^p)^(1/p)
		return pow( result , 1/(double)p );
	}
	else
	{
		//p == -1
		return lan.compute_maximum();
	}
}

double compute_max_norm_discance_of_landscapes( const Persistence_landscape& first, const Persistence_landscape& second )
{
    return std::max( compute_maximal_distance_non_symmetric(first,second) , compute_maximal_distance_non_symmetric(second,first) );
}


bool comparePairsForMerging( std::pair< double , unsigned > first , std::pair< double , unsigned > second )
{
    return (first.first < second.first);
}




double compute_inner_product( const Persistence_landscape& l1 , const Persistence_landscape& l2 )
{
    bool dbg = false;
    double result = 0;

    for ( size_t level = 0 ; level != std::min( l1.size() , l2.size() ) ; ++level )
    {
        if ( dbg ){cerr << "Computing inner product for a level : " << level << endl;getchar();}
        if ( l1.land[level].size() * l2.land[level].size() == 0 )continue;

        //endpoints of the interval on which we will compute the inner product of two locally linear functions:
        double x1 = -std::numeric_limits<int>::max();
        double x2;
        if ( l1.land[level][1].first < l2.land[level][1].first )
        {
            x2 = l1.land[level][1].first;
        }
        else
        {
            x2 = l2.land[level][1].first;
        }

        //iterators for the landscapes l1 and l2
        size_t l1It = 0;
        size_t l2It = 0;

        while ( (l1It < l1.land[level].size()-1) && (l2It < l2.land[level].size()-1) )
        {
            //compute the value of a inner product on a interval [x1,x2]

            double a,b,c,d;

            a = (l1.land[level][l1It+1].second - l1.land[level][l1It].second)/(l1.land[level][l1It+1].first - l1.land[level][l1It].first);
            b = l1.land[level][l1It].second - a*l1.land[level][l1It].first;
            c = (l2.land[level][l2It+1].second - l2.land[level][l2It].second)/(l2.land[level][l2It+1].first - l2.land[level][l2It].first);
            d = l2.land[level][l2It].second - c*l2.land[level][l2It].first;

            double contributionFromThisPart
            =
            (a*c*x2*x2*x2/3 + (a*d+b*c)*x2*x2/2 + b*d*x2) - (a*c*x1*x1*x1/3 + (a*d+b*c)*x1*x1/2 + b*d*x1);

            result += contributionFromThisPart;

            if ( dbg )
            {
                cerr << "[l1.land[level][l1It].first,l1.land[level][l1It+1].first] : " << l1.land[level][l1It].first << " , " << l1.land[level][l1It+1].first << endl;
                cerr << "[l2.land[level][l2It].first,l2.land[level][l2It+1].first] : " << l2.land[level][l2It].first << " , " << l2.land[level][l2It+1].first << endl;
                cerr << "a : " << a << ", b : " << b << " , c: " << c << ", d : " << d << endl;
                cerr << "x1 : " << x1 << " , x2 : " << x2 << endl;
                cerr << "contributionFromThisPart : " << contributionFromThisPart << endl;
                cerr << "result : " << result << endl;
                getchar();
            }

            //we have two intervals in which functions are constant:
            //[l1.land[level][l1It].first , l1.land[level][l1It+1].first]
            //and
            //[l2.land[level][l2It].first , l2.land[level][l2It+1].first]
            //We also have an interval [x1,x2]. Since the intervals in the landscapes cover the whole R, then it is clear that x2
            //is either l1.land[level][l1It+1].first of l2.land[level][l2It+1].first or both. Lets test it.
            if ( x2 == l1.land[level][l1It+1].first )
            {
                if ( x2 == l2.land[level][l2It+1].first )
                {
                    //in this case, we increment both:
                    ++l2It;
                    if ( dbg ){cerr << "Incrementing both \n";}
                }
                else
                {
                    if ( dbg ){cerr << "Incrementing first \n";}
                }
                ++l1It;
            }
            else
            {
                //in this case we increment l2It
                ++l2It;
                if ( dbg ){cerr << "Incrementing second \n";}
            }
            //Now, we shift x1 and x2:
            x1 = x2;
            if ( l1.land[level][l1It+1].first < l2.land[level][l2It+1].first )
            {
                x2 = l1.land[level][l1It+1].first;
            }
            else
            {
                x2 = l2.land[level][l2It+1].first;
            }

        }

    }
    return result;
}


void Persistence_landscape::plot( const char* filename , int from, int to , double xRangeBegin , double xRangeEnd , double yRangeBegin , double yRangeEnd )
{
    //this program create a gnuplot script file that allows to plot persistence diagram.
    ofstream out;

    std::ostringstream nameSS;
    nameSS << filename << "_GnuplotScript";
    std::string nameStr = nameSS.str();
    out.open( (char*)nameStr.c_str() );

    if ( (xRangeBegin != -1) || (xRangeEnd != -1) || (yRangeBegin != -1) || (yRangeEnd != -1)  )
    {
        out << "set xrange [" << xRangeBegin << " : " << xRangeEnd << "]" << endl;
        out << "set yrange [" << yRangeBegin << " : " << yRangeEnd << "]" << endl;
    }

    if ( from == -1 ){from = 0;}
    if ( to == -1 ){to = this->land.size();}

    out << "plot ";
    for ( size_t lambda= std::min((size_t)from,this->land.size()) ; lambda != std::min((size_t)to,this->land.size()) ; ++lambda )
    {
        out << "     '-' using 1:2 title 'l" << lambda << "' with lp";
        if ( lambda+1 != std::min((size_t)to,this->land.size()) )
        {
            out << ", \\";
        }
        out << endl;
    }

    for ( size_t lambda= std::min((size_t)from,this->land.size()) ; lambda != std::min((size_t)to,this->land.size()) ; ++lambda )
    {
        for ( size_t i = 1 ; i != this->land[lambda].size()-1 ; ++i )
        {
            out << this->land[lambda][i].first << " " << this->land[lambda][i].second << endl;
        }
        out << "EOF" << endl;
    }
    cout << "Gnuplot script to visualize persistence diagram written to the file: " << nameStr << ". Type load '" << nameStr << "' in gnuplot to visualize." << endl;
}




}//namespace gudhi stat
}//namespace gudhi


#endif

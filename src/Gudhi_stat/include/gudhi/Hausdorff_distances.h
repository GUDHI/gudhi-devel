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

#ifndef HAUSDORFF_DISTANCES_H
#define HAUSDORFF_DISTANCES_H

#include <cmath>
#include <limits>
#include <vector>
#include <cstdlib>
#include <iostream>


/**
 * This file contains various implementations of Hausrodff distances between collections of points. It contains various implementations that can be used for specific case.
**/ 

/**
 * The implementation below works for a case of a metric space given by a distance matrix, and a subspace of a metric space. Our task is to find a Hausdorff distance between the subspace and the whole space. 
 * The input is a distance matrix (lower triangular part) and a vector of bools indicating a subspace (elementss set to true belongs to the space, the elements set to false do not). 
**/ 
class Hausdorff_distance_between_subspace_and_the_whole_metric_space
{
public:
	Hausdorff_distance_between_subspace_and_the_whole_metric_space( const std::vector< std::vector<double> >& distance_matrix ):distance_matrix(distance_matrix){}	
	double operator()( const std::vector< bool >& is_subspace  )
	{
		double maximal_distance = -std::numeric_limits<double>::max();
		for ( size_t j = 0 ; j != this->distance_matrix.size() ; ++j )  
		{			
			double minimal_distance = std::numeric_limits<double>::max();
			for ( size_t i = 0 ; i != this->distance_matrix.size() ; ++i )
			{
				if ( !is_subspace[i] )continue;
				double distance = 0;
				if ( i < j )
				{
					distance = this->distance_matrix[i][j];
				}
				else
				{
					if ( i > j )distance = this->distance_matrix[j][i];
				}
				if ( distance > minimal_distance )minimal_distance = distance;
			}		
			if ( maximal_distance < minimal_distance )maximal_distance = minimal_distance;
		}
		return maximal_distance;	
	}//Hausdorff_distance_between_subspace_and_the_whole_metric_space
	
	double operator()( const std::vector< size_t >& subspace , const std::vector< size_t >& space = std::vector< size_t >()  )
	{
		bool dbg = false;
		if ( dbg )
		{
			std::cerr << "Calling double operator()( const std::vector< size_t >& subspace , const std::vector< size_t >& space = std::vector< size_t >()  ) method \n";
			std::cerr << "subspace.size() : " << subspace.size() << std::endl;
		}
		double maximal_distance = -std::numeric_limits<double>::max();
		for ( size_t j = 0 ; j != this->distance_matrix.size() ; ++j )
		{			
			double minimal_distance = std::numeric_limits<double>::max();
			for ( size_t i = 0 ; i != subspace.size() ; ++i )  
			{				
				double distance = 0;
				if ( subspace[i] < j )
				{
					distance = this->distance_matrix[ j ][ subspace[i] ];
				}
				else
				{
					if ( subspace[i] > j )distance = this->distance_matrix[ subspace[i] ][ j ];
				}
				if ( distance < minimal_distance )minimal_distance = distance;
			}
			if ( maximal_distance < minimal_distance )maximal_distance = minimal_distance;
		}
		if ( dbg )
		{
			std::cerr << "maximal_distance : " << maximal_distance << std::endl;
			getchar();
		}
		return maximal_distance;	
	}//Hausdorff_distance_between_subspace_and_the_whole_metric_space
private:
	const std::vector< std::vector<double> >& distance_matrix;
};

template <typename PointType , typename distanceFunction >
std::vector< std::vector<double> > compute_all_to_all_distance_matrix_between_points( const std::vector< PointType >& points )
{
	bool dbg = false;
	std::vector< std::vector<double> > result;
	result.reserve( points.size() );
	distanceFunction f;
	for ( size_t i = 0 ; i != points.size() ; ++i )
	{
		std::vector<double> this_row;
		this_row.reserve( i );
		for ( size_t j = 0  ; j != i ; ++j )
		{
			double distance = f( points[i] , points[j] );
			if ( dbg ){std::cerr << "The distance between point : " << i << " and the point : " << j << " is :" << distance << std::endl;}
			this_row.push_back( distance );
		}
		result.push_back( this_row );
	}
	return result;
}//compute_all_to_all_distance_matrix_between_points

template <typename T>
class identity
{
public:
	T& operator()( T& org )
	{
		return org;
	}
};



#endif

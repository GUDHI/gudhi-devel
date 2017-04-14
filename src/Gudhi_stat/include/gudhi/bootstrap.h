/*    Thicharacteristic_of_all_pointss file is part of the Gudhi Library. The Gudhi library
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

#ifndef BOOTSTRAP_H
#define BOOTSTRAP_H


#ifdef GUDHI_USE_TBB
#include <tbb/parallel_sort.h>
#include <tbb/task_scheduler_init.h>
#endif

#include <vector>
#include <algorithm>
#include <random>
#include <ctime>

namespace Gudhi
{
namespace Gudhi_stat
{

/**
* This is a generic function to perform bootstrap.
* In this function we assume that there is a class to compute characteristic of collection of points (PointCloudCharacteristics) and that it stores coordinates of all points. It allows to compute the characteristic
* of the whole point cloud (by using CharacteristicFunction) or of it proper subset of the whole point cloud (given the list of numers of points in the subset). 
* Both functionalities will be used in this implementation.
* The characteristic of point cloud, may be the poit cloud itself, its persistence diagram in a fixed dimension, or anything else. We only assume that space of points characteristics is a metric space
* and that we can compute a distance between two characteristics of collections of points by using DistanceBetweenPointsCharacteristics function. 
**/



template < typename PointCloudCharacteristics , typename CharacteristicFunction , typename DistanceBetweenPointsCharacteristics >
double bootstrap( size_t number_of_points , CharacteristicFunction f , DistanceBetweenPointsCharacteristics distance , size_t number_of_repetitions , size_t size_of_subsample , double quantile = 0.95 , size_t maximal_number_of_threads_in_TBB = std::numeric_limits<size_t>::max() )
{
	bool dbg = false;
	
	#ifdef GUDHI_USE_TBB
	tbb::task_scheduler_init init(maximal_number_of_threads_in_TBB == std::numeric_limits<size_t>::max() ? tbb::task_scheduler_init::automatic : maximal_number_of_threads_in_TBB);
	#endif
	
	if ( size_of_subsample >= number_of_points )
	{
		std::cerr << "Size of subsample is greater or equal to the number of points. The bootstrap procedure do not make sense in this case. \n";
		return 0;
	}
	
	//initialization of a random number generator:
    std::srand ( unsigned ( std::time(0) ) );
	
	//we will shuffle the vector of numbers 0,1,2,...,points.size()-1 in order to pick a subset of a size size_of_subsample		
	std::vector<size_t> numbers_to_sample_(number_of_points) ; //create vector of size_t of a size number_of_points
	std::iota (std::begin(numbers_to_sample_), std::end(numbers_to_sample_), 0);//populate it with 1 2 3 ... number_of_points.
	
	//now we compute the characteristic od all the points:	
	PointCloudCharacteristics characteristic_of_all_points = f( numbers_to_sample_ );
	
	//vector to keep the distances between characteristic_of_points and characteristic_of_subsample:
	std::vector< double > vector_of_distances( number_of_repetitions , 0 );
	


//TODO- at the moment, the operations I am doing over here do not seems to be threat safe. When using TBB, I am getting wrong results.
//It is quite likelly because I am not using a method to compute persistence which is threat safe. VERIFY this as soon as I merge with 
//the new metod to compute persistence. 

//	#ifdef GUDHI_USE_TBB
//    tbb::parallel_for ( tbb::blocked_range<size_t>(0, number_of_repetitions), [&](const tbb::blocked_range<size_t>& range) 
//    {
//    for  ( size_t it_no = range.begin() ;  it_no != range.end() ; ++it_no )
//	#else
	for ( size_t it_no = 0 ;  it_no < number_of_repetitions ; ++it_no )
//	#endif	
	{
		if ( dbg )
		{
			std::cout << "Still : " << number_of_repetitions-it_no << " tests to go. \n The subsampled vector consist of points number : ";
			std::cout << "it_no : " << it_no << std::endl;
			std::cout << "number_of_points : " << number_of_points << std::endl;
		}
		//do a random shuffle of vector_of_characteristics_of_poits
		std::vector<size_t> numbers_to_sample(number_of_points) ; //create vector of size_t of a size number_of_points
	    std::iota (std::begin(numbers_to_sample), std::end(numbers_to_sample), 0);//populate it with 1 2 3 ... number_of_points.	
		//TODO: consider doing it in a smarter/faster way.
		std::random_shuffle( numbers_to_sample.begin() , numbers_to_sample.end() );
		
		//construct a vector< PointType > of a size size_of_subsample:
		std::vector< size_t > subsampled_points;
		subsampled_points.reserve( size_of_subsample );
		for ( size_t i = 0 ; i != size_of_subsample ; ++i )
		{
			subsampled_points.push_back( numbers_to_sample[i] );
			if ( dbg )std::cout << numbers_to_sample[i] << " , ";
		}
		
		
		//now we can compute characteristic of subsampled_points:
		PointCloudCharacteristics characteristic_of_subsampled_points = f( subsampled_points );
		if ( dbg )std::cout << std::endl << "Characteristic of subsampled points computed.\n";					
		
		//and now we compute distance between characteristic_of_points and characteristic_of_subsample. Note that subsampled points go first, and this is neded, since sometimes all points are not needed.
		double dist = distance( characteristic_of_subsampled_points , characteristic_of_all_points );
		
		if ( dbg )
		{
			 std::cout << "The distance between characteristic of all points and the characteristic of subsample is : " << dist << std::endl;
			 getchar();
		 }
		
		vector_of_distances[it_no] = dist;		
	}
//	#ifdef GUDHI_USE_TBB
//    }
//    );
//	#endif
			
	size_t position_of_quantile = floor(quantile*vector_of_distances.size());
	if ( position_of_quantile ) --position_of_quantile;
	
	
	
	if ( dbg )
	{
		std::cerr << "quantile : " << quantile << std::endl;
		std::cerr << "position_of_quantile : " << position_of_quantile << std::endl;
		
		std::sort( vector_of_distances.begin() , vector_of_distances.end() );
		//std::cout << "position_of_quantile : " << position_of_quantile << ", and here is the array : " << std::endl;
		for ( size_t i = 0 ; i != vector_of_distances.size() ; ++i )
		{
			std::cout << vector_of_distances[i] << " " ;
		}
		std::cout << std::endl;
	}
	  
	//now we need to sort the vector_of_distances and find the quantile:
	std::nth_element (vector_of_distances.begin(), vector_of_distances.begin()+position_of_quantile, vector_of_distances.end());
	
	
	//for Hausdorff bootrstra I have to multily it by 2.
	//In case of other bootsraps, I do not have to do it. We need a special variable saying if Ineed this multiplication or not.//This should be done outside the bootstrap, since the fact hat we need it do not come from bootstrab, but from geometry of bottleneck distance
	
	if ( dbg )std::cout << "Result : " << vector_of_distances[ position_of_quantile ] << std::endl;
	
	return vector_of_distances[ position_of_quantile ];
	
}//bootstrap



}//namespace Gudhi_stat
}//namespace Gudhi

#endif

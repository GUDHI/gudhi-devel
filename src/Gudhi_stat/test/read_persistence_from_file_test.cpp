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



#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "gudhi_stat"
#include <boost/test/unit_test.hpp>
#include <gudhi/reader_utils.h>
#include <gudhi/read_persistence_from_file.h>

#include <iostream>



using namespace Gudhi;
using namespace Gudhi::Gudhi_stat;



	

BOOST_AUTO_TEST_CASE(test_read_file_with_four_elements_per_line) 
{	
	std::vector< std::pair< double,double > > what_we_should_get;
	what_we_should_get.push_back( std::make_pair( 0 , 2 ) );
	what_we_should_get.push_back( std::make_pair( 10 , 1000 ) );
	what_we_should_get.push_back( std::make_pair( 10 , 90 ) );
	what_we_should_get.push_back( std::make_pair( 4 , 4 ) );		
	std::vector<std::pair<double,double> > what_we_get = read_persistence_intervals_in_one_dimension_from_file("data/persistence_file_with_four_entries_per_line" , 1 , 1000 );
	
	
	//for ( size_t i = 0 ; i != what_we_get.size() ; ++i )
	//{
	//	std::cerr << what_we_get[i].first << " , " << what_we_get[i].second << std::endl;
	//}
	
	BOOST_CHECK( what_we_should_get.size() == what_we_get.size() );
	
	for ( size_t i = 0 ; i != what_we_get.size() ; ++i )
	{
		BOOST_CHECK( what_we_should_get[i] == what_we_get[i] );
	}
}


BOOST_AUTO_TEST_CASE(test_read_file_with_three_elements_per_line) 
{	
	std::vector< std::pair< double,double > > what_we_should_get;	
	what_we_should_get.push_back( std::make_pair( 4 , 9999 ) );
	what_we_should_get.push_back( std::make_pair( 0 , 1 ) );
	what_we_should_get.push_back( std::make_pair( 1 , 9999 ) );
	what_we_should_get.push_back( std::make_pair( 10 , 90 ) );
	what_we_should_get.push_back( std::make_pair( 4 , 4 ) );
	
	std::vector<std::pair<double,double> > what_we_get = read_persistence_intervals_in_one_dimension_from_file("data/persistence_file_with_three_entries_per_line" , 1 , 9999 );
	
	
	//for ( size_t i = 0 ; i != what_we_get.size() ; ++i )
	//{
	//	std::cerr << what_we_get[i].first << " , " << what_we_get[i].second << std::endl;
	//}
	
	BOOST_CHECK( what_we_should_get.size() == what_we_get.size() );
	
	for ( size_t i = 0 ; i != what_we_get.size() ; ++i )
	{
		BOOST_CHECK( what_we_should_get[i] == what_we_get[i] );
	}
}


BOOST_AUTO_TEST_CASE(test_read_file_with_two_elements_per_line) 
{	
	std::vector< std::pair< double,double > > what_we_should_get;	
	what_we_should_get.push_back( std::make_pair( 4 , 10 ) );
	what_we_should_get.push_back( std::make_pair( 4 , 9999 ) );
	what_we_should_get.push_back( std::make_pair( 0 , 1 ) );
	what_we_should_get.push_back( std::make_pair( 1 , 4 ) );	
	
	std::vector<std::pair<double,double> > what_we_get = read_persistence_intervals_in_one_dimension_from_file("data/persistence_file_with_two_entries_per_line" , 1 , 9999 );
	
	
	//for ( size_t i = 0 ; i != what_we_get.size() ; ++i )
	//{
	//	std::cerr << what_we_get[i].first << " , " << what_we_get[i].second << std::endl;
	//}
	
	BOOST_CHECK( what_we_should_get.size() == what_we_get.size() );
	
	for ( size_t i = 0 ; i != what_we_get.size() ; ++i )
	{
		BOOST_CHECK( what_we_should_get[i] == what_we_get[i] );
	}
}


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
#include <gudhi/abstract_classes/Abs_Topological_data.h>
#include <gudhi/concretizations/Persistence_landscape.h>

#include <iostream>



using namespace Gudhi;
using namespace Gudhi::Gudhi_stat;


double epsilon = 0.0000005;

using namespace std;

	

BOOST_AUTO_TEST_CASE(check_construction_of_landscape) 
{	
	Persistence_landscape p( "data/file_with_diagram" );
	
	Persistence_landscape q;
	q.load_landscape_from_file( "data/file_with_landscape_from_file_with_diagram" );
		
	BOOST_CHECK( p == q );
}


BOOST_AUTO_TEST_CASE(check_computations_of_integrals) 
{
	Persistence_landscape p( "data/file_with_diagram" );	
	double integral = p.compute_integral_of_landscape();
	//cerr << integral  << " " << 2.34992 << endl;
	BOOST_CHECK( fabs( integral - 2.34992 ) <= 0.00001 );
}


BOOST_AUTO_TEST_CASE(check_computations_of_integrals_for_each_level_separatelly) 
{
	Persistence_landscape p( "data/file_with_diagram" );	
	
	std::vector< double > integrals_fir_different_levels;
	integrals_fir_different_levels.push_back(	0.216432	);
	integrals_fir_different_levels.push_back(	0.204763	);
	integrals_fir_different_levels.push_back(	0.188793	);
	integrals_fir_different_levels.push_back(	0.178856	);
	integrals_fir_different_levels.push_back(	0.163142	);
	integrals_fir_different_levels.push_back(	0.155015	);
	integrals_fir_different_levels.push_back(	0.143046	);
	integrals_fir_different_levels.push_back(	0.133765	);
	integrals_fir_different_levels.push_back(	0.123531	);
	integrals_fir_different_levels.push_back(	0.117393	);
	integrals_fir_different_levels.push_back(	0.111269	);
	integrals_fir_different_levels.push_back(	0.104283	);
	integrals_fir_different_levels.push_back(	0.0941308	);
	integrals_fir_different_levels.push_back(	0.0811208	);
	integrals_fir_different_levels.push_back(	0.0679001	);
	integrals_fir_different_levels.push_back(	0.0580801	);
	integrals_fir_different_levels.push_back(	0.0489647	);
	integrals_fir_different_levels.push_back(	0.0407936	);
	integrals_fir_different_levels.push_back(	0.0342599	);
	integrals_fir_different_levels.push_back(	0.02896	);
	integrals_fir_different_levels.push_back(	0.0239881	);
	integrals_fir_different_levels.push_back(	0.0171792	);
	integrals_fir_different_levels.push_back(	0.0071511	);
	integrals_fir_different_levels.push_back(	0.00462067	);
	integrals_fir_different_levels.push_back(	0.00229033	);
	integrals_fir_different_levels.push_back(	0.000195296	);


	
	
	for ( size_t level = 0 ; level != p.size() ; ++level )
	{
		double integral = p.compute_integral_of_landscape( level );
		BOOST_CHECK( fabs( integral - integrals_fir_different_levels[level] ) <= 0.00001 );
	}
	
}

BOOST_AUTO_TEST_CASE(check_computations_of_integrals_of_powers_of_landscape) 
{
	Persistence_landscape p( "data/file_with_diagram" );	
	
	std::vector<double> integrals_fir_different_powers;
 	integrals_fir_different_powers.push_back(	0.216432	);
	integrals_fir_different_powers.push_back(	0.204763	);
	integrals_fir_different_powers.push_back(	0.188793	);
	integrals_fir_different_powers.push_back(	0.178856	);
	integrals_fir_different_powers.push_back(	0.163142	);
	
	for ( size_t power = 0 ; power != 5 ; ++power ) 
	{
		double integral = p.compute_integral_of_landscape( power );
		BOOST_CHECK( fabs( integral - integrals_fir_different_powers[power] ) <= 0.00001 );
	}
}

BOOST_AUTO_TEST_CASE(check_computations_of_values_on_different_points) 
{
	Persistence_landscape p( "data/file_with_diagram" );	
	
	
	BOOST_CHECK( fabs( p.compute_value_at_a_given_point(1,0.0)  ) <= 0.00001 );
	BOOST_CHECK( fabs(  p.compute_value_at_a_given_point(1,0.1) - 0.0692324 ) <= 0.00001 );
	BOOST_CHECK( fabs( p.compute_value_at_a_given_point(1,0.2) - 0.163369 ) <= 0.00001 );
	BOOST_CHECK( fabs( p.compute_value_at_a_given_point(1,0.3) - 0.217115 ) <= 0.00001 );
	BOOST_CHECK( fabs( p.compute_value_at_a_given_point(2,0.0) ) <= 0.00001 );
	BOOST_CHECK( fabs( p.compute_value_at_a_given_point(2,0.1) - 0.0633688 ) <= 0.00001 );
	BOOST_CHECK( fabs( p.compute_value_at_a_given_point(2,0.2) - 0.122361 ) <= 0.00001 );
	BOOST_CHECK( fabs( p.compute_value_at_a_given_point(2,0.3) - 0.195401 ) <= 0.00001 );			
	BOOST_CHECK( fabs( p.compute_value_at_a_given_point(3,0.0)  ) <= 0.00001 );
	BOOST_CHECK( fabs( p.compute_value_at_a_given_point(3,0.1) - 0.0455386 ) <= 0.00001 );
	BOOST_CHECK( fabs( p.compute_value_at_a_given_point(3,0.2) - 0.0954012 ) <= 0.00001 );
	BOOST_CHECK( fabs( p.compute_value_at_a_given_point(3,0.3) - 0.185282 ) <= 0.00001 );
}


BOOST_AUTO_TEST_CASE(check_computations_sum_differences_and_multiplications) 
{	
	Persistence_landscape p( "data/file_with_diagram" );	
	Persistence_landscape second;
	second.load_landscape_from_file("data/file_with_landscape_from_file_with_diagram_1" );
	
	Persistence_landscape sum = p + second;
	Persistence_landscape difference = p - second;
	Persistence_landscape multiply_by_scalar = 10*p;

	
	Persistence_landscape template_sum;
	template_sum.load_landscape_from_file( "data/sum" );
	
	Persistence_landscape template_difference;
	template_difference.load_landscape_from_file( "data/difference" );
	
	Persistence_landscape template_multiply_by_scalar;
	template_multiply_by_scalar.load_landscape_from_file( "data/multiply_by_scalar" );
	
	BOOST_CHECK( sum == template_sum );
	BOOST_CHECK( difference == template_difference );
	BOOST_CHECK( multiply_by_scalar == template_multiply_by_scalar );		
}



BOOST_AUTO_TEST_CASE(check_computations_of_maxima_and_norms) 
{	
	Persistence_landscape p( "data/file_with_diagram" );	
	Persistence_landscape second;
	second.load_landscape_from_file("data/file_with_landscape_from_file_with_diagram_1" );	
	Persistence_landscape sum = p + second;
	
	BOOST_CHECK( fabs( p.compute_maximum() - 0.431313 ) <= 0.00001 );
	BOOST_CHECK( fabs( p.compute_norm_of_landscape(1) - 2.34992 ) <= 0.00001 );
	BOOST_CHECK( fabs( p.compute_norm_of_landscape(2) - 0.706095 ) <= 0.00001 );
	BOOST_CHECK( fabs( p.compute_norm_of_landscape(3) - 0.501867 ) <= 0.00001 );
	BOOST_CHECK( fabs( compute_discance_of_landscapes(p,sum,1) - 27.9323 ) <= 0.00005 );	
	BOOST_CHECK( fabs( compute_discance_of_landscapes(p,sum,2) - 2.35199 ) <= 0.00001 );
	BOOST_CHECK( fabs(compute_discance_of_landscapes(p,sum,-1) - 0.464478 ) <= 0.00001 );		
}

BOOST_AUTO_TEST_CASE(check_computations_of_averages) 
{
	Persistence_landscape p( "data/file_with_diagram" );
	Persistence_landscape q( "data/file_with_diagram_1" );
	std::vector< Persistence_landscape* > to_average;
	to_average.push_back( &p );
	to_average.push_back( &q );
	Persistence_landscape av;
	av.compute_average( to_average );

	Persistence_landscape template_average;
	template_average.load_landscape_from_file( "data/average" );
	BOOST_CHECK ( template_average == av );
}




BOOST_AUTO_TEST_CASE(check_computations_of_distances)
{
	Persistence_landscape p( "data/file_with_diagram" );
	Persistence_landscape q( "data/file_with_diagram_1" );
	BOOST_CHECK( fabs( p.distance( q )- 25.5824) <= 0.00005 );	
	BOOST_CHECK( fabs( p.distance( q , 2) - 2.12636 ) <= 0.00001 );	
	BOOST_CHECK( fabs( p.distance( q , -1 )-0.359068 ) <= 0.00001 );	
}
	

BOOST_AUTO_TEST_CASE(check_computations_of_scalar_product)
{
	Persistence_landscape p( "data/file_with_diagram" );
	Persistence_landscape q( "data/file_with_diagram_1" );
	BOOST_CHECK( fabs(  p.compute_scalar_product( q ) - 0.754498 ) <= 0.00001 );	
}


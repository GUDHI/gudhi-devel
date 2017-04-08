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
#include <gudhi/persistence_representations/Persistence_landscape_on_grid.h>

#include <iostream>



using namespace Gudhi;
using namespace Gudhi::Gudhi_stat;


double epsilon = 0.0000005;



	

BOOST_AUTO_TEST_CASE(check_construction_of_landscape) 
{	
	Persistence_landscape_on_grid l( "data/file_with_diagram_1" , 100 );
	l.print_to_file( "landscape_from_file_with_diagram_1" );
	
	Persistence_landscape_on_grid g;
	g.load_landscape_from_file( "landscape_from_file_with_diagram_1" );
	
	BOOST_CHECK( l == g );
}

BOOST_AUTO_TEST_CASE(check_construction_of_landscape_using_only_ten_levels) 
{	
	//TODO
	size_t number = 10;
	Persistence_landscape_on_grid l( "data/file_with_diagram_1" , 100  ,number );	
	Persistence_landscape_on_grid g( "data/file_with_diagram_1" , 100 );
	//cut all the elements of order > 10 in g. 
	
	for ( size_t level = 0 ; level != number ; ++level )
	{
		std::vector<double> v1 = l.vectorize(level);
		std::vector<double> v2 = g.vectorize(level);				
		BOOST_CHECK( v1.size() == v2.size() );
		for ( size_t i = 0 ; i != v1.size() ; ++i )
		{
			BOOST_CHECK( v1[i] == v2[i] );
		}
	}	
}



BOOST_AUTO_TEST_CASE(check_computations_of_integrals) 
{
	Persistence_landscape_on_grid p( "data/file_with_diagram_1" , 100 );	
	double integral = p.compute_integral_of_landscape();
	//cerr << "integral : " << integral << endl;getchar();
	BOOST_CHECK( fabs( integral - 27.343 ) <= 0.00005 );
}


BOOST_AUTO_TEST_CASE(check_computations_of_integrals_for_each_level_separatelly) 
{
	Persistence_landscape_on_grid p( "data/file_with_diagram_1" , 100 );	
	
	std::vector< double > integrals_fir_different_levels;	
	//integrals_fir_different_levels.push_back();
	integrals_fir_different_levels.push_back(0.241168);
	integrals_fir_different_levels.push_back(0.239276);
	integrals_fir_different_levels.push_back(0.237882);
	integrals_fir_different_levels.push_back(0.235193);
	integrals_fir_different_levels.push_back(0.230115);
	integrals_fir_different_levels.push_back(0.227626);
	integrals_fir_different_levels.push_back(0.226132);
	integrals_fir_different_levels.push_back(0.223643);
	integrals_fir_different_levels.push_back(0.221651);
	integrals_fir_different_levels.push_back(0.220556);
	integrals_fir_different_levels.push_back(0.21727);
	integrals_fir_different_levels.push_back(0.215976);
	integrals_fir_different_levels.push_back(0.213685);
	integrals_fir_different_levels.push_back(0.211993);
	integrals_fir_different_levels.push_back(0.2102);
	integrals_fir_different_levels.push_back(0.208707);
	integrals_fir_different_levels.push_back(0.207014);
	integrals_fir_different_levels.push_back(0.205122);
	integrals_fir_different_levels.push_back(0.204226);
	integrals_fir_different_levels.push_back(0.202633);

	
	for ( size_t level = 0 ; level != integrals_fir_different_levels.size() ; ++level )
	{
		double integral = p.compute_integral_of_landscape( level );		
		//cerr << integral << endl;
		BOOST_CHECK( fabs( integral - integrals_fir_different_levels[level] ) <= 0.00005 );
	}	
	
}



BOOST_AUTO_TEST_CASE(check_computations_of_integrals_of_powers_of_landscape) 
{
	Persistence_landscape_on_grid p( "data/file_with_diagram_1" , 100 );	
	
	std::vector<double> integrals_fir_different_powers;
	integrals_fir_different_powers.push_back(	0.241168);
	integrals_fir_different_powers.push_back(	0.239276);
	integrals_fir_different_powers.push_back(	0.237882);
	integrals_fir_different_powers.push_back(	0.235193);
	integrals_fir_different_powers.push_back(	0.23011);
 
	for ( size_t power = 0 ; power != 5 ; ++power ) 
	{
		double integral = p.compute_integral_of_landscape( power );
		//cerr << integral << endl;
		BOOST_CHECK( fabs( integral - integrals_fir_different_powers[power] ) <= 0.00001 );
	}
}


BOOST_AUTO_TEST_CASE(check_computations_of_values_on_different_points) 
{
	Persistence_landscape_on_grid p( "data/file_with_diagram_1" , 100 );	

	std::vector< double > results_level_0;	
	results_level_0.push_back(0.00997867);
	results_level_0.push_back(0.0521921);
	results_level_0.push_back(0.104312);
	results_level_0.push_back(0.156432);
	results_level_0.push_back(0.208552);
	results_level_0.push_back(0.260672);
	results_level_0.push_back(0.312792);
	results_level_0.push_back(0.364912);
	results_level_0.push_back(0.417032);
	results_level_0.push_back(0.429237);
	
	std::vector< double > results_level_10;
	results_level_10.push_back(7.21433e-05);
	results_level_10.push_back(0.0422135);
	results_level_10.push_back(0.0943335);
	results_level_10.push_back(0.146453);
	results_level_10.push_back(0.198573);
	results_level_10.push_back(0.240715);
	results_level_10.push_back(0.272877);
	results_level_10.push_back(0.324997);
	results_level_10.push_back(0.359232);
	results_level_10.push_back(0.379344);
	
	double x = 0.0012321;
	 double dx = 0.05212;
	 for ( size_t i = 0 ; i != 10 ; ++i )
	 {  		
		BOOST_CHECK( almost_equal( p.compute_value_at_a_given_point(0,x) , results_level_0[i] ) );
		BOOST_CHECK( almost_equal( p.compute_value_at_a_given_point(10,x) , results_level_10[i] ) );
		x += dx;
	 }
}


BOOST_AUTO_TEST_CASE(check_computations_sum_differences_and_multiplications) 
{	
	Persistence_landscape_on_grid p( "data/file_with_diagram_1" ,100 );	
	Persistence_landscape_on_grid second("data/file_with_diagram_1" , 100 );	
	
	Persistence_landscape_on_grid sum = p + second;
	Persistence_landscape_on_grid difference = p - second;
	Persistence_landscape_on_grid multiply_by_scalar = 10*p;	;

	
	Persistence_landscape_on_grid template_sum;
	template_sum.load_landscape_from_file( "data/sum_on_grid_test" );
	
	Persistence_landscape_on_grid template_difference;
	template_difference.load_landscape_from_file( "data/difference_on_grid_test" );
	
	Persistence_landscape_on_grid template_multiply_by_scalar;
	template_multiply_by_scalar.load_landscape_from_file( "data/multiply_by_scalar_on_grid_test" );
	
	BOOST_CHECK( sum == template_sum );
	BOOST_CHECK( difference == template_difference );
	BOOST_CHECK( multiply_by_scalar == template_multiply_by_scalar );		
}



BOOST_AUTO_TEST_CASE(check_computations_of_maxima_and_norms) 
{	
	Persistence_landscape_on_grid p( "data/file_with_diagram_1" , 0 , 1 , 100 );	
	Persistence_landscape_on_grid second("data/file_with_diagram_2" , 0 , 1 , 100 );	
	Persistence_landscape_on_grid sum = p + second;
	
	//cerr << p.compute_maximum() << endl;
	//cerr <<  p.compute_norm_of_landscape(1) << endl;
	//cerr << p.compute_norm_of_landscape(2) << endl;
	//cerr <<  p.compute_norm_of_landscape(3) << endl;
	//cerr <<  compute_distance_of_landscapes_on_grid(p,sum,1) << endl;
	//cerr <<  compute_distance_of_landscapes_on_grid(p,sum,2) << endl;
	//cerr <<  compute_distance_of_landscapes_on_grid(p,sum,std::numeric_limits<double>::max())  << endl;
	
	BOOST_CHECK( fabs( p.compute_maximum() - 0.46 ) <= 0.00001 );
	BOOST_CHECK( fabs( p.compute_norm_of_landscape(1) - 27.3373 ) <= 0.00001 );
	BOOST_CHECK( fabs( p.compute_norm_of_landscape(2) - 1.84143 ) <= 0.00001 );
	BOOST_CHECK( fabs( p.compute_norm_of_landscape(3) - 0.927067 ) <= 0.00001 );
	BOOST_CHECK( fabs( compute_distance_of_landscapes_on_grid(p,sum,1) - 16.8519 ) <= 0.00005 );	
	BOOST_CHECK( fabs( compute_distance_of_landscapes_on_grid(p,sum,2) - 1.44542 ) <= 0.00001 );	
	BOOST_CHECK( fabs(compute_distance_of_landscapes_on_grid(p,sum,std::numeric_limits<double>::max()) - 0.45 ) <= 0.00001 );		
}




BOOST_AUTO_TEST_CASE(check_computations_of_averages) 
{
	Persistence_landscape_on_grid p( "data/file_with_diagram", 0,1,100 );
	Persistence_landscape_on_grid q( "data/file_with_diagram_1", 0,1,100 );
	Persistence_landscape_on_grid av;	
	av.compute_average( {&p,&q} );
	
	Persistence_landscape_on_grid template_average;
	template_average.load_landscape_from_file( "data/average_on_a_grid" );
	BOOST_CHECK ( template_average == av );
}




BOOST_AUTO_TEST_CASE(check_computations_of_distances)
{
	Persistence_landscape_on_grid p( "data/file_with_diagram", 0,1,10000 );
	Persistence_landscape_on_grid q( "data/file_with_diagram_1", 0,1,10000 );
	BOOST_CHECK( fabs( p.distance( q )- 25.5779) <= 0.00005 );	
	BOOST_CHECK( fabs( p.distance( q , 2) - 2.04891) <= 0.00001 );	
	BOOST_CHECK( fabs( p.distance( q , std::numeric_limits<double>::max() )-0.359		 ) <= 0.00001 );
}
	

BOOST_AUTO_TEST_CASE(check_computations_of_scalar_product)
{
	Persistence_landscape_on_grid p( "data/file_with_diagram" , 0,1,10000);
	Persistence_landscape_on_grid q( "data/file_with_diagram_1", 0,1,10000 );
	//std::cerr << p.compute_scalar_product( q ) << std::endl;
	BOOST_CHECK( almost_equal(  p.compute_scalar_product( q ) , 0.754367 ) );	
}




//Below I am storing the code used to generate tests for that functionality.
/*
	Persistence_landscape_on_grid l( "file_with_diagram_1" , 100 );
	l.print_to_file( "landscape_from_file_with_diagram_1" );
	
	Persistence_landscape_on_grid g;
	g.load_landscape_from_file( "landscape_from_file_with_diagram_1" );
	
	cerr << ( l == g );
	*/
	
	/*
	Persistence_landscape_on_grid l( "file_with_diagram_1" , 100 );
	cerr << l << endl;
	cerr << l.compute_integral_of_landscape() << endl;
	*/
	 
	 /*
	 Persistence_landscape_on_grid p( "file_with_diagram_1" , 100 );	
	 for ( size_t level = 0 ; level != 30 ; ++level )
	 {
		 double integral = p.compute_integral_of_landscape( level );
		 cerr << integral << endl;
	 }
	 */
	 
	 /*
	 Persistence_landscape_on_grid p( "file_with_diagram_1" , 100 );	
	 for ( size_t power = 0 ; power != 5 ; ++power )
	 {
		 double integral = p.compute_integral_of_landscape( (double)power );
		 cerr << integral << endl;
	 }
	 */
	 
	 /*
	 Persistence_landscape_on_grid p( "file_with_diagram_1" , 100 );	
	 double x = 0.0012321;
	 double dx = 0.05212;
	 for ( size_t i = 0 ; i != 10 ; ++i )
	 {  		
		cerr << p.compute_value_at_a_given_point(10,x) << endl;
		x += dx;
	 }
	 */
	 
	 /*	
	Persistence_landscape_on_grid p( "file_with_diagram_1",100 );	
	Persistence_landscape_on_grid second("file_with_diagram_1",100 );		
	Persistence_landscape_on_grid sum = p + second;
	Persistence_landscape_on_grid difference = p - second;
	Persistence_landscape_on_grid multiply_by_scalar = 10*p;
	sum.print_to_file( "sum_on_grid_test" );
	difference.print_to_file( "difference_on_grid_test" );
	multiply_by_scalar.print_to_file( "multiply_by_scalar_on_grid_test" );
	*/
	
	
	/*	
	Persistence_landscape_on_grid p( "file_with_diagram_1" , 0 , 1 , 100 );
	Persistence_landscape_on_grid second("file_with_diagram_1", 0 , 1 , 100 );	
	Persistence_landscape_on_grid sum = p + second;
	
	cerr << "max : " << p.compute_maximum() << endl;
	cerr << "1-norm : " << p.compute_norm_of_landscape(1) << endl;
	cerr << "2-norm : " << p.compute_norm_of_landscape(2) << endl;
	cerr << "3-norm : " << p.compute_norm_of_landscape(3) << endl;
	
	cerr <<  compute_distance_of_landscapes_on_grid(p,sum,1) << endl;
	cerr <<  compute_distance_of_landscapes_on_grid(p,sum,2) << endl;
	cerr <<  compute_distance_of_landscapes_on_grid(p,sum,-1)  << endl;
	*/
	
	/*
	Persistence_landscape_on_grid p( "file_with_diagram", 0,1,100 );
	Persistence_landscape_on_grid q( "file_with_diagram_1", 0,1,100 );		
	Persistence_landscape_on_grid av;	
	av.compute_average( {&p,&q} );
	av.print_to_file("average_on_a_grid");

	Persistence_landscape_on_grid template_average;
	template_average.load_landscape_from_file( "average_on_a_grid" );
	if ( template_average == av )
	{
		cerr << "OK OK \n";
	}*/
	
	/*
	Persistence_landscape_on_grid p( "file_with_diagram" , 0,1,10000);
	Persistence_landscape_on_grid q( "file_with_diagram_1" , 0,1,10000);
	cerr <<  p.distance( &q )<< endl;
	cerr <<  p.distance( &q , 2 ) << endl;
	cerr <<  p.distance( &q , std::numeric_limits<double>::max() ) << endl;
	*/
	
/*
	Persistence_landscape_on_grid p( "file_with_diagram", 0,1,10000 );
	Persistence_landscape_on_grid q( "file_with_diagram_1", 0,1,10000 );
	
	//std::vector< std::pair< double,double > > aa;
	//aa.push_back( std::make_pair( 0,1 ) );
	//Persistence_landscape_on_grid p( aa, 0,1,10 );
	//Persistence_landscape_on_grid q( aa, 0,1,10 );
	cerr <<  p.compute_scalar_product( &q ) << endl;
*/

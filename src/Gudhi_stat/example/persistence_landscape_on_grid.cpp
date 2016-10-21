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



#include <gudhi/abstract_classes/Abs_Topological_data.h>
#include <gudhi/concretizations/Persistence_landscape_on_grid.h>



using namespace Gudhi;
using namespace Gudhi::Gudhi_stat;

#include <iostream>


int main( int argc , char** argv )
{
	//create two simple vectors with birth--death pairs:
	
	std::vector< std::pair< double , double > > persistence1;
	std::vector< std::pair< double , double > > persistence2;
	
	persistence1.push_back( std::make_pair(1,2) );
	persistence1.push_back( std::make_pair(6,8) );
	persistence1.push_back( std::make_pair(0,4) );
	persistence1.push_back( std::make_pair(3,8) );
	
	persistence2.push_back( std::make_pair(2,9) );
	persistence2.push_back( std::make_pair(1,6) );
	persistence2.push_back( std::make_pair(3,5) );
	persistence2.push_back( std::make_pair(6,10) );
	
	//create two persistence landscapes based on persistence1 and persistence2:
	Persistence_landscape_on_grid l1( persistence1 , 0 , 11 , 20 );
	Persistence_landscape_on_grid l2( persistence2 , 0 , 11 , 20 );
	
	//This is how to compute integral of landscapes:
	std::cout << "Integral of the first landscape : " << l1.compute_integral_of_landscape() << endl;
	std::cout << "Integral of the second landscape : " << l2.compute_integral_of_landscape() << endl;
	
	//And here how to write landscapes to stream:
	std::cout << "l1 : " << l1 << std::endl;
	std::cout << "l2 : " << l2 << std::endl;
	
	//here are the maxima of the functions:
	std::cout << "Maximum of l1 : " << l1.compute_maximum() << std::endl;
	std::cout << "Maximum of l2 : " << l2.compute_maximum() << std::endl;
	
	//here are the norms of landscapes:
	std::cout << "L^1 Norm of l1 : " << l1.compute_norm_of_landscape( 1. ) << std::endl;
	std::cout << "L^1 Norm of l2 : " << l2.compute_norm_of_landscape( 1. ) << std::endl;
	
	//here is the average of landscapes:
	Persistence_landscape_on_grid average;
	std::vector< Persistence_landscape_on_grid* > to_average;
	to_average.push_back( &l1 );
	to_average.push_back( &l2 );
	average.compute_average( to_average );
	std::cout << "average : " << average << std::endl;
	
	//here is the distance of landscapes:
	std::cout << "Distance : " << l1.distance( l2 ) << std::endl;
	
	//here is the scalar product of landscapes:
	std::cout << "Scalar product : " << l1.compute_scalar_product( l2 ) << std::endl;
	
	//here is how to create a file which is suitable for vizualization via gnuplot:
	average.plot( "average_landscape" );
		
		
	return 0;
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
	
	cerr <<  compute_discance_of_landscapes_on_grid(p,sum,1) << endl;
	cerr <<  compute_discance_of_landscapes_on_grid(p,sum,2) << endl;
	cerr <<  compute_discance_of_landscapes_on_grid(p,sum,-1)  << endl;
	*/
	
	/*
	Persistence_landscape_on_grid p( "file_with_diagram", 0,1,100 );
	Persistence_landscape_on_grid q( "file_with_diagram_1", 0,1,100 );	
	std::vector< Abs_Topological_data_with_averages* > to_average;
	to_average.push_back( &p );
	to_average.push_back( &q );
	Persistence_landscape_on_grid av;	
	av.compute_average( to_average );
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
	cerr <<  p.distance( &q , -1 ) << endl;
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

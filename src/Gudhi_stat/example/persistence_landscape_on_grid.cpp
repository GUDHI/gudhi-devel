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
	

	Persistence_landscape_on_grid p( "file_with_diagram", 0,1,10000 );
	Persistence_landscape_on_grid q( "file_with_diagram_1", 0,1,10000 );
	
	//std::vector< std::pair< double,double > > aa;
	//aa.push_back( std::make_pair( 0,1 ) );
	//Persistence_landscape_on_grid p( aa, 0,1,10 );
	//Persistence_landscape_on_grid q( aa, 0,1,10 );
	cerr <<  p.compute_scalar_product( &q ) << endl;

		
	return 0;
}

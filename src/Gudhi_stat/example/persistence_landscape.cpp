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



#include <gudhi/abstract_classes/Abs_Topological_data.h>
#include <gudhi/concretizations/Persistence_landscapes.h>



using namespace Gudhi;
using namespace Gudhi::Gudhi_stat;

#include <iostream>


int main( int argc , char** argv )
{
	
	if ( argc != 2 )
	{
		std::cerr << "To run this program, please provide a name of a file with persistence landscape \n";
		//return 1;
	}		
	Persistence_landscape p("../test/data/file_with_diagram");
	
	Persistence_landscape q;
	q.load_landscape_from_file( "file_with_landscape_from_file_with_diagram" );
	
	if ( p != q )
	{
		cout << "Not equal \n";
	}
	
	double integral = p.compute_integral_of_landscape();
	cout << "integral : " << integral <<endl;
	
	//compute integral for each level separatelly
	for ( size_t level = 0 ; level != p.size() ; ++level )
	{
		cout << p.compute_integral_of_landscape( level ) << endl;
	}
	
	//compute integral of p-th power of landscspe
	for ( size_t power = 0 ; power != 5 ; ++power ) 
	{
		cout << p.compute_integral_of_landscape( power ) << endl;
	}
	
	cout << "Value of level 1 at 0 : " <<  p.compute_value_at_a_given_point(1,0.0) << endl;
	cout << "Value of level 1  at 1 : " <<  p.compute_value_at_a_given_point(1,0.1) << endl;
	cout << "Value of level 1  at 2 : " <<  p.compute_value_at_a_given_point(1,0.2) << endl;
	cout << "Value of level 1  at 3 : " <<  p.compute_value_at_a_given_point(1,0.3) << endl;
	
		
	cout << "Value of level 2 at 0 : " <<  p.compute_value_at_a_given_point(2,0.0) << endl;
	cout << "Value of level 2  at 1 : " <<  p.compute_value_at_a_given_point(2,0.1) << endl;
	cout << "Value of level 2  at 2 : " <<  p.compute_value_at_a_given_point(2,0.2) << endl;
	cout << "Value of level 2  at 3 : " <<  p.compute_value_at_a_given_point(2,0.3) << endl;
	
		
	cout << "Value of level 3 at 0 : " <<  p.compute_value_at_a_given_point(3,0.0) << endl;
	cout << "Value of level 3  at 1 : " <<  p.compute_value_at_a_given_point(3,0.1) << endl;
	cout << "Value of level 3  at 2 : " <<  p.compute_value_at_a_given_point(3,0.2) << endl;
	cout << "Value of level 3  at 3 : " <<  p.compute_value_at_a_given_point(3,0.3) << endl;
	
	
	
	Persistence_landscape second;
	second.load_landscape_from_file("file_with_landscape_from_file_with_diagram_1" );
	
	Persistence_landscape sum = p + second;
	Persistence_landscape difference = p - second;
	Persistence_landscape multiply_by_scalar = 10*p;
	
	//sum.print_to_file("sum");
	//difference.print_to_file("difference");
	//multiply_by_scalar.print_to_file("multiply_by_scalar");
	
	Persistence_landscape template_sum;
	template_sum.load_landscape_from_file( "sum" );
	Persistence_landscape template_difference;
	template_difference.load_landscape_from_file( "difference" );
	Persistence_landscape template_multiply_by_scalar;
	template_multiply_by_scalar.load_landscape_from_file( "multiply_by_scalar" );
	
	if ( sum != template_sum )
	{
		cerr << "Problem with sums \n";
	}
	if ( difference != template_difference )
	{
		cerr << "Problem with differences \n";
	}
	if ( multiply_by_scalar != template_multiply_by_scalar )
	{
		cerr << "Problem with multiplication by scalar \n";
	}
	
	
	
	cout << "Maximum : " << p.compute_maximum() << endl;
	
	cout << "L^1 norm : " << p.compute_norm_of_landscape(1) << endl;
	cout << "L^2 norm : " << p.compute_norm_of_landscape(2) << endl;
	cout << "L^3 norm : " << p.compute_norm_of_landscape(3) << endl;	
	
	
	cout << "L^1 distance : " << compute_discance_of_landscapes(p,sum,1) << endl;
	cout << "L^2 distance : " << compute_discance_of_landscapes(p,sum,2) << endl;
	cout << "L^infty distance : " << compute_discance_of_landscapes(p,sum,-1) << endl;
	
	{
		Persistence_landscape p( "../test/data/file_with_diagram" );
		Persistence_landscape q( "../test/data/file_with_diagram_1" );
		std::vector< Abs_Topological_data_with_averages* > to_average;
		to_average.push_back( &p );
		to_average.push_back( &q );
		Persistence_landscape av;
		av.compute_average( to_average );
	
		Persistence_landscape template_average;
		template_average.load_landscape_from_file( "average" );
		if ( template_average != av )
		{
			cerr << "We have a problem with average \n";
		}
	}
	
	
	{
		Persistence_landscape p( "../test/data/file_with_diagram" );
		Persistence_landscape q( "../test/data/file_with_diagram_1" );
		cout << "L^1 distance : " <<  p.distance( &q ) << endl;
		cout << "L^2 distance : " <<  p.distance( &q , 2) << endl;
		cout << "L^infty distance : " <<  p.distance( &q , -1 ) << endl;
	}
	
	
	{
		Persistence_landscape p( "../test/data/file_with_diagram" );
		Persistence_landscape q( "../test/data/file_with_diagram_1" );
		cout << "Scalar product : " <<  p.compute_scalar_product( &q ) << endl;
		
	}
	
	return 0;
}

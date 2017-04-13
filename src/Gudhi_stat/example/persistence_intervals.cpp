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



#include <gudhi/reader_utils.h>
#include <gudhi/persistence_representations/Persistence_intervals.h>

#include <iostream>



using namespace Gudhi;
using namespace Gudhi::Gudhi_stat;

	
int main( int argc , char** argv )
{
	if ( argc != 2 )
	{
		std::cout << "To run this program, please provide the name of a file with persistence diagram \n";
		return 1;
	}
	
	Persistence_intervals p( argv[1] );
	std::pair<double,double> min_max_ = p.get_x_range();	
	std::cout << "Birth-death range : " <<  min_max_.first << " " << min_max_.second << std::endl;
	
	
	std::vector<double> dominant_ten_intervals_length = p.length_of_dominant_intervals(10);
	std::cout << "Lendth of ten dominant intervals : " <<std::endl;	
	for ( size_t i = 0 ; i != dominant_ten_intervals_length.size() ; ++i )
	{
		std::cout << dominant_ten_intervals_length[i] <<std::endl;
	}
	
	std::vector< std::pair<double,double> > ten_dominant_intervals = p.dominant_intervals( 10 );	
	std::cout << "Here are the dominant intervals : " <<std::endl;
	for ( size_t i = 0 ; i != ten_dominant_intervals.size() ; ++i )
	{
		std::cout << "( " << ten_dominant_intervals[i].first<< "," << ten_dominant_intervals[i].second <<std::endl;
	}
	
	std::vector< size_t > histogram = p.histogram_of_lengths( 10  );
	std::cout << "Here is the histogram of barcode's length : " <<std::endl;
	for ( size_t i = 0 ; i != histogram.size() ; ++i )
	{
		std::cout << histogram[i]  << " ";
	}		
	std::cout <<std::endl;

	
	std::vector< size_t > cumulative_histogram = p.cumulative_histogram_of_lengths( 10  );
	std::cout<< "Cumuative histogram : " <<std::endl;	
	for ( size_t i = 0 ; i != cumulative_histogram.size() ; ++i )
	{
		std::cout << cumulative_histogram[i] << " ";
	}		
	std::cout <<std::endl;

	std::vector< double  > char_funct_diag = p.characteristic_function_of_diagram( min_max_.first , min_max_.second );	
	std::cout << "Characteristic function of diagram : " <<std::endl;
	for ( size_t i = 0 ; i != char_funct_diag.size() ; ++i )
	{
		std::cout << char_funct_diag[i] << " ";
	}		
	std::cout <<std::endl;

	std::vector< double  > cumul_char_funct_diag = p.cumulative_characteristic_function_of_diagram( min_max_.first , min_max_.second );
	std::cout << "Cumulative characteristic function of diagram : " <<std::endl;
	for ( size_t i = 0 ; i != cumul_char_funct_diag.size() ; ++i )
	{		
		std::cout << cumul_char_funct_diag[i] << " ";
	}	
	std::cout <<std::endl;

	std::cout << "Persistence Betti numbers \n";
	std::vector< std::pair< double , size_t > > pbns = p.compute_persistent_betti_numbers();	
	for ( size_t i = 0 ; i != pbns.size() ; ++i )
	{
		std::cout <<  pbns[i].first << " " << pbns[i].second <<std::endl;	
	}
	
	return 0;
}








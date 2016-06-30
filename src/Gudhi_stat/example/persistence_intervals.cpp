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



#include <gudhi/reader_utils.h>
#include <gudhi/abstract_classes/Abs_Topological_data.h>
#include "gudhi/concretizations/persistence_intervals.h"

#include <iostream>



using namespace Gudhi;
using namespace Gudhi::Gudhi_stat;


double epsilon = 0.0000005;

using namespace std;

	
int main( int argc , char** argv )
{
	if ( argc != 2 )
	{
		cout << "To run this program, please provide the name of a file with persistence diagram \n";
		return 1;
	}
	
	Persistence_intervals p( argv[1] );
	std::pair<double,double> min_max_ = p.min_max();
	
	cout << "Birth-death range : " <<  min_max_.first << " " << min_max_.second << endl;
	
	
	std::vector<double> dominant_ten_intervals_length = p.length_of_dominant_intervals(10);
	cout << "Lendth of ten dominant intervals : " << endl;	
	for ( size_t i = 0 ; i != dominant_ten_intervals_length.size() ; ++i )
	{
		cout << dominant_ten_intervals_length[i] << endl;
	}
	std::vector< std::pair<double,double> > ten_dominant_intervals = p.dominant_intervals( 10 );
	
	cout << "Here are the dominant intervals : " << endl;
	for ( size_t i = 0 ; i != ten_dominant_intervals.size() ; ++i )
	{
		cout << "( " << ten_dominant_intervals[i].first<< "," << ten_dominant_intervals[i].second << endl;
	}
	
	std::vector< size_t > histogram = p.histograms_of_lengths( 10  );
	cout << "Here is the histogram of barcode's length : " << endl;
	for ( size_t i = 0 ; i != histogram.size() ; ++i )
	{
		cout << histogram[i]  << " ";
	}		
	cout << endl;

	
	std::vector< size_t > cumulative_histogram = p.cumulative_histograms_of_lengths( 10  );
	cout<< "Cumuative histogram : " << endl;	
	for ( size_t i = 0 ; i != cumulative_histogram.size() ; ++i )
	{
		cout << cumulative_histogram[i] << " ";
	}		
	cout << endl;

	std::vector< double  > char_funct_diag = p.characteristic_function_of_diagram( min_max_.first , min_max_.second );	
	cout << "Characteristic function of diagram : " << endl;
	for ( size_t i = 0 ; i != char_funct_diag.size() ; ++i )
	{
		cout << char_funct_diag[i] << " ";
	}		
	cout << endl;

	std::vector< double  > cumul_char_funct_diag = p.cumulative_characteristic_function_of_diagram( min_max_.first , min_max_.second );
	cout << "Cumulative characteristic function of diagram : " << endl;
	for ( size_t i = 0 ; i != cumul_char_funct_diag.size() ; ++i )
	{		
		cout << cumul_char_funct_diag[i] << " ";
	}	
	cout << endl;

	cout << "Persistence Betti numbers \n";
	std::vector< std::pair< double , size_t > > pbns = p.compute_persistent_betti_numbers();	
	for ( size_t i = 0 ; i != pbns.size() ; ++i )
	{
		cout <<  pbns[i].first << " " << pbns[i].second << endl;	
	}
	
	return 0;
}








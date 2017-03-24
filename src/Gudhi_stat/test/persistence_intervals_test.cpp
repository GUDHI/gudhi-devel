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
#include "gudhi/persistence_representations/Persistence_intervals.h"

#include <iostream>



using namespace Gudhi;
using namespace Gudhi::Gudhi_stat;


double epsilon = 0.0000005;


	


//cout << "Left most end of the interval : " << min_max_.first << std::endl;
//cout << "Right most end of the interval : " << min_max_.second << std::endl;
BOOST_AUTO_TEST_CASE(check_min_max_function) 
{
	Persistence_intervals p( "data/file_with_diagram" );
	std::pair<double,double> min_max_ = p.min_max();
	
	//cout << min_max_.first << " " << min_max_.second << std::endl;getchar();
	
	BOOST_CHECK( fabs ( min_max_.first - 0.0290362 ) <= epsilon );
	BOOST_CHECK( fabs ( min_max_.second - 0.994537 ) <= epsilon );
}





BOOST_AUTO_TEST_CASE(check_length_of_dominant_intervals) 
{
	Persistence_intervals p( "data/file_with_diagram" );
	std::vector<double> dominant_ten_intervals_length = p.length_of_dominant_intervals(10);
	std::vector< double > dominant_intervals_length;
	dominant_intervals_length.push_back(0.862625);
	dominant_intervals_length.push_back(0.800893);
	dominant_intervals_length.push_back(0.762061);
	dominant_intervals_length.push_back(0.756501);
	dominant_intervals_length.push_back(0.729367);
	dominant_intervals_length.push_back(0.718177);
	dominant_intervals_length.push_back(0.708395);
	dominant_intervals_length.push_back(0.702844);
	dominant_intervals_length.push_back(0.700468);
	dominant_intervals_length.push_back(0.622177);
	for ( size_t i = 0 ; i != dominant_ten_intervals_length.size() ; ++i )
	{
		BOOST_CHECK( fabs ( dominant_ten_intervals_length[i] - dominant_intervals_length[i] ) <= epsilon );
	}
}
BOOST_AUTO_TEST_CASE(check_dominant_intervals) 
{
	Persistence_intervals p( "data/file_with_diagram" );
	std::vector< std::pair<double,double> > ten_dominant_intervals = p.dominant_intervals( 10 );
	
	std::vector< std::pair<double,double> > templ;
	templ.push_back( std::pair<double,double>( 0.114718,0.977343) );
	templ.push_back( std::pair<double,double>( 0.133638,0.93453) );
	templ.push_back( std::pair<double,double>( 0.104599,0.866659) );
	templ.push_back( std::pair<double,double>( 0.149798,0.906299) );
	templ.push_back( std::pair<double,double>( 0.247352,0.976719) );
	templ.push_back( std::pair<double,double>( 0.192675,0.910852) );
	templ.push_back( std::pair<double,double>( 0.191836,0.900231) );
	templ.push_back( std::pair<double,double>( 0.284998,0.987842) );
	templ.push_back( std::pair<double,double>( 0.294069,0.994537) );
	templ.push_back( std::pair<double,double>( 0.267421,0.889597) );
		
	for ( size_t i = 0 ; i != ten_dominant_intervals.size() ; ++i )
	{
		BOOST_CHECK( fabs ( ten_dominant_intervals[i].first - templ[i].first ) <= epsilon );
		BOOST_CHECK( fabs ( ten_dominant_intervals[i].second - templ[i].second ) <= epsilon );
	}
}

BOOST_AUTO_TEST_CASE(check_histograms_of_lengths) 
{
	Persistence_intervals p( "data/file_with_diagram" );
	std::vector< size_t > histogram = p.histograms_of_lengths( 10  );
	std::vector< size_t > template_histogram;
	template_histogram.push_back(10); 
	template_histogram.push_back(5); 
	template_histogram.push_back(3); 
	template_histogram.push_back(4); 
	template_histogram.push_back(4); 
	template_histogram.push_back(3); 
	template_histogram.push_back(6);
	template_histogram.push_back(1); 
	template_histogram.push_back(7);
	template_histogram.push_back(1); 
	template_histogram.push_back(1);		
	for ( size_t i = 0 ; i != histogram.size() ; ++i )
	{
		BOOST_CHECK( fabs ( histogram[i] - template_histogram[i] ) <= epsilon );
	}		
}

BOOST_AUTO_TEST_CASE(check_cumulative_histograms_of_lengths) 
{
	Persistence_intervals p( "data/file_with_diagram" );
	std::vector< size_t > cumulative_histogram = p.cumulative_histograms_of_lengths( 10  );
	std::vector< size_t > template_cumulative_histogram;
	template_cumulative_histogram.push_back(10); 
	template_cumulative_histogram.push_back(15); 
	template_cumulative_histogram.push_back(18);
	template_cumulative_histogram.push_back(22); 
	template_cumulative_histogram.push_back(26); 
	template_cumulative_histogram.push_back(29); 
	template_cumulative_histogram.push_back(35); 
	template_cumulative_histogram.push_back(36); 
	template_cumulative_histogram.push_back(43); 
	template_cumulative_histogram.push_back(44); 
	template_cumulative_histogram.push_back(45);
	
	for ( size_t i = 0 ; i != cumulative_histogram.size() ; ++i )
	{
		BOOST_CHECK( fabs ( cumulative_histogram[i] - template_cumulative_histogram[i] ) <= epsilon );
	}		
}	
BOOST_AUTO_TEST_CASE(check_characteristic_function_of_diagram) 
{
	Persistence_intervals p( "data/file_with_diagram" );
	std::pair<double,double> min_max_ = p.min_max();
	std::vector< double  > char_funct_diag = p.characteristic_function_of_diagram( min_max_.first , min_max_.second );	
	std::vector< double > template_char_funct_diag;
	template_char_funct_diag.push_back(	0.370665	);
	template_char_funct_diag.push_back(	0.84058	);
	template_char_funct_diag.push_back(	1.24649	);
	template_char_funct_diag.push_back(	1.3664	);
	template_char_funct_diag.push_back(	1.34032	);
	template_char_funct_diag.push_back(	1.31904	);
	template_char_funct_diag.push_back(	1.14076	);
	template_char_funct_diag.push_back(	0.991259	);
	template_char_funct_diag.push_back(	0.800714	);
	template_char_funct_diag.push_back(	0.0676303	);

	for ( size_t i = 0 ; i != char_funct_diag.size() ; ++i )
	{
		//cout << char_funct_diag[i] << std::endl;
		if ( fabs ( char_funct_diag[i] - template_char_funct_diag[i] ) >= 0.0001  )
		{
			std::cout << "Boost test fail  check_characteristic_function_of_diagram : " << std::endl;
			std::cerr << char_funct_diag[i] << " " << template_char_funct_diag[i] << std::endl;
			std::cerr << fabs ( char_funct_diag[i] - template_char_funct_diag[i] ) << std::endl;
			std::cerr << 0.0001 << std::endl;
			getchar();
		}
		BOOST_CHECK( fabs ( char_funct_diag[i] - template_char_funct_diag[i] ) <= 0.0001 );
	}		
}	


BOOST_AUTO_TEST_CASE(check_cumulative_characteristic_function_of_diagram) 
{
	Persistence_intervals p( "data/file_with_diagram" );
	std::pair<double,double> min_max_ = p.min_max();
	std::vector< double  > cumul_char_funct_diag = p.cumulative_characteristic_function_of_diagram( min_max_.first , min_max_.second );
	std::vector< double > template_char_funct_diag_cumul;		
	
	template_char_funct_diag_cumul.push_back(0.370665);
	template_char_funct_diag_cumul.push_back(1.21125);
	template_char_funct_diag_cumul.push_back(2.45774);
	template_char_funct_diag_cumul.push_back(3.82414);
	template_char_funct_diag_cumul.push_back(5.16446);
	template_char_funct_diag_cumul.push_back(6.4835);
	template_char_funct_diag_cumul.push_back(7.62426);
	template_char_funct_diag_cumul.push_back(8.61552);
	template_char_funct_diag_cumul.push_back(9.41623);
	template_char_funct_diag_cumul.push_back(9.48386);

	for ( size_t i = 0 ; i != cumul_char_funct_diag.size() ; ++i )
	{		
		//cout << cumul_char_funct_diag[i] << std::endl;
		BOOST_CHECK( fabs (cumul_char_funct_diag[i] - template_char_funct_diag_cumul[i] ) <= 0.0001 );
	}	
}


BOOST_AUTO_TEST_CASE(check_compute_persistent_betti_numbers) 
{
	Persistence_intervals p( "data/file_with_diagram" );
	std::vector< std::pair< double , double > > pbns;
	pbns.push_back(std::pair<double,double>(	0.0290362,1	));
	pbns.push_back(std::pair<double,double>(	0.0307676,2	));
	pbns.push_back(std::pair<double,double>(	0.0366312,3	));
	pbns.push_back(std::pair<double,double>(	0.0544614,4	));
	pbns.push_back(std::pair<double,double>(	0.0920033,5	));
	pbns.push_back(std::pair<double,double>(	0.104599,6	));
	pbns.push_back(std::pair<double,double>(	0.114718,7	));
	pbns.push_back(std::pair<double,double>(	0.117379,8	));
	pbns.push_back(std::pair<double,double>(	0.123493,9	));
	pbns.push_back(std::pair<double,double>(	0.133638,10	));
	pbns.push_back(std::pair<double,double>(	0.137798,9	));
	pbns.push_back(std::pair<double,double>(	0.149798,10	));
	pbns.push_back(std::pair<double,double>(	0.155421,11	));
	pbns.push_back(std::pair<double,double>(	0.158443,12	));
	pbns.push_back(std::pair<double,double>(	0.176956,13	));
	pbns.push_back(std::pair<double,double>(	0.183234,12	));
	pbns.push_back(std::pair<double,double>(	0.191069,13	));
	pbns.push_back(std::pair<double,double>(	0.191333,14	));
	pbns.push_back(std::pair<double,double>(	0.191836,15	));
	pbns.push_back(std::pair<double,double>(	0.192675,16	));
	pbns.push_back(std::pair<double,double>(	0.208564,17	));
	pbns.push_back(std::pair<double,double>(	0.218425,18	));
	pbns.push_back(std::pair<double,double>(	0.219902,17	));
	pbns.push_back(std::pair<double,double>(	0.23233,16	));
	pbns.push_back(std::pair<double,double>(	0.234558,17	));
	pbns.push_back(std::pair<double,double>(	0.237166,16	));
	pbns.push_back(std::pair<double,double>(	0.247352,17	));
	pbns.push_back(std::pair<double,double>(	0.267421,18	));
	pbns.push_back(std::pair<double,double>(	0.268093,19	));
	pbns.push_back(std::pair<double,double>(	0.278734,18	));
	pbns.push_back(std::pair<double,double>(	0.284722,19	));
	pbns.push_back(std::pair<double,double>(	0.284998,20	));
	pbns.push_back(std::pair<double,double>(	0.294069,21	));
	pbns.push_back(std::pair<double,double>(	0.306293,22	));
	pbns.push_back(std::pair<double,double>(	0.322361,21	));
	pbns.push_back(std::pair<double,double>(	0.323152,22	));
	pbns.push_back(std::pair<double,double>(	0.371021,23	));
	pbns.push_back(std::pair<double,double>(	0.372395,24	));
	pbns.push_back(std::pair<double,double>(	0.387744,25	));
	pbns.push_back(std::pair<double,double>(	0.435537,26	));
	pbns.push_back(std::pair<double,double>(	0.462911,25	));
	pbns.push_back(std::pair<double,double>(	0.483569,26	));
	pbns.push_back(std::pair<double,double>(	0.489209,25	));
	pbns.push_back(std::pair<double,double>(	0.517115,24	));
	pbns.push_back(std::pair<double,double>(	0.522197,23	));
	pbns.push_back(std::pair<double,double>(	0.532665,22	));
	pbns.push_back(std::pair<double,double>(	0.545262,23	));
	pbns.push_back(std::pair<double,double>(	0.587227,22	));
	pbns.push_back(std::pair<double,double>(	0.593036,23	));
	pbns.push_back(std::pair<double,double>(	0.602647,24	));
	pbns.push_back(std::pair<double,double>(	0.605044,25	));
	pbns.push_back(std::pair<double,double>(	0.621962,24	));
	pbns.push_back(std::pair<double,double>(	0.629449,23	));
	pbns.push_back(std::pair<double,double>(	0.636719,22	));
	pbns.push_back(std::pair<double,double>(	0.64957,21	));
	pbns.push_back(std::pair<double,double>(	0.650781,22	));
	pbns.push_back(std::pair<double,double>(	0.654951,23	));
	pbns.push_back(std::pair<double,double>(	0.683489,24	));
	pbns.push_back(std::pair<double,double>(	0.687172,23	));
	pbns.push_back(std::pair<double,double>(	0.69703,22	));
	pbns.push_back(std::pair<double,double>(	0.701174,21	));
	pbns.push_back(std::pair<double,double>(	0.717623,22	));
	pbns.push_back(std::pair<double,double>(	0.722023,21	));
	pbns.push_back(std::pair<double,double>(	0.722298,20	));
	pbns.push_back(std::pair<double,double>(	0.725347,19	));
	pbns.push_back(std::pair<double,double>(	0.73071,18	));
	pbns.push_back(std::pair<double,double>(	0.758355,17	));
	pbns.push_back(std::pair<double,double>(	0.770913,18	));
	pbns.push_back(std::pair<double,double>(	0.790833,17	));
	pbns.push_back(std::pair<double,double>(	0.821211,16	));
	pbns.push_back(std::pair<double,double>(	0.849305,17	));
	pbns.push_back(std::pair<double,double>(	0.853669,16	));
	pbns.push_back(std::pair<double,double>(	0.866659,15	));
	pbns.push_back(std::pair<double,double>(	0.872896,16	));
	pbns.push_back(std::pair<double,double>(	0.889597,15	));
	pbns.push_back(std::pair<double,double>(	0.900231,14	));
	pbns.push_back(std::pair<double,double>(	0.903847,13	));
	pbns.push_back(std::pair<double,double>(	0.906299,12	));
	pbns.push_back(std::pair<double,double>(	0.910852,11	));
	pbns.push_back(std::pair<double,double>(	0.93453,10	));
	pbns.push_back(std::pair<double,double>(	0.944757,9	));
	pbns.push_back(std::pair<double,double>(	0.947812,8	));
	pbns.push_back(std::pair<double,double>(	0.959154,7	));
	pbns.push_back(std::pair<double,double>(	0.975654,6	));
	pbns.push_back(std::pair<double,double>(	0.976719,5	));
	pbns.push_back(std::pair<double,double>(	0.977343,4	));
	pbns.push_back(std::pair<double,double>(	0.980129,3	));
	pbns.push_back(std::pair<double,double>(	0.987842,2	));
	pbns.push_back(std::pair<double,double>(	0.990127,1	));
	pbns.push_back(std::pair<double,double>(	0.994537,0	));




	std::vector< std::pair< double , size_t > > pbns_new = p.compute_persistent_betti_numbers();	
	for ( size_t i = 0 ; i != pbns.size() ; ++i )
	{
		//cout << pbns_new[i].first << "," << pbns_new[i].second << std::endl;
		BOOST_CHECK( fabs ( pbns[i].first - pbns_new[i].first ) <= epsilon );		
		BOOST_CHECK( fabs ( pbns[i].second - pbns_new[i].second ) <= epsilon );
	}
}



BOOST_AUTO_TEST_CASE(check_k_n_n) 
{
	Persistence_intervals p( "data/file_with_diagram" );
	std::vector< double > knn = p.k_n_n( 5 );	
	std::vector< double > knn_template;
	knn_template.push_back(	1.04208	);
	knn_template.push_back(	1.00344	);
	knn_template.push_back(	0.979395	);
	knn_template.push_back(	0.890643	);
	knn_template.push_back(	0.874769	);
	knn_template.push_back(	0.845787	);
	knn_template.push_back(	0.819713	);
	knn_template.push_back(	0.803984	);
	knn_template.push_back(	0.799864	);
	knn_template.push_back(	0.786945	);
		
	for ( size_t i = 0 ; i != knn.size() ; ++i )
	{
		//cout << knn[i] << std::endl;
		BOOST_CHECK( fabs ( knn[i] - knn_template[i] ) <= 0.000005 );
	}
}





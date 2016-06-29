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



#include <gudhi/concretizations/Vector_distances_in_diagram.h>
#include <iostream>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "gudhi_stat"
#include <boost/test/unit_test.hpp>
#include <gudhi/reader_utils.h>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace Gudhi;
using namespace Gudhi::Gudhi_stat;

using namespace std;


double epsilon = 0.000005;




BOOST_AUTO_TEST_CASE(check_sortev_vector_distances_template) 
{	
	Vector_distances_in_diagram< euclidean_distance<double> > p( "data/file_with_diagram" , 100 );
	std::vector< double > sortev_vector_distances_template;
	sortev_vector_distances_template.push_back(0.609968);
	sortev_vector_distances_template.push_back(0.566317);
	sortev_vector_distances_template.push_back(0.538858);
	sortev_vector_distances_template.push_back(0.534927);
	sortev_vector_distances_template.push_back(0.534927);
	sortev_vector_distances_template.push_back(0.534927);
	sortev_vector_distances_template.push_back(0.534927);
	sortev_vector_distances_template.push_back(0.515741);
	sortev_vector_distances_template.push_back(0.515741);
	sortev_vector_distances_template.push_back(0.515741);
	sortev_vector_distances_template.push_back(0.515741);
	sortev_vector_distances_template.push_back(0.507828);
	sortev_vector_distances_template.push_back(0.507828);
	sortev_vector_distances_template.push_back(0.507828);
	sortev_vector_distances_template.push_back(0.500911);
	sortev_vector_distances_template.push_back(0.500911);
	sortev_vector_distances_template.push_back(0.500911);
	sortev_vector_distances_template.push_back(0.500911);
	sortev_vector_distances_template.push_back(0.500911);
	sortev_vector_distances_template.push_back(0.496986);
	sortev_vector_distances_template.push_back(0.496986);
	sortev_vector_distances_template.push_back(0.496986);
	sortev_vector_distances_template.push_back(0.496986);
	sortev_vector_distances_template.push_back(0.495306);
	sortev_vector_distances_template.push_back(0.495306);
	sortev_vector_distances_template.push_back(0.495306);
	sortev_vector_distances_template.push_back(0.495306);
	sortev_vector_distances_template.push_back(0.439945);
	sortev_vector_distances_template.push_back(0.439945);
	sortev_vector_distances_template.push_back(0.439945);
	sortev_vector_distances_template.push_back(0.439945);
	sortev_vector_distances_template.push_back(0.439945);
	sortev_vector_distances_template.push_back(0.439945);
	sortev_vector_distances_template.push_back(0.439945);
	sortev_vector_distances_template.push_back(0.424097);
	sortev_vector_distances_template.push_back(0.424097);
	sortev_vector_distances_template.push_back(0.424097);
	sortev_vector_distances_template.push_back(0.424097);
	sortev_vector_distances_template.push_back(0.424097);
	sortev_vector_distances_template.push_back(0.413891);
	sortev_vector_distances_template.push_back(0.413891);
	sortev_vector_distances_template.push_back(0.413891);
	sortev_vector_distances_template.push_back(0.413891);
	sortev_vector_distances_template.push_back(0.413891);
	sortev_vector_distances_template.push_back(0.413891);
	sortev_vector_distances_template.push_back(0.413891);
	sortev_vector_distances_template.push_back(0.410613);
	sortev_vector_distances_template.push_back(0.410613);
	sortev_vector_distances_template.push_back(0.410613);
	sortev_vector_distances_template.push_back(0.410613);
	sortev_vector_distances_template.push_back(0.410613);
	sortev_vector_distances_template.push_back(0.407853);
	sortev_vector_distances_template.push_back(0.407853);
	sortev_vector_distances_template.push_back(0.407853);
	sortev_vector_distances_template.push_back(0.407853);
	sortev_vector_distances_template.push_back(0.407853);
	sortev_vector_distances_template.push_back(0.407853);
	sortev_vector_distances_template.push_back(0.402306);
	sortev_vector_distances_template.push_back(0.402306);
	sortev_vector_distances_template.push_back(0.402306);
	sortev_vector_distances_template.push_back(0.402306);
	sortev_vector_distances_template.push_back(0.402306);
	sortev_vector_distances_template.push_back(0.401937);
	sortev_vector_distances_template.push_back(0.377605);
	sortev_vector_distances_template.push_back(0.377605);
	sortev_vector_distances_template.push_back(0.377605);
	sortev_vector_distances_template.push_back(0.377605);
	sortev_vector_distances_template.push_back(0.377605);
	sortev_vector_distances_template.push_back(0.363859);
	sortev_vector_distances_template.push_back(0.357765);
	sortev_vector_distances_template.push_back(0.357765);
	sortev_vector_distances_template.push_back(0.357765);
	sortev_vector_distances_template.push_back(0.357765);
	sortev_vector_distances_template.push_back(0.357765);
	sortev_vector_distances_template.push_back(0.357765);
	sortev_vector_distances_template.push_back(0.357765);
	sortev_vector_distances_template.push_back(0.357765);
	sortev_vector_distances_template.push_back(0.357765);
	sortev_vector_distances_template.push_back(0.357765);
	sortev_vector_distances_template.push_back(0.357765);
	sortev_vector_distances_template.push_back(0.357765);
	sortev_vector_distances_template.push_back(0.357765);
	sortev_vector_distances_template.push_back(0.357765);
	sortev_vector_distances_template.push_back(0.345124);
	sortev_vector_distances_template.push_back(0.345124);
	sortev_vector_distances_template.push_back(0.345124);
	sortev_vector_distances_template.push_back(0.345124);
	sortev_vector_distances_template.push_back(0.345124);
	sortev_vector_distances_template.push_back(0.345124);
	sortev_vector_distances_template.push_back(0.345124);
	sortev_vector_distances_template.push_back(0.345124);
	sortev_vector_distances_template.push_back(0.345124);
	sortev_vector_distances_template.push_back(0.345124);
	sortev_vector_distances_template.push_back(0.345124);
	sortev_vector_distances_template.push_back(0.345124);
	sortev_vector_distances_template.push_back(0.345124);
	sortev_vector_distances_template.push_back(0.345124);
	sortev_vector_distances_template.push_back(0.345124);
	sortev_vector_distances_template.push_back(0.34469);
	sortev_vector_distances_template.push_back(0.34469);  
	size_t proj_no = p.number_of_vectorize_functions();
	std::vector< double > aa = p.vectorize(proj_no);
	
	for ( size_t i = 0 ; i != aa.size() ; ++i )	
	{
		BOOST_CHECK( fabs ( sortev_vector_distances_template[i] - aa[i] ) <= epsilon );
	}
}	



	
	
	
	
	
	
	
BOOST_AUTO_TEST_CASE(check_projections_to_R) 
{
	Vector_distances_in_diagram< euclidean_distance<double> > p( "data/file_with_diagram" , 100 );
	std::vector< double > proj;
	proj.push_back(	0	);
	proj.push_back(	0.6099679993	);
	proj.push_back(	1.176284775	);
	proj.push_back(	1.715142954	);
	proj.push_back(	2.25006986	);
	proj.push_back(	2.784996767	);
	proj.push_back(	3.319923673	);
	proj.push_back(	3.854850579	);
	proj.push_back(	4.370591225	);
	proj.push_back(	4.886331872	);
	proj.push_back(	5.402072518	);
	proj.push_back(	5.917813164	);
	proj.push_back(	6.425641089	);
	proj.push_back(	6.933469014	);
	proj.push_back(	7.441296939	);
	proj.push_back(	7.942207817	);
	proj.push_back(	8.443118694	);
	proj.push_back(	8.944029572	);
	proj.push_back(	9.44494045	);
	proj.push_back(	9.945851328	);
	proj.push_back(	10.44283706	);
	proj.push_back(	10.93982279	);
	proj.push_back(	11.43680853	);
	proj.push_back(	11.93379426	);
	proj.push_back(	12.42910009	);
	proj.push_back(	12.92440592	);
	proj.push_back(	13.41971176	);
	proj.push_back(	13.91501759	);
	proj.push_back(	14.35496286	);
	proj.push_back(	14.79490814	);
	proj.push_back(	15.23485341	);
	proj.push_back(	15.67479869	);
	proj.push_back(	16.11474396	);
	proj.push_back(	16.55468923	);
	proj.push_back(	16.99463451	);
	proj.push_back(	17.41873131	);
	proj.push_back(	17.84282811	);
	proj.push_back(	18.26692491	);
	proj.push_back(	18.69102171	);
	proj.push_back(	19.11511851	);
	proj.push_back(	19.52900989	);
	proj.push_back(	19.94290127	);
	proj.push_back(	20.35679265	);
	proj.push_back(	20.77068403	);
	proj.push_back(	21.18457541	);
	proj.push_back(	21.59846679	);
	proj.push_back(	22.01235817	);
	proj.push_back(	22.42297162	);
	proj.push_back(	22.83358506	);
	proj.push_back(	23.2441985	);
	proj.push_back(	23.65481194	);
	proj.push_back(	24.06542538	);
	proj.push_back(	24.47327795	);
	proj.push_back(	24.88113051	);
	proj.push_back(	25.28898308	);
	proj.push_back(	25.69683564	);
	proj.push_back(	26.1046882	);
	proj.push_back(	26.51254077	);
	proj.push_back(	26.91484703	);
	proj.push_back(	27.3171533	);
	proj.push_back(	27.71945957	);
	proj.push_back(	28.12176583	);
	proj.push_back(	28.5240721	);
	proj.push_back(	28.92600955	);
	proj.push_back(	29.30361425	);
	proj.push_back(	29.68121895	);
	proj.push_back(	30.05882365	);
	proj.push_back(	30.43642834	);
	proj.push_back(	30.81403304	);
	proj.push_back(	31.17789162	);
	proj.push_back(	31.53565711	);
	proj.push_back(	31.89342259	);
	proj.push_back(	32.25118808	);
	proj.push_back(	32.60895356	);
	proj.push_back(	32.96671905	);
	proj.push_back(	33.32448454	);
	proj.push_back(	33.68225002	);
	proj.push_back(	34.04001551	);
	proj.push_back(	34.397781	);
	proj.push_back(	34.75554648	);
	proj.push_back(	35.11331197	);
	proj.push_back(	35.47107745	);
	proj.push_back(	35.82884294	);
	proj.push_back(	36.18660843	);
	proj.push_back(	36.53173245	);
	proj.push_back(	36.87685646	);
	proj.push_back(	37.22198048	);
	proj.push_back(	37.5671045	);
	proj.push_back(	37.91222852	);
	proj.push_back(	38.25735254	);
	proj.push_back(	38.60247656	);
	proj.push_back(	38.94760057	);
	proj.push_back(	39.29272459	);
	proj.push_back(	39.63784861	);
	proj.push_back(	39.98297263	);
	proj.push_back(	40.32809665	);
	proj.push_back(	40.67322067	);
	proj.push_back(	41.01834468	);
	proj.push_back(	41.3634687	);
	proj.push_back(	41.70815824	);



	for ( size_t proj_no = 0 ; proj_no != p.number_of_projections_to_R() ; ++proj_no )	 
	{		
		//cout << std::setprecision(10)  << 	p.project_to_R(proj_no) << endl;
		BOOST_CHECK( fabs ( p.project_to_R(proj_no) - proj[proj_no] ) <= epsilon );
	}
}





BOOST_AUTO_TEST_CASE(check_distance_computations) 
{
	Vector_distances_in_diagram< euclidean_distance<double> > p( "data/file_with_diagram" , 100 );
	Vector_distances_in_diagram< euclidean_distance<double> > p_prime( "data/file_with_diagram" , 10 );	
	std::vector< std::pair<double,double> > intervals(10);
	intervals[0] = std::pair<double,double>( 1,2 );
	intervals[1] = std::pair<double,double>( 3,4 );
	intervals[2] = std::pair<double,double>( 5,6 );
	intervals[3] = std::pair<double,double>( 7,8 );
	intervals[4] = std::pair<double,double>( 9,10 );
	intervals[5] = std::pair<double,double>( 11,12 );
	intervals[6] = std::pair<double,double>( 13,14 );
	intervals[7] = std::pair<double,double>( 15,16 );
	intervals[8] = std::pair<double,double>( 17,18 );
	intervals[9] = std::pair<double,double>( 19,20 );		
	Vector_distances_in_diagram< euclidean_distance<double> > p_bis( intervals , 10 );		
	BOOST_CHECK( fabs ( p_prime.distance( (Abs_Topological_data_with_distances*)(&p_bis) , 1 ) - 1.669 ) <= epsilon );		
}



BOOST_AUTO_TEST_CASE(check_compute_average) 
{
	Vector_distances_in_diagram< euclidean_distance<double> > p( "data/file_with_diagram" , 100 );
	//compute average 
	std::vector< std::pair<double,double> > i1(3);
	i1[0] = std::pair<double,double>( 1,2 );
	i1[1] = std::pair<double,double>( 3,4 );
	i1[2] = std::pair<double,double>( 5,6 );
	
	std::vector< std::pair<double,double> > i2(3);
	i2[0] = std::pair<double,double>( 2,3 );
	i2[1] = std::pair<double,double>( 4,5);
	i2[2] = std::pair<double,double>( 6,7 );
	
	Vector_distances_in_diagram< euclidean_distance<double> > A( i1 , 3 );
	Vector_distances_in_diagram< euclidean_distance<double> > B( i1 , 3 );
	
	std::vector< Abs_Topological_data_with_averages* > to_average;
	to_average.push_back( (Abs_Topological_data_with_averages*)(&A) );
	to_average.push_back( (Abs_Topological_data_with_averages*)(&B) );
	
	Vector_distances_in_diagram< euclidean_distance<double> > average; 
	average.compute_average( to_average );
	
	
	std::vector< double > vect = average.output_for_visualization();
	
}



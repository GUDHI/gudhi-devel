/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016  INRIA (France)
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
#define BOOST_TEST_MODULE "Persistence_landscapes_test"
#include <boost/test/unit_test.hpp>
#include <gudhi/reader_utils.h>
#include <gudhi/Persistence_landscape.h>
#include <gudhi/Unitary_tests_utils.h>

#include <iostream>
#include <limits>

using namespace Gudhi;
using namespace Gudhi::Persistence_representations;

double epsilon = 0.0005;

BOOST_AUTO_TEST_CASE(check_construction_of_landscape) {	
  std::vector<std::pair<double, double> > diag =
  read_persistence_intervals_in_one_dimension_from_file("data/file_with_diagram");  
  Persistence_landscape p(diag);
  Persistence_landscape q;  
  q.load_landscape_from_file("data/file_with_landscape_from_file_with_diagram");   
  BOOST_CHECK(p == q);
}

BOOST_AUTO_TEST_CASE(check_construction_of_landscape_form_gudhi_style_file) {
  Persistence_landscape p("data/persistence_file_with_four_entries_per_line", 1);
  // p.print_to_file("persistence_file_with_four_entries_per_line_landscape");
  Persistence_landscape q;
  q.load_landscape_from_file("data/persistence_file_with_four_entries_per_line_landscape");  
  BOOST_CHECK(p == q);
}



BOOST_AUTO_TEST_CASE(check_computations_of_integrals) {
  std::vector<std::pair<double, double> > diag =
      read_persistence_intervals_in_one_dimension_from_file("data/file_with_diagram");
  Persistence_landscape p(diag);
  //double integral = p.compute_integral_of_landscape();
  //BOOST_CHECK(fabs(integral - 2.34992) <= 0.00001);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.compute_integral_of_landscape(), 2.34992, epsilon);
}

BOOST_AUTO_TEST_CASE(check_computations_of_integrals_for_each_level_separatelly) {
  std::vector<std::pair<double, double> > diag =
      read_persistence_intervals_in_one_dimension_from_file("data/file_with_diagram");
  Persistence_landscape p(diag);

  std::vector<double> integrals_for_different_levels;
  integrals_for_different_levels.push_back(0.216432);
  integrals_for_different_levels.push_back(0.204763);
  integrals_for_different_levels.push_back(0.188793);
  integrals_for_different_levels.push_back(0.178856);
  integrals_for_different_levels.push_back(0.163142);
  integrals_for_different_levels.push_back(0.155015);
  integrals_for_different_levels.push_back(0.143046);
  integrals_for_different_levels.push_back(0.133765);
  integrals_for_different_levels.push_back(0.123531);
  integrals_for_different_levels.push_back(0.117393);
  integrals_for_different_levels.push_back(0.111269);
  integrals_for_different_levels.push_back(0.104283);
  integrals_for_different_levels.push_back(0.0941308);
  integrals_for_different_levels.push_back(0.0811208);
  integrals_for_different_levels.push_back(0.0679001);
  integrals_for_different_levels.push_back(0.0580801);
  integrals_for_different_levels.push_back(0.0489647);
  integrals_for_different_levels.push_back(0.0407936);
  integrals_for_different_levels.push_back(0.0342599);
  integrals_for_different_levels.push_back(0.02896);
  integrals_for_different_levels.push_back(0.0239881);
  integrals_for_different_levels.push_back(0.0171792);
  integrals_for_different_levels.push_back(0.0071511);
  integrals_for_different_levels.push_back(0.00462067);
  integrals_for_different_levels.push_back(0.00229033);
  integrals_for_different_levels.push_back(0.000195296);

  for (size_t level = 0; level != p.size(); ++level) {
    //double integral = p.compute_integral_of_a_level_of_a_landscape(level);
    GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.compute_integral_of_a_level_of_a_landscape(level),
                                    integrals_for_different_levels[level], epsilon);
    //BOOST_CHECK(fabs(integral - integrals_for_different_levels[level]) <= 0.00001);
  }
}

BOOST_AUTO_TEST_CASE(check_computations_of_integrals_of_powers_of_landscape) {
  std::vector<std::pair<double, double> > diag =
      read_persistence_intervals_in_one_dimension_from_file("data/file_with_diagram");
  Persistence_landscape p(diag);

  std::vector<double> integrals_for_different_powers;
  integrals_for_different_powers.push_back(17.1692);
  integrals_for_different_powers.push_back(2.34992);
  integrals_for_different_powers.push_back(0.49857);
  integrals_for_different_powers.push_back(0.126405);
  integrals_for_different_powers.push_back(0.0355235);

  for (size_t power = 0; power != 5; ++power) {
    GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.compute_integral_of_landscape((double)power),
                                    integrals_for_different_powers[power], epsilon);
  }
}

BOOST_AUTO_TEST_CASE(check_computations_of_values_on_different_points) {
  std::vector<std::pair<double, double> > diag =
      read_persistence_intervals_in_one_dimension_from_file("data/file_with_diagram");
  Persistence_landscape p(diag);

  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.compute_value_at_a_given_point(1, 0.0), 0.       , epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.compute_value_at_a_given_point(1, 0.1), 0.0692324, epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.compute_value_at_a_given_point(1, 0.2), 0.163369 , epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.compute_value_at_a_given_point(1, 0.3), 0.217115 , epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.compute_value_at_a_given_point(2, 0.0), 0.       , epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.compute_value_at_a_given_point(2, 0.1), 0.0633688, epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.compute_value_at_a_given_point(2, 0.2), 0.122361 , epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.compute_value_at_a_given_point(2, 0.3), 0.195401 , epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.compute_value_at_a_given_point(3, 0.0), 0.       , epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.compute_value_at_a_given_point(3, 0.1), 0.0455386, epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.compute_value_at_a_given_point(3, 0.2), 0.0954012, epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.compute_value_at_a_given_point(3, 0.3), 0.185282 , epsilon);
}

BOOST_AUTO_TEST_CASE(check_computations_sum_differences_and_multiplications) {
  std::vector<std::pair<double, double> > diag =
      read_persistence_intervals_in_one_dimension_from_file("data/file_with_diagram");
  Persistence_landscape p(diag);
  Persistence_landscape second;
  second.load_landscape_from_file("data/file_with_landscape_from_file_with_diagram_1");

  Persistence_landscape sum = p + second;
  Persistence_landscape difference = p - second;
  Persistence_landscape multiply_by_scalar = 10 * p;

  Persistence_landscape template_sum;
  template_sum.load_landscape_from_file("data/sum");

  Persistence_landscape template_difference;
  template_difference.load_landscape_from_file("data/difference");

  Persistence_landscape template_multiply_by_scalar;
  template_multiply_by_scalar.load_landscape_from_file("data/multiply_by_scalar");

  BOOST_CHECK(sum == template_sum);
  BOOST_CHECK(difference == template_difference);
  BOOST_CHECK(multiply_by_scalar == template_multiply_by_scalar);
}

BOOST_AUTO_TEST_CASE(check_computations_of_maxima_and_norms) {
  std::vector<std::pair<double, double> > diag =
      read_persistence_intervals_in_one_dimension_from_file("data/file_with_diagram");
  Persistence_landscape p(diag);
  Persistence_landscape second;
  second.load_landscape_from_file("data/file_with_landscape_from_file_with_diagram_1");
  Persistence_landscape sum = p + second;

  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.compute_maximum()           , 0.431313, epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.compute_norm_of_landscape(1), 2.34992 , epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.compute_norm_of_landscape(2), 0.706095, epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.compute_norm_of_landscape(3), 0.501867, epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(compute_distance_of_landscapes(p, sum, 1), 27.9323, epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(compute_distance_of_landscapes(p, sum, 2), 2.35199, epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(compute_distance_of_landscapes(p, sum, std::numeric_limits<double>::max()),
                                  0.464478, epsilon);
}

BOOST_AUTO_TEST_CASE(check_default_parameters_of_distances) {
  std::vector<std::pair<double, double> > diag =
      read_persistence_intervals_in_one_dimension_from_file("data/file_with_diagram");
  Persistence_landscape p(diag);

  std::vector<std::pair<double, double> > diag1 =
      read_persistence_intervals_in_one_dimension_from_file("data/file_with_diagram_1");
  Persistence_landscape q(diag1);

  double dist_numeric_limit_max = p.distance(q, std::numeric_limits<double>::max());
  double dist_infinity = p.distance(q, std::numeric_limits<double>::infinity());

  GUDHI_TEST_FLOAT_EQUALITY_CHECK(dist_numeric_limit_max, dist_infinity);
}

BOOST_AUTO_TEST_CASE(check_computations_of_averages) {
  std::vector<std::pair<double, double> > diag =
      read_persistence_intervals_in_one_dimension_from_file("data/file_with_diagram");
  Persistence_landscape p(diag);
  std::vector<std::pair<double, double> > diag2 =
      read_persistence_intervals_in_one_dimension_from_file("data/file_with_diagram_1");
  Persistence_landscape q(diag2);
  Persistence_landscape av;
  av.compute_average({&p, &q});

  Persistence_landscape template_average;
  template_average.load_landscape_from_file("data/average");
  BOOST_CHECK(template_average == av);
}

BOOST_AUTO_TEST_CASE(check_computations_of_distances) {
  std::vector<std::pair<double, double> > diag =
      read_persistence_intervals_in_one_dimension_from_file("data/file_with_diagram");
  Persistence_landscape p(diag);
  std::vector<std::pair<double, double> > diag2 =
      read_persistence_intervals_in_one_dimension_from_file("data/file_with_diagram_1");
  Persistence_landscape q(diag2);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.distance(q), 25.5824, epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.distance(q, 2), 2.1264, epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.distance(q, std::numeric_limits<double>::max()), 0.359068, epsilon);
}

BOOST_AUTO_TEST_CASE(check_computations_of_scalar_product) {
  std::vector<std::pair<double, double> > diag =
      read_persistence_intervals_in_one_dimension_from_file("data/file_with_diagram");
  Persistence_landscape p(diag);
  std::vector<std::pair<double, double> > diag2 =
      read_persistence_intervals_in_one_dimension_from_file("data/file_with_diagram_1");
  Persistence_landscape q(diag2);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.compute_scalar_product(q), 0.754498, epsilon);
}

























// Below I am storing the code used to generate tests for that functionality.
/*
if ( argc != 2 )
        {
                std::cerr << "To run this program, please provide a name of a file with persistence landscape \n";
                //return 1;
        }
        Persistence_landscape p("data/file_with_diagram");

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


        cout << "L^1 distance : " << compute_distance_of_landscapes(p,sum,1) << endl;
        cout << "L^2 distance : " << compute_distance_of_landscapes(p,sum,2) << endl;
        cout << "L^infty distance : " << compute_distance_of_landscapes(p,sum,std::numeric_limits<double>::max() ) <<
endl;

        {
                Persistence_landscape p( "data/file_with_diagram" );
                Persistence_landscape q( "data/file_with_diagram_1" );
                Persistence_landscape av;
                av.compute_average( {&p,&q} );

                Persistence_landscape template_average;
                template_average.load_landscape_from_file( "average" );
                if ( template_average != av )
                {
                        cerr << "We have a problem with average \n";
                }
        }


        {
                Persistence_landscape p( "data/file_with_diagram" );
                Persistence_landscape q( "data/file_with_diagram_1" );
                cout << "L^1 distance : " <<  p.distance( &q ) << endl;
                cout << "L^2 distance : " <<  p.distance( &q , 2) << endl;
                cout << "L^infty distance : " <<  p.distance( &q , std::numeric_limits<double>::max() ) << endl;
        }


        {
                Persistence_landscape p( "data/file_with_diagram" );
                Persistence_landscape q( "data/file_with_diagram_1" );
                cout << "Scalar product : " <<  p.compute_scalar_product( &q ) << endl;
        }
*/

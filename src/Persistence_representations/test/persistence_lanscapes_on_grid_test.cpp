/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "Persistence_landscapes_on_grid_test"
#include <boost/test/unit_test.hpp>
#include <gudhi/reader_utils.h>
#include <gudhi/Persistence_landscape_on_grid.h>
#include <gudhi/Unitary_tests_utils.h>

#include <iostream>

using namespace Gudhi;
using namespace Gudhi::Persistence_representations;

double epsilon = 0.0005;

BOOST_AUTO_TEST_CASE(check_construction_of_landscape) {
  Persistence_landscape_on_grid l("data/file_with_diagram_1", 100, std::numeric_limits<unsigned short>::max());
  l.print_to_file("landscape_from_file_with_diagram_1");

  Persistence_landscape_on_grid g;
  g.load_landscape_from_file("landscape_from_file_with_diagram_1");

  BOOST_CHECK(l == g);
}

BOOST_AUTO_TEST_CASE(check_construction_of_landscape_using_only_ten_levels) {
  // TODO
  unsigned number = 10;
  Persistence_landscape_on_grid l("data/file_with_diagram_1", 100, number);
  Persistence_landscape_on_grid g("data/file_with_diagram_1", 100, std::numeric_limits<unsigned short>::max());
  // cut all the elements of order > 10 in g.

  for (size_t level = 0; level != number; ++level) {
    std::vector<double> v1 = l.vectorize(level);
    std::vector<double> v2 = g.vectorize(level);
    BOOST_CHECK(v1.size() == v2.size());
    for (size_t i = 0; i != v1.size(); ++i) {
      GUDHI_TEST_FLOAT_EQUALITY_CHECK(v1[i], v2[i]);
    }
  }
}

BOOST_AUTO_TEST_CASE(check_computations_of_integrals) {
  Persistence_landscape_on_grid p("data/file_with_diagram_1", 100, std::numeric_limits<unsigned short>::max());
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.compute_integral_of_landscape(), 27.343, epsilon);
}

BOOST_AUTO_TEST_CASE(check_computations_of_integrals_for_each_level_separatelly) {
  Persistence_landscape_on_grid p("data/file_with_diagram_1", 100, std::numeric_limits<unsigned short>::max());

  std::vector<double> integrals_fir_different_levels;
  // integrals_fir_different_levels.push_back();
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

  for (size_t level = 0; level != integrals_fir_different_levels.size(); ++level) {
    GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.compute_integral_of_landscape(level), integrals_fir_different_levels[level],
                                    epsilon);
  }
}

BOOST_AUTO_TEST_CASE(check_computations_of_integrals_of_powers_of_landscape) {
  Persistence_landscape_on_grid p("data/file_with_diagram_1", 100, std::numeric_limits<unsigned short>::max());

  std::vector<double> integrals_fir_different_powers;
  integrals_fir_different_powers.push_back(0.241168);
  integrals_fir_different_powers.push_back(0.239276);
  integrals_fir_different_powers.push_back(0.237882);
  integrals_fir_different_powers.push_back(0.235193);
  integrals_fir_different_powers.push_back(0.23011);

  for (size_t power = 0; power != 5; ++power) {
    GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.compute_integral_of_landscape(power), integrals_fir_different_powers[power],
                                    epsilon);
  }
}

BOOST_AUTO_TEST_CASE(check_computations_of_values_on_different_points) {
  Persistence_landscape_on_grid p("data/file_with_diagram_1", 100, std::numeric_limits<unsigned short>::max());

  std::vector<double> results_level_0;
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

  std::vector<double> results_level_10;
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
  for (size_t i = 0; i != 10; ++i) {
    GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.compute_value_at_a_given_point(0, x), results_level_0[i], epsilon);
    GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.compute_value_at_a_given_point(10, x), results_level_10[i], epsilon);
    x += dx;
  }
}

BOOST_AUTO_TEST_CASE(check_computations_sum_differences_and_multiplications) {
  Persistence_landscape_on_grid p("data/file_with_diagram_1", 100, std::numeric_limits<unsigned short>::max());
  Persistence_landscape_on_grid second("data/file_with_diagram_1", 100, std::numeric_limits<unsigned short>::max());

  Persistence_landscape_on_grid sum = p + second;
  Persistence_landscape_on_grid difference = p - second;
  Persistence_landscape_on_grid multiply_by_scalar = 10 * p;
  ;

  Persistence_landscape_on_grid template_sum;
  template_sum.load_landscape_from_file("data/sum_on_grid_test");

  Persistence_landscape_on_grid template_difference;
  template_difference.load_landscape_from_file("data/difference_on_grid_test");

  Persistence_landscape_on_grid template_multiply_by_scalar;
  template_multiply_by_scalar.load_landscape_from_file("data/multiply_by_scalar_on_grid_test");

  BOOST_CHECK(sum == template_sum);
  BOOST_CHECK(difference == template_difference);
  BOOST_CHECK(multiply_by_scalar == template_multiply_by_scalar);
}

BOOST_AUTO_TEST_CASE(check_computations_of_maxima_and_norms) {
  Persistence_landscape_on_grid p("data/file_with_diagram_1", 0., 1., 100);
  Persistence_landscape_on_grid second("data/file_with_diagram_2", 0., 1., 100);
  Persistence_landscape_on_grid sum = p + second;

  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.compute_maximum(), 0.46, epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.compute_norm_of_landscape(1), 27.3373, epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.compute_norm_of_landscape(2), 1.84143, epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.compute_norm_of_landscape(3), 0.927067, epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(compute_distance_of_landscapes_on_grid(p, sum, 1), 16.8519, epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(compute_distance_of_landscapes_on_grid(p, sum, 2), 1.44542, epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(compute_distance_of_landscapes_on_grid(p, sum, std::numeric_limits<double>::max()),
                                  0.45, epsilon);
}

BOOST_AUTO_TEST_CASE(check_default_parameters_of_distances) {
  std::vector<std::pair<double, double> > diag = read_persistence_intervals_in_dimension("data/file_with_diagram");
  Persistence_landscape_on_grid p(diag, 0., 1., 100);

  std::vector<std::pair<double, double> > diag1 = read_persistence_intervals_in_dimension("data/file_with_diagram_1");
  Persistence_landscape_on_grid q(diag1, 0., 1., 100);

  double dist_numeric_limit_max = p.distance(q, std::numeric_limits<double>::max());
  double dist_infinity = p.distance(q, std::numeric_limits<double>::infinity());

  GUDHI_TEST_FLOAT_EQUALITY_CHECK(dist_numeric_limit_max, dist_infinity);
}

BOOST_AUTO_TEST_CASE(check_computations_of_averages) {
  Persistence_landscape_on_grid p("data/file_with_diagram", 0., 1., 100);
  Persistence_landscape_on_grid q("data/file_with_diagram_1", 0., 1., 100);
  Persistence_landscape_on_grid av;
  av.compute_average({&p, &q});

  Persistence_landscape_on_grid template_average;
  template_average.load_landscape_from_file("data/average_on_a_grid");
  BOOST_CHECK(template_average == av);
}

BOOST_AUTO_TEST_CASE(check_computations_of_distances) {
  Persistence_landscape_on_grid p("data/file_with_diagram", 0., 1., 10000);
  Persistence_landscape_on_grid q("data/file_with_diagram_1", 0., 1., 10000);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.distance(q), 25.5779, epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.distance(q, 2), 2.04891, epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.distance(q, std::numeric_limits<double>::max()), 0.359, epsilon);
}

BOOST_AUTO_TEST_CASE(check_computations_of_scalar_product) {
  Persistence_landscape_on_grid p("data/file_with_diagram", 0., 1., 10000);
  Persistence_landscape_on_grid q("data/file_with_diagram_1", 0., 1., 10000);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.compute_scalar_product(q), 0.754367, epsilon);
}

// Below I am storing the code used to generate tests for that functionality.
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

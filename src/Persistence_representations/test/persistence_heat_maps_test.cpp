/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <iostream>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "Persistence_heat_maps_test"
#include <boost/test/unit_test.hpp>

#include <gudhi/reader_utils.h>
#include <gudhi/Persistence_heat_maps.h>
#include <gudhi/Unitary_tests_utils.h>

using namespace Gudhi;
using namespace Gudhi::Persistence_representations;

double epsilon = 0.0005;

BOOST_AUTO_TEST_CASE(check_construction_of_heat_maps)
{
  std::vector<std::vector<double> > filter = create_Gaussian_filter(100, 1);
  Persistence_heat_maps<constant_scaling_function> p("data/file_with_diagram", filter, false, 1000, 0, 1);
  p.print_to_file("data/persistence_heat_map_from_file_with_diagram");

  Persistence_heat_maps<constant_scaling_function> q;
  q.load_from_file("data/persistence_heat_map_from_file_with_diagram");

  BOOST_CHECK(p == q);
}

BOOST_AUTO_TEST_CASE(check_averages_of_heat_maps)
{
  std::vector<std::vector<double> > filter = create_Gaussian_filter(30, 1);
  Persistence_heat_maps<constant_scaling_function> p("data/file_with_diagram", filter, false, 1000, 0, 10);
  Persistence_heat_maps<constant_scaling_function> q("data/file_with_diagram_1", filter, false, 1000, 0, 10);
  Persistence_heat_maps<constant_scaling_function> r("data/file_with_diagram_2", filter, false, 1000, 0, 10);

  Persistence_heat_maps<constant_scaling_function> av;
  av.compute_average({&p, &q, &r});

  Persistence_heat_maps<constant_scaling_function> template_average;
  template_average.load_from_file("data/template_average_of_heat_maps");

  BOOST_CHECK(av == template_average);
}

BOOST_AUTO_TEST_CASE(check_median_of_heat_maps)
{
  std::vector<std::vector<double> > filter = create_Gaussian_filter(30, 1);
  Persistence_heat_maps<constant_scaling_function> p("data/file_with_diagram", filter, false, 1000, 0, 1);
  Persistence_heat_maps<constant_scaling_function> q("data/file_with_diagram_1", filter, false, 1000, 0, 1);
  Persistence_heat_maps<constant_scaling_function> r("data/file_with_diagram_2", filter, false, 1000, 0, 1);

  std::vector<Persistence_heat_maps<constant_scaling_function>*> to_compute_median;
  to_compute_median.push_back(&p);
  to_compute_median.push_back(&q);
  to_compute_median.push_back(&r);
  Persistence_heat_maps<constant_scaling_function> median;
  median.compute_median(to_compute_median);

  Persistence_heat_maps<constant_scaling_function> template_median;
  template_median.load_from_file("data/template_median_of_heat_maps");

  BOOST_CHECK(median == template_median);
}

BOOST_AUTO_TEST_CASE(check_compute_percentage_of_active_of_heat_maps)
{
  std::vector<std::vector<double> > filter = create_Gaussian_filter(30, 1);
  Persistence_heat_maps<constant_scaling_function> p("data/file_with_diagram", filter, false, 1000, 0, 1);
  Persistence_heat_maps<constant_scaling_function> q("data/file_with_diagram_1", filter, false, 1000, 0, 1);
  Persistence_heat_maps<constant_scaling_function> r("data/file_with_diagram_2", filter, false, 1000, 0, 1);

  std::vector<Persistence_heat_maps<constant_scaling_function>*> to_compute_percentage_of_active;
  to_compute_percentage_of_active.push_back(&p);
  to_compute_percentage_of_active.push_back(&q);
  to_compute_percentage_of_active.push_back(&r);
  Persistence_heat_maps<constant_scaling_function> percentage_of_active;
  percentage_of_active.compute_percentage_of_active(to_compute_percentage_of_active, 0);

  Persistence_heat_maps<constant_scaling_function> template_percentage_of_active;
  template_percentage_of_active.load_from_file("data/template_percentage_of_active_of_heat_maps");

  BOOST_CHECK(percentage_of_active == template_percentage_of_active);
}

BOOST_AUTO_TEST_CASE(check_vectorize_for_heat_maps)
{
  std::vector<std::vector<double> > filter = create_Gaussian_filter(30, 1);
  Persistence_heat_maps<constant_scaling_function> p("data/file_with_diagram", filter, false, 5, 0, 1);

  std::vector<double> p_vect_template;

  p_vect_template.push_back(0.0606728);
  p_vect_template.push_back(0.0610023);
  p_vect_template.push_back(0.0607978);
  p_vect_template.push_back(0.0600647);
  p_vect_template.push_back(0.0588224);
  p_vect_template.push_back(0.0619829);
  p_vect_template.push_back(0.0623218);
  p_vect_template.push_back(0.0621152);
  p_vect_template.push_back(0.0613686);
  p_vect_template.push_back(0.0601016);
  p_vect_template.push_back(0.0627679);
  p_vect_template.push_back(0.0631134);
  p_vect_template.push_back(0.0629066);
  p_vect_template.push_back(0.0621528);
  p_vect_template.push_back(0.0608719);
  p_vect_template.push_back(0.0630073);
  p_vect_template.push_back(0.0633564);
  p_vect_template.push_back(0.0631511);
  p_vect_template.push_back(0.0623968);
  p_vect_template.push_back(0.0611132);
  p_vect_template.push_back(0.0626947);
  p_vect_template.push_back(0.0630445);
  p_vect_template.push_back(0.0628425);
  p_vect_template.push_back(0.0620941);
  p_vect_template.push_back(0.060819);

  std::vector<double> p_vect = p.vectorize(0);
  for (size_t i = 0; i != p_vect.size(); ++i) {
    BOOST_CHECK(almost_equal(p_vect_template[i], p_vect[i]));
  }
}

BOOST_AUTO_TEST_CASE(check_distance_for_heat_maps)
{
  std::vector<std::vector<double> > filter = create_Gaussian_filter(30, 1);
  Persistence_heat_maps<constant_scaling_function> p("data/file_with_diagram", filter, false, 1000, 0, 1);
  Persistence_heat_maps<constant_scaling_function> q("data/file_with_diagram_1", filter, false, 1000, 0, 1);
  Persistence_heat_maps<constant_scaling_function> r("data/file_with_diagram_2", filter, false, 1000, 0, 1);

  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.distance(p), 0., epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.distance(q), 624.183, epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.distance(r), 415.815, epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(q.distance(p), 624.183, epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(q.distance(q), 0., epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(q.distance(r), 528.066, epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(r.distance(p), 415.815, epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(r.distance(q), 528.066, epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(r.distance(r), 0., epsilon);
}

BOOST_AUTO_TEST_CASE(check_projections_to_R_for_heat_maps)
{
  std::vector<std::vector<double> > filter = create_Gaussian_filter(30, 1);
  Persistence_heat_maps<constant_scaling_function> p("data/file_with_diagram", filter, false, 1000, 0, 1);
  Persistence_heat_maps<constant_scaling_function> q("data/file_with_diagram_1", filter, false, 1000, 0, 1);
  Persistence_heat_maps<constant_scaling_function> r("data/file_with_diagram_2", filter, false, 1000, 0, 1);

  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.project_to_R(0), 44.3308, epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(q.project_to_R(0), 650.568, epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(r.project_to_R(0), 429.287, epsilon);
}

BOOST_AUTO_TEST_CASE(check_scalar_products_for_heat_maps)
{
  std::vector<std::vector<double> > filter = create_Gaussian_filter(30, 1);
  Persistence_heat_maps<constant_scaling_function> p("data/file_with_diagram", filter, false, 1000, 0, 1);
  Persistence_heat_maps<constant_scaling_function> q("data/file_with_diagram_1", filter, false, 1000, 0, 1);
  Persistence_heat_maps<constant_scaling_function> r("data/file_with_diagram_2", filter, false, 1000, 0, 1);

  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.compute_scalar_product(p), 0.0345687, epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.compute_scalar_product(q), 0.0509357, epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(p.compute_scalar_product(r), 0.0375608, epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(q.compute_scalar_product(p), 0.0509357, epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(q.compute_scalar_product(q), 1.31293, epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(q.compute_scalar_product(r), 0.536799, epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(r.compute_scalar_product(p), 0.0375608, epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(r.compute_scalar_product(q), 0.536799, epsilon);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(r.compute_scalar_product(r), 0.672907, epsilon);
}

BOOST_AUTO_TEST_CASE(check_arythmetic_operations_for_heat_maps)
{
  std::vector<std::vector<double> > filter = create_Gaussian_filter(30, 1);
  Persistence_heat_maps<constant_scaling_function> p("data/file_with_diagram", filter, false, 1000, 0, 1);
  Persistence_heat_maps<constant_scaling_function> q("data/file_with_diagram_1", filter, false, 1000, 0, 1);

  Persistence_heat_maps<constant_scaling_function> sum = p + q;
  Persistence_heat_maps<constant_scaling_function> difference = p - q;
  Persistence_heat_maps<constant_scaling_function> multiply_by_scalar = 2 * p;

  // sum.print_to_file( "sum" );
  // difference.print_to_file( "difference" );
  // multiply_by_scalar.print_to_file( "multiply_by_scalar" );

  Persistence_heat_maps<constant_scaling_function> sum_template;
  sum_template.load_from_file("data/heat_map_sum");
  Persistence_heat_maps<constant_scaling_function> difference_template;
  difference_template.load_from_file("data/heat_map_difference");
  Persistence_heat_maps<constant_scaling_function> multiply_by_scalar_template;
  multiply_by_scalar_template.load_from_file("data/heat_map_multiply_by_scalar");

  BOOST_CHECK(sum == sum_template);
}

BOOST_AUTO_TEST_CASE(check_distance_of_heat_maps_infinite_power_parameters)
{
  std::vector<std::vector<double> > filter = create_Gaussian_filter(100, 1);
  Persistence_heat_maps<constant_scaling_function> p("data/file_with_diagram", filter, false, 1000, 0, 1);

  std::vector<std::vector<double> > filter_2 = create_Gaussian_filter(150, 1);
  Persistence_heat_maps<constant_scaling_function> q("data/file_with_diagram", filter_2, true, 1000, 0, 1);

  double distance_max_double_parameter = p.distance(q, std::numeric_limits<double>::max());
  double distance_inf_double_parameter = p.distance(q, std::numeric_limits<double>::infinity());

  // std::clog << "distance_max_double_parameter: " << distance_max_double_parameter << std::endl;
  // std::clog << "distance_inf_double_parameter: " << distance_inf_double_parameter << std::endl;

  BOOST_CHECK(distance_max_double_parameter == distance_inf_double_parameter);
}

// Below I am storing the code used to generate tests for that functionality.
/*
        std::vector< std::pair< double,double > > intervals;
        intervals.push_back( std::make_pair(0.5,0.5) );
        std::vector< std::vector<double> > filter = create_Gaussian_filter(5,1);
        Persistence_heat_maps p( intervals , filter ,  constant_function, false , 100 , 0 , 1 );
        p.plot( "heat_map_1" );


        std::vector< std::pair< double,double > > intervals2;
        intervals2.push_back( std::make_pair(7,12) );
        Persistence_heat_maps q( intervals2 , filter ,  constant_function, false , 100 , 0 , 10 );
        q.plot( "heat_map_2" );
*/
/*
        std::vector< std::pair< double,double > > intervals;
        intervals.push_back( std::make_pair(0.5,0.5) );
        std::vector< std::vector<double> > filter = create_Gaussian_filter(5,1);
        Persistence_heat_maps p( intervals , filter ,  constant_function, false , 10 , 0 , 1 );
        p.write_to_file( "aaa" );

        Persistence_heat_maps q;
        q.load_from_file( "aaa" );

        cerr << ( p == q ) << endl;
        */

/*
        std::vector< std::vector<double> > filter = create_Gaussian_filter(30,1);
        Persistence_heat_maps p( "file_with_diagram" , filter ,  constant_function, false , 100 , 0 , 1 );
        p.plot( "heat_map_1" );
*/

/*
   //test to construct persistence heat map:
    std::vector< std::vector<double> > filter = create_Gaussian_filter(100,1);
        Persistence_heat_maps p( "file_with_diagram" , filter ,  constant_function, false , 1000 , 0 , 1 );
        p.print_to_file( "persistence_heat_map_from_file_with_diagram" );

        Persistence_heat_maps q;
        q.load_from_file( "persistence_heat_map_from_file_with_diagram" );

        cerr << (p == q) << endl;
*/
/*
        //test of computations of a mean:
        std::vector< std::pair< double,double > > intervals;
        intervals.push_back( std::make_pair(5,5) );
        std::vector< std::vector<double> > filter = create_Gaussian_filter(5,1);
        Persistence_heat_maps p( intervals , filter ,  constant_function, false , 100 , 0 , 10 );
        p.plot( "heat_map_1" );


        std::vector< std::pair< double,double > > intervals2;
        intervals2.push_back( std::make_pair(7,7) );
        Persistence_heat_maps q( intervals2 , filter ,  constant_function, false , 100 , 0 , 10 );
        q.plot( "heat_map_2" );


        Persistence_heat_maps av;
        av.compute_average( { &P , &q } );
        av.plot( "average" );
*/

/*
        std::vector< std::vector<double> > filter = create_Gaussian_filter(30,1);
        Persistence_heat_maps p( "file_with_diagram" , filter ,  constant_function, false , 1000 , 0 , 10 );
        Persistence_heat_maps q( "file_with_diagram_1" , filter ,  constant_function, false , 1000 , 0 , 10 );
        Persistence_heat_maps r( "file_with_diagram_2" , filter ,  constant_function, false , 1000 , 0 , 10 );
        Persistence_heat_maps av;
        av.compute_average( {&p,&q,&r} );

        av.print_to_file( "template_average_of_heat_maps" );
*/

/*
        std::vector< std::pair< double,double > > intervals;
        intervals.push_back( std::make_pair(5,5) );
        std::vector< std::vector<double> > filter = create_Gaussian_filter(5,1);
        Persistence_heat_maps p( intervals , filter ,  constant_function, false , 10 , 0 , 10 );
        p.plot( "heat_map_1" );

        std::vector< std::pair< double,double > > intervals2;
        intervals2.push_back( std::make_pair(7,7) );
        Persistence_heat_maps q( intervals2 , filter ,  constant_function, false , 10 , 0 , 10 );
        q.plot( "heat_map_2" );

        Persistence_heat_maps median;
        median.compute_median( {&p,&q} );
        median.plot( "median" );
*/

/*
        std::vector< std::vector<double> > filter = create_Gaussian_filter(30,1);
        Persistence_heat_maps p( "file_with_diagram" , filter ,  constant_function, false , 1000 , 0 , 1 );
        Persistence_heat_maps q( "file_with_diagram_1" , filter ,  constant_function, false , 1000 , 0 , 1 );
        Persistence_heat_maps r( "file_with_diagram_2" , filter ,  constant_function, false , 1000 , 0 , 1 );
        Persistence_heat_maps median;
        median.compute_median( {&p,&q,&r} );
        median.print_to_file( "template_median_of_heat_maps" );
*/

/*
        std::vector< std::vector<double> > filter = create_Gaussian_filter(30,1);
        Persistence_heat_maps p( "file_with_diagram" , filter ,  constant_function, false , 1000 , 0 , 1 );
        Persistence_heat_maps q( "file_with_diagram_1" , filter ,  constant_function, false , 1000 , 0 , 1 );
        Persistence_heat_maps r( "file_with_diagram_2" , filter ,  constant_function, false , 1000 , 0 , 1 );

        Persistence_heat_maps percentage_of_active;
        percentage_of_active.compute_percentage_of_active( {&p,&q,&r} , 0.1 );

        percentage_of_active.print_to_file( "template_percentage_of_active_of_heat_maps" );
        //percentage_of_active.plot( "template_percentage_of_active_of_heat_maps" );
*/

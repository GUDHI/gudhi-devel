/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016 Inria
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

#include <gudhi/Persistence_vectors.h>
#include <gudhi/common_persistence_representations.h>
#include <gudhi/read_persistence_from_file.h>
#include <iostream>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "vector_representation_test"
#include <boost/test/unit_test.hpp>
#include <gudhi/reader_utils.h>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace Gudhi;
using namespace Gudhi::Persistence_representations;

BOOST_AUTO_TEST_CASE(check_read_write_to_files) {
  std::vector<std::pair<double, double> > intervals;
  intervals.push_back(std::make_pair(2, 3));
  intervals.push_back(std::make_pair(4, 7));
  intervals.push_back(std::make_pair(9, 10));
  intervals.push_back(std::make_pair(3, 11));
  Vector_distances_in_diagram<Euclidean_distance> p(intervals, -1);
  p.write_to_file("test_vector_representation_write_read");

  Vector_distances_in_diagram<Euclidean_distance> q;
  q.load_from_file("test_vector_representation_write_read");

  BOOST_CHECK(p == q);
}

BOOST_AUTO_TEST_CASE(check_sortev_vector_distances_template) {
  Vector_distances_in_diagram<Euclidean_distance> p("data/file_with_diagram", 100);
  std::vector<double> sortev_vector_distances_template;

  sortev_vector_distances_template.push_back(0.609968);
  sortev_vector_distances_template.push_back(0.566317);
  sortev_vector_distances_template.push_back(0.538858);
  sortev_vector_distances_template.push_back(0.534927);
  sortev_vector_distances_template.push_back(0.515741);
  sortev_vector_distances_template.push_back(0.507828);
  sortev_vector_distances_template.push_back(0.500911);
  sortev_vector_distances_template.push_back(0.496986);
  sortev_vector_distances_template.push_back(0.495306);
  sortev_vector_distances_template.push_back(0.439945);
  sortev_vector_distances_template.push_back(0.424097);
  sortev_vector_distances_template.push_back(0.413891);
  sortev_vector_distances_template.push_back(0.413891);
  sortev_vector_distances_template.push_back(0.413891);
  sortev_vector_distances_template.push_back(0.412621);
  sortev_vector_distances_template.push_back(0.410613);
  sortev_vector_distances_template.push_back(0.407853);
  sortev_vector_distances_template.push_back(0.407853);
  sortev_vector_distances_template.push_back(0.402306);
  sortev_vector_distances_template.push_back(0.401937);
  sortev_vector_distances_template.push_back(0.377605);
  sortev_vector_distances_template.push_back(0.363859);
  sortev_vector_distances_template.push_back(0.357765);
  sortev_vector_distances_template.push_back(0.357765);
  sortev_vector_distances_template.push_back(0.357765);
  sortev_vector_distances_template.push_back(0.357765);
  sortev_vector_distances_template.push_back(0.357765);
  sortev_vector_distances_template.push_back(0.353401);
  sortev_vector_distances_template.push_back(0.348004);
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
  sortev_vector_distances_template.push_back(0.339466);
  sortev_vector_distances_template.push_back(0.33935);
  sortev_vector_distances_template.push_back(0.32834);
  sortev_vector_distances_template.push_back(0.327276);
  sortev_vector_distances_template.push_back(0.318626);
  sortev_vector_distances_template.push_back(0.318082);
  sortev_vector_distances_template.push_back(0.30603);
  sortev_vector_distances_template.push_back(0.30525);
  sortev_vector_distances_template.push_back(0.297308);
  sortev_vector_distances_template.push_back(0.296333);
  sortev_vector_distances_template.push_back(0.296333);
  sortev_vector_distances_template.push_back(0.296333);
  sortev_vector_distances_template.push_back(0.296333);
  sortev_vector_distances_template.push_back(0.293372);
  sortev_vector_distances_template.push_back(0.292666);
  sortev_vector_distances_template.push_back(0.292666);
  sortev_vector_distances_template.push_back(0.292666);
  sortev_vector_distances_template.push_back(0.292666);
  sortev_vector_distances_template.push_back(0.292666);
  sortev_vector_distances_template.push_back(0.292666);
  sortev_vector_distances_template.push_back(0.292666);
  sortev_vector_distances_template.push_back(0.292666);
  sortev_vector_distances_template.push_back(0.292666);
  sortev_vector_distances_template.push_back(0.292666);
  sortev_vector_distances_template.push_back(0.292666);
  sortev_vector_distances_template.push_back(0.292666);
  sortev_vector_distances_template.push_back(0.292666);
  sortev_vector_distances_template.push_back(0.292666);
  sortev_vector_distances_template.push_back(0.292666);
  sortev_vector_distances_template.push_back(0.29029);
  sortev_vector_distances_template.push_back(0.290218);
  sortev_vector_distances_template.push_back(0.289782);
  sortev_vector_distances_template.push_back(0.288128);
  sortev_vector_distances_template.push_back(0.286416);
  sortev_vector_distances_template.push_back(0.285969);
  sortev_vector_distances_template.push_back(0.282046);
  sortev_vector_distances_template.push_back(0.28154);
  sortev_vector_distances_template.push_back(0.281085);
  sortev_vector_distances_template.push_back(0.280227);
  sortev_vector_distances_template.push_back(0.279273);
  sortev_vector_distances_template.push_back(0.278936);
  sortev_vector_distances_template.push_back(0.278706);
  sortev_vector_distances_template.push_back(0.278507);
  sortev_vector_distances_template.push_back(0.278097);
  sortev_vector_distances_template.push_back(0.276293);
  sortev_vector_distances_template.push_back(0.276293);
  sortev_vector_distances_template.push_back(0.276293);
  sortev_vector_distances_template.push_back(0.276293);
  sortev_vector_distances_template.push_back(0.276293);
  sortev_vector_distances_template.push_back(0.276293);
  sortev_vector_distances_template.push_back(0.276293);
  sortev_vector_distances_template.push_back(0.276293);
  sortev_vector_distances_template.push_back(0.276293);
  sortev_vector_distances_template.push_back(0.276169);
  sortev_vector_distances_template.push_back(0.270563);
  sortev_vector_distances_template.push_back(0.264009);

  size_t proj_no = p.number_of_vectorize_functions();
  std::vector<double> aa = p.vectorize(proj_no);

  for (size_t i = 0; i != aa.size(); ++i) {
    BOOST_CHECK(almost_equal(sortev_vector_distances_template[i], aa[i]));
  }
}

BOOST_AUTO_TEST_CASE(check_projections_to_R) {
  Vector_distances_in_diagram<Euclidean_distance> p("data/file_with_diagram", 100);
  std::vector<double> proj;
  proj.push_back(0);
  proj.push_back(0.6099679993);
  proj.push_back(1.176284775);
  proj.push_back(1.715142954);
  proj.push_back(2.25006986);
  proj.push_back(2.765810506);
  proj.push_back(3.273638431);
  proj.push_back(3.774549309);
  proj.push_back(4.271535042);
  proj.push_back(4.766840875);
  proj.push_back(5.206786149);
  proj.push_back(5.63088295);
  proj.push_back(6.04477433);
  proj.push_back(6.45866571);
  proj.push_back(6.87255709);
  proj.push_back(7.285177939);
  proj.push_back(7.695791381);
  proj.push_back(8.103643945);
  proj.push_back(8.511496508);
  proj.push_back(8.913802775);
  proj.push_back(9.315740229);
  proj.push_back(9.693344927);
  proj.push_back(10.0572035);
  proj.push_back(10.41496899);
  proj.push_back(10.77273448);
  proj.push_back(11.13049996);
  proj.push_back(11.48826545);
  proj.push_back(11.84603094);
  proj.push_back(12.19943233);
  proj.push_back(12.5474364);
  proj.push_back(12.89256042);
  proj.push_back(13.23768444);
  proj.push_back(13.58280846);
  proj.push_back(13.92793248);
  proj.push_back(14.2730565);
  proj.push_back(14.61818051);
  proj.push_back(14.96330453);
  proj.push_back(15.30842855);
  proj.push_back(15.65355257);
  proj.push_back(15.99867659);
  proj.push_back(16.34380061);
  proj.push_back(16.68892462);
  proj.push_back(17.03404864);
  proj.push_back(17.37917266);
  proj.push_back(17.7238622);
  proj.push_back(18.06332781);
  proj.push_back(18.40267754);
  proj.push_back(18.73101759);
  proj.push_back(19.05829313);
  proj.push_back(19.3769189);
  proj.push_back(19.69500045);
  proj.push_back(20.0010306);
  proj.push_back(20.30628026);
  proj.push_back(20.60358868);
  proj.push_back(20.89992192);
  proj.push_back(21.19625516);
  proj.push_back(21.4925884);
  proj.push_back(21.78892164);
  proj.push_back(22.08229394);
  proj.push_back(22.37495987);
  proj.push_back(22.66762581);
  proj.push_back(22.96029174);
  proj.push_back(23.25295768);
  proj.push_back(23.54562361);
  proj.push_back(23.83828955);
  proj.push_back(24.13095549);
  proj.push_back(24.42362142);
  proj.push_back(24.71628736);
  proj.push_back(25.00895329);
  proj.push_back(25.30161923);
  proj.push_back(25.59428516);
  proj.push_back(25.8869511);
  proj.push_back(26.17961703);
  proj.push_back(26.47228297);
  proj.push_back(26.76257262);
  proj.push_back(27.05279049);
  proj.push_back(27.34257265);
  proj.push_back(27.63070097);
  proj.push_back(27.91711687);
  proj.push_back(28.20308566);
  proj.push_back(28.48513176);
  proj.push_back(28.76667161);
  proj.push_back(29.04775635);
  proj.push_back(29.32798359);
  proj.push_back(29.60725702);
  proj.push_back(29.88619335);
  proj.push_back(30.16489915);
  proj.push_back(30.44340655);
  proj.push_back(30.72150329);
  proj.push_back(30.99779604);
  proj.push_back(31.27408878);
  proj.push_back(31.55038153);
  proj.push_back(31.82667427);
  proj.push_back(32.10296702);
  proj.push_back(32.37925976);
  proj.push_back(32.6555525);
  proj.push_back(32.93184525);
  proj.push_back(33.20813799);
  proj.push_back(33.48430662);
  proj.push_back(33.7548692);

  for (size_t proj_no = 0; proj_no != p.number_of_projections_to_R(); ++proj_no) {
    // cout << std::setprecision(10)  << 	p.project_to_R(proj_no) << endl;
    BOOST_CHECK(almost_equal(p.project_to_R(proj_no), proj[proj_no]));
  }
}

BOOST_AUTO_TEST_CASE(check_distance_computations) {
  Vector_distances_in_diagram<Euclidean_distance> p("data/file_with_diagram", 100);
  Vector_distances_in_diagram<Euclidean_distance> p_prime("data/file_with_diagram", 10);
  std::vector<std::pair<double, double> > intervals(10);
  intervals[0] = std::pair<double, double>(1, 2);
  intervals[1] = std::pair<double, double>(3, 4);
  intervals[2] = std::pair<double, double>(5, 6);
  intervals[3] = std::pair<double, double>(7, 8);
  intervals[4] = std::pair<double, double>(9, 10);
  intervals[5] = std::pair<double, double>(11, 12);
  intervals[6] = std::pair<double, double>(13, 14);
  intervals[7] = std::pair<double, double>(15, 16);
  intervals[8] = std::pair<double, double>(17, 18);
  intervals[9] = std::pair<double, double>(19, 20);
  Vector_distances_in_diagram<Euclidean_distance> p_bis(intervals, 10);
  // cerr << "p_prime.distance( (Abs_Topological_data_with_distances*)(&p_bis) , 1 ) : " << p_prime.distance(
  // (Abs_Topological_data_with_distances*)(&p_bis) , 1 ) << endl;
  BOOST_CHECK(almost_equal(p_prime.distance(p_bis, 1), 1.86428));
}

BOOST_AUTO_TEST_CASE(check_default_parameters_of_distances) {
  std::vector<std::pair<double, double> > diag = read_persistence_intervals_in_dimension("data/file_with_diagram");
  Vector_distances_in_diagram<Euclidean_distance> p(diag, 100);

  std::vector<std::pair<double, double> > diag1 = read_persistence_intervals_in_dimension("data/file_with_diagram_1");
  Vector_distances_in_diagram<Euclidean_distance> q(diag1, 100);

  double dist_numeric_limit_max = p.distance(q, std::numeric_limits<double>::max());
  double dist_infinity = p.distance(q, std::numeric_limits<double>::infinity());

  BOOST_CHECK(dist_numeric_limit_max == dist_infinity);
}

BOOST_AUTO_TEST_CASE(check_compute_average) {
  Vector_distances_in_diagram<Euclidean_distance> p("data/file_with_diagram", 100);
  // compute average
  std::vector<std::pair<double, double> > i1(3);
  i1[0] = std::pair<double, double>(1, 2);
  i1[1] = std::pair<double, double>(3, 8);
  i1[2] = std::pair<double, double>(1, 6);

  std::vector<std::pair<double, double> > i2(3);
  i2[0] = std::pair<double, double>(2, 9);
  i2[1] = std::pair<double, double>(2, 15);
  i2[2] = std::pair<double, double>(6, 17);

  Vector_distances_in_diagram<Euclidean_distance> A(i1, -1);
  Vector_distances_in_diagram<Euclidean_distance> B(i1, -1);

  Vector_distances_in_diagram<Euclidean_distance> average;
  average.compute_average({&A, &B});

  Vector_distances_in_diagram<Euclidean_distance> template_average;
  template_average.load_from_file("data/average_of_persistence_vectors");

  BOOST_CHECK(template_average == average);
}

BOOST_AUTO_TEST_CASE(check_arythmetic_operations) {
  std::vector<std::pair<double, double> > i1(3);
  i1[0] = std::pair<double, double>(1, 2);
  i1[1] = std::pair<double, double>(3, 8);
  i1[2] = std::pair<double, double>(1, 6);

  std::vector<std::pair<double, double> > i2(3);
  i2[0] = std::pair<double, double>(2, 9);
  i2[1] = std::pair<double, double>(2, 15);
  i2[2] = std::pair<double, double>(6, 17);

  Vector_distances_in_diagram<Euclidean_distance> A(i1, -1);
  Vector_distances_in_diagram<Euclidean_distance> B(i1, -1);

  Vector_distances_in_diagram<Euclidean_distance> sum = A + B;
  Vector_distances_in_diagram<Euclidean_distance> difference = A - B;

  BOOST_CHECK(almost_equal(sum.vector_in_position(0), 7.07107));
  BOOST_CHECK(almost_equal(sum.vector_in_position(1), 7.07107));
  BOOST_CHECK(almost_equal(sum.vector_in_position(2), 5.65685));
  BOOST_CHECK(almost_equal(sum.vector_in_position(3), 1.41421));
  BOOST_CHECK(almost_equal(sum.vector_in_position(4), 1.41421));
  BOOST_CHECK(almost_equal(sum.vector_in_position(5), 1.41421));

  BOOST_CHECK(almost_equal(difference.vector_in_position(0), 0));
  BOOST_CHECK(almost_equal(difference.vector_in_position(1), 0));
  BOOST_CHECK(almost_equal(difference.vector_in_position(2), 0));
  BOOST_CHECK(almost_equal(difference.vector_in_position(3), 0));
  BOOST_CHECK(almost_equal(difference.vector_in_position(4), 0));
  BOOST_CHECK(almost_equal(difference.vector_in_position(5), 0));

  Vector_distances_in_diagram<Euclidean_distance> prod = 2. * A;
  BOOST_CHECK(almost_equal(prod.vector_in_position(0), 7.07107));
  BOOST_CHECK(almost_equal(prod.vector_in_position(1), 7.07107));
  BOOST_CHECK(almost_equal(prod.vector_in_position(2), 5.65685));
  BOOST_CHECK(almost_equal(prod.vector_in_position(3), 1.41421));
  BOOST_CHECK(almost_equal(prod.vector_in_position(4), 1.41421));
  BOOST_CHECK(almost_equal(prod.vector_in_position(5), 1.41421));

  Vector_distances_in_diagram<Euclidean_distance> prod1 = A * 2;
  BOOST_CHECK(almost_equal(prod1.vector_in_position(0), 7.07107));
  BOOST_CHECK(almost_equal(prod1.vector_in_position(1), 7.07107));
  BOOST_CHECK(almost_equal(prod1.vector_in_position(2), 5.65685));
  BOOST_CHECK(almost_equal(prod1.vector_in_position(3), 1.41421));
  BOOST_CHECK(almost_equal(prod1.vector_in_position(4), 1.41421));
  BOOST_CHECK(almost_equal(prod1.vector_in_position(5), 1.41421));
}

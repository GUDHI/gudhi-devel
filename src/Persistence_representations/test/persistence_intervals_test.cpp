/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 *      - 2019/12 Vincent Rouvreau: Fix #118 - Make histogram_of_lengths and cumulative_histogram_of_lengths
 *          return the exact number_of_bins (was failing on x86)
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "Persistence_intervals_test"
#include <boost/test/unit_test.hpp>
#include <gudhi/reader_utils.h>
#include "gudhi/Persistence_intervals.h"
#include <gudhi/common_persistence_representations.h>
#include <gudhi/Unitary_tests_utils.h>

#include <iostream>

using namespace Gudhi;
using namespace Gudhi::Persistence_representations;

BOOST_AUTO_TEST_CASE(check_min_max_function) {
  Persistence_intervals p("data/file_with_diagram");
  std::pair<double, double> min_max_ = p.get_x_range();

  GUDHI_TEST_FLOAT_EQUALITY_CHECK(min_max_.first, 0.0290362, Gudhi::Persistence_representations::epsi);
  GUDHI_TEST_FLOAT_EQUALITY_CHECK(min_max_.second, 0.994537, Gudhi::Persistence_representations::epsi);
}

BOOST_AUTO_TEST_CASE(check_length_of_dominant_intervals) {
  Persistence_intervals p("data/file_with_diagram");
  std::vector<double> dominant_ten_intervals_length = p.length_of_dominant_intervals(10);
  std::vector<double> dominant_intervals_length{0.862625, 0.800893, 0.762061, 0.756501, 0.729367,
                                                0.718177, 0.708395, 0.702844, 0.700468, 0.622177};
  for (size_t i = 0; i != dominant_ten_intervals_length.size(); ++i) {
    GUDHI_TEST_FLOAT_EQUALITY_CHECK(dominant_ten_intervals_length[i], dominant_intervals_length[i],
                                    Gudhi::Persistence_representations::epsi);
  }
}
BOOST_AUTO_TEST_CASE(check_dominant_intervals) {
  Persistence_intervals p("data/file_with_diagram");
  std::vector<std::pair<double, double> > ten_dominant_intervals = p.dominant_intervals(10);

  std::vector<std::pair<double, double> > templ{ {0.114718, 0.977343}, {0.133638, 0.93453},
                                                 {0.104599, 0.866659}, {0.149798, 0.906299},
                                                 {0.247352, 0.976719}, {0.192675, 0.910852},
                                                 {0.191836, 0.900231}, {0.284998, 0.987842},
                                                 {0.294069, 0.994537}, {0.267421, 0.889597} };

  for (size_t i = 0; i != ten_dominant_intervals.size(); ++i) {
    GUDHI_TEST_FLOAT_EQUALITY_CHECK(ten_dominant_intervals[i].first, templ[i].first,
                                    Gudhi::Persistence_representations::epsi);
    GUDHI_TEST_FLOAT_EQUALITY_CHECK(ten_dominant_intervals[i].second, templ[i].second,
                                    Gudhi::Persistence_representations::epsi);
  }
}

BOOST_AUTO_TEST_CASE(check_histogram_of_lengths) {
  Persistence_intervals p("data/file_with_diagram");
  std::vector<size_t> histogram = p.histogram_of_lengths(10);
  std::vector<size_t> template_histogram{10, 5, 3, 4, 4, 3, 6, 1, 7, 2};
  for (size_t i = 0; i != histogram.size(); ++i) {
    BOOST_CHECK(histogram[i] == template_histogram[i]);
  }
}

BOOST_AUTO_TEST_CASE(check_cumulative_histograms_of_lengths) {
  Persistence_intervals p("data/file_with_diagram");
  std::vector<size_t> cumulative_histogram = p.cumulative_histogram_of_lengths(10);
  std::vector<size_t> template_cumulative_histogram{10, 15, 18, 22, 26, 29, 35, 36, 43, 45};

  for (size_t i = 0; i != cumulative_histogram.size(); ++i) {
    BOOST_CHECK(cumulative_histogram[i] == template_cumulative_histogram[i]);
  }
}
BOOST_AUTO_TEST_CASE(check_characteristic_function_of_diagram) {
  Persistence_intervals p("data/file_with_diagram");
  std::pair<double, double> min_max_ = p.get_x_range();
  std::vector<double> char_funct_diag = p.characteristic_function_of_diagram(min_max_.first, min_max_.second);
  std::vector<double> template_char_funct_diag{0.370665, 0.84058, 1.24649, 1.3664, 1.34032,
                                               1.31904, 1.14076, 0.991259, 0.800714, 0.0676303};

  for (size_t i = 0; i != char_funct_diag.size(); ++i) {
    GUDHI_TEST_FLOAT_EQUALITY_CHECK(char_funct_diag[i], template_char_funct_diag[i],
                                    Gudhi::Persistence_representations::epsi);
  }
}

BOOST_AUTO_TEST_CASE(check_cumulative_characteristic_function_of_diagram) {
  Persistence_intervals p("data/file_with_diagram");
  std::pair<double, double> min_max_ = p.get_x_range();
  std::vector<double> cumul_char_funct_diag =
      p.cumulative_characteristic_function_of_diagram(min_max_.first, min_max_.second);
  std::vector<double> template_char_funct_diag_cumul{0.370665, 1.21125, 2.45774, 3.82414, 5.16446,
                                                     6.4835, 7.62426, 8.61552, 9.41623, 9.48386};

  for (size_t i = 0; i != cumul_char_funct_diag.size(); ++i) {
    GUDHI_TEST_FLOAT_EQUALITY_CHECK(cumul_char_funct_diag[i], template_char_funct_diag_cumul[i],
                                    Gudhi::Persistence_representations::epsi);
  }
}

BOOST_AUTO_TEST_CASE(check_compute_persistent_betti_numbers) {
  Persistence_intervals p("data/file_with_diagram");
  std::vector<std::pair<double, size_t> > pbns{ {0.0290362, 1}, {0.0307676, 2}, {0.0366312, 3}, {0.0544614, 4},
                                                {0.0920033, 5}, {0.104599, 6}, {0.114718, 7}, {0.117379, 8},
                                                {0.123493, 9}, {0.133638, 10}, {0.137798, 9}, {0.149798, 10},
                                                {0.155421, 11}, {0.158443, 12}, {0.176956, 13}, {0.183234, 12},
                                                {0.191069, 13}, {0.191333, 14}, {0.191836, 15}, {0.192675, 16},
                                                {0.208564, 17}, {0.218425, 18}, {0.219902, 17}, {0.23233, 16},
                                                {0.234558, 17}, {0.237166, 16}, {0.247352, 17}, {0.267421, 18},
                                                {0.268093, 19}, {0.278734, 18}, {0.284722, 19}, {0.284998, 20},
                                                {0.294069, 21}, {0.306293, 22}, {0.322361, 21}, {0.323152, 22},
                                                {0.371021, 23}, {0.372395, 24}, {0.387744, 25}, {0.435537, 26},
                                                {0.462911, 25}, {0.483569, 26}, {0.489209, 25}, {0.517115, 24},
                                                {0.522197, 23}, {0.532665, 22}, {0.545262, 23}, {0.587227, 22},
                                                {0.593036, 23}, {0.602647, 24}, {0.605044, 25}, {0.621962, 24},
                                                {0.629449, 23}, {0.636719, 22}, {0.64957, 21}, {0.650781, 22},
                                                {0.654951, 23}, {0.683489, 24}, {0.687172, 23}, {0.69703, 22},
                                                {0.701174, 21}, {0.717623, 22}, {0.722023, 21}, {0.722298, 20},
                                                {0.725347, 19}, {0.73071, 18}, {0.758355, 17}, {0.770913, 18},
                                                {0.790833, 17}, {0.821211, 16}, {0.849305, 17}, {0.853669, 16},
                                                {0.866659, 15}, {0.872896, 16}, {0.889597, 15}, {0.900231, 14},
                                                {0.903847, 13}, {0.906299, 12}, {0.910852, 11}, {0.93453, 10},
                                                {0.944757, 9}, {0.947812, 8}, {0.959154, 7}, {0.975654, 6},
                                                {0.976719, 5}, {0.977343, 4}, {0.980129, 3}, {0.987842, 2},
                                                {0.990127, 1}, {0.994537, 0} };

  std::vector<std::pair<double, size_t> > pbns_new = p.compute_persistent_betti_numbers();
  for (size_t i = 0; i != pbns.size(); ++i) {
    GUDHI_TEST_FLOAT_EQUALITY_CHECK(pbns[i].first, pbns_new[i].first, Gudhi::Persistence_representations::epsi);
    BOOST_CHECK(pbns[i].second == pbns_new[i].second);
  }
}

BOOST_AUTO_TEST_CASE(check_k_n_n) {
  Persistence_intervals p("data/file_with_diagram");
  std::vector<double> knn = p.k_n_n(5);
  std::vector<double> knn_template{1.04208, 1.00344, 0.979395, 0.890643, 0.874769,
                                   0.845787, 0.819713, 0.803984, 0.799864, 0.786945};

  for (size_t i = 0; i != knn.size(); ++i) {
    GUDHI_TEST_FLOAT_EQUALITY_CHECK(knn[i], knn_template[i], Gudhi::Persistence_representations::epsi);
  }
}

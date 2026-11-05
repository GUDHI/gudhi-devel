/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author:       Francois Godi
 *
 *    Copyright (C) 2015 Inria
 *
 *    Modification(s):
 *      - 2026/04 Vincent Rouvreau: Use Gudhi::random in place of c++ custom use
 *      - YYYY/MM Author: Description of the modification
 */


#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "bottleneck distance"
#include <boost/test/unit_test.hpp>

#include <gudhi/Bottleneck.h>
#include <gudhi/Unitary_tests_utils.h>
#include <gudhi/random.h>

using namespace Gudhi::persistence_diagram;

int n1 = 81;  // a natural number >0
int n2 = 180;  // a natural number >0
double upper_bound = 406.43;  // any real >0


std::vector< std::pair<double, double> > v1, v2;

BOOST_AUTO_TEST_CASE(persistence_graph) {
  // Random construction
  for (int i = 0; i < n1; i++) {
    double a = Gudhi::random::get_uniform<double>(0., upper_bound);;
    double b = Gudhi::random::get_uniform<double>(0., upper_bound);;
    v1.emplace_back(std::min(a, b), std::max(a, b));
  }
  for (int i = 0; i < n2; i++) {
    double a = Gudhi::random::get_uniform<double>(0., upper_bound);;
    double b = Gudhi::random::get_uniform<double>(0., upper_bound);;
    v2.emplace_back(std::min(a, b), std::max(a, b));
  }
  Persistence_graph g(v1, v2, 0.);
  std::vector<double> d(g.sorted_distances());
  //
  BOOST_CHECK(!g.on_the_u_diagonal(n1 - 1));
  BOOST_CHECK(!g.on_the_u_diagonal(n1));
  BOOST_CHECK(!g.on_the_u_diagonal(n2 - 1));
  BOOST_CHECK(g.on_the_u_diagonal(n2));
  BOOST_CHECK(!g.on_the_v_diagonal(n1 - 1));
  BOOST_CHECK(g.on_the_v_diagonal(n1));
  BOOST_CHECK(g.on_the_v_diagonal(n2 - 1));
  BOOST_CHECK(g.on_the_v_diagonal(n2));
  //
  BOOST_CHECK(g.corresponding_point_in_u(0) == n2);
  BOOST_CHECK(g.corresponding_point_in_u(n1) == 0);
  BOOST_CHECK(g.corresponding_point_in_v(0) == n1);
  BOOST_CHECK(g.corresponding_point_in_v(n2) == 0);
  //
  BOOST_CHECK(g.size() == (n1 + n2));
  //
  BOOST_CHECK((int) d.size() == (n1 + n2)*(n1 + n2) + n1 + n2 + 1);
  BOOST_CHECK(std::count(d.begin(), d.end(), GUDHI_PROTECT_FLOAT(g.distance(0, 0))) > 0);
  BOOST_CHECK(std::count(d.begin(), d.end(), GUDHI_PROTECT_FLOAT(g.distance(0, n1 - 1))) > 0);
  BOOST_CHECK(std::count(d.begin(), d.end(), GUDHI_PROTECT_FLOAT(g.distance(0, n1))) > 0);
  BOOST_CHECK(std::count(d.begin(), d.end(), GUDHI_PROTECT_FLOAT(g.distance(0, n2 - 1))) > 0);
  BOOST_CHECK(std::count(d.begin(), d.end(), GUDHI_PROTECT_FLOAT(g.distance(0, n2))) > 0);
  BOOST_CHECK(std::count(d.begin(), d.end(), GUDHI_PROTECT_FLOAT(g.distance(0, (n1 + n2) - 1))) > 0);
  BOOST_CHECK(std::count(d.begin(), d.end(), GUDHI_PROTECT_FLOAT(g.distance(n1, 0))) > 0);
  BOOST_CHECK(std::count(d.begin(), d.end(), GUDHI_PROTECT_FLOAT(g.distance(n1, n1 - 1))) > 0);
  BOOST_CHECK(std::count(d.begin(), d.end(), GUDHI_PROTECT_FLOAT(g.distance(n1, n1))) > 0);
  BOOST_CHECK(std::count(d.begin(), d.end(), GUDHI_PROTECT_FLOAT(g.distance(n1, n2 - 1))) > 0);
  BOOST_CHECK(std::count(d.begin(), d.end(), GUDHI_PROTECT_FLOAT(g.distance(n1, n2))) > 0);
  BOOST_CHECK(std::count(d.begin(), d.end(), GUDHI_PROTECT_FLOAT(g.distance(n1, (n1 + n2) - 1))) > 0);
  BOOST_CHECK(std::count(d.begin(), d.end(), GUDHI_PROTECT_FLOAT(g.distance((n1 + n2) - 1, 0))) > 0);
  BOOST_CHECK(std::count(d.begin(), d.end(), GUDHI_PROTECT_FLOAT(g.distance((n1 + n2) - 1, n1 - 1))) > 0);
  BOOST_CHECK(std::count(d.begin(), d.end(), GUDHI_PROTECT_FLOAT(g.distance((n1 + n2) - 1, n1))) > 0);
  BOOST_CHECK(std::count(d.begin(), d.end(), GUDHI_PROTECT_FLOAT(g.distance((n1 + n2) - 1, n2 - 1))) > 0);
  BOOST_CHECK(std::count(d.begin(), d.end(), GUDHI_PROTECT_FLOAT(g.distance((n1 + n2) - 1, n2))) > 0);
  BOOST_CHECK(std::count(d.begin(), d.end(), GUDHI_PROTECT_FLOAT(g.distance((n1 + n2) - 1, (n1 + n2) - 1))) > 0);
}

BOOST_AUTO_TEST_CASE(neighbors_finder) {
  Persistence_graph g(v1, v2, 0.);
  Neighbors_finder nf(g, 1.);
  for (int v_point_index = 1; v_point_index < ((n2 + n1)*9 / 10); v_point_index += 2)
    nf.add(v_point_index);
  //
  int v_point_index_1 = nf.pull_near(n2 / 2);
  BOOST_CHECK((v_point_index_1 == -1) || (g.distance(n2 / 2, v_point_index_1) <= 1.));
  std::vector<int> l = nf.pull_all_near(n2 / 2);
  bool v = true;
  for (auto it = l.cbegin(); it != l.cend(); ++it)
    v = v && (g.distance(n2 / 2, *it) > 1.);
  BOOST_CHECK(v);
  int v_point_index_2 = nf.pull_near(n2 / 2);
  BOOST_CHECK(v_point_index_2 == -1);
}

BOOST_AUTO_TEST_CASE(layered_neighbors_finder) {
  Persistence_graph g(v1, v2, 0.);
  Layered_neighbors_finder lnf(g, 1.);
  for (int v_point_index = 1; v_point_index < ((n2 + n1)*9 / 10); v_point_index += 2)
    lnf.add(v_point_index, v_point_index % 7);
  //
  int v_point_index_1 = lnf.pull_near(n2 / 2, 6);
  BOOST_CHECK((v_point_index_1 == -1) || (g.distance(n2 / 2, v_point_index_1) <= 1.));
  int v_point_index_2 = lnf.pull_near(n2 / 2, 6);
  BOOST_CHECK(v_point_index_2 == -1);
  v_point_index_1 = lnf.pull_near(n2 / 2, 0);
  BOOST_CHECK((v_point_index_1 == -1) || (g.distance(n2 / 2, v_point_index_1) <= 1.));
  v_point_index_2 = lnf.pull_near(n2 / 2, 0);
  BOOST_CHECK(v_point_index_2 == -1);
}

BOOST_AUTO_TEST_CASE(graph_matching) {
  Persistence_graph g(v1, v2, 0.);
  Graph_matching m1(g);
  m1.set_r(0.);
  int e = 0;
  while (m1.multi_augment())
    ++e;
  BOOST_CHECK(e > 0);
  BOOST_CHECK(e <= 2 * sqrt(2 * (n1 + n2)));
  Graph_matching m2 = m1;
  BOOST_CHECK(!m2.multi_augment());
  m2.set_r(upper_bound);
  e = 0;
  while (m2.multi_augment())
    ++e;
  BOOST_CHECK(e <= 2 * sqrt(2 * (n1 + n2)));
  BOOST_CHECK(m2.perfect());
  BOOST_CHECK(!m1.perfect());
}

void check_bottleneck_counterexample(const std::vector<std::pair<double, double>>& a,
                                     const std::vector<std::pair<double, double>>& b,
                                     double expected) {
  double exact_ab = bottleneck_distance(a, b, 0.);
  double exact_ba = bottleneck_distance(b, a, 0.);
  // Python's e=None maps to this overload's default value, the smallest positive double.
  double default_ab = bottleneck_distance(a, b);
  double default_ba = bottleneck_distance(b, a);

  BOOST_CHECK_EQUAL(exact_ab, exact_ba);
  BOOST_CHECK_CLOSE_FRACTION(exact_ab, expected, 1e-14);
  BOOST_CHECK_EQUAL(default_ab, default_ba);
  BOOST_CHECK_CLOSE_FRACTION(default_ab, expected, 1e-14);
}

// With the old shortcut, MSVC reaches n=6, k=3, l=3: k>sqrt(n), l>=sqrt(n), and k(l+1)=2n.
BOOST_AUTO_TEST_CASE(bottleneck_windows_counterexample) {
  std::vector<std::pair<double, double>> a = {
      {0.8657832199712898, 1.1256902536912314},
      {0.16921569257803723, 0.6031897594038342},
      {-0.04675641745672246, 0.5202109447215970},
  };
  std::vector<std::pair<double, double>> b = {
      {0.15043295756523575, 0.18414854990571342},
      {0.15043295756523575, 0.18414854990571342},
      {0.6796253715781101, 1.3001211960521455},
  };

  check_bottleneck_counterexample(a, b, 0.28348368108915978);
}

// A non-affine perturbation of the previous diagram reaches the same (n,k,l) on MSVC.
BOOST_AUTO_TEST_CASE(bottleneck_windows_perturbed_counterexample) {
  std::vector<std::pair<double, double>> a = {
      {0.86588321997128981, 1.1257902536912314},
      {0.16901569257803722, 0.60298975940383426},
      {-0.046456417456722458, 0.52051094472159698},
  };
  std::vector<std::pair<double, double>> b = {
      {0.15033295756523576, 0.18404854990571343},
      {0.15033295756523576, 0.18404854990571343},
      {0.67982537157811007, 1.3003211960521455},
  };

  check_bottleneck_counterexample(a, b, 0.28348368108915972);
}

// With the old shortcut, GCC reaches n=8, k=3, l=3: k>sqrt(n), l>=sqrt(n), and k(l+1)<2n.
BOOST_AUTO_TEST_CASE(bottleneck_linux_counterexample) {
  std::vector<std::pair<double, double>> a = {
      {-0.7734086850435968, -0.7734065850126042},
      {0.4057382451550824, 0.40574118001625126},
      {-0.8218366177063472, -0.8136106266959778},
      {0.4057382451550824, 0.40574118001625126},
  };
  std::vector<std::pair<double, double>> b = {
      {-0.6759744467664093, 0.17199054436205374},
      {-0.4787571504589536, -0.47692522493813944},
      {-0.6759744467664093, 0.17199054436205374},
      {-0.6759744467664093, 0.17199054436205374},
  };

  check_bottleneck_counterexample(a, b, 0.4239824955642315);
}

// A non-affine perturbation of the previous diagram reaches the same (n,k,l) on GCC.
BOOST_AUTO_TEST_CASE(bottleneck_linux_perturbed_counterexample) {
  std::vector<std::pair<double, double>> a = {
      {-0.77340768504359680, -0.77340558501260415},
      {0.40573624515508239, 0.40573918001625126},
      {-0.82183361770634722, -0.81360762669597786},
      {0.40573624515508239, 0.40573918001625126},
  };
  std::vector<std::pair<double, double>> b = {
      {-0.67597544676640930, 0.17198954436205374},
      {-0.47875515045895362, -0.47692322493813943},
      {-0.67597544676640930, 0.17198954436205374},
      {-0.67597544676640930, 0.17198954436205374},
  };

  check_bottleneck_counterexample(a, b, 0.4239824955642315);
}

BOOST_AUTO_TEST_CASE(global) {
  double delta_min = upper_bound / 1000.;
  double delta_max = upper_bound / 100.;

  std::vector< std::pair<double, double> > v1, v2;
  for (int i = 0; i < n1; i++) {
    double a = Gudhi::random::get_uniform<double>(0., upper_bound);
    double b = Gudhi::random::get_uniform<double>(0., upper_bound);
    double x = Gudhi::random::get_uniform<double>(delta_min, delta_max);
    double y = Gudhi::random::get_uniform<double>(delta_min, delta_max);
    v1.emplace_back(std::min(a, b), std::max(a, b));
    v2.emplace_back(std::min(a, b) + std::min(x, y), std::max(a, b) + std::max(x, y));
    if (i % 5 == 0)
      v1.emplace_back(std::min(a, b), std::min(a, b) + x);
    if (i % 3 == 0)
      v2.emplace_back(std::max(a, b), std::max(a, b) + y);
  }
  BOOST_CHECK(bottleneck_distance(v1, v2, 0.) <= upper_bound / 100.);
  BOOST_CHECK(bottleneck_distance(v1, v2, upper_bound / 10000.) <= upper_bound / 100. + upper_bound / 10000.);
  BOOST_CHECK(std::abs(bottleneck_distance(v1, v2, 0.) - bottleneck_distance(v1, v2, upper_bound / 10000.)) <= upper_bound / 10000.);

  std::vector< std::pair<double, double> > empty;
  std::vector< std::pair<double, double> > one = {{8, 10}};
  BOOST_CHECK(bottleneck_distance(empty, empty) == 0);
  BOOST_CHECK(bottleneck_distance(empty, one) == 1);
}

BOOST_AUTO_TEST_CASE(neg_global) {
  double delta_min = upper_bound / 1000.;
  double delta_max = upper_bound / 100.;

  std::vector< std::pair<double, double> > v1, v2;
  for (int i = 0; i < n1; i++) {
    double a = std::log(Gudhi::random::get_uniform<double>(0., upper_bound));
    double b = std::log(Gudhi::random::get_uniform<double>(0., upper_bound));
    double x = std::log(Gudhi::random::get_uniform<double>(delta_min, delta_max));
    double y = std::log(Gudhi::random::get_uniform<double>(delta_min, delta_max));
    v1.emplace_back(std::min(a, b), std::max(a, b));
    v2.emplace_back(std::min(a, b) + std::min(x, y), std::max(a, b) + std::max(x, y));
    if (i % 5 == 0)
      v1.emplace_back(std::min(a, b), std::min(a, b) + x);
    if (i % 3 == 0)
      v2.emplace_back(std::max(a, b), std::max(a, b) + y);
  }
  BOOST_CHECK(bottleneck_distance(v1, v2, 0.) <= upper_bound / 100.);
  BOOST_CHECK(bottleneck_distance(v1, v2, upper_bound / 10000.) <= upper_bound / 100. + upper_bound / 10000.);
  BOOST_CHECK(std::abs(bottleneck_distance(v1, v2, 0.) - bottleneck_distance(v1, v2, upper_bound / 10000.)) <= upper_bound / 10000.);

  std::vector< std::pair<double, double> > empty;
  std::vector< std::pair<double, double> > one = {{8, 10}};
  BOOST_CHECK(bottleneck_distance(empty, empty) == 0);
  BOOST_CHECK(bottleneck_distance(empty, one) == 1);
}

BOOST_AUTO_TEST_CASE(bottleneck_simple_test) {
  std::vector< std::pair<double, double> > v1, v2;
  double inf = std::numeric_limits<double>::infinity();
  double neginf = -inf;
  double b;

  v1.emplace_back(9.6, 14.);
  v2.emplace_back(9.5, 14.1);

  b = Gudhi::persistence_diagram::bottleneck_distance(v1, v2, 0.);
  BOOST_CHECK(b > 0.09 && b < 0.11);

  v1.emplace_back(-34.974, -34.2);

  b = Gudhi::persistence_diagram::bottleneck_distance(v1, v2, 0.);
  BOOST_CHECK(b > 0.386 && b < 0.388);

  v1.emplace_back(neginf, 3.7);

  b = Gudhi::persistence_diagram::bottleneck_distance(v1, v2, 0.);
  BOOST_CHECK_EQUAL(b, inf);

  v2.emplace_back(neginf, 4.45);

  b = Gudhi::persistence_diagram::bottleneck_distance(v1, v2, 0.);
  BOOST_CHECK(b > 0.74 && b < 0.76);

  v1.emplace_back(-60.6, 52.1);
  v2.emplace_back(-61.5, 53.);

  b = Gudhi::persistence_diagram::bottleneck_distance(v1, v2, 0.);
  BOOST_CHECK(b > 0.89 && b < 0.91);

  v1.emplace_back(3., inf);
  v2.emplace_back(3.2, inf);

  b = Gudhi::persistence_diagram::bottleneck_distance(v1, v2, 0.);
  BOOST_CHECK(b > 0.89 && b < 0.91);

  v1.emplace_back(neginf, inf);
  v2.emplace_back(neginf, inf);

  b = Gudhi::persistence_diagram::bottleneck_distance(v1, v2, 0.);
  BOOST_CHECK(b > 0.89 && b < 0.91);

  v2.emplace_back(6, inf);

  b = Gudhi::persistence_diagram::bottleneck_distance(v1, v2, 0.);
  BOOST_CHECK_EQUAL(b, inf);
}

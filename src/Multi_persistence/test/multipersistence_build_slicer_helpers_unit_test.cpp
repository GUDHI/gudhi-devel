/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <initializer_list>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "multi_persistence"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/Simplex_tree.h>
#include <gudhi/Multi_parameter_filtered_complex.h>
#include <gudhi/Multi_parameter_filtration.h>
#include <gudhi/Dynamic_multi_parameter_filtration.h>
#include <gudhi/slicer_helpers.h>
#include <gudhi/Slicer.h>
#include <gudhi/Multi_persistence/Persistence_interface_cohomology.h>
#include <gudhi/multi_simplex_tree_helpers.h>

using Gudhi::multi_filtration::Dynamic_multi_parameter_filtration;
using Gudhi::multi_filtration::Multi_parameter_filtration;
using Gudhi::multi_persistence::build_complex_from_bitmap;
using Gudhi::multi_persistence::build_complex_from_scc_file;
using Gudhi::multi_persistence::build_complex_from_simplex_tree;
using Gudhi::multi_persistence::build_slicer_from_bitmap;
using Gudhi::multi_persistence::build_slicer_from_scc_file;
using Gudhi::multi_persistence::build_slicer_from_simplex_tree;
using Gudhi::multi_persistence::Multi_parameter_filtered_complex;
using Gudhi::multi_persistence::Persistence_interface_cohomology;
using Gudhi::multi_persistence::Slicer;
using Gudhi::multi_persistence::write_complex_to_scc_file;
using Gudhi::multi_persistence::Simplex_tree_options_multidimensional_filtration;
using Gudhi::multi_persistence::Persistence_interface_cohomology;

typedef boost::mpl::list<Multi_parameter_filtration<double>, Dynamic_multi_parameter_filtration<double> >
    list_of_tested_variants;

template <class Fil>
Multi_parameter_filtered_complex<Fil> build_input_complex()
{
  using Complex = Multi_parameter_filtered_complex<Fil>;
  using FC = typename Complex::Filtration_value_container;
  using BC = typename Complex::Boundary_container;
  using DC = typename Complex::Dimension_container;
  using ini = std::initializer_list<typename Fil::value_type>;

  BC bc = {{}, {}, {}, {}, {}, {0, 4}, {1, 4}, {0, 1}, {2, 4}, {3, 4}, {2, 3}, {1, 3}, {5, 6, 7}, {8, 9, 10}};

  DC dc = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2};

  FC fc = {ini{0.636869, 0.0729406},
           ini{1.27032, 0.0729406},
           ini{0.686347, 0.188994},
           ini{1.27611, 0.188994},
           ini{0.144742, 0.2130470},
           ini{0.287974, 0.2130470},
           ini{0.151118, 0.2130471},
           ini{0.29683, 0.2130471},
           ini{0.232049, 0.226328},
           ini{0.416981, 0.226328},
           ini{0.454557, 0.226328},
           ini{0.635158, 0.226328},
           ini{0.639283, 0.226328},
           ini{0.82801, 0.226328}};

  return Complex(bc, dc, fc);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(slicer_scc_file_io, Fil, list_of_tested_variants)
{
  using Complex = Multi_parameter_filtered_complex<Fil>;
  using FC = typename Complex::Filtration_value_container;
  using BC = typename Complex::Boundary_container;
  using DC = typename Complex::Dimension_container;
  using ini = std::initializer_list<typename Fil::value_type>;
  using Slicer = Slicer<Fil, Persistence_interface_cohomology<Fil> >;

  const std::string& filePath = "scc_test.txt";
  int degree = 1;
  bool rivetCompatible = false;
  bool ignoreLastGenerators = false;
  bool reverse = false;
  int dimShift = 0;

  Complex in_cpx = build_input_complex<Fil>();
  write_complex_to_scc_file(filePath, in_cpx, degree, rivetCompatible, ignoreLastGenerators, false, reverse);
  auto out_cpx = build_complex_from_scc_file<Fil>(filePath, rivetCompatible, reverse, dimShift);

  write_slicer_to_scc_file(filePath, Slicer(in_cpx), degree, rivetCompatible, ignoreLastGenerators, false, reverse);
  Slicer out_slicer = build_slicer_from_scc_file<Slicer>(filePath, rivetCompatible, reverse, dimShift);

  BOOST_CHECK(out_slicer.get_boundaries() == out_cpx.get_boundaries());
  BOOST_CHECK(out_slicer.get_dimensions() == out_cpx.get_dimensions());
  BOOST_CHECK(out_slicer.get_filtration_values() == out_cpx.get_filtration_values());
  BOOST_CHECK_EQUAL(out_slicer.get_max_dimension(), out_cpx.get_max_dimension());

  // columns within a same dimension are reversed in order when reverse == false
  BC bc = {{}, {}, {}, {}, {}, {1, 3}, {1, 2}, {0, 1}, {0, 2}, {3, 4}, {0, 3}, {0, 4}, {6, 7, 8}, {9, 10, 11}};
  DC dc = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2};
  FC fc = {ini{0.144742, 0.2130470},
           ini{1.27611, 0.188994},
           ini{0.686347, 0.188994},
           ini{1.27032, 0.0729406},
           ini{0.636869, 0.0729406},
           ini{0.635158, 0.226328},
           ini{0.454557, 0.226328},
           ini{0.416981, 0.226328},
           ini{0.232049, 0.226328},
           ini{0.29683, 0.2130471},
           ini{0.151118, 0.2130471},
           ini{0.287974, 0.2130470},
           ini{0.82801, 0.226328},
           ini{0.639283, 0.226328}};

  BOOST_CHECK(bc == out_cpx.get_boundaries());
  BOOST_CHECK(dc == out_cpx.get_dimensions());
  BOOST_CHECK(fc == out_cpx.get_filtration_values());
  BOOST_CHECK_EQUAL(in_cpx.get_max_dimension(), out_cpx.get_max_dimension());

  rivetCompatible = true;

  write_complex_to_scc_file(filePath, in_cpx, degree, rivetCompatible, ignoreLastGenerators, false, reverse);
  out_cpx = build_complex_from_scc_file<Fil>(filePath, rivetCompatible, reverse, dimShift);

  BOOST_CHECK(bc == out_cpx.get_boundaries());
  BOOST_CHECK(dc == out_cpx.get_dimensions());
  BOOST_CHECK(fc == out_cpx.get_filtration_values());
  BOOST_CHECK_EQUAL(in_cpx.get_max_dimension(), out_cpx.get_max_dimension());

  rivetCompatible = false;
  reverse = true;

  write_complex_to_scc_file(filePath, in_cpx, degree, rivetCompatible, ignoreLastGenerators, false, reverse);
  out_cpx = build_complex_from_scc_file<Fil>(filePath, rivetCompatible, reverse, dimShift);

  BOOST_CHECK(in_cpx.get_boundaries() == out_cpx.get_boundaries());
  BOOST_CHECK(in_cpx.get_dimensions() == out_cpx.get_dimensions());
  BOOST_CHECK(in_cpx.get_filtration_values() == out_cpx.get_filtration_values());
  BOOST_CHECK_EQUAL(in_cpx.get_max_dimension(), out_cpx.get_max_dimension());

  ignoreLastGenerators = true;
  reverse = false;

  bc = {{}, {}, {}, {}, {}, {}, {}, {1, 2, 3}, {4, 5, 6}};
  dc = {1, 1, 1, 1, 1, 1, 1, 2, 2};
  fc = {ini{0.635158, 0.226328},
        ini{0.454557, 0.226328},
        ini{0.416981, 0.226328},
        ini{0.232049, 0.226328},
        ini{0.29683, 0.2130471},
        ini{0.151118, 0.2130471},
        ini{0.287974, 0.2130470},
        ini{0.82801, 0.226328},
        ini{0.639283, 0.226328}};

  write_complex_to_scc_file(filePath, in_cpx, degree, rivetCompatible, ignoreLastGenerators, false, reverse);
  out_cpx = build_complex_from_scc_file<Fil>(filePath, rivetCompatible, reverse, dimShift);

  BOOST_CHECK(bc == out_cpx.get_boundaries());
  BOOST_CHECK(dc == out_cpx.get_dimensions());
  BOOST_CHECK(fc == out_cpx.get_filtration_values());
  BOOST_CHECK_EQUAL(2, out_cpx.get_max_dimension());

  reverse = true;

  bc = {{}, {}, {}, {}, {}, {0, 4}, {1, 4}, {0, 1}, {2, 4}, {3, 4}, {2, 3}, {1, 3}};
  dc = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1};
  fc = {ini{0.636869, 0.0729406},
        ini{1.27032, 0.0729406},
        ini{0.686347, 0.188994},
        ini{1.27611, 0.188994},
        ini{0.144742, 0.2130470},
        ini{0.287974, 0.2130470},
        ini{0.151118, 0.2130471},
        ini{0.29683, 0.2130471},
        ini{0.232049, 0.226328},
        ini{0.416981, 0.226328},
        ini{0.454557, 0.226328},
        ini{0.635158, 0.226328}};

  write_complex_to_scc_file(filePath, in_cpx, degree, rivetCompatible, ignoreLastGenerators, false, reverse);
  out_cpx = build_complex_from_scc_file<Fil>(filePath, rivetCompatible, reverse, dimShift);

  BOOST_CHECK(bc == out_cpx.get_boundaries());
  BOOST_CHECK(dc == out_cpx.get_dimensions());
  BOOST_CHECK(fc == out_cpx.get_filtration_values());
  BOOST_CHECK_EQUAL(1, out_cpx.get_max_dimension());

  ignoreLastGenerators = false;
  reverse = false;
  degree = 0;
  dimShift = -1;

  bc = {{}, {}, {}, {}, {}, {1, 3}, {1, 2}, {0, 1}, {0, 2}, {3, 4}, {0, 3}, {0, 4}, {6, 7, 8}, {9, 10, 11}};
  dc = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2};
  fc = {ini{0.144742, 0.2130470},
        ini{1.27611, 0.188994},
        ini{0.686347, 0.188994},
        ini{1.27032, 0.0729406},
        ini{0.636869, 0.0729406},
        ini{0.635158, 0.226328},
        ini{0.454557, 0.226328},
        ini{0.416981, 0.226328},
        ini{0.232049, 0.226328},
        ini{0.29683, 0.2130471},
        ini{0.151118, 0.2130471},
        ini{0.287974, 0.2130470},
        ini{0.82801, 0.226328},
        ini{0.639283, 0.226328}};

  write_complex_to_scc_file(filePath, in_cpx, degree, rivetCompatible, ignoreLastGenerators, false, reverse);
  out_cpx = build_complex_from_scc_file<Fil>(filePath, rivetCompatible, reverse, dimShift);

  BOOST_CHECK(bc == out_cpx.get_boundaries());
  BOOST_CHECK(dc == out_cpx.get_dimensions());
  BOOST_CHECK(fc == out_cpx.get_filtration_values());
  BOOST_CHECK_EQUAL(in_cpx.get_max_dimension(), out_cpx.get_max_dimension());

  reverse = true;

  write_complex_to_scc_file(filePath, in_cpx, degree, rivetCompatible, ignoreLastGenerators, false, reverse);
  out_cpx = build_complex_from_scc_file<Fil>(filePath, rivetCompatible, reverse, dimShift);

  BOOST_CHECK(in_cpx.get_boundaries() == out_cpx.get_boundaries());
  BOOST_CHECK(in_cpx.get_dimensions() == out_cpx.get_dimensions());
  BOOST_CHECK(in_cpx.get_filtration_values() == out_cpx.get_filtration_values());
  BOOST_CHECK_EQUAL(in_cpx.get_max_dimension(), out_cpx.get_max_dimension());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(slicer_bitmap_io, Fil, list_of_tested_variants)
{
  using Complex = Multi_parameter_filtered_complex<Fil>;
  using FC = typename Complex::Filtration_value_container;
  using BC = typename Complex::Boundary_container;
  using DC = typename Complex::Dimension_container;
  using ini = std::initializer_list<typename Fil::value_type>;
  using Slicer = Slicer<Fil, Persistence_interface_cohomology<Fil> >;

  std::vector<Fil> vertexValues = {ini{0.636869, 0.0729406},
                                   ini{1.27032, 0.0729406},
                                   ini{0.686347, 0.188994},
                                   ini{1.27611, 0.188994},
                                   ini{0.144742, 0.2130470},
                                   ini{0.287974, 0.2130470},
                                   ini{0.151118, 0.2130471},
                                   ini{0.29683, 0.2130471},
                                   ini{0.232049, 0.226328}};
  std::vector<unsigned int> shape = {3, 3};

  auto cpx = build_complex_from_bitmap(vertexValues, shape);

  BC bc = {{},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {0, 1},
           {1, 2},
           {0, 3},
           {1, 4},
           {2, 5},
           {3, 4},
           {4, 5},
           {3, 6},
           {4, 7},
           {5, 8},
           {6, 7},
           {7, 8},
           {9, 11, 12, 14},
           {10, 12, 13, 15},
           {14, 16, 17, 19},
           {15, 17, 18, 20}};
  DC dc = {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2};
  FC fc = {ini{0.636869, 0.0729406},  // 0
           ini{1.27032, 0.0729406},   // 1
           ini{0.686347, 0.188994},   // 2
           ini{1.27611, 0.188994},    // 3
           ini{0.144742, 0.2130470},  // 4
           ini{0.287974, 0.2130470},  // 5
           ini{0.151118, 0.2130471},  // 6
           ini{0.29683, 0.2130471},   // 7
           ini{0.232049, 0.226328},   // 8
           ini{1.27032, 0.0729406},   // 0 1
           ini{1.27032, 0.188994},    // 1 2
           ini{1.27611, 0.188994},    // 0 3
           ini{1.27032, 0.2130470},   // 1 4
           ini{0.686347, 0.2130470},  // 2 5
           ini{1.27611, 0.2130470},   // 3 4
           ini{0.287974, 0.2130470},  // 4 5
           ini{1.27611, 0.2130471},   // 3 6
           ini{0.29683, 0.2130471},   // 4 7
           ini{0.287974, 0.226328},   // 5 8
           ini{0.29683, 0.2130471},   // 6 7
           ini{0.29683, 0.226328},    // 7 8
           ini{1.27611, 0.2130470},  ini{1.27032, 0.2130470}, ini{1.27611, 0.2130471}, ini{0.29683, 0.226328}};

  BOOST_CHECK(bc == cpx.get_boundaries());
  BOOST_CHECK(dc == cpx.get_dimensions());
  BOOST_CHECK(fc == cpx.get_filtration_values());
  BOOST_CHECK_EQUAL(2, cpx.get_max_dimension());

  auto slicer = build_slicer_from_bitmap<Slicer>(vertexValues, shape);

  BOOST_CHECK(bc == slicer.get_boundaries());
  BOOST_CHECK(dc == slicer.get_dimensions());
  BOOST_CHECK(fc == slicer.get_filtration_values());
  BOOST_CHECK_EQUAL(2, slicer.get_max_dimension());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(slicer_simplex_tree_io, Fil, list_of_tested_variants)
{
  using ST = Gudhi::Simplex_tree<Simplex_tree_options_multidimensional_filtration<Fil>>;
  using Complex = Multi_parameter_filtered_complex<Fil>;
  using FC = typename Complex::Filtration_value_container;
  using BC = typename Complex::Boundary_container;
  using DC = typename Complex::Dimension_container;
  using ini = std::initializer_list<typename Fil::value_type>;
  using Slicer1 = Slicer<Fil, Persistence_interface_cohomology<Fil> >;
  using Fil2 = Multi_parameter_filtration<int>; // Will always be different from ST::Filtration_value
  using Slicer2 = Slicer<Fil2, Persistence_interface_cohomology<Fil2> >;

  ST simplexTree;

  simplexTree.insert_simplex_and_subfaces({0,1,2}, ini{0, 3});
  simplexTree.insert_simplex_and_subfaces({1,3}, ini{0, 4});
  simplexTree.insert_simplex_and_subfaces({4,5}, ini{0, 6});
  simplexTree.insert_simplex_and_subfaces({3,4,5,6}, ini{0, 5});
  simplexTree.insert_simplex_and_subfaces({2,6}, ini{0, 7});
  simplexTree.insert_simplex_and_subfaces({3,4}, ini{0, 8});
  simplexTree.insert_simplex_and_subfaces({0,1,2}, ini{0, 2});
  simplexTree.insert_simplex_and_subfaces({4,5,6}, ini{0, 4});
  simplexTree.insert_simplex_and_subfaces({2,6}, ini{0, 1});
  simplexTree.insert_simplex_and_subfaces({1,3}, ini{0, 8});

  auto cpx = build_complex_from_simplex_tree<Fil>(simplexTree);

  BC bc = {{},     {},     {},        {},           {},           {},           {},           {0, 1},
           {0, 2}, {1, 2}, {1, 3},    {2, 6},       {3, 4},       {3, 5},       {3, 6},       {4, 5},
           {4, 6}, {5, 6}, {7, 8, 9}, {12, 13, 15}, {12, 14, 16}, {13, 14, 17}, {15, 16, 17}, {19, 20, 21, 22}};
  DC dc = {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3};
  FC fc = {ini{0, 2}, ini{0, 2}, ini{0, 1}, ini{0, 4}, ini{0, 4}, ini{0, 4}, ini{0, 1}, ini{0, 2},
           ini{0, 2}, ini{0, 2}, ini{0, 4}, ini{0, 1}, ini{0, 5}, ini{0, 5}, ini{0, 5}, ini{0, 4},
           ini{0, 4}, ini{0, 4}, ini{0, 2}, ini{0, 5}, ini{0, 5}, ini{0, 5}, ini{0, 4}, ini{0, 5}};

  BOOST_CHECK(bc == cpx.get_boundaries());
  BOOST_CHECK(dc == cpx.get_dimensions());
  BOOST_CHECK(fc == cpx.get_filtration_values());
  BOOST_CHECK_EQUAL(3, cpx.get_max_dimension());

  auto slicer = build_slicer_from_simplex_tree<Slicer1>(simplexTree);

  BOOST_CHECK(bc == slicer.get_boundaries());
  BOOST_CHECK(dc == slicer.get_dimensions());
  BOOST_CHECK(fc == slicer.get_filtration_values());
  BOOST_CHECK_EQUAL(3, slicer.get_max_dimension());

  auto slicer2 = build_slicer_from_simplex_tree<Slicer2>(simplexTree);

  BOOST_CHECK(bc == slicer2.get_boundaries());
  BOOST_CHECK(dc == slicer2.get_dimensions());
  BOOST_CHECK_EQUAL(slicer.get_filtration_values().size(), slicer2.get_filtration_values().size());
  BOOST_CHECK_EQUAL(3, slicer2.get_max_dimension());
}


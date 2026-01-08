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

#include <gudhi/Multi_parameter_filtered_complex.h>
#include <gudhi/Multi_parameter_filtration.h>
#include <gudhi/Dynamic_multi_parameter_filtration.h>
#include <gudhi/Projective_cover_kernel.h>

using Gudhi::multi_filtration::Dynamic_multi_parameter_filtration;
using Gudhi::multi_filtration::Multi_parameter_filtration;
using Gudhi::multi_persistence::Multi_parameter_filtered_complex;
using Gudhi::multi_persistence::Projective_cover_kernel;

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

  BC bc = {{1},      {1},      {3},  {3},      {5},  {5},      {5, 7}, {7},  {10}, {13},         {10},
           {1},      {1, 3},   {13}, {18},     {19}, {13},     {1, 3}, {18}, {19}, {18},         {22},
           {22},     {10, 24}, {24}, {24, 26}, {26}, {26, 28}, {28},   {36}, {38}, {28, 36, 39}, {24, 26},
           {28, 40}, {22},     {37}, {36},     {37}, {38},     {39},   {40}};

  DC dc = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

  FC fc = {ini{0.636869, 0.0729406}, ini{1.27032, 0.0729406},  ini{0.686347, 0.188994},  ini{1.27611, 0.188994},
           ini{0.144742, 0.2130470}, ini{0.287974, 0.2130470}, ini{0.151118, 0.2130471}, ini{0.29683, 0.2130471},
           ini{0.232049, 0.226328},  ini{0.416981, 0.226328},  ini{0.454557, 0.226328},  ini{0.635158, 0.226328},
           ini{0.639283, 0.226328},  ini{0.82801, 0.226328},   ini{0.381793, 0.2263281}, ini{0.408723, 0.2263281},
           ini{0.414005, 0.2263281}, ini{0.638054, 0.2263281}, ini{0.755329, 0.2263281}, ini{0.809926, 0.2263281},
           ini{0.377664, 0.228439},  ini{0.495411, 0.228439},  ini{0.989342, 0.228439},  ini{0.227305, 0.2284391},
           ini{0.399651, 0.2284391}, ini{0.461108, 0.2284391}, ini{0.911316, 0.2284391}, ini{0.456372, 0.2284392},
           ini{0.868285, 0.2284392}, ini{0.281964, 0.2284393}, ini{0.361102, 0.2284393}, ini{0.434142, 0.2284393},
           ini{0.455947, 0.2284393}, ini{0.48357, 0.2284393},  ini{0.494671, 0.2284393}, ini{0.511928, 0.2284393},
           ini{0.524031, 0.2284393}, ini{0.67094, 0.2284393},  ini{0.711337, 0.2284393}, ini{0.857271, 0.2284393},
           ini{0.967129, 0.2284393}};

  return Complex(bc, dc, fc);
}

template <class Fil>
Multi_parameter_filtered_complex<Fil> build_output_complex()
{
  using Complex = Multi_parameter_filtered_complex<Fil>;
  using FC = typename Complex::Filtration_value_container;
  using BC = typename Complex::Boundary_container;
  using DC = typename Complex::Dimension_container;
  using ini = std::initializer_list<typename Fil::value_type>;

  BC bc = {{},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {},
           {4, 5},
           {4, 6, 7},
           {14, 20},
           {8, 23, 24},
           {9, 16},
           {8, 10},
           {25, 32},
           {21, 34},
           {29, 36},
           {0, 11},
           {12, 17},
           {35, 37},
           {0, 2, 12},
           {30, 38},
           {14, 18},
           {15, 19},
           {9, 13},
           {8, 23, 25, 27, 29, 31, 39},
           {8, 23, 25, 27, 28},
           {8, 23, 25, 26},
           {8, 23, 25, 27, 33, 40},
           {21, 22},
           {0, 1},
           {2, 3}};

  DC dc = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};

  FC fc = {
      ini{0.636869, 0.0729406}, ini{1.27032, 0.0729406},  ini{0.686347, 0.188994},  ini{1.27611, 0.188994},
      ini{0.144742, 0.213047},  ini{0.287974, 0.213047},  ini{0.151118, 0.2130471}, ini{0.29683, 0.2130471},
      ini{0.232049, 0.226328},  ini{0.416981, 0.226328},  ini{0.454557, 0.226328},  ini{0.635158, 0.226328},
      ini{0.639283, 0.226328},  ini{0.82801, 0.226328},   ini{0.381793, 0.2263281}, ini{0.408723, 0.2263281},
      ini{0.414005, 0.2263281}, ini{0.638054, 0.2263281}, ini{0.755329, 0.2263281}, ini{0.809926, 0.2263281},
      ini{0.377664, 0.228439},  ini{0.495411, 0.228439},  ini{0.989342, 0.228439},  ini{0.227305, 0.2284391},
      ini{0.399651, 0.2284391}, ini{0.461108, 0.2284391}, ini{0.911316, 0.2284391}, ini{0.456372, 0.2284392},
      ini{0.868285, 0.2284392}, ini{0.281964, 0.2284393}, ini{0.361102, 0.2284393}, ini{0.434142, 0.2284393},
      ini{0.455947, 0.2284393}, ini{0.48357, 0.2284393},  ini{0.494671, 0.2284393}, ini{0.511928, 0.2284393},
      ini{0.524031, 0.2284393}, ini{0.67094, 0.2284393},  ini{0.711337, 0.2284393}, ini{0.857271, 0.2284393},
      ini{0.967129, 0.2284393}, ini{0.287974, 0.213047},  ini{0.29683, 0.2130471},  ini{0.381793, 0.228439},
      ini{0.399651, 0.2284391}, ini{0.416981, 0.2263281}, ini{0.454557, 0.226328},  ini{0.461108, 0.2284393},
      ini{0.495411, 0.2284393}, ini{0.524031, 0.2284393}, ini{0.636869, 0.226328},  ini{0.639283, 0.2263281},
      ini{0.67094, 0.2284393},  ini{0.686347, 0.226328},  ini{0.711337, 0.2284393}, ini{0.755329, 0.2263281},
      ini{0.809926, 0.2263281}, ini{0.82801, 0.226328},   ini{0.857271, 0.2284393}, ini{0.868285, 0.2284392},
      ini{0.911316, 0.2284391}, ini{0.967129, 0.2284393}, ini{0.989342, 0.228439},  ini{1.27032, 0.0729406},
      ini{1.27611, 0.188994},
  };

  return Complex(bc, dc, fc);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(projective_cover_kernel, Fil, list_of_tested_variants)
{
  using Complex = Multi_parameter_filtered_complex<Fil>;
  using BC = typename Complex::Boundary_container;

  Complex input = build_input_complex<Fil>();
  Complex output = build_output_complex<Fil>();

  Projective_cover_kernel kernel(input, 0);

  BC gen = kernel.build_generators();
  BC realGen(output.get_boundaries().begin() + 41, output.get_boundaries().end());
  BOOST_CHECK(gen == realGen);

  Complex res = kernel.create_complex();
  BOOST_CHECK(res.get_boundaries() == output.get_boundaries());
  BOOST_CHECK(res.get_dimensions() == output.get_dimensions());
  BOOST_CHECK(res.get_filtration_values() == output.get_filtration_values());
  BOOST_CHECK_EQUAL(res.get_max_dimension(), output.get_max_dimension());
}

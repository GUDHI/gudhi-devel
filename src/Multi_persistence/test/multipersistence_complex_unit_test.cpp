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
#include <utility>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "multi_persistence"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/Multi_parameter_filtered_complex.h>
#include <gudhi/Multi_parameter_filtration.h>
#include <gudhi/Dynamic_multi_parameter_filtration.h>

using Gudhi::multi_filtration::Dynamic_multi_parameter_filtration;
using Gudhi::multi_filtration::Multi_parameter_filtration;
using Gudhi::multi_persistence::Multi_parameter_filtered_complex;

using I = std::uint32_t;
using D = int;

using list_of_tested_variants = boost::mpl::list<Multi_parameter_filtration<double>,
                                                 Dynamic_multi_parameter_filtration<double>,
                                                 Multi_parameter_filtration<int>,
                                                 Dynamic_multi_parameter_filtration<int> >;

BOOST_AUTO_TEST_CASE_TEMPLATE(multi_complex_constructors, Fil, list_of_tested_variants)
{
  using Complex = Multi_parameter_filtered_complex<Fil, I, D>;
  using FC = typename Complex::Filtration_value_container;
  using BC = typename Complex::Boundary_container;
  using DC = typename Complex::Dimension_container;
  using ini = std::initializer_list<typename Fil::value_type>;

  Multi_parameter_filtered_complex<Fil, I, D> emptyC;
  BOOST_CHECK_EQUAL(emptyC.get_number_of_cycle_generators(), 0);
  BOOST_CHECK_EQUAL(emptyC.get_number_of_parameters(), 0);
  BOOST_CHECK(emptyC.is_ordered_by_dimension());
  BOOST_CHECK_EQUAL(emptyC.get_filtration_values().size(), 0);
  BOOST_CHECK_EQUAL(emptyC.get_dimensions().size(), 0);
  BOOST_CHECK_EQUAL(emptyC.get_boundaries().size(), 0);
  BOOST_CHECK_EQUAL(emptyC.get_max_dimension(), -1);

  BC bc = {{}, {}, {}, {0, 1}, {0, 2}, {3, 4}};
  DC dc = {0, 0, 0, 1, 1, 2};
  FC fc = {ini{0, 1, 2}, ini{0, 1, 2}, ini{0, 1, 2}, ini{3, 4, 5}, ini{3, 4, 5}, ini{6, 7, 8}};

  Multi_parameter_filtered_complex<Fil, I, D> copyC(bc, dc, fc);
  BOOST_CHECK_EQUAL(copyC.get_number_of_cycle_generators(), 6);
  BOOST_CHECK_EQUAL(copyC.get_number_of_parameters(), 3);
  BOOST_CHECK(copyC.is_ordered_by_dimension());
  BOOST_CHECK_EQUAL(copyC.get_filtration_values().size(), 6);
  BOOST_CHECK_EQUAL(copyC.get_dimensions().size(), 6);
  BOOST_CHECK_EQUAL(copyC.get_boundaries().size(), 6);
  BOOST_CHECK_EQUAL(copyC.get_max_dimension(), 2);
  // BOOST_CHECK_EQUAL(copyC.get_dimension(0), 0);
  // BOOST_CHECK_EQUAL(copyC.get_dimension(1), 0);
  // BOOST_CHECK_EQUAL(copyC.get_dimension(2), 0);
  // BOOST_CHECK_EQUAL(copyC.get_dimension(3), 1);
  // BOOST_CHECK_EQUAL(copyC.get_dimension(4), 1);
  // BOOST_CHECK_EQUAL(copyC.get_dimension(5), 2);

  bc = {{}, {}, {0, 1}, {}, {0, 3}, {2, 4}};
  dc = {0, 0, 1, 0, 1, 2};
  fc = {ini{0, 1, 2}, ini{0, 1, 2}, ini{3, 4, 5}, ini{0, 1, 2}, ini{3, 4, 5}, ini{6, 7, 8}};

  Multi_parameter_filtered_complex<Fil, I, D> moveC(std::move(bc), std::move(dc), std::move(fc));
  BOOST_CHECK_EQUAL(moveC.get_number_of_cycle_generators(), 6);
  BOOST_CHECK_EQUAL(moveC.get_number_of_parameters(), 3);
  BOOST_CHECK(!moveC.is_ordered_by_dimension());
  BOOST_CHECK_EQUAL(moveC.get_filtration_values().size(), 6);
  BOOST_CHECK_EQUAL(moveC.get_dimensions().size(), 6);
  BOOST_CHECK_EQUAL(moveC.get_boundaries().size(), 6);
  BOOST_CHECK_EQUAL(moveC.get_max_dimension(), 2);
  // BOOST_CHECK_EQUAL(moveC.get_dimension(0), 0);
  // BOOST_CHECK_EQUAL(moveC.get_dimension(1), 0);
  // BOOST_CHECK_EQUAL(moveC.get_dimension(2), 1);
  // BOOST_CHECK_EQUAL(moveC.get_dimension(3), 0);
  // BOOST_CHECK_EQUAL(moveC.get_dimension(4), 1);
  // BOOST_CHECK_EQUAL(moveC.get_dimension(5), 2);

  Multi_parameter_filtered_complex<Multi_parameter_filtration<long int>, I, D> copyCC(copyC);
  BOOST_CHECK_EQUAL(copyCC.get_number_of_cycle_generators(), 6);
  BOOST_CHECK_EQUAL(copyCC.get_number_of_parameters(), 3);
  BOOST_CHECK(copyCC.is_ordered_by_dimension());
  BOOST_CHECK_EQUAL(copyCC.get_filtration_values().size(), 6);
  BOOST_CHECK_EQUAL(copyCC.get_dimensions().size(), 6);
  BOOST_CHECK_EQUAL(copyCC.get_boundaries().size(), 6);
  BOOST_CHECK_EQUAL(copyCC.get_max_dimension(), 2);
  // BOOST_CHECK_EQUAL(copyCC.get_dimension(0), 0);
  // BOOST_CHECK_EQUAL(copyCC.get_dimension(1), 0);
  // BOOST_CHECK_EQUAL(copyCC.get_dimension(2), 0);
  // BOOST_CHECK_EQUAL(copyCC.get_dimension(3), 1);
  // BOOST_CHECK_EQUAL(copyCC.get_dimension(4), 1);
  // BOOST_CHECK_EQUAL(copyCC.get_dimension(5), 2);

  Multi_parameter_filtered_complex<Multi_parameter_filtration<long int>, I, D> copyCC2 = copyC;
  BOOST_CHECK_EQUAL(copyCC2.get_number_of_cycle_generators(), 6);
  BOOST_CHECK_EQUAL(copyCC2.get_number_of_parameters(), 3);
  BOOST_CHECK(copyCC2.is_ordered_by_dimension());
  BOOST_CHECK_EQUAL(copyCC2.get_filtration_values().size(), 6);
  BOOST_CHECK_EQUAL(copyCC2.get_dimensions().size(), 6);
  BOOST_CHECK_EQUAL(copyCC2.get_boundaries().size(), 6);
  BOOST_CHECK_EQUAL(copyCC2.get_max_dimension(), 2);
  // BOOST_CHECK_EQUAL(copyCC2.get_dimension(0), 0);
  // BOOST_CHECK_EQUAL(copyCC2.get_dimension(1), 0);
  // BOOST_CHECK_EQUAL(copyCC2.get_dimension(2), 0);
  // BOOST_CHECK_EQUAL(copyCC2.get_dimension(3), 1);
  // BOOST_CHECK_EQUAL(copyCC2.get_dimension(4), 1);
  // BOOST_CHECK_EQUAL(copyCC2.get_dimension(5), 2);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(multi_complex_sorts, Fil, list_of_tested_variants)
{
  using Complex = Multi_parameter_filtered_complex<Fil, I, D>;
  using FC = typename Complex::Filtration_value_container;
  using BC = typename Complex::Boundary_container;
  using DC = typename Complex::Dimension_container;
  using Index = typename Complex::Index;
  using ini = std::initializer_list<typename Fil::value_type>;

  BC bc = {{2, 6}, {}, {1, 3}, {}, {}, {0}, {1, 4}};
  DC dc = {2, 0, 1, 0, 0, 3, 1};
  FC fc = {ini{3, 4, 5}, ini{0, 1, 2}, ini{0, 2, 2}, ini{0, 2, 1}, ini{1, 2, 3}, ini{6, 4, 5}, ini{1, 3, 3}};

  Multi_parameter_filtered_complex<Fil, I, D> cpx(bc, dc, fc);
  BOOST_CHECK(!cpx.is_ordered_by_dimension());

  bc = {{}, {}, {}, {0, 1}, {1, 2}, {3, 4}, {5}};
  dc = {0, 0, 0, 1, 1, 2, 3};
  fc = {ini{0, 2, 1}, ini{0, 1, 2}, ini{1, 2, 3}, ini{0, 2, 2}, ini{1, 3, 3}, ini{3, 4, 5}, ini{6, 4, 5}};

  Multi_parameter_filtered_complex<Fil, I, D> byDimC(cpx);
  byDimC.sort_by_dimension_co_lexicographically();
  BOOST_CHECK(byDimC.is_ordered_by_dimension());
  BOOST_CHECK(byDimC.get_boundaries() == bc);
  BOOST_CHECK(byDimC.get_dimensions() == dc);
  BOOST_CHECK(byDimC.get_filtration_values() == fc);

  bc = {{}, {}, {0, 1}, {}, {0, 3}, {2, 4}, {5}};
  dc = {0, 0, 1, 0, 1, 2, 3};
  fc = {ini{0, 1, 2}, ini{0, 2, 1}, ini{0, 2, 2}, ini{1, 2, 3}, ini{1, 3, 3}, ini{3, 4, 5}, ini{6, 4, 5}};

  Multi_parameter_filtered_complex<Fil, I, D> lexC(cpx);
  lexC.sort([&lexC](Index i, Index j) -> bool {
    const auto& filtrationValues = lexC.get_filtration_values();
    return is_strict_less_than_lexicographically(filtrationValues[i], filtrationValues[j]);
  });
  BOOST_CHECK(!lexC.is_ordered_by_dimension());
  BOOST_CHECK(lexC.get_boundaries() == bc);
  BOOST_CHECK(lexC.get_dimensions() == dc);
  BOOST_CHECK(lexC.get_filtration_values() == fc);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(multi_complex_other, Fil, list_of_tested_variants)
{
  using Complex = Multi_parameter_filtered_complex<Fil, I, D>;
  using FC = typename Complex::Filtration_value_container;
  using BC = typename Complex::Boundary_container;
  using DC = typename Complex::Dimension_container;
  using Index = typename Complex::Index;
  using ini = std::initializer_list<typename Fil::value_type>;

  BC bc = {{}, {}, {}, {0, 1}, {1, 2}, {3, 4}, {5}};
  DC dc = {0, 0, 0, 1, 1, 2, 3};
  FC fc = {ini{0, 2, 1}, ini{0, 1, 2}, ini{1, 2, 3}, ini{0, 2, 2}, ini{1, 3, 3}, ini{3, 4, 5}, ini{6, 4, 5}};

  Multi_parameter_filtered_complex<Fil, I, D> cpx(bc, dc, fc);
  Multi_parameter_filtered_complex<Fil, I, D> cpxToGrid(bc, dc, fc);
  Multi_parameter_filtered_complex<Fil, I, D> cpxToGridCoord(bc, dc, fc);

  std::vector<std::vector<typename Fil::value_type> > grid = {{0, 2, 4, 8}, {0, 3, 6, 9}, {0, 4, 8, 16}};
  fc = {ini{0, 3, 0}, ini{0, 0, 4}, ini{2, 3, 4}, ini{0, 3, 4}, ini{2, 3, 4}, ini{4, 3, 4}, ini{8, 3, 4}};

  cpxToGrid.coarsen_on_grid(grid, false);
  BOOST_CHECK(cpxToGrid.is_ordered_by_dimension());
  BOOST_CHECK(cpxToGrid.get_boundaries() == bc);
  BOOST_CHECK(cpxToGrid.get_dimensions() == dc);
  BOOST_CHECK(cpxToGrid.get_filtration_values() == fc);
  BOOST_CHECK_EQUAL(cpxToGrid.get_max_dimension(), 3);

  fc = {ini{0, 1, 0}, ini{0, 0, 1}, ini{1, 1, 1}, ini{0, 1, 1}, ini{1, 1, 1}, ini{2, 1, 1}, ini{3, 1, 1}};

  cpxToGridCoord.coarsen_on_grid(grid);
  BOOST_CHECK(cpxToGridCoord.is_ordered_by_dimension());
  BOOST_CHECK(cpxToGridCoord.get_boundaries() == bc);
  BOOST_CHECK(cpxToGridCoord.get_dimensions() == dc);
  BOOST_CHECK(cpxToGridCoord.get_filtration_values() == fc);
  BOOST_CHECK_EQUAL(cpxToGridCoord.get_max_dimension(), 3);

  bc = {{}, {}, {}, {0, 1}, {1, 2}};
  dc = {0, 0, 0, 1, 1};
  fc = {ini{0, 2, 1}, ini{0, 1, 2}, ini{1, 2, 3}, ini{0, 2, 2}, ini{1, 3, 3}};

  BOOST_CHECK(cpx.is_ordered_by_dimension());
  Index idx = cpx.prune_above_dimension(1);
  BOOST_CHECK(cpx.is_ordered_by_dimension());
  BOOST_CHECK(cpx.get_boundaries() == bc);
  BOOST_CHECK(cpx.get_dimensions() == dc);
  BOOST_CHECK(cpx.get_filtration_values() == fc);
  BOOST_CHECK_EQUAL(cpx.get_max_dimension(), 1);
  BOOST_CHECK_EQUAL(idx, 5);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(multi_complex_friend, Fil, list_of_tested_variants)
{
  using Complex = Multi_parameter_filtered_complex<Fil, I, D>;
  using FC = typename Complex::Filtration_value_container;
  using BC = typename Complex::Boundary_container;
  using DC = typename Complex::Dimension_container;
  using ini = std::initializer_list<typename Fil::value_type>;

  BC bc = {{2, 6}, {}, {1, 3}, {}, {}, {0}, {1, 4}};
  DC dc = {2, 0, 1, 0, 0, 3, 1};
  FC fc = {ini{3, 4, 5}, ini{0, 1, 2}, ini{0, 2, 2}, ini{0, 2, 1}, ini{1, 2, 3}, ini{6, 4, 5}, ini{1, 3, 3}};

  Multi_parameter_filtered_complex<Fil, I, D> cpx(bc, dc, fc);
  BOOST_CHECK(!cpx.is_ordered_by_dimension());

  bc = {{}, {}, {}, {0, 1}, {1, 2}, {3, 4}, {5}};
  dc = {0, 0, 0, 1, 1, 2, 3};
  fc = {ini{0, 2, 1}, ini{0, 1, 2}, ini{1, 2, 3}, ini{0, 2, 2}, ini{1, 3, 3}, ini{3, 4, 5}, ini{6, 4, 5}};

  auto [byDimC, perm] = build_permuted_complex(cpx);
  BOOST_CHECK(byDimC.is_ordered_by_dimension());
  BOOST_CHECK(byDimC.get_boundaries() == bc);
  BOOST_CHECK(byDimC.get_dimensions() == dc);
  BOOST_CHECK(byDimC.get_filtration_values() == fc);

  Multi_parameter_filtered_complex<Fil, I, D> byDimC2 = build_permuted_complex(cpx, perm);
  BOOST_CHECK(byDimC2.is_ordered_by_dimension());
  BOOST_CHECK(byDimC2.get_boundaries() == bc);
  BOOST_CHECK(byDimC2.get_dimensions() == dc);
  BOOST_CHECK(byDimC2.get_filtration_values() == fc);

  std::vector<std::vector<typename Fil::value_type> > grid = {{0, 2, 4, 8}, {0, 3, 6, 9}, {0, 4, 8, 16}};

  auto cpxToGridCoord = build_complex_coarsen_on_grid(byDimC, grid);

  using Complex2 = decltype(cpxToGridCoord);
  using FC2 = typename Complex2::Filtration_value_container;
  using ini2 = std::initializer_list<std::int32_t>;

  FC2 fc2 = {ini2{0, 1, 0}, ini2{0, 0, 1}, ini2{1, 1, 1}, ini2{0, 1, 1}, ini2{1, 1, 1}, ini2{2, 1, 1}, ini2{3, 1, 1}};

  BOOST_CHECK(cpxToGridCoord.is_ordered_by_dimension());
  BOOST_CHECK(cpxToGridCoord.get_boundaries() == bc);
  BOOST_CHECK(cpxToGridCoord.get_dimensions() == dc);
  BOOST_CHECK(cpxToGridCoord.get_filtration_values() == fc2);
  BOOST_CHECK_EQUAL(cpxToGridCoord.get_max_dimension(), 3);
}

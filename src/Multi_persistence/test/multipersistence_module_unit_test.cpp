/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

// #include <array>
// #include <initializer_list>
// #include <limits>
// #include <utility>
// #include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "multi_persistence"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/Multi_persistence/Module.h>

using Gudhi::multi_persistence::Module;

using I = std::uint32_t;
using D = int;
using T = double;
// using Bar = Gudhi::persistence_matrix::Persistence_interval<int, T>;
// using Barcode = std::vector<Bar>;
// using Flat_barcode = std::vector<std::array<T, 2>>;
// using Multi_dimensional_barcode = std::vector<Barcode>;
// using Multi_dimensional_flat_barcode = std::vector<Flat_barcode>;
// using Test_barcode = std::vector<std::pair<T, T>>;
// using Test_multi_dimensional_barcode = std::vector<Test_barcode>;

using list_of_tested_variants = boost::mpl::list<>;

BOOST_AUTO_TEST_CASE_TEMPLATE(module_constructors, Module_t, list_of_tested_variants)
{
}


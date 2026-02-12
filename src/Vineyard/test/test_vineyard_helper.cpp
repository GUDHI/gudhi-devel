/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "vineyard"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/vineyard_helper.h>

#include "vy_test_utilities.h"

using option_list = boost::mpl::list<Chain_vineyard_options, RU_vineyard_options>;

BOOST_AUTO_TEST_CASE_TEMPLATE(vy_initialization, Option, option_list) {
  
}

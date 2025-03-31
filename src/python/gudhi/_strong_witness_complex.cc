/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Simplex_tree.h>
#include <gudhi/Strong_witness_complex.h>

#include "Simplex_tree_interface.h"
#include "Simplex_tree_interface_cython.h"

#include <utility> // std::pair
#include <vector>

namespace Gudhi {

class Simplex_tree_interface;

namespace witness_complex {

// /////////////////////////////////////////////////////////////////////////////
//  Strong_witness_complex_interface declaration
// /////////////////////////////////////////////////////////////////////////////

class Strong_witness_complex_interface
{
public:
    using Nearest_landmark_range = std::vector<std::pair<std::size_t, double>>;
    using Nearest_landmark_table = std::vector<Nearest_landmark_range>;

public:
    Strong_witness_complex_interface(const Nearest_landmark_table& nlt);
    ~Strong_witness_complex_interface();

    void create_simplex_tree(Simplex_tree_interface* simplex_tree, double  max_alpha_square, std::size_t limit_dimension);
    void create_simplex_tree(Simplex_tree_interface* simplex_tree, double  max_alpha_square);

private:
    Strong_witness_complex<Nearest_landmark_table>* witness_complex_;
};

// /////////////////////////////////////////////////////////////////////////////
//  Strong_witness_complex_interface definition
// /////////////////////////////////////////////////////////////////////////////

Strong_witness_complex_interface::Strong_witness_complex_interface(const Nearest_landmark_table& nlt)
{
    witness_complex_ = new Strong_witness_complex<Nearest_landmark_table>(nlt);
}

Strong_witness_complex_interface::~Strong_witness_complex_interface()
{
    delete witness_complex_;
}

void Strong_witness_complex_interface::create_simplex_tree(Simplex_tree_interface* simplex_tree,
                                                           double  max_alpha_square,
                                                           std::size_t limit_dimension)
{
    witness_complex_->create_complex(*simplex_tree, max_alpha_square, limit_dimension);
}

void Strong_witness_complex_interface::create_simplex_tree(Simplex_tree_interface* simplex_tree, double  max_alpha_square)
{
    witness_complex_->create_complex(*simplex_tree, max_alpha_square);
}

}  // namespace witness_complex

}  // namespace Gudhi

// /////////////////////////////////////////////////////////////////////////////
// Strong_witness_complex_interface wrapping
// /////////////////////////////////////////////////////////////////////////////

#include <nanobind/nanobind.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;

namespace gwc = Gudhi::witness_complex;
using gwci = gwc::Strong_witness_complex_interface;


NB_MODULE(_strong_witness_complex_ext, m) {
  m.attr("__license__") = "GPL v3";

  nb::class_<gwci>(m, "Strong_witness_complex_interface")
      .def(nb::init<const gwci::Nearest_landmark_table&>(), "Constructor")
      .def("create_simplex_tree", nb::overload_cast<Gudhi::Simplex_tree_interface*, double, std::size_t>(&gwci::create_simplex_tree), "")
      .def("create_simplex_tree", nb::overload_cast<Gudhi::Simplex_tree_interface*, double>(&gwci::create_simplex_tree), "");

}

//
// _strong_witness_complex.cc ends here

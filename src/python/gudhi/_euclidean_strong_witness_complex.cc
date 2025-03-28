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
#include <gudhi/Euclidean_strong_witness_complex.h>

#include "Simplex_tree_interface.h"

#include <CGAL/Epick_d.h>

#include <vector>
#include <utility>  // std::pair
#include <iostream>
#include <cstddef>

namespace Gudhi {

namespace witness_complex {

// /////////////////////////////////////////////////////////////////////////////
// Euclidean_strong_witness_complex_interface declaration
// /////////////////////////////////////////////////////////////////////////////

class Euclidean_strong_witness_complex_interface
{
    using Dynamic_kernel = CGAL::Epick_d< CGAL::Dynamic_dimension_tag >;
    using Point_d = Dynamic_kernel::Point_d;
    using Simplex_key = typename Simplex_tree<>::Simplex_key;

public:
     Euclidean_strong_witness_complex_interface(const std::vector<std::vector<double>>& landmarks,
                                                const std::vector<std::vector<double>>& witnesses);
    ~Euclidean_strong_witness_complex_interface();

    void create_simplex_tree(Simplex_tree_interface* simplex_tree, double max_alpha_square, std::size_t limit_dimension);
    void create_simplex_tree(Simplex_tree_interface* simplex_tree, double max_alpha_square);
    std::vector<double> get_point(unsigned vh);

private:
    std::vector<Point_d> landmarks_;
    Euclidean_strong_witness_complex<Dynamic_kernel>* witness_complex_;
};

// /////////////////////////////////////////////////////////////////////////////
// Euclidean_strong_witness_complex_interface definition
// /////////////////////////////////////////////////////////////////////////////

Euclidean_strong_witness_complex_interface::Euclidean_strong_witness_complex_interface(const std::vector<std::vector<double>>& landmarks,
                                                                                       const std::vector<std::vector<double>>& witnesses)
{
    landmarks_.reserve(landmarks.size());
    for (auto& landmark : landmarks)
        landmarks_.emplace_back(landmark.begin(), landmark.end());
    witness_complex_ = new Euclidean_strong_witness_complex<Dynamic_kernel>(landmarks_, witnesses);
}

Euclidean_strong_witness_complex_interface::~Euclidean_strong_witness_complex_interface()
{
    delete witness_complex_;
}

void Euclidean_strong_witness_complex_interface::create_simplex_tree(Simplex_tree_interface* simplex_tree, double max_alpha_square, std::size_t limit_dimension)
{
    witness_complex_->create_complex(*simplex_tree, max_alpha_square, limit_dimension);
}

void Euclidean_strong_witness_complex_interface::create_simplex_tree(Simplex_tree_interface* simplex_tree, double max_alpha_square)
{
    witness_complex_->create_complex(*simplex_tree, max_alpha_square);
}

std::vector<double> Euclidean_strong_witness_complex_interface::get_point(unsigned vh)
{
    std::vector<double> vd;
    if (vh < landmarks_.size()) {
        Point_d ph = witness_complex_->get_point(vh);
        for (auto coord = ph.cartesian_begin(); coord < ph.cartesian_end(); ++coord) {
            vd.push_back(*coord);
        }
    }
    return vd;
}

}  // namespace witness_complex

}  // namespace Gudhi

// /////////////////////////////////////////////////////////////////////////////
// Euclidean_strong_witness_complex_interface wrapping
// /////////////////////////////////////////////////////////////////////////////

#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;
namespace gwc = Gudhi::witness_complex;
using gwci = gwc::Euclidean_strong_witness_complex_interface;

NB_MODULE(_euclidean_strong_witness_complex_ext, m) {
  m.attr("__license__") = "GPL v3";

  nb::class_<gwci>(m, "Euclidean_strong_witness_complex_interface")
      .def(nb::init<const std::vector<std::vector<double>>&, const std::vector<std::vector<double>>&>(),
          "Constructs strong witness complex for given sets of witnesses and landmarks in Euclidean space.")
      .def("create_simplex_tree", nb::overload_cast<Gudhi::Simplex_tree_interface*, double>(&gwci::create_simplex_tree), "")
      .def("create_simplex_tree", nb::overload_cast<Gudhi::Simplex_tree_interface*, double, std::size_t>(&gwci::create_simplex_tree), "")
      .def("get_point", &gwci::get_point,
R"doc(
    This function returns the point corresponding to a given vertex.
    Args:
        arg (int): The vertex id.
    Returns:
        list of float: The coordinates of the vertex.
)doc")
      ;
}

//
// _euclidean_strong_witness_complex.cc ends here

//
// _euclidean_strong_witness_complex.cc ends here

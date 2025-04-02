/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <Simplex_tree_interface.h>

#include <gudhi/Simplex_tree.h>
#include <gudhi/Rips_complex.h>
#include <gudhi/Sparse_rips_complex.h>
#include <gudhi/distance_functions.h>

#include <optional>

#include <iostream>
#include <vector>
#include <utility>  // std::pair
#include <string>

namespace Gudhi {

namespace rips_complex {

// /////////////////////////////////////////////////////////////////////////////
// Rips_complex_interface declaration
// /////////////////////////////////////////////////////////////////////////////

class Rips_complex_interface
{
    using Point_d = std::vector<double>;
    using Distance_matrix = std::vector<std::vector<Simplex_tree_interface::Filtration_value>>;

public:
     Rips_complex_interface() = default;
    ~Rips_complex_interface() = default;

    void init_points(const std::vector<std::vector<double>>& points, double threshold);
    void init_matrix(const std::vector<std::vector<double>>& matrix, double threshold);

    void init_points_sparse(const std::vector<std::vector<double>>& points, double threshold, double epsilon);
    void init_matrix_sparse(const std::vector<std::vector<double>>& matrix, double threshold, double epsilon);

    void create_simplex_tree(Simplex_tree_interface* simplex_tree, int dim_max);

private:
    // std::variant would work, but we don't require C++17 yet, and boost::variant is not super convenient.
    // Anyway, storing a graph would make more sense. Or changing the interface completely so there is no such storage.
    std::optional<Rips_complex<Simplex_tree_interface::Filtration_value>> rips_complex_;
    std::optional<Sparse_rips_complex<Simplex_tree_interface::Filtration_value>> sparse_rips_complex_;
};

// /////////////////////////////////////////////////////////////////////////////
// Rips_complex_interface definition
// /////////////////////////////////////////////////////////////////////////////

void Rips_complex_interface::init_points(const std::vector<std::vector<double>>& points, double threshold)
{
    rips_complex_.emplace(points, threshold, Gudhi::Euclidean_distance());
}

void Rips_complex_interface::init_matrix(const std::vector<std::vector<double>>& matrix, double threshold)
{
    rips_complex_.emplace(matrix, threshold);
}

void Rips_complex_interface::init_points_sparse(const std::vector<std::vector<double>>& points, double threshold, double epsilon)
{
    sparse_rips_complex_.emplace(points, Gudhi::Euclidean_distance(), epsilon, -std::numeric_limits<double>::infinity(), threshold);
}

void Rips_complex_interface::init_matrix_sparse(const std::vector<std::vector<double>>& matrix, double threshold, double epsilon)
{
    sparse_rips_complex_.emplace(matrix, epsilon, -std::numeric_limits<double>::infinity(), threshold);
}

void Rips_complex_interface::create_simplex_tree(Simplex_tree_interface* simplex_tree, int dim_max)
{
    if (rips_complex_) {
        rips_complex_->create_complex(*simplex_tree, dim_max);
    } else {
        sparse_rips_complex_->create_complex(*simplex_tree, dim_max);
    }
}


}  // namespace rips_complex

}  // namespace Gudhi

// /////////////////////////////////////////////////////////////////////////////
// Rips_complex_interface wrapping
// /////////////////////////////////////////////////////////////////////////////

#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;
namespace grc = Gudhi::rips_complex;
using grci = grc::Rips_complex_interface;

NB_MODULE(_rips_complex_ext, m) {
    m.attr("__license__") = "GPL v3";

    nb::class_<grci>(m, "Rips_complex_interface")
            .def(nb::init<>(), "Constructor")
            .def("init_points", &grci::init_points, "")
            .def("init_matrix", &grci::init_matrix, "")
            .def("init_points_sparse", &grci::init_points_sparse, "")
            .def("init_matrix_sparse", &grci::init_matrix_sparse, "")
            .def("create_simplex_tree", &grci::create_simplex_tree, "");
}


//
// _rips_complex.cc ends here

/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - 2025/03 Thibaud Kloczko: Use nanobind instead of Cython for python bindings.
 *      - 2025/04 Hannah Schreiber: Re-add possibility of tensors (numpy, torch etc.) as input.
 *      - YYYY/MM Author: Description of the modification
 */

#include <optional>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/ndarray.h>

#include <gudhi/Simplex_tree.h>
#include <gudhi/Rips_complex.h>
#include <gudhi/Sparse_rips_complex.h>
#include <gudhi/distance_functions.h>
#include <python_interfaces/Simplex_tree_interface.h>
#include <python_interfaces/points_utils.h>

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

  void init_points(const Sequence& points, double threshold)
  {
    rips_complex_.emplace(points, threshold, Gudhi::Euclidean_distance());
  }

  void init_points(const Tensor& points, double threshold)
  {
    init_points(_get_sequence_from_tensor(points), threshold);
  }

  void init_matrix(const Sequence& matrix, double threshold) { rips_complex_.emplace(matrix, threshold); }

  void init_matrix(const Tensor& matrix, double threshold)
  {
    init_matrix(_get_sequence_from_tensor(matrix), threshold);
  }

  void init_points_sparse(const Sequence& points, double threshold, double epsilon)
  {
    sparse_rips_complex_.emplace(
        points, Gudhi::Euclidean_distance(), epsilon, -std::numeric_limits<double>::infinity(), threshold);
  }

  void init_points_sparse(const Tensor& points, double threshold, double epsilon)
  {
    init_points_sparse(_get_sequence_from_tensor(points), threshold, epsilon);
  }

  void init_matrix_sparse(const Sequence& matrix, double threshold, double epsilon)
  {
    sparse_rips_complex_.emplace(matrix, epsilon, -std::numeric_limits<double>::infinity(), threshold);
  }

  void init_matrix_sparse(const Tensor& matrix, double threshold, double epsilon)
  {
    init_matrix_sparse(_get_sequence_from_tensor(matrix), threshold, epsilon);
  }

  void create_simplex_tree(Simplex_tree_interface* simplex_tree, int dim_max)
  {
    if (rips_complex_) {
      rips_complex_->create_complex(*simplex_tree, dim_max);
    } else {
      sparse_rips_complex_->create_complex(*simplex_tree, dim_max);
    }
  }

 private:
  // std::variant would work, but we don't require C++17 yet, and boost::variant is not super convenient.
  // Anyway, storing a graph would make more sense. Or changing the interface completely so there is no such storage.
  std::optional<Rips_complex<Simplex_tree_interface::Filtration_value>> rips_complex_;
  std::optional<Sparse_rips_complex<Simplex_tree_interface::Filtration_value>> sparse_rips_complex_;
};

}  // namespace rips_complex
}  // namespace Gudhi

// /////////////////////////////////////////////////////////////////////////////
// Rips_complex_interface wrapping
// /////////////////////////////////////////////////////////////////////////////

namespace nb = nanobind;
namespace grc = Gudhi::rips_complex;
using grci = grc::Rips_complex_interface;

NB_MODULE(_rips_complex_ext, m)
{
  m.attr("__license__") = "MIT";

  nb::class_<grci>(m, "Rips_complex_interface")
      .def(nb::init<>(), "Constructor")
      .def("init_points",
           nb::overload_cast<const Sequence&, double>(&grci::init_points),
           nb::arg("points"),
           nb::arg("threshold"))
      .def("init_matrix",
           nb::overload_cast<const Sequence&, double>(&grci::init_matrix),
           nb::arg("matrix"),
           nb::arg("threshold"))
      .def("init_points_sparse",
           nb::overload_cast<const Sequence&, double, double>(&grci::init_points_sparse),
           nb::arg("points"),
           nb::arg("threshold"),
           nb::arg("epsilon"))
      .def("init_matrix_sparse",
           nb::overload_cast<const Sequence&, double, double>(&grci::init_matrix_sparse),
           nb::arg("points"),
           nb::arg("threshold"),
           nb::arg("epsilon"))
      .def("init_points_with_tensor",
           nb::overload_cast<const Tensor&, double>(&grci::init_points),
           nb::arg("points"),
           nb::arg("threshold"))
      .def("init_matrix_with_tensor",
           nb::overload_cast<const Tensor&, double>(&grci::init_matrix),
           nb::arg("matrix"),
           nb::arg("threshold"))
      .def("init_points_sparse_with_tensor",
           nb::overload_cast<const Tensor&, double, double>(&grci::init_points_sparse),
           nb::arg("points"),
           nb::arg("threshold"),
           nb::arg("epsilon"))
      .def("init_matrix_sparse_with_tensor",
           nb::overload_cast<const Tensor&, double, double>(&grci::init_matrix_sparse),
           nb::arg("points"),
           nb::arg("threshold"),
           nb::arg("epsilon"))
      .def("create_simplex_tree", &grci::create_simplex_tree, "");
}

//
// _rips_complex.cc ends here

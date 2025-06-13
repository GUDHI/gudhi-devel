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

#include <boost/range/iterator_range_core.hpp>
#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>

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
  using Tensor_point_iterator = Numpy_array_element_iterator<double>;
  using Tensor_point = boost::iterator_range<Tensor_point_iterator>;

 public:
  // All init methods (see old cythonization strategy) were merged into the constructors, such that nanobind
  // is the one choosing what is a Sequence and what is a Tensor. Testing this in Python is uselessly tricky.

  Rips_complex_interface(const Sequence2D& array, double threshold, bool isPoints)
  {
    if (isPoints)
      rips_complex_.emplace(array, threshold, Gudhi::Euclidean_distance());
    else
      rips_complex_.emplace(array, threshold);
  }

  Rips_complex_interface(const Tensor2D& array, double threshold, bool isPoints)
  {
    if (isPoints) {
      auto view = array.view();
      auto distance = [&view](const double& p1, const double& p2) -> Simplex_tree_interface::Filtration_value {
        auto end1 = &p1 + view.shape(1) * view.stride(1);
        auto r1 = Tensor_point(Tensor_point_iterator(&p1, end1, view.stride(1)), Tensor_point_iterator(end1));
        auto end2 = &p2 + view.shape(1) * view.stride(1);
        auto r2 = Tensor_point(Tensor_point_iterator(&p2, end2, view.stride(1)), Tensor_point_iterator(end2));
        return Gudhi::Euclidean_distance()(r1, r2);
      };
      rips_complex_.emplace(Numpy_2d_span<double>(array), threshold, distance);
    } else {
      rips_complex_.emplace(_get_sequence_from_tensor(array), threshold);
    }
  }

  Rips_complex_interface(const Sequence2D& array, double threshold, double epsilon, bool isPoints)
  {
    if (isPoints)
      sparse_rips_complex_.emplace(
          array, Gudhi::Euclidean_distance(), epsilon, -std::numeric_limits<double>::infinity(), threshold);
    else
      sparse_rips_complex_.emplace(array, epsilon, -std::numeric_limits<double>::infinity(), threshold);
  }

  Rips_complex_interface(const Tensor2D& array, double threshold, double epsilon, bool isPoints)
  {
    if (isPoints) {
      auto view = array.view();
      auto distance = [&view](const double* p1, const double* p2) -> Simplex_tree_interface::Filtration_value {
        auto end1 = p1 + view.shape(1) * view.stride(1);
        auto r1 = Tensor_point(Tensor_point_iterator(p1, end1, view.stride(1)), Tensor_point_iterator(end1));
        auto end2 = p2 + view.shape(1) * view.stride(1);
        auto r2 = Tensor_point(Tensor_point_iterator(p2, end2, view.stride(1)), Tensor_point_iterator(end2));
        return Gudhi::Euclidean_distance()(r1, r2);
      };
      sparse_rips_complex_.emplace(
          Numpy_2d_span<double>(array), distance, epsilon, -std::numeric_limits<double>::infinity(), threshold);
    } else {
      sparse_rips_complex_.emplace(
          _get_sequence_from_tensor(array), epsilon, -std::numeric_limits<double>::infinity(), threshold);
    }
  }

  ~Rips_complex_interface() = default;

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
      .def(nb::init<const Sequence2D&, double, bool>(), nb::call_guard<nb::gil_scoped_release>())
      .def(nb::init<const Tensor2D&, double, bool>(), nb::call_guard<nb::gil_scoped_release>())
      .def(nb::init<const Sequence2D&, double, double, bool>(), nb::call_guard<nb::gil_scoped_release>())
      .def(nb::init<const Tensor2D&, double, double, bool>(), nb::call_guard<nb::gil_scoped_release>())
      .def("create_simplex_tree", &grci::create_simplex_tree, nb::call_guard<nb::gil_scoped_release>());
}

//
// _rips_complex.cc ends here

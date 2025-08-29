/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - 2024/03 Vincent Rouvreau: Renamed Alpha_complex_factory as Delaunay_complex_factory for DelaunayCechComplex.
 *                                  Factorize create_complex
 *      - 2024/10 Vincent Rouvreau: Add square root filtration values interface
 *      - 2025/03 Thibaud Kloczko: Use nanobind instead of Cython for python bindings
 *      - 2025/04 Hannah Schreiber: Re-add possibility of tensors (numpy, torch etc.) as input
 *      - YYYY/MM Author: Description of the modification
 */

#include <vector>
#include <memory>   // for std::unique_ptr
#include <cstddef>  // for std::size_t

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <python_interfaces/Simplex_tree_interface.h>
#include <python_interfaces/points_utils.h>
#include <python_interfaces/delaunay_complex_interface.h>

namespace Gudhi {
namespace delaunay_complex {

std::unique_ptr<Abstract_delaunay_complex> make_2d_epick_complex_ptr(const Sequence2D&, const Sequence1D&, bool);
std::unique_ptr<Abstract_delaunay_complex> make_2d_epeck_complex_ptr(const Sequence2D&, const Sequence1D&, bool);
std::unique_ptr<Abstract_delaunay_complex> make_2d_epick_complex_ptr(const Sequence2D&, bool);
std::unique_ptr<Abstract_delaunay_complex> make_2d_epeck_complex_ptr(const Sequence2D&, bool);
std::unique_ptr<Abstract_delaunay_complex> make_3d_epick_complex_ptr(const Sequence2D&, const Sequence1D&, bool);
std::unique_ptr<Abstract_delaunay_complex> make_3d_epeck_complex_ptr(const Sequence2D&, const Sequence1D&, bool);
std::unique_ptr<Abstract_delaunay_complex> make_3d_epick_complex_ptr(const Sequence2D&, bool);
std::unique_ptr<Abstract_delaunay_complex> make_3d_epeck_complex_ptr(const Sequence2D&, bool);
std::unique_ptr<Abstract_delaunay_complex> make_dynamic_d_epick_complex_ptr(const Sequence2D&, const Sequence1D&, bool);
std::unique_ptr<Abstract_delaunay_complex> make_dynamic_d_epeck_complex_ptr(const Sequence2D&, const Sequence1D&, bool);
std::unique_ptr<Abstract_delaunay_complex> make_dynamic_d_epick_complex_ptr(const Sequence2D&, bool);
std::unique_ptr<Abstract_delaunay_complex> make_dynamic_d_epeck_complex_ptr(const Sequence2D&, bool);

class Delaunay_complex_interface
{
 public:
  Delaunay_complex_interface(const Sequence2D& points, const Sequence1D& weights, bool fast_version, bool exact_version)
  {
    // Specific cases for dimensions 2 and 3
    const std::size_t dimension = ((points.size() > 0) ? points[0].size() : 0);
    if (fast_version) {
      if (dimension == 2) {
        delaunay_ptr_ = make_2d_epick_complex_ptr(points, weights, exact_version);
      } else if (dimension == 3) {
        delaunay_ptr_ = make_3d_epick_complex_ptr(points, weights, exact_version);
      } else {
        delaunay_ptr_ = make_dynamic_d_epick_complex_ptr(points, weights, exact_version);
      }
    } else {
      if (dimension == 2) {
        delaunay_ptr_ = make_2d_epeck_complex_ptr(points, weights, exact_version);
      } else if (dimension == 3) {
        delaunay_ptr_ = make_3d_epeck_complex_ptr(points, weights, exact_version);
      } else {
        delaunay_ptr_ = make_dynamic_d_epeck_complex_ptr(points, weights, exact_version);
      }
    }
  }

  Delaunay_complex_interface(const Tensor2D& points, const Tensor1D& weights, bool fast_version, bool exact_version)
      : Delaunay_complex_interface(_get_sequence_from_tensor(points),
                                   _get_sequence_from_tensor(weights),
                                   fast_version,
                                   exact_version)
  {}

  Delaunay_complex_interface(const Sequence2D& points, bool fast_version, bool exact_version)
  {
    // Specific cases for dimensions 2 and 3
    const std::size_t dimension = ((points.size() > 0) ? points[0].size() : 0);
    if (fast_version) {
      if (dimension == 2) {
        delaunay_ptr_ = make_2d_epick_complex_ptr(points, exact_version);
      } else if (dimension == 3) {
        delaunay_ptr_ = make_3d_epick_complex_ptr(points, exact_version);
      } else {
        delaunay_ptr_ = make_dynamic_d_epick_complex_ptr(points, exact_version);
      }
    } else {
      if (dimension == 2) {
        delaunay_ptr_ = make_2d_epeck_complex_ptr(points, exact_version);
      } else if (dimension == 3) {
        delaunay_ptr_ = make_3d_epeck_complex_ptr(points, exact_version);
      } else {
        delaunay_ptr_ = make_dynamic_d_epeck_complex_ptr(points, exact_version);
      }
    }
  }

  Delaunay_complex_interface(const Tensor2D& points, bool fast_version, bool exact_version)
      : Delaunay_complex_interface(_get_sequence_from_tensor(points), fast_version, exact_version)
  {}

  std::vector<double> get_point(int vh) { return delaunay_ptr_->get_point(vh); }

  void create_simplex_tree(Simplex_tree_interface* simplex_tree,
                           double max_alpha_square,
                           Delaunay_filtration filtration,
                           bool output_squared_values)
  {
    // Nothing to be done in case of an empty point set
    if (delaunay_ptr_->num_vertices() > 0) {
      delaunay_ptr_->create_simplex_tree(simplex_tree, max_alpha_square, filtration, output_squared_values);
    }
  }

  static void set_float_relative_precision(double precision)
  {
    // cf. CGAL::Epeck_d kernel type in Delaunay_complex_interface
    if (precision <= 0 || precision >= 1) {
      throw std::invalid_argument("Precision must be strictly greater than 0 and lower than 0");
    }
    CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>::FT::set_relative_precision_of_to_double(precision);
  }

  static double get_float_relative_precision()
  {
    // cf. CGAL::Epeck_d kernel type in Delaunay_complex_interface
    return CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>::FT::get_relative_precision_of_to_double();
  }

 private:
  std::unique_ptr<Abstract_delaunay_complex> delaunay_ptr_;
};

}  // namespace delaunay_complex
}  // namespace Gudhi

// /////////////////////////////////////////////////////////////////////////////
// Delaunay_complex_interface wrapping
// /////////////////////////////////////////////////////////////////////////////

namespace nb = nanobind;
namespace gdc = Gudhi::delaunay_complex;
using gdci = gdc::Delaunay_complex_interface;

NB_MODULE(_delaunay_complex_ext, m)
{
  m.attr("__license__") = "GPL v3";

  nb::enum_<gdc::Delaunay_filtration>(m, "Filtration")
      .value("NONE", gdc::Delaunay_filtration::NONE, "Default Delaunay Complex")
      .value("CECH", gdc::Delaunay_filtration::CECH, "Delaunay Cech Complex")
      .value("ALPHA", gdc::Delaunay_filtration::ALPHA, "Alpha Complex");

  nb::class_<gdci>(m, "Delaunay_complex_interface")
      .def(nb::init<const Sequence2D&, const Sequence1D&, bool, bool>(), nb::call_guard<nb::gil_scoped_release>())
      .def(nb::init<const Tensor2D&, const Tensor1D&, bool, bool>(), nb::call_guard<nb::gil_scoped_release>())
      .def(nb::init<const Sequence2D&, bool, bool>(), nb::call_guard<nb::gil_scoped_release>())
      .def(nb::init<const Tensor2D&, bool, bool>(), nb::call_guard<nb::gil_scoped_release>())
      .def("create_simplex_tree", &gdci::create_simplex_tree, nb::call_guard<nb::gil_scoped_release>())
      .def("get_point", &gdci::get_point, R"doc(
This function returns the point corresponding to a given vertex from the :class:`~gudhi.SimplexTree` (the
same as the k-th input point, where `k=vertex`)

Args:
    vertex: The vertex.
Returns:
    the point.

:raises IndexError: In case the point has no associated vertex in the diagram (because of weights or because it
    is a duplicate).
        )doc")
      .def_static("get_float_relative_precision", &gdci::get_float_relative_precision, R"doc(
Get the float relative precision of filtration values computation when constructing with :code:`precision = 'safe'`
(the default).

Returns:
    The float relative precision.
        )doc")
      .def_static("set_float_relative_precision", &gdci::set_float_relative_precision, R"doc(
Set the float relative precision of filtration values computation when constructing with :code:`precision = 'safe'`
(the default).

Args:
    precision: When constructing :func:`~gudhi.delaunay_cech_complex`, :func:`~gudhi.alpha_complex`, or
        :func:`~gudhi.weighted_alpha_complex` with :code:`precision = 'safe'` (the default), one can
        set the float relative precision of filtration values computed. Default is :code:`1e-5` (cf.
        :func:`~gudhi.DelaunayComplex.get_float_relative_precision`). For more details, please refer to
        `CGAL::Lazy_exact_nt<NT>::set_relative_precision_of_to_double <https://doc.cgal.org/latest/Number_types/classCGAL_1_1Lazy__exact__nt.html>`_

:raises ValueError: If precision is not in (0, 1).
        )doc");
}

//
// _delaunay_complex.cc ends here

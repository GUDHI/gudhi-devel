/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - 2025/03 Hannah Schreiber: Use nanobind instead of Cython for python bindings.
 *      - YYYY/MM Author: Description of the modification
 */

#include <string>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/pair.h>
#include <nanobind/ndarray.h>

#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Bitmap_cubical_complex_base.h>
#include <gudhi/Bitmap_cubical_complex_periodic_boundary_conditions_base.h>
#include <python_interfaces/Persistent_cohomology_interface.h>

namespace Gudhi {
namespace cubical_complex {

class Cubical_complex_interface : public Bitmap_cubical_complex<Bitmap_cubical_complex_base<double>>
{
  using Base = Bitmap_cubical_complex<Bitmap_cubical_complex_base<double>>;

 public:
  using Base::Base;  // inheriting constructors

  explicit Cubical_complex_interface(const std::string& perseus_style_file) : Base(perseus_style_file.c_str()) {}

  // TODO: nanobind is probably making a copy here (to verify), as it is only used privately we could think
  // at another strategy?
  // But as the vector is probably very small (number of dimensions), it is perhaps not worth it.
  const std::vector<unsigned>& shape() { return this->sizes; };

  nanobind::ndarray<double, nanobind::numpy> get_numpy_array()
  {
    return nanobind::ndarray<double, nanobind::numpy>(Base::data.data(), {Base::data.size()});
  }
};

class Periodic_cubical_complex_interface
    : public Bitmap_cubical_complex<Bitmap_cubical_complex_periodic_boundary_conditions_base<double>>
{
  using Base = Bitmap_cubical_complex<Bitmap_cubical_complex_periodic_boundary_conditions_base<double>>;

 public:
  using Base::Base;  // inheriting constructors

  explicit Periodic_cubical_complex_interface(const std::string& perseus_style_file) : Base(perseus_style_file.c_str())
  {}

  // TODO: nanobind is probably making a copy here (to verify), as it is only used privately we could think
  // of another strategy?
  // But as the vector is probably very small (number of dimensions), it is perhaps not worth it.
  const std::vector<unsigned>& shape() { return this->sizes; };

  // TODO: nanobind is probably making a copy here (to verify), as it is only used privately we could think
  // of another strategy?
  // But as the vector is probably very small (number of dimensions), it is perhaps not worth it.
  const std::vector<bool>& periodicities() { return this->directions_in_which_periodic_b_cond_are_to_be_imposed; }

  nanobind::ndarray<double, nanobind::numpy> get_numpy_array()
  {
    return nanobind::ndarray<double, nanobind::numpy>(Base::data.data(), {Base::data.size()});
  }
};

}  // namespace cubical_complex
}  // namespace Gudhi

namespace nb = nanobind;

using CC = Gudhi::cubical_complex::Cubical_complex_interface;
using CPers = Gudhi::Persistent_cohomology_interface<CC>;

using PCC = Gudhi::cubical_complex::Periodic_cubical_complex_interface;
using PCPers = Gudhi::Persistent_cohomology_interface<PCC>;

NB_MODULE(_cubical_complex_ext, m)
{
  m.attr("__license__") = "MIT";

  nb::class_<CC>(m, "_Bitmap_cubical_complex_interface")
      .def(nb::init<const std::vector<unsigned int>&, const std::vector<double>&, bool>(),
           nb::call_guard<nb::gil_scoped_release>())
      .def(nb::init<const std::string&>(), nb::call_guard<nb::gil_scoped_release>())
      .def("num_simplices", &CC::num_simplices, nb::call_guard<nb::gil_scoped_release>(), R"doc(
This function returns the number of all cubes in the complex.

:returns:  int -- the number of all cubes in the complex.
           )doc")
      .def("dimension",
           nb::overload_cast<>(&CC::dimension, nb::const_),
           nb::call_guard<nb::gil_scoped_release>(),
           R"doc(
This function returns the dimension of the complex.

:returns:  int -- the complex dimension.
           )doc")
      .def("shape", &CC::shape)
      .def("_get_numpy_array", &CC::get_numpy_array, nb::rv_policy::reference_internal);

  nb::class_<CPers>(m, "_Cubical_complex_persistence_interface")
      .def(nb::init<CC&, bool>(), nb::call_guard<nb::gil_scoped_release>())
      .def("_compute_persistence", &CPers::compute_persistence, nb::call_guard<nb::gil_scoped_release>())
      .def("_get_persistence", &CPers::get_persistence)
      .def("_cofaces_of_cubical_persistence_pairs", &CPers::cofaces_of_cubical_persistence_pairs)
      .def("_vertices_of_cubical_persistence_pairs", &CPers::vertices_of_cubical_persistence_pairs)
      .def("_betti_numbers", &CPers::betti_numbers)
      .def("_persistent_betti_numbers", &CPers::persistent_betti_numbers)
      .def("_intervals_in_dimension", &CPers::intervals_in_dimension);

  nb::class_<PCC>(m, "_Periodic_cubical_complex_interface")
      .def(nb::init<const std::vector<unsigned int>&, const std::vector<double>&, const std::vector<bool>&, bool>(),
           nb::call_guard<nb::gil_scoped_release>())
      .def(nb::init<const std::string&>(), nb::call_guard<nb::gil_scoped_release>())
      .def("num_simplices", &PCC::num_simplices, nb::call_guard<nb::gil_scoped_release>(), R"doc(
This function returns the number of all cubes in the complex.

:returns:  int -- the number of all cubes in the complex.
           )doc")
      .def("dimension",
           nb::overload_cast<>(&PCC::dimension, nb::const_),
           nb::call_guard<nb::gil_scoped_release>(),
           R"doc(
This function returns the dimension of the complex.

:returns:  int -- the complex dimension.
           )doc")
      .def("shape", &PCC::shape)
      .def("periodicities", &PCC::periodicities)
      .def("_get_numpy_array", &PCC::get_numpy_array, nb::rv_policy::reference_internal);

  nb::class_<PCPers>(m, "_Periodic_cubical_complex_persistence_interface")
      .def(nb::init<PCC&, bool>(), nb::call_guard<nb::gil_scoped_release>())
      .def("_compute_persistence", &PCPers::compute_persistence, nb::call_guard<nb::gil_scoped_release>())
      .def("_get_persistence", &PCPers::get_persistence)
      .def("_cofaces_of_cubical_persistence_pairs", &PCPers::cofaces_of_cubical_persistence_pairs)
      .def("_vertices_of_cubical_persistence_pairs", &PCPers::vertices_of_cubical_persistence_pairs)
      .def("_betti_numbers", &PCPers::betti_numbers)
      .def("_persistent_betti_numbers", &PCPers::persistent_betti_numbers)
      .def("_intervals_in_dimension", &PCPers::intervals_in_dimension);
}

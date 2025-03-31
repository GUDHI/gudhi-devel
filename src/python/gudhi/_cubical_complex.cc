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
#include <Persistent_cohomology_interface.h>

namespace Gudhi {
namespace cubical_complex {

class Cubical_complex_interface : public Bitmap_cubical_complex<Bitmap_cubical_complex_base<double>>
{
  using Base = Bitmap_cubical_complex<Bitmap_cubical_complex_base<double>>;

 public:
  using Base::Base;  // inheriting constructors
  using Base::data;

  explicit Cubical_complex_interface(const std::string& perseus_style_file) : Base(perseus_style_file.c_str()) {}

  // not const because cython does not handle const very well
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
  using Base::data;

  explicit Periodic_cubical_complex_interface(const std::string& perseus_style_file) : Base(perseus_style_file.c_str())
  {}

  // not const because cython does not handle const very well
  const std::vector<unsigned>& shape() { return this->sizes; };

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
      .def(nb::init<const std::vector<unsigned int>&, const std::vector<double>&, bool>())
      .def(nb::init<const std::string&>())
      .def("num_simplices", &CC::num_simplices, R"pbdoc(
        """This function returns the number of all cubes in the complex.

        :returns:  int -- the number of all cubes in the complex.
        """
        )pbdoc")
      .def("dimension", nb::overload_cast<>(&CC::dimension, nb::const_), R"pbdoc(
        """This function returns the dimension of the complex.

        :returns:  int -- the complex dimension.
        """
        )pbdoc")
      .def("shape", &CC::shape)
      .def("get_numpy_array", &CC::get_numpy_array, nb::rv_policy::reference_internal)
      .def_rw("data", &CC::data);

  nb::class_<CPers>(m, "_Cubical_complex_persistence_interface")
      .def(nb::init<CC&, bool>())
      .def("compute_persistence",
           &CPers::compute_persistence,
           nb::arg("homology_coeff_field"),
           nb::arg("double min_persistence"))
      .def("get_persistence", &CPers::get_persistence)
      .def("cofaces_of_cubical_persistence_pairs", &CPers::cofaces_of_cubical_persistence_pairs)
      .def("vertices_of_cubical_persistence_pairs", &CPers::vertices_of_cubical_persistence_pairs)
      .def("betti_numbers", &CPers::betti_numbers)
      .def("persistent_betti_numbers", &CPers::persistent_betti_numbers, nb::arg("from_value"), nb::arg("to_value"))
      .def("intervals_in_dimension", &CPers::intervals_in_dimension, nb::arg("dimension"));

  nb::class_<PCC>(m, "_Periodic_cubical_complex_interface")
      .def(nb::init<const std::vector<unsigned int>&, const std::vector<double>&, const std::vector<bool>&, bool>())
      .def(nb::init<const std::string&>())
      .def("num_simplices", &PCC::num_simplices, R"pbdoc(
        """This function returns the number of all cubes in the complex.

        :returns:  int -- the number of all cubes in the complex.
        """
        )pbdoc")
      .def("dimension", nb::overload_cast<>(&PCC::dimension, nb::const_), R"pbdoc(
        """This function returns the dimension of the complex.

        :returns:  int -- the complex dimension.
        """
        )pbdoc")
      .def("shape", &PCC::shape)
      .def("periodicities", &PCC::periodicities)
      .def("get_numpy_array", &PCC::get_numpy_array, nb::rv_policy::reference_internal)
      .def_rw("data", &PCC::data);

  nb::class_<PCPers>(m, "_Periodic_cubical_complex_persistence_interface")
      .def(nb::init<PCC&, bool>())
      .def("compute_persistence",
           &PCPers::compute_persistence,
           nb::arg("homology_coeff_field"),
           nb::arg("double min_persistence"))
      .def("get_persistence", &PCPers::get_persistence)
      .def("cofaces_of_cubical_persistence_pairs", &PCPers::cofaces_of_cubical_persistence_pairs)
      .def("vertices_of_cubical_persistence_pairs", &PCPers::vertices_of_cubical_persistence_pairs)
      .def("betti_numbers", &PCPers::betti_numbers)
      .def("persistent_betti_numbers", &PCPers::persistent_betti_numbers, nb::arg("from_value"), nb::arg("to_value"))
      .def("intervals_in_dimension", &PCPers::intervals_in_dimension, nb::arg("dimension"));
}

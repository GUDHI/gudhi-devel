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

#include <boost/mpl/for_each.hpp>

#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Bitmap_cubical_complex_base.h>
#include <gudhi/Bitmap_cubical_complex_periodic_boundary_conditions_base.h>
#include <python_interfaces/Persistent_cohomology_interface.h>
#include <python_interfaces/filtration_value_types_support.h>

namespace nb = nanobind;

namespace Gudhi {
namespace cubical_complex {

template <class Filtration_value>
class Cubical_complex_interface : public Bitmap_cubical_complex<Bitmap_cubical_complex_base<Filtration_value>>
{
  using Base = Bitmap_cubical_complex<Bitmap_cubical_complex_base<Filtration_value>>;

 public:
  using Base::Base;  // inheriting constructors

  explicit Cubical_complex_interface(const std::string& perseus_style_file) : Base(perseus_style_file.c_str()) {}

  // TODO: nanobind is probably making a copy here (to verify), as it is only used privately we could think
  // at another strategy?
  // But as the vector is probably very small (number of dimensions), it is perhaps not worth it.
  const std::vector<unsigned>& shape() { return this->sizes; };

  nb::ndarray<Filtration_value, nb::numpy> get_numpy_array()
  {
    return nb::ndarray<Filtration_value, nb::numpy>(Base::data.data(), {Base::data.size()});
  }
};

template <class Filtration_value>
class Periodic_cubical_complex_interface
    : public Bitmap_cubical_complex<Bitmap_cubical_complex_periodic_boundary_conditions_base<Filtration_value>>
{
  using Base = Bitmap_cubical_complex<Bitmap_cubical_complex_periodic_boundary_conditions_base<Filtration_value>>;

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

  nb::ndarray<Filtration_value, nb::numpy> get_numpy_array()
  {
    return nb::ndarray<Filtration_value, nb::numpy>(Base::data.data(), {Base::data.size()});
  }
};

}  // namespace cubical_complex
}  // namespace Gudhi

NB_MODULE(_cubical_complex_ext, m)
{
  m.attr("__license__") = "MIT";
  
  // Cubical_complex
  auto add_class_for_cubical_complex_interface = [&m](auto identity_tag) {
      using CC = typename decltype(identity_tag)::type;
      using Filtration_value = typename CC::Filtration_value;

      std::string class_name =
        Gudhi::get_class_name_for_filtration_type<Filtration_value>("_Bitmap_cubical_complex_interface");
      std::string persistence_class_name =
        Gudhi::get_class_name_for_filtration_type<Filtration_value>("_Cubical_complex_persistence_interface");

      nb::class_<CC>(m, class_name.c_str())
          .def(nb::init<const std::vector<unsigned int>&, const std::vector<Filtration_value>&, bool>(),
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
          .def("_shape", &CC::shape)
          .def("_get_numpy_array", &CC::get_numpy_array, nb::rv_policy::reference_internal);

      using CPers = Gudhi::Persistent_cohomology_interface<CC>;
      nb::class_<CPers>(m, persistence_class_name.c_str())
          .def(nb::init<CC&, bool>(), nb::call_guard<nb::gil_scoped_release>())
          .def("_compute_persistence", &CPers::compute_persistence, nb::call_guard<nb::gil_scoped_release>())
          .def("_get_persistence", &CPers::get_persistence)
          .def("_cofaces_of_cubical_persistence_pairs",
               &CPers::cofaces_of_cubical_persistence_pairs,
               nb::call_guard<nb::gil_scoped_release>())
          .def("_vertices_of_cubical_persistence_pairs",
               &CPers::vertices_of_cubical_persistence_pairs,
               nb::call_guard<nb::gil_scoped_release>())
          .def("_betti_numbers", &CPers::betti_numbers)
          .def("_persistent_betti_numbers", &CPers::persistent_betti_numbers)
          .def("_intervals_in_dimension", &CPers::intervals_in_dimension);

  };

  Gudhi::for_each_filtration_value_type(
    Gudhi::Supported_for_data_structure<Gudhi::cubical_complex::Cubical_complex_interface>::List{},
    add_class_for_cubical_complex_interface
  );
  
  // Periodic_cubical_complex
  auto add_class_for_periodic_cubical_complex_interface = [&m](auto identity_tag) {
      using PCC = typename decltype(identity_tag)::type;
      using Filtration_value = typename PCC::Filtration_value;

      std::string class_name =
        Gudhi::get_class_name_for_filtration_type<Filtration_value>("_Periodic_cubical_complex_interface");
      std::string persistence_class_name = 
        Gudhi::get_class_name_for_filtration_type<Filtration_value>("_Periodic_cubical_complex_persistence_interface");
      
      nb::class_<PCC>(m, class_name.c_str())
          .def(nb::init<const std::vector<unsigned int>&, const std::vector<Filtration_value>&, const std::vector<bool>&, bool>(),
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
          .def("_shape", &PCC::shape)
          .def("_periodicities", &PCC::periodicities)
          .def("_get_numpy_array", &PCC::get_numpy_array, nb::rv_policy::reference_internal);
    
      using PCPers = Gudhi::Persistent_cohomology_interface<PCC>;
      nb::class_<PCPers>(m, persistence_class_name.c_str())
          .def(nb::init<PCC&, bool>(), nb::call_guard<nb::gil_scoped_release>())
          .def("_compute_persistence", &PCPers::compute_persistence, nb::call_guard<nb::gil_scoped_release>())
          .def("_get_persistence", &PCPers::get_persistence)
          .def("_cofaces_of_cubical_persistence_pairs",
               &PCPers::cofaces_of_cubical_persistence_pairs,
               nb::call_guard<nb::gil_scoped_release>())
          .def("_vertices_of_cubical_persistence_pairs",
               &PCPers::vertices_of_cubical_persistence_pairs,
               nb::call_guard<nb::gil_scoped_release>())
          .def("_betti_numbers", &PCPers::betti_numbers)
          .def("_persistent_betti_numbers", &PCPers::persistent_betti_numbers)
          .def("_intervals_in_dimension", &PCPers::intervals_in_dimension);


  };
  
  Gudhi::for_each_filtration_value_type(
    Gudhi::Supported_for_data_structure<Gudhi::cubical_complex::Periodic_cubical_complex_interface>::List{},
    add_class_for_periodic_cubical_complex_interface
  );

}

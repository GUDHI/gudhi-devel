/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Marc Glisse
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <vector>
#include <array>
#include <limits>
#include <stdexcept>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

#include <boost/range/counting_range.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include <gudhi/Persistence_on_a_line.h>
#include <gudhi/Persistence_on_rectangle.h>
#include <gudhi/Debug_utils.h>

namespace py = pybind11;
typedef std::vector<std::array< float, 2>> Vf;
typedef std::vector<std::array<double, 2>> Vd;
PYBIND11_MAKE_OPAQUE(Vf);
PYBIND11_MAKE_OPAQUE(Vd);

template<class T>
py::array wrap_persistence_1d(py::array_t<T> data) {
  py::buffer_info buf = data.request();
  if(buf.ndim!=1)
    throw std::runtime_error("Data must be a 1-dimensional array");
  auto cnt = boost::counting_range<py::ssize_t>(0, buf.shape[0]);
  char* p = static_cast<char*>(buf.ptr);
  auto h = buf.strides[0];
  auto proj = [=](py::ssize_t i){ return *reinterpret_cast<T*>(p + i * h); };
  auto r = boost::adaptors::transform(cnt, proj);
  std::vector<std::array<T, 2>> dgm;
  {
    py::gil_scoped_release release;
    Gudhi::persistent_cohomology::compute_persistence_of_function_on_line(r, [&](T b, T d){ dgm.push_back({b, d}); });
  }
  return py::array(py::cast(std::move(dgm)));
}

py::list wrap_persistence_2d(py::array_t<double, py::array::c_style | py::array::forcecast> data, double min_persistence) {
  py::buffer_info buf = data.request();
  if(buf.ndim!=2)
    throw std::runtime_error("Data must be a 2-dimensional array");
  // If we make this function public, it should probably be enhanced to handle these cases.
  if(buf.shape[0] < 2 || buf.shape[1] < 2)
    throw std::runtime_error("The Python caller is supposed to ensure that shape[i]>=2");
  std::vector<std::array<double, 2>> dgm0, dgm1;
  {
    py::gil_scoped_release release;
    double mini = Gudhi::cubical_complex::persistence_on_rectangle_from_top_cells(
        static_cast<double const*>(buf.ptr),
        static_cast<unsigned>(buf.shape[0]),
        static_cast<unsigned>(buf.shape[1]),
        [&](double b, double d){ if (d - b > min_persistence) dgm0.push_back({b, d}); },
        [&](double b, double d){ if (d - b > min_persistence) dgm1.push_back({b, d}); });
    dgm0.push_back({mini, std::numeric_limits<double>::infinity()});
  }
  py::list ret;
  ret.append(py::array(py::cast(std::move(dgm0))));
  ret.append(py::array(py::cast(std::move(dgm1))));
  return ret;
}

PYBIND11_MODULE(_pers_cub_low_dim, m) {
  py::bind_vector<Vf>(m, "VectorPairFloat" , py::buffer_protocol());
  py::bind_vector<Vd>(m, "VectorPairDouble", py::buffer_protocol());
  m.def("_persistence_on_a_line", wrap_persistence_1d<float>, py::arg().noconvert());
  m.def("_persistence_on_a_line", wrap_persistence_1d<double>);
  m.def("_persistence_on_rectangle_from_top_cells", wrap_persistence_2d);
}

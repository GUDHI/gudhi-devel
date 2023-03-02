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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

#include <boost/range/counting_range.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include <gudhi/Persistence_on_a_line.h>

namespace py = pybind11;
typedef std::vector<std::array< float, 2>> Vf;
typedef std::vector<std::array<double, 2>> Vd;
PYBIND11_MAKE_OPAQUE(Vf);
PYBIND11_MAKE_OPAQUE(Vd);

template<class T>
py::array fun(py::array_t<T> data) {
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

PYBIND11_MODULE(_persline, m) {
  py::bind_vector<Vf>(m, "VectorPairFloat" , py::buffer_protocol());
  py::bind_vector<Vd>(m, "VectorPairDouble", py::buffer_protocol());
  m.def("persistence_on_a_line", fun<float>, py::arg().noconvert());
  m.def("persistence_on_a_line", fun<double>);
}

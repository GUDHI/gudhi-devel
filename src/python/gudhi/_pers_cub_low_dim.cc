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

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/bind_vector.h>

#include <boost/range/counting_range.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include <gudhi/Persistence_on_a_line.h>
#include <gudhi/Persistence_on_rectangle.h>

namespace py = nanobind;
typedef std::vector<std::array< float, 2>> Vf;
typedef std::vector<std::array<double, 2>> Vd;
NB_MAKE_OPAQUE(Vf);
NB_MAKE_OPAQUE(Vd);

template<class T>
auto wrap_persistence_1d(py::ndarray<T, py::ndim<1>> data) {
  auto cnt = boost::counting_range<py::ssize_t>(0, data.shape(0));
  auto h = data.stride(0);
  auto proj = [=](py::ssize_t i){ return *reinterpret_cast<T*>(data.data() + i * h); };
  auto r = boost::adaptors::transform(cnt, proj);
  std::vector<std::array<T, 2>> dgm;
  {
    py::gil_scoped_release release;
    Gudhi::persistent_cohomology::compute_persistence_of_function_on_line(r, [&](T b, T d){ dgm.push_back({b, d}); });
  }
  return py::ndarray<py::numpy, py::ndim<2>, T>(dgm.data(), {dgm.size(), 2}).cast();
}

py::list wrap_persistence_2d(py::ndarray<double, py::ndim<2>, py::c_contig> data, double min_persistence) {
  std::vector<std::array<double, 2>> dgm0, dgm1;
  {
    py::gil_scoped_release release;
    double mini = Gudhi::cubical_complex::persistence_on_rectangle_from_top_cells(
        static_cast<double const*>(data.data()),
        static_cast<unsigned>(data.shape(0)),
        static_cast<unsigned>(data.shape(1)),
        [&](double b, double d){ if (d - b > min_persistence) dgm0.push_back({b, d}); },
        [&](double b, double d){ if (d - b > min_persistence) dgm1.push_back({b, d}); });
    dgm0.push_back({mini, std::numeric_limits<double>::infinity()});
  }
  py::list ret;
  //ret.append(py::cast(std::move(dgm0)));
  //ret.append(py::cast(std::move(dgm1)));
  ret.append(py::ndarray<py::numpy, py::ndim<2>, double>(dgm0.data(), {dgm0.size(), 2}).cast());
  ret.append(py::ndarray<py::numpy, py::ndim<2>, double>(dgm1.data(), {dgm1.size(), 2}).cast());
  return ret;
}

NB_MODULE(_pers_cub_low_dim, m) {
  py::bind_vector<Vf>(m, "VectorPairFloat");
  py::bind_vector<Vd>(m, "VectorPairDouble");
  m.def("_persistence_on_a_line", wrap_persistence_1d<float>, py::arg().noconvert());
  m.def("_persistence_on_a_line", wrap_persistence_1d<double>);
  m.def("_persistence_on_rectangle_from_top_cells", wrap_persistence_2d);
}

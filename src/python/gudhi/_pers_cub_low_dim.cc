/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Marc Glisse
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - 2025/01 Vincent Rouvreau: Use nanobind instead of PyBind11 for python bindings
 *      - YYYY/MM Author: Description of the modification
 */

#include <vector>
#include <array>
#include <limits>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <boost/range/counting_range.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include <gudhi/Persistence_on_a_line.h>
#include <gudhi/Persistence_on_rectangle.h>

namespace nb = nanobind;

template <class T>
auto wrap_persistence_1d(nb::ndarray<T, nb::ndim<1>> data)
{
  auto cnt = boost::counting_range<nb::ssize_t>(0, data.shape(0));
  auto h = data.stride(0);
  auto proj = [=](nb::ssize_t i) { return *reinterpret_cast<T*>(data.data() + i * h); };
  auto r = boost::adaptors::transform(cnt, proj);
  std::vector<std::array<T, 2>> dgm;
  {
    nb::gil_scoped_release release;
    Gudhi::persistent_cohomology::compute_persistence_of_function_on_line(r, [&](T b, T d) { dgm.push_back({b, d}); });
  }
  return nb::ndarray<nb::numpy, nb::ndim<2>, T>(dgm.data(), {dgm.size(), 2}).cast();
}

nb::list wrap_persistence_2d(nb::ndarray<double, nb::ndim<2>, nb::c_contig> data, double min_persistence)
{
  std::vector<std::array<double, 2>> dgm0, dgm1;
  {
    nb::gil_scoped_release release;
    double mini = Gudhi::cubical_complex::persistence_on_rectangle_from_top_cells(
        static_cast<double const*>(data.data()),
        static_cast<unsigned>(data.shape(0)),
        static_cast<unsigned>(data.shape(1)),
        [&](double b, double d) {
          if (d - b > min_persistence) dgm0.push_back({b, d});
        },
        [&](double b, double d) {
          if (d - b > min_persistence) dgm1.push_back({b, d});
        });
    dgm0.push_back({mini, std::numeric_limits<double>::infinity()});
  }
  nb::list ret;
  ret.append(nb::ndarray<nb::numpy, nb::ndim<2>, double>(dgm0.data(), {dgm0.size(), 2}).cast());
  ret.append(nb::ndarray<nb::numpy, nb::ndim<2>, double>(dgm1.data(), {dgm1.size(), 2}).cast());
  return ret;
}

NB_MODULE(_pers_cub_low_dim_ext, m)
{
  m.def("_persistence_on_a_line", &wrap_persistence_1d<float>, nb::arg().noconvert());
  m.def("_persistence_on_a_line", &wrap_persistence_1d<double>);
  m.def("_persistence_on_rectangle_from_top_cells", &wrap_persistence_2d);
}

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

#include <array>
#include <vector>
#include <limits>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <boost/range/counting_range.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include <gudhi/Persistence_on_a_line.h>
#include <gudhi/Persistence_on_rectangle.h>
#include <python_interfaces/numpy_utils.h>

namespace nb = nanobind;

template <class T>
auto wrap_persistence_1d(nb::ndarray<const T, nb::ndim<1> > data)
{
  auto data_view = data.view();
  auto cnt = boost::counting_range<nb::ssize_t>(0, data_view.shape(0));
  auto proj = [&data_view](nb::ssize_t i) { return data_view(i); };
  auto r = boost::adaptors::transform(cnt, proj);
  std::vector<std::array<T, 2> >* dgm = new std::vector<std::array<T, 2> >();
  dgm->reserve(data_view.shape(0) + 1);  // rough upper bound
  {
    nb::gil_scoped_release release;
    Gudhi::persistent_cohomology::compute_persistence_of_function_on_line(
        r, [&](T b, T d) { dgm->emplace_back(std::array<T, 2>{b, d}); });
  }
  if (dgm->size() < data_view.shape(0) / 2) dgm->shrink_to_fit();
  return _wrap_as_numpy_array(dgm);
}

nb::list wrap_persistence_2d(nb::ndarray<const double, nb::ndim<2>, nb::c_contig> data, double min_persistence)
{
  std::vector<std::array<double, 2> >* dgm0 = new std::vector<std::array<double, 2> >();
  std::vector<std::array<double, 2> >* dgm1 = new std::vector<std::array<double, 2> >();
  dgm0->reserve(data.shape(0) * data.shape(1) + data.shape(0) + data.shape(1) + 1);  // rough upper bound
  dgm1->reserve(data.shape(0) * data.shape(1) + 1);  // rough upper bound
  {
    nb::gil_scoped_release release;
    double mini = Gudhi::cubical_complex::persistence_on_rectangle_from_top_cells(
        static_cast<double const*>(data.data()),
        static_cast<unsigned int>(data.shape(0)),
        static_cast<unsigned int>(data.shape(1)),
        [&](double b, double d) {
          if (d - b > min_persistence) {
            dgm0->emplace_back(std::array<double, 2>{b, d});
          }
        },
        [&](double b, double d) {
          if (d - b > min_persistence) {
            dgm1->emplace_back(std::array<double, 2>{b, d});
          }
        });
    dgm0->emplace_back(std::array<double, 2>{mini, std::numeric_limits<double>::infinity()});
  }
  if (dgm0->size() < dgm0->capacity() / 2) dgm0->shrink_to_fit();
  if (dgm1->size() < dgm1->capacity() / 2) dgm1->shrink_to_fit();
  nb::list ret;
  ret.append(_wrap_as_numpy_array(dgm0));
  ret.append(_wrap_as_numpy_array(dgm1));
  return ret;
}

NB_MODULE(_pers_cub_low_dim_ext, m)
{
  m.attr("__license__") = "MIT";
  m.def("_persistence_on_a_line", &wrap_persistence_1d<float>, nb::arg().noconvert());
  m.def("_persistence_on_a_line", &wrap_persistence_1d<double>);
  m.def("_persistence_on_rectangle_from_top_cells", &wrap_persistence_2d);
}

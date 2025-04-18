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
  auto proj = [&](nb::ssize_t i) { return data(i); };
  auto r = boost::adaptors::transform(cnt, proj);
  std::vector<T>* dgm = new std::vector<T>();
  dgm->reserve(data.shape(0) * 2);  // rough upper bound, but makes it a difference in performance
  {
    nb::gil_scoped_release release;
    Gudhi::persistent_cohomology::compute_persistence_of_function_on_line(r, [&](T b, T d) {
      dgm->push_back(b);
      dgm->push_back(d);
    });
  }
  return nb::ndarray<T, nb::numpy>(dgm->data(), {dgm->size() / 2, 2}, nb::capsule(dgm, [](void* p) noexcept {
                                     delete reinterpret_cast<std::vector<T>*>(p);
                                   }));
}

nb::list wrap_persistence_2d(nb::ndarray<double, nb::ndim<2>, nb::c_contig> data, double min_persistence)
{
  std::vector<double>* dgm0 = new std::vector<double>();
  std::vector<double>* dgm1 = new std::vector<double>();
  dgm0->reserve(data.shape(0) * data.shape(1) * 2);  // rough upper bound
  dgm1->reserve(data.shape(0) * data.shape(1) * 2);  // rough upper bound
  {
    nb::gil_scoped_release release;
    double mini = Gudhi::cubical_complex::persistence_on_rectangle_from_top_cells(
        static_cast<double const*>(data.data()),
        static_cast<unsigned>(data.shape(0)),
        static_cast<unsigned>(data.shape(1)),
        [&](double b, double d) {
          if (d - b > min_persistence) {
            dgm0->push_back(b);
            dgm0->push_back(d);
          }
        },
        [&](double b, double d) {
          if (d - b > min_persistence) {
            dgm1->push_back(b);
            dgm1->push_back(d);
          }
        });
    dgm0->push_back(mini);
    dgm0->push_back(std::numeric_limits<double>::infinity());
  }
  nb::list ret;
  ret.append(
      nb::ndarray<nb::numpy, double>(dgm0->data(), {dgm0->size() / 2, 2}, nb::capsule(dgm0, [](void* p) noexcept {
                                       delete reinterpret_cast<std::vector<double>*>(p);
                                     })));
  ret.append(
      nb::ndarray<nb::numpy, double>(dgm1->data(), {dgm1->size() / 2, 2}, nb::capsule(dgm1, [](void* p) noexcept {
                                       delete reinterpret_cast<std::vector<double>*>(p);
                                     })));
  return ret;
}

NB_MODULE(_pers_cub_low_dim_ext, m)
{
  m.attr("__license__") = "MIT";
  m.def("_persistence_on_a_line", &wrap_persistence_1d<float>, nb::arg().noconvert());
  m.def("_persistence_on_a_line", &wrap_persistence_1d<double>);
  m.def("_persistence_on_rectangle_from_top_cells", &wrap_persistence_2d);
}

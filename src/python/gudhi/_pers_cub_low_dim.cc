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

// #include <array>
#include <cstddef>
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
auto wrap_persistence_1d(nb::ndarray<T, nb::ndim<1>> data)
{
  auto data_view = data.view();
  auto cnt = boost::counting_range<nb::ssize_t>(0, data_view.shape(0));
  auto proj = [&](nb::ssize_t i) { return data_view(i); };
  auto r = boost::adaptors::transform(cnt, proj);
  auto dgm = new T[data_view.shape(0) * 2];
  std::size_t count = 0;
  {
    nb::gil_scoped_release release;
    Gudhi::persistent_cohomology::compute_persistence_of_function_on_line(r, [&](T b, T d) {
      dgm[count] = b;
      dgm[count + 1] = d;
      count += 2;
    });
  }
  return _wrap_as_numpy_array(dgm, count / 2, 2);
  // std::vector<std::array<T, 2> >* dgm = new std::vector<std::array<T, 2> >();
  // dgm->reserve(data_view.shape(0) * 2);
  // {
  //   nb::gil_scoped_release release;
  //   Gudhi::persistent_cohomology::compute_persistence_of_function_on_line(r, [&](T b, T d) {
  //     dgm->emplace_back(std::array<T, 2>{b,d});
  //   });
  // }
  // dgm->shrink_to_fit();
  // return nb::ndarray<T, nb::numpy>(dgm->data(), {dgm->size(), 2}, nb::capsule(dgm, [](void* p) noexcept {
  //   delete reinterpret_cast<std::vector<std::array<T, 2> >*>(p);
  // }));
  // std::vector<T>* dgm = new std::vector<T>();
  // // dgm->reserve(data_view.shape(0) * 2);  // rough upper bound, but it makes a difference in performance
  // {
  //   nb::gil_scoped_release release;
  //   Gudhi::persistent_cohomology::compute_persistence_of_function_on_line(r, [&](T b, T d) {
  //     dgm->push_back(b);
  //     dgm->push_back(d);
  //   });
  // }
  // // dgm->shrink_to_fit();
  // return _wrap_as_numpy_array(dgm, dgm->size() / 2, 2);
}

nb::list wrap_persistence_2d(nb::ndarray<double, nb::ndim<2>, nb::c_contig> data, double min_persistence)
{
  std::vector<double>* dgm0 = new std::vector<double>();
  std::vector<double>* dgm1 = new std::vector<double>();
  // dgm0->reserve(data.shape(0) * data.shape(1) * 2);  // rough upper bound
  // dgm1->reserve(data.shape(0) * data.shape(1) * 2);  // rough upper bound
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
  // if (dgm0->size() < dgm0->capacity() / 2) dgm0->shrink_to_fit();
  // if (dgm1->size() < dgm1->capacity() / 2) dgm1->shrink_to_fit();
  nb::list ret;
  ret.append(_wrap_as_numpy_array(dgm0, dgm0->size() / 2, 2));
  ret.append(_wrap_as_numpy_array(dgm1, dgm1->size() / 2, 2));
  return ret;
}

NB_MODULE(_pers_cub_low_dim_ext, m)
{
  m.attr("__license__") = "MIT";
  m.def("_persistence_on_a_line", &wrap_persistence_1d<float>, nb::arg().noconvert());
  m.def("_persistence_on_a_line", &wrap_persistence_1d<double>);
  m.def("_persistence_on_rectangle_from_top_cells", &wrap_persistence_2d);
}

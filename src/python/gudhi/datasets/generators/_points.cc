/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hind Montassif
 *
 *    Copyright (C) 2021 Inria
 *
 *    Modification(s):
 *      - 2025/01 Vincent Rouvreau: Use nanobind instead of PyBind11 for python bindings
 *      - YYYY/MM Author: Description of the modification
 */

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>

#include <gudhi/random_point_generators.h>
#include <gudhi/Debug_utils.h>

#include <CGAL/Epick_d.h>

#include <vector>

namespace nb = nanobind;

typedef CGAL::Epick_d< CGAL::Dynamic_dimension_tag > Kern;

auto generate_points_on_sphere(const size_t n_samples, const int ambient_dim, double radius, std::string sample) {
  if (sample != "random") {
    throw nb::value_error("This sample type is not supported");
  }

  std::vector<typename Kern::Point_d> points_generated;
  {
    nb::gil_scoped_release release;
    points_generated = Gudhi::generate_points_on_sphere_d<Kern>(n_samples, ambient_dim, radius);
  }

  // Reserve sufficient memory space to copy data
  std::vector<double> points(n_samples * ambient_dim);
  double *ptr = points.data();
  for (size_t i = 0; i < n_samples; i++)
    for (int j = 0; j < ambient_dim; j++)
      ptr[i*ambient_dim+j] = points_generated[i][j];

  return nb::ndarray<nb::numpy, double>(points.data(), { n_samples, (size_t)ambient_dim }).cast();
}

auto generate_points_on_torus(size_t n_samples, int dim, std::string sample) {
  if ( (sample != "random") && (sample != "grid")) {
    throw nb::value_error("This sample type is not supported");
  }

  std::vector<typename Kern::Point_d> points_generated;
  {
    nb::gil_scoped_release release;
    points_generated = Gudhi::generate_points_on_torus_d<Kern>(n_samples, dim, sample);
  }

  size_t npoints = points_generated.size();
  GUDHI_CHECK(2*dim == points_generated[0].size(), "Py array second dimension not matching the double torus dimension");

  // Reserve sufficient memory space to copy data
  std::vector<double> points(npoints * (size_t)2*dim);
  double *ptr = points.data();
  for (size_t i = 0; i < npoints; i++)
    for (int j = 0; j < 2*dim; j++)
      ptr[i*(2*dim)+j] = points_generated[i][j];

  return nb::ndarray<nb::numpy, double>(points.data(), { npoints, (size_t)2*dim }).cast();
}

NB_MODULE(_points_ext, m) {
    m.attr("__license__") = "LGPL v3";

    m.def("sphere", &generate_points_on_sphere,
          nb::arg("n_samples"), nb::arg("ambient_dim"),
          nb::arg("radius") = 1., nb::arg("sample") = "random",
          R"pbdoc(
          Generate random i.i.d. points uniformly on a (d-1)-sphere in R^d

          :param n_samples: The number of points to be generated.
          :type n_samples: integer
          :param ambient_dim: The ambient dimension d.
          :type ambient_dim: integer
          :param radius: The radius. Default value is `1.`.
          :type radius: float
          :param sample: The sample type. Default and only available value is `"random"`.
          :type sample: string
          :returns: the generated points on a sphere.
          )pbdoc");

    m.def("ctorus", &generate_points_on_torus,
          nb::arg("n_samples"), nb::arg("dim"), nb::arg("sample") = "random",
          R"pbdoc(
          Generate random i.i.d. points on a d-torus in R^2d or as a grid

          :param n_samples: The number of points to be generated.
          :type n_samples: integer
          :param dim: The dimension of the torus on which points would be generated in R^2*dim.
          :type dim: integer
          :param sample: The sample type. Available values are: `"random"` and `"grid"`. Default value is `"random"`.
          :type sample: string
          :returns: the generated points on a torus.

          The shape of returned numpy array is:

          If sample is 'random': (n_samples, 2*dim).

          If sample is 'grid': (⌊n_samples**(1./dim)⌋**dim, 2*dim), where shape[0] is rounded down to the closest perfect 'dim'th power.
          )pbdoc");
}

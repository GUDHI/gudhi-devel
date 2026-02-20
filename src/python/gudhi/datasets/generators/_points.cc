/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hind Montassif
 *
 *    Copyright (C) 2021 Inria
 *
 *    Modification(s):
 *      - 2025/01 Vincent Rouvreau: Use nanobind instead of PyBind11 for python bindings
 *      - 2026/01 Vincent Rouvreau: Add Gudhi::random::Random_generator support to set the seed
 *      - YYYY/MM Author: Description of the modification
 */

#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include <CGAL/Epick_d.h>

// For Windows, where points is a client of random here (random.dll is the provider)
#define RANDOM_DLL_IMPORT
#include <gudhi/random_point_generators.h>

#include <gudhi/Debug_utils.h>
#include <python_interfaces/numpy_utils.h>
#include <python_interfaces/random_utils.h>

namespace nb = nanobind;

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> Kern;

auto generate_points_on_sphere(const size_t n_samples, const int ambient_dim, double radius, std::string sample,
                               Gudhi::random::Random_generator& rng) {
  if (sample != "random") {
    throw nb::value_error("This sample type is not supported");
  }

  std::vector<typename Kern::Point_d> points_generated;
  {
    nb::gil_scoped_release release;
    points_generated = Gudhi::generate_points_on_sphere_d<Kern>(n_samples, ambient_dim, radius, 0.,
                                                                rng.get_default_random());
  }

  // Reserve sufficient memory space to copy data
  auto points = new double[n_samples * ambient_dim];
  for (size_t i = 0; i < n_samples; i++)
    for (int j = 0; j < ambient_dim; j++) points[i * ambient_dim + j] = points_generated[i][j];

  return _wrap_as_numpy_array(points, n_samples, ambient_dim);
}

auto generate_points_on_torus(size_t n_samples, int dim, std::string sample, Gudhi::random::Random_generator& rng)
{
  if ((sample != "random") && (sample != "grid")) {
    throw nb::value_error("This sample type is not supported");
  }

  std::vector<typename Kern::Point_d> points_generated;
  {
    nb::gil_scoped_release release;
    points_generated = Gudhi::generate_points_on_torus_d<Kern>(n_samples, dim, sample, 0., rng.get_default_random());
  }

  size_t npoints = points_generated.size();
  GUDHI_CHECK(2 * dim == points_generated[0].size(),
              "Py array second dimension not matching the double torus dimension");

  // Reserve sufficient memory space to copy data
  auto points = new double[npoints * 2 * dim];
  for (size_t i = 0; i < npoints; i++)
    for (int j = 0; j < 2 * dim; j++) points[i * (2 * dim) + j] = points_generated[i][j];

  return _wrap_as_numpy_array(points, npoints, 2 * dim);
}

NB_MODULE(_points_ext, m)
{
  m.attr("__license__") = "LGPL v3";

  m.def("_sphere",
        &generate_points_on_sphere,
        nb::arg("n_samples"),
        nb::arg("ambient_dim"),
        nb::arg("radius"),
        nb::arg("sample"),
        nb::arg("rng"));

  m.def("_ctorus",
        &generate_points_on_torus,
        nb::arg("n_samples"),
        nb::arg("dim"),
        nb::arg("sample"),
        nb::arg("rng"));
}

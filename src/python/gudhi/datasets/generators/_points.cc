/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hind Montassif
 *
 *    Copyright (C) 2021 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <gudhi/random_point_generators.h>
#include <gudhi/Debug_utils.h>

#include <CGAL/Epick_d.h>

namespace py = pybind11;


typedef CGAL::Epick_d< CGAL::Dynamic_dimension_tag > Kern;

py::array_t<double> generate_points_on_sphere(size_t n_samples, int ambient_dim, double radius, std::string sample) {

    if (sample != "random") {
        throw pybind11::value_error("This sample type is not supported");
    }

    py::array_t<double> points({n_samples, (size_t)ambient_dim});

    py::buffer_info buf = points.request();
    double *ptr = static_cast<double *>(buf.ptr);

    GUDHI_CHECK(n_samples == buf.shape[0], "Py array first dimension not matching n_samples on sphere");
    GUDHI_CHECK(ambient_dim == buf.shape[1], "Py array second dimension not matching the ambient space dimension");


    std::vector<typename Kern::Point_d> points_generated;

    {
        py::gil_scoped_release release;
        points_generated = Gudhi::generate_points_on_sphere_d<Kern>(n_samples, ambient_dim, radius);
    }

    for (size_t i = 0; i < n_samples; i++)
        for (int j = 0; j < ambient_dim; j++)
            ptr[i*ambient_dim+j] = points_generated[i][j];

    return points;
}

py::array_t<double> generate_points_on_torus(size_t n_samples, int dim, std::string sample) {

    if ( (sample != "random") && (sample != "grid")) {
        throw pybind11::value_error("This sample type is not supported");
    }

    std::vector<typename Kern::Point_d> points_generated;

    {
        py::gil_scoped_release release;
        points_generated = Gudhi::generate_points_on_torus_d<Kern>(n_samples, dim, sample);
    }

    size_t npoints = points_generated.size();

    GUDHI_CHECK(2*dim == points_generated[0].size(), "Py array second dimension not matching the double torus dimension");

    py::array_t<double> points({npoints, (size_t)2*dim});

    py::buffer_info buf = points.request();
    double *ptr = static_cast<double *>(buf.ptr);

    for (size_t i = 0; i < npoints; i++)
        for (int j = 0; j < 2*dim; j++)
            ptr[i*(2*dim)+j] = points_generated[i][j];

    return points;
}

PYBIND11_MODULE(_points, m) {
    m.attr("__license__") = "LGPL v3";

    m.def("sphere", &generate_points_on_sphere,
          py::arg("n_samples"), py::arg("ambient_dim"),
          py::arg("radius") = 1., py::arg("sample") = "random",
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
          py::arg("n_samples"), py::arg("dim"), py::arg("sample") = "random",
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

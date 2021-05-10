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


py::array_t<double> generate_points_on_torus(size_t n_samples, int dim, bool uniform) {

    std::vector<typename Kern::Point_d> points_generated;

    {
        py::gil_scoped_release release;
        points_generated = Gudhi::generate_points_on_torus_d<Kern>(n_samples, dim, uniform);
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
      m.def("torus", &generate_points_on_torus,
          py::arg("n_samples"), py::arg("dim"), py::arg("uniform") = false,
          R"pbdoc(
    Generate random i.i.d. points on a d-torus in R^2d

    :param n_samples: The number of points to be generated.
    :type n_samples: integer
    :param dim: The dimension.
    :type dim: integer
    :param uniform: A flag to define if the points generation is uniform (i.e generated as a grid).
    :type uniform: bool
    :rtype: numpy array of float
    :returns: the generated points on a torus.
    )pbdoc");
}

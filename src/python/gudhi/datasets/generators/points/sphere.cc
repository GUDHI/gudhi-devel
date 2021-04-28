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

py::array_t<double> generate_points_on_sphere(size_t num_points, int dim, double radius) {

    py::array_t<double> points({num_points, (size_t)dim});

    py::buffer_info buf = points.request();
    double *ptr = static_cast<double *>(buf.ptr);

    GUDHI_CHECK(num_points == buf.shape[0], "Py array first dimension not matching num_points on sphere");
    GUDHI_CHECK(dim == buf.shape[1], "Py array second dimension not matching the ambient space dimension");


    py::gil_scoped_release release;
    auto points_generated = Gudhi::generate_points_on_sphere_d<Kern>(num_points, dim, radius);

    for (size_t i = 0; i < num_points; i++)
        for (int j = 0; j < dim; j++)
            ptr[i*dim+j] = points_generated[i][j];

    return points;
}

PYBIND11_MODULE(sphere, m) {
      m.attr("__license__") = "LGPL v3";
      m.def("generate_random_points", &generate_points_on_sphere,
          py::arg("num_points"), py::arg("dim"), py::arg("radius") = 1,
          R"pbdoc(
    Generate random i.i.d. points uniformly on a (d-1)-sphere in R^d

    :param num_points: The number of points to be generated.
    :type num_points: unsigned integer
    :param dim: The dimension.
    :type dim: integer
    :param radius: The radius.
    :type radius: float
    :rtype: numpy array of float
    :returns: the generated points on a sphere.
    )pbdoc");
}

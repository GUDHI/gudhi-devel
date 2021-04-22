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

#include <CGAL/Epick_d.h>

namespace py = pybind11;


typedef CGAL::Epick_d< CGAL::Dynamic_dimension_tag > Kern;

template <typename Kernel>
py::array_t<double> generate_points_on_sphere(py::object num_points, py::object dim, py::object radius) {
    int npoints = num_points.cast<int>();
    int d = dim.cast<int>();
    double rad  = radius.cast<double>();
    
    py::gil_scoped_release release;

    auto points_generated = Gudhi::generate_points_on_sphere_d<Kernel>(npoints, d, rad);
    
    py::gil_scoped_acquire acquire;
    
    py::array_t<double> points({npoints, d});
 
    py::buffer_info buf = points.request();

    double *ptr = static_cast<double *>(buf.ptr);

    assert(npoints == buf.shape[0]);
    assert(d == buf.shape[1]);
    

    for (size_t i = 0; i < (size_t)npoints; i++)
        for (size_t j = 0; j < (size_t)d; j++)
            ptr[i*d+j] = points_generated.at(i).at(j);

    return points;
}

PYBIND11_MODULE(random_point_generators, m) {
      m.attr("__license__") = "LGPL v3";
      m.def("generate_points_on_sphere_d", &generate_points_on_sphere<Kern>,
          py::arg("num_points"), py::arg("dim"), py::arg("radius"),
          R"pbdoc(
    Generate points on a sphere

    :param num_points: The number of points to be generated.
    :type num_points: integer
    :param dim: The sphere dimension.
    :type dim: integer
    :param radius: The sphere radius.
    :type radius: float
    :rtype: numpy array of points
    :returns: the generated points on a sphere.
    )pbdoc");
}

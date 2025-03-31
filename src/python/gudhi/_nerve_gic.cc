/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2018 Inria
 *
 *    Modification(s):
 *      - 2025/03 Hannah Schreiber: Use nanobind instead of Cython for python bindings
 *      - YYYY/MM Author: Description of the modification
 */

#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>

#include <gudhi/distance_functions.h>
#include <gudhi/GIC.h>

#include "Simplex_tree_interface.h"

namespace Gudhi {

namespace cover_complex {

class Nerve_gic_interface : public Cover_complex<std::vector<double>>
{
 public:
  void create_simplex_tree(Simplex_tree_interface* simplex_tree) { create_complex(*simplex_tree); }

  void set_cover_from_euclidean_voronoi(int m) { set_cover_from_Voronoi(Gudhi::Euclidean_distance(), m); }

  double set_graph_from_automatic_euclidean_rips(int N)
  {
    return set_graph_from_automatic_rips(Gudhi::Euclidean_distance(), N);
  }

  void set_graph_from_euclidean_rips(double threshold) { set_graph_from_rips(threshold, Gudhi::Euclidean_distance()); }
};

}  // namespace cover_complex

}  // namespace Gudhi

namespace nb = nanobind;

using NGI = Gudhi::cover_complex::Nerve_gic_interface;

NB_MODULE(_nerve_gic_ext, m)
{
  m.attr("__license__") = "GPL v3";

  nb::class_<NGI>(m, "_Nerve_gic_interface")
      .def(nb::init<>())
      .def("compute_confidence_level_from_distance",
           &NGI::compute_confidence_level_from_distance,
           nb::arg("distance"),
           R"pbdoc(
        """Computes the confidence level of a specific bottleneck distance
        threshold.

        :param distance: Bottleneck distance.
        :type distance: double

        :rtype: double
        :returns: Confidence level.
        """
        )pbdoc")
      .def("compute_distance_from_confidence_level",
           &NGI::compute_distance_from_confidence_level,
           nb::arg("alpha"),
           R"pbdoc(
        """Computes the bottleneck distance threshold corresponding to a
        specific confidence level.

        :param alpha: Confidence level.
        :type alpha: double

        :rtype: double
        :returns: Bottleneck distance.
        """
        )pbdoc")
      .def("compute_distribution", &NGI::compute_distribution, nb::arg("N") = 100, R"pbdoc(
        """Computes bootstrapped distances distribution.

        :param N: Loop number (default value is 100).
        :type alpha: int
        """
        )pbdoc")
      .def("compute_p_value", &NGI::compute_p_value, R"pbdoc(
        """Computes the p-value, i.e. the opposite of the confidence level of
        the largest bottleneck distance preserving the points in the
        persistence diagram of the output simplicial complex.

        :rtype: double
        :returns: p-value.
        """
        )pbdoc")
      .def("compute_PD", &NGI::compute_PD, R"pbdoc(
        """Computes the extended persistence diagram of the complex.
        """
        )pbdoc")
      .def("find_simplices", &NGI::find_simplices, R"pbdoc(
        """Computes the simplices of the simplicial complex.
        """
        )pbdoc")
      .def("create_simplex_tree", &NGI::create_simplex_tree, nb::arg("simplex_tree"))
      .def("read_point_cloud", &NGI::read_point_cloud, nb::arg("off_file_name"))
      .def("set_automatic_resolution", &NGI::set_automatic_resolution, R"pbdoc(
        """Computes the optimal length of intervals (i.e. the smallest interval
        length avoiding discretization artifacts - see :cite:`Carriere17c`) for a
        functional cover.

        :rtype: double
        :returns: resolution interval length used to compute the cover.
        """
        )pbdoc")
      .def("set_color_from_coordinate", &NGI::set_color_from_coordinate, nb::arg("k") = 0, R"pbdoc(
        """Computes the function used to color the nodes of the simplicial
        complex from the k-th coordinate.

        :param k: Coordinate to use (start at 0). Default value is 0.
        :type k: int
        """
        )pbdoc")
      .def("set_color_from_file", &NGI::set_color_from_file, nb::arg("color_file_name"))
      .def("set_color_from_range", &NGI::set_color_from_range, nb::arg("color"), R"pbdoc(
        """Computes the function used to color the nodes of the simplicial
        complex from a vector stored in memory.

        :param color: Input vector of values.
        :type color: vector[double]
        """
        )pbdoc")
      .def("set_cover_from_file", &NGI::set_cover_from_file, nb::arg("cover_file_name"))
      .def("set_cover_from_range",
           &NGI::set_cover_from_range<std::vector<std::vector<int>>>,
           nb::arg("assignments"),
           R"pbdoc(
        """Creates a cover C from a vector stored in memory.

        :param assignments: Vector containing the assignments of the points to their corresponding cover elements. For instance, if the i-th point belongs to the 1st and 3rd cover elements, then assignments[i] = [1,3].
        :type assignments: List[List[int]]
        """
        )pbdoc")
      .def("set_cover_from_function", &NGI::set_cover_from_function, R"pbdoc(
        """Creates a cover C from the preimages of the function f.
        """
        )pbdoc")
      .def("set_cover_from_Voronoi", &NGI::set_cover_from_euclidean_voronoi, nb::arg("m") = 100, R"pbdoc(
        """Creates the cover C from the Voronoï cells of a subsampling of the
        point cloud.

        :param m: Number of points in the subsample. Default value is 100.
        :type m: int
        """
        )pbdoc")
      .def("set_function_from_coordinate", &NGI::set_function_from_coordinate, nb::arg("k"), R"pbdoc(
        """Creates the function f from the k-th coordinate of the point cloud.

        :param k: Coordinate to use (start at 0).
        :type k: int
        """
        )pbdoc")
      .def("set_function_from_file", &NGI::set_function_from_file, nb::arg("func_file_name"))
      .def("set_function_from_range",
           &NGI::set_function_from_range<std::vector<double>>,
           nb::arg("function"),
           R"pbdoc(
        """Creates the function f from a vector stored in memory.

        :param function: Input vector of values.
        :type function: vector[double]
        """
        )pbdoc")
      .def("set_gain", &NGI::set_gain, nb::arg("g") = 0.3, R"pbdoc(
        """Sets a gain from a value stored in memory.

        :param g: Gain (default value is 0.3).
        :type g: double
        """
        )pbdoc")
      .def("set_graph_from_automatic_rips",
           &NGI::set_graph_from_automatic_euclidean_rips,
           nb::arg("N") = 100,
           R"pbdoc(
        """Creates a graph G from a Rips complex whose threshold value is
        automatically tuned with subsampling - see :cite:`Carriere17c`.

        :param N: Number of subsampling iteration (the default reasonable value is 100, but there is no guarantee on how to choose it).
        :type N: int
        :rtype: double
        :returns: Delta threshold used for computing the Rips complex.
        """
        )pbdoc")
      .def("set_graph_from_file", &NGI::set_graph_from_file, nb::arg("graph_file_name"))
      .def("set_graph_from_OFF", &NGI::set_graph_from_OFF, R"pbdoc(
        """Creates a graph G from the triangulation given by the input OFF
        file.
        """
        )pbdoc")
      .def("set_graph_from_euclidean_rips",
           &NGI::set_graph_from_euclidean_rips,
           nb::arg("threshold"),
           R"pbdoc(
        """Creates a graph G from a Rips complex.

        :param threshold: Threshold value for the Rips complex.
        :type threshold: double
        """
        )pbdoc")
      .def("set_mask", &NGI::set_mask, nb::arg("nodemask"), R"pbdoc(
        """Sets the mask, which is a threshold integer such that nodes in the
        complex that contain a number of data points which is less than or
        equal to this threshold are not displayed.

        :param nodemask: Threshold.
        :type nodemask: int
        """
        )pbdoc")
      .def("set_resolution_with_interval_length",
           &NGI::set_resolution_with_interval_length,
           nb::arg("resolution"),
           R"pbdoc(
        """Sets a length of intervals from a value stored in memory.

        :param resolution: Length of intervals.
        :type resolution: double
        """
        )pbdoc")
      .def("set_resolution_with_interval_number",
           &NGI::set_resolution_with_interval_number,
           nb::arg("resolution"),
           R"pbdoc(
        """Sets a number of intervals from a value stored in memory.

        :param resolution: Number of intervals.
        :type resolution: int
        """
        )pbdoc")
      .def("set_subsampling", &NGI::set_subsampling, nb::arg("constant"), nb::arg("power"), R"pbdoc(
        """Sets the constants used to subsample the data set. These constants
        are explained in :cite:`Carriere17c`.

        :param constant: Constant.
        :type constant: double

        :param power: Power.
        :type resolution: double
        """
        )pbdoc")
      .def("set_type", &NGI::set_type, nb::arg("type"), R"pbdoc(
        """Specifies whether the type of the output simplicial complex.

        :param type: either "GIC" or "Nerve".
        :type type: string
        """
        )pbdoc")
      .def("set_verbose", &NGI::set_verbose, nb::arg("verbose") = false, R"pbdoc(
        """Specifies whether the program should display information or not.

        :param verbose: true = display info, false = do not display info.
        :type verbose: boolean
        """
        )pbdoc")
      .def("subpopulation", &NGI::subpopulation, nb::arg("c"), R"pbdoc(
        """Returns the data subset corresponding to a specific node of the
        created complex.

        :param c: ID of the node.
        :type c: int

        :rtype: vector[int]
        :returns: Vector of IDs of data points.
        """
        )pbdoc")
      .def("subcolor", &NGI::subcolor, nb::arg("c"), R"pbdoc(
        """Returns the mean color value corresponding to a specific node of the
        created complex.

        :param c: ID of the node.
        :type c: int

        :rtype: float
        :returns: Mean color value of data points.
        """
        )pbdoc")
      .def("write_info", &NGI::write_info, R"pbdoc(
        """Creates a .txt file called SC.txt describing the 1-skeleton, which can
        then be plotted with e.g. KeplerMapper.
        """
        )pbdoc")
      .def("plot_dot", &NGI::plot_DOT, R"pbdoc(
        """Creates a .dot file called SC.dot for neato (part of the graphviz
        package) once the simplicial complex is computed to get a visualization of
        its 1-skeleton in a .pdf file.
        """
        )pbdoc")
      .def("plot_off", &NGI::plot_OFF, R"pbdoc(
        """Creates a .off file called SC.off for 3D visualization, which contains
        the 2-skeleton of the GIC. This function assumes that the cover has been
        computed with Voronoi. If data points are in 1D or 2D, the remaining
        coordinates of the points embedded in 3D are set to 0.
        """
        )pbdoc")
      .def("set_point_cloud_from_range", &NGI::set_point_cloud_from_range, nb::arg("cloud"), R"pbdoc(
        """ Reads and stores the input point cloud from a vector stored in
        memory.

        :param cloud: Input vector containing the point cloud.
        :type cloud: vector[vector[double]]
        """
        )pbdoc")
      .def("set_distances_from_range", &NGI::set_distances_from_range, nb::arg("distance_matrix"), R"pbdoc(
        """ Reads and stores the input distance matrix from a vector stored in
        memory.

        :param distance_matrix: Input vector containing the distance matrix.
        :type distance_matrix: vector[vector[double]]
        """
        )pbdoc");
}

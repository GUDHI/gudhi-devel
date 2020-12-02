/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Points_off_io.h>
#include <gudhi/reader_utils.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <vector>
#include <utility>  // for std::pair
#include <string>

namespace py = pybind11;

std::vector<std::vector<double>> read_points_from_OFF_file(const std::string& off_file) {
  Gudhi::Points_off_reader<std::vector<double>> off_reader(off_file);
  return off_reader.get_point_cloud();
}


PYBIND11_MODULE(reader_utility, m) {
      m.attr("__license__") = "MIT";
      PYBIND11_NUMPY_DTYPE(std::pair<double, double>, first, second);
      // read_points_from_off_file
      m.def("read_points_from_OFF_file", &read_points_from_OFF_file,
          py::arg("off_file"),
          R"pbdoc(
    Read points from OFF file.

    :param off_file: An OFF file style name.
    :type off_file: string

    :returns:  The point set.
    :rtype: List[List[float]]
    )pbdoc");
      // read_lower_triangular_matrix_from_csv_file
      m.def("read_matrix_from_csv_file", &Gudhi::read_lower_triangular_matrix_from_csv_file<double>,
          py::arg("filename"), py::arg("separator") = static_cast<const char>(';'),
          R"pbdoc(
    Read lower triangular matrix from a CSV style file.

    :param csv_file: A CSV file style name.
    :type csv_file: string
    :param separator: The value separator in the CSV file. Default value is ';'
    :type separator: char

    :returns:  The lower triangular matrix.
    :rtype: List[List[float]]
    )pbdoc");
      // read_persistence_intervals_grouped_by_dimension
      m.def("read_pers_intervals_grouped_by_dimension", &Gudhi::read_persistence_intervals_grouped_by_dimension,
          py::arg("persistence_file"),
          R"pbdoc(
    Reads a file containing persistence intervals.
    Each line might contain 2, 3 or 4 values: [[field] dimension] birth death
    The return value is a `dict(dim, list(tuple(birth, death)))`
    where `dim` is an `int`, `birth` a `float`, and `death` a `float`.
    Note: the function does not check that birth <= death.

    :param persistence_file: A persistence file style name.
    :type persistence_file: string

    :returns:  The persistence pairs grouped by dimension.
    :rtype: Dict[int, List[Tuple[float, float]]]
    )pbdoc");
      // read_persistence_intervals_in_dimension
      m.def("read_pers_intervals_in_dimension", [](const std::string& persistence_file, int only_this_dim) {
            std::vector<std::pair<double, double>> pers = Gudhi::read_persistence_intervals_in_dimension(persistence_file, only_this_dim);
            return py::array(pers.size(), pers.data());
          },
          py::arg("persistence_file"), py::arg("only_this_dim") = static_cast<int>(-1),
          R"pbdoc(
    Reads a file containing persistence intervals.
    Each line of persistence_file might contain 2, 3 or 4 values:
    [[field] dimension] birth death
    Note: the function does not check that birth <= death.

    :param persistence_file: A persistence file style name.
    :type persistence_file: string
    :param only_this_dim: The specific dimension. Default value is -1.
        If `only_this_dim` = -1, dimension is ignored and all lines are returned.
        If `only_this_dim` is >= 0, only the lines where dimension =
        `only_this_dim` (or where dimension is not specified) are returned.
    :type only_this_dim: int.

    :returns:  The persistence intervals.
    :rtype: numpy array of dimension 2
    )pbdoc");
}

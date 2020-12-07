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

#include <boost/function_output_iterator.hpp>  // for boost::make_function_output_iterator

#include <vector>
#include <array>
#include <string>

namespace py = pybind11;


PYBIND11_MODULE(reader_utils, m) {
      m.attr("__license__") = "MIT";
      m.def("read_points_from_off_file", [](const std::string& off_file) {
            Gudhi::Points_off_reader<std::vector<double>> off_reader(off_file);
            return off_reader.get_point_cloud();
          },
        py::arg("off_file"),
        R"pbdoc(
    Read points from OFF file.

    :param off_file: Path to an `OFF <fileformats.html#off-file-format>`_ file.
    :type off_file: string

    :returns:  The point set.
    :rtype: List[List[float]]
        )pbdoc"
      );
      m.def("read_lower_triangular_matrix_from_csv_file", &Gudhi::read_lower_triangular_matrix_from_csv_file<double>,
        py::arg("csv_file"), py::arg("separator") = ';',
        R"pbdoc(
    Read lower triangular matrix from a CSV style file.

    :param csv_file: Path to a CSV file.
    :type csv_file: string
    :param separator: The value separator in the CSV file. Default value is ';'
    :type separator: char

    :returns:  The lower triangular matrix.
    :rtype: List[List[float]]
        )pbdoc"
      );
      m.def("read_persistence_intervals_grouped_by_dimension", &Gudhi::read_persistence_intervals_grouped_by_dimension,
        py::arg("persistence_file"),
        R"pbdoc(
    Reads a file containing persistence intervals.
    Each line might contain 2, 3 or 4 values: [[field] dimension] birth death
    The return value is a `dict(dim, list(tuple(birth, death)))`
    where `dim` is an `int`, `birth` a `float`, and `death` a `float`.
    Note: the function does not check that birth <= death.

    :param persistence_file: Path to a `persistence diagram <fileformats.html#persistence-diagram>`_ file.
    :type persistence_file: string

    :returns:  The persistence pairs grouped by dimension.
    :rtype: Dict[int, List[Tuple[float, float]]]
        )pbdoc"
      );
      m.def("read_persistence_intervals_in_dimension", [](const std::string& persistence_file, int only_this_dim) {
            // std::vector<std::pair<double, double>> as in read_persistence_intervals_in_dimension
            // was not giving the good results because of memory alignment
            // std::vector is mandatory as we don't know in advance the number of persistence intervals to be read
            std::vector<std::array<double,2>> ret;
            Gudhi::read_persistence_intervals_and_dimension(persistence_file,
              boost::make_function_output_iterator([only_this_dim, &ret](std::tuple<int, double, double> t) {
                if (only_this_dim == get<0>(t) || only_this_dim == -1) ret.push_back({get<1>(t), get<2>(t)});
              }));
            return py::array_t<double>(std::vector<std::ptrdiff_t>{static_cast<std::ptrdiff_t>(ret.size()), 2},
              &ret[0][0]);
          },
        py::arg("persistence_file"), py::arg("only_this_dim") = static_cast<int>(-1),
        R"pbdoc(
    Reads a file containing persistence intervals.
    Each line of persistence_file might contain 2, 3 or 4 values:
    [[field] dimension] birth death
    Note: the function does not check that birth <= death.

    :param persistence_file: Path to a `persistence diagram <fileformats.html#persistence-diagram>`_ file.
    :type persistence_file: string
    :param only_this_dim: The specific dimension. Default value is -1.
        If `only_this_dim` = -1, dimension is ignored and all lines are returned.
        If `only_this_dim` is >= 0, only the lines where dimension =
        `only_this_dim` (or where dimension is not specified) are returned.
    :type only_this_dim: int.

    :returns:  The persistence intervals.
    :rtype: numpy array of dimension 2
        )pbdoc"
      );
}

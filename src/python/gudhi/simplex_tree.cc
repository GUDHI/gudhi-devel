/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 */

#include "Simplex_tree_interface.h"

#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;

using gsti = Gudhi::Simplex_tree_interface;

NB_MODULE(simplex_tree_ext, m) {
    m.attr("__license__") = "GPL v3";

    nb::class_<gsti>(m, "Simplex_tree_interface")
            .def(nb::init<>())
            .def("simplex_filtration",
                    &gsti::simplex_filtration,
                    nb::arg("simplex"),
                    R"pbdoc(TODO)pbdoc")
            .def("assign_simplex_filtration",
                    &gsti::assign_simplex_filtration,
                    nb::arg("simplex"), nb::arg("filtration"),
                    R"pbdoc(TODO)pbdoc")
            .def("initialize_filtration",
                    &gsti::initialize_filtration,
                    R"pbdoc(TODO)pbdoc")
            .def("num_vertices",
                    &gsti::num_vertices,
                    R"pbdoc(TODO)pbdoc")
            .def("num_simplices",
                    nb::overload_cast<>(&gsti::num_simplices, nb::const_),
                    R"pbdoc(TODO)pbdoc")
            .def("is_empty",
                    &gsti::is_empty,
                    R"pbdoc(TODO)pbdoc")
            .def("set_dimension",
                    &gsti::set_dimension,
                    nb::arg("dimension"), nb::arg("exact")=true,
                    R"pbdoc(TODO)pbdoc")
            .def("dimension",
                    nb::overload_cast<>(&gsti::dimension, nb::const_),
                    R"pbdoc(TODO)pbdoc")
            .def("upper_bound_dimension",
                    &gsti::upper_bound_dimension,
                    R"pbdoc(TODO)pbdoc")
            .def("find_simplex",
                    &gsti::find_simplex,
                    nb::arg("simplex"),
                    R"pbdoc(TODO)pbdoc")
            .def("insert",
                    &gsti::insert,
                    nb::arg("simplex"), nb::arg("filtration"),
                    R"pbdoc(TODO)pbdoc")
            .def("insert_matrix",
                    &gsti::insert_matrix,
                    nb::arg("filtrations"), nb::arg("n"), nb::arg("stride0"), nb::arg("stride1"), nb::arg("max_filtration"),
                    R"pbdoc(TODO)pbdoc")
            .def("insert_batch_vertices",
                    &gsti::insert_batch_vertices<std::vector<int>>,
                    R"pbdoc(TODO)pbdoc")
            .def("get_star",
                    &gsti::get_star,
                    nb::arg("simplex"),
                    R"pbdoc(TODO)pbdoc")
            .def("get_cofaces",
                    &gsti::get_cofaces,
                    nb::arg("simplex"), nb::arg("dimension"),
                    R"pbdoc(TODO)pbdoc")
            .def("expansion",
                    &gsti::expansion,
                    nb::arg("max_dim"),
                    R"pbdoc(TODO)pbdoc")
            .def("remove_maximal_simplex",
                    &gsti::remove_maximal_simplex,
                    nb::arg("simplex"),
                    R"pbdoc(TODO)pbdoc")
            .def("prune_above_filtration",
                    &gsti::prune_above_filtration,
                    nb::arg("filtration"),
                    R"pbdoc(TODO)pbdoc")
            .def("prune_above_dimension",
                    &gsti::prune_above_dimension,
                    nb::arg("dimension"),
                    R"pbdoc(TODO)pbdoc")
            .def("make_filtration_non_decreasing",
                    &gsti::make_filtration_non_decreasing,
                    R"pbdoc(TODO)pbdoc")
            .def("compute_extended_filtration",
                    &gsti::compute_extended_filtration,
                    R"pbdoc(TODO)pbdoc")
            .def("collapse_edges",
                    &gsti::collapse_edges,
                    nb::arg("nb_collapse_iteration"),
                    R"pbdoc(TODO)pbdoc")
            .def("reset_filtration",
                    &gsti::reset_filtration,
                    nb::arg("filtration"), nb::arg("dimension")=0,
                    R"pbdoc(TODO)pbdoc")
            .def("get_simplex_and_filtration",
                    &gsti::get_simplex_and_filtration,
                    nb::arg("f_simplex"),
                    R"pbdoc(TODO)pbdoc")
            .def("get_simplices_iterator_begin",
                    &gsti::get_simplices_iterator_begin,
                    R"pbdoc(TODO)pbdoc")
            .def("get_simplices_iterator_end",
                    &gsti::get_simplices_iterator_end,
                    R"pbdoc(TODO)pbdoc")
            .def("get_filtration_iterator_begin",
                    &gsti::get_filtration_iterator_begin,
                    R"pbdoc(TODO)pbdoc")
            .def("get_filtration_iterator_end",
                    &gsti::get_filtration_iterator_end,
                    R"pbdoc(TODO)pbdoc")
            .def("get_skeleton_iterator_begin",
                    &gsti::get_skeleton_iterator_begin,
                    nb::arg("dimension"),
                    R"pbdoc(TODO)pbdoc")
            .def("get_skeleton_iterator_end",
                    &gsti::get_skeleton_iterator_end,
                    nb::arg("dimension"),
                    R"pbdoc(TODO)pbdoc")
            .def("get_boundary_iterators",
                    &gsti::get_boundary_iterators,
                    nb::arg("simplex"),
                    R"pbdoc(TODO)pbdoc")
            ;
}

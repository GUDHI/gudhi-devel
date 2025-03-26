/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 */

#include <../include/Simplex_tree_interface.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>

namespace nb=nanobind;

NB_MODULE(simplex_tree_ext, m) {
      m.attr("__license__") = "GPL v3";

      nb::class_<Gudhi::Simplex_tree_interface>(m, "Simplex_tree_interface")
          .def(nb::init<>())
          .def("simplex_filtration",
               &Gudhi::Simplex_tree_interface::simplex_filtration,
               nb::arg("simplex"),
               R"pbdoc(TODO)pbdoc")
          .def("assign_simplex_filtration",
               &Gudhi::Simplex_tree_interface::assign_simplex_filtration,
               nb::arg("simplex"), nb::arg("filtration"),
               R"pbdoc(TODO)pbdoc")
          // .def("initialize_filtration",
          //      &Gudhi::Simplex_tree_interface::initialize_filtration,
          //      nb::arg("ignore_infinite_values")=false,
          //      R"pbdoc(TODO)pbdoc")
          .def("num_vertices",
               &Gudhi::Simplex_tree_interface::num_vertices,
               R"pbdoc(TODO)pbdoc")
          // .def("num_simplices",
          //      &Gudhi::Simplex_tree_interface::num_simplices,
          //      R"pbdoc(TODO)pbdoc")
          .def("is_empty",
               &Gudhi::Simplex_tree_interface::is_empty,
               R"pbdoc(TODO)pbdoc")
          .def("set_dimension",
               &Gudhi::Simplex_tree_interface::set_dimension,
               nb::arg("dimension"), nb::arg("exact")=true,
               R"pbdoc(TODO)pbdoc")
          // .def("dimension",
          //      &Gudhi::Simplex_tree_interface::dimension,
          //      R"pbdoc(TODO)pbdoc")
          .def("upper_bound_dimension",
               &Gudhi::Simplex_tree_interface::upper_bound_dimension,
               R"pbdoc(TODO)pbdoc")
          .def("find_simplex",
               &Gudhi::Simplex_tree_interface::find_simplex,
               nb::arg("simplex"),
               R"pbdoc(TODO)pbdoc")
          .def("insert",
               &Gudhi::Simplex_tree_interface::insert,
               nb::arg("simplex"), nb::arg("filtration"),
               R"pbdoc(TODO)pbdoc")
          .def("insert_matrix",
               &Gudhi::Simplex_tree_interface::insert_matrix,
               nb::arg("filtrations"), nb::arg("n"), nb::arg("stride0"), nb::arg("stride1"), nb::arg("max_filtration"),
               R"pbdoc(TODO)pbdoc")
          // .def("insert_batch_vertices",
          //      &Gudhi::Simplex_tree_interface::insert_batch_vertices,
          //      //nb::arg("v"), nb::arg("f"), nb::arg("filt")=Filtration_value(),
          //      nb::arg("v"), nb::arg("f"), nb::arg("filt"),
          //      R"pbdoc(TODO)pbdoc")
          .def("get_star",
               &Gudhi::Simplex_tree_interface::get_star,
               nb::arg("simplex"),
               R"pbdoc(TODO)pbdoc")
          .def("get_cofaces",
               &Gudhi::Simplex_tree_interface::get_cofaces,
               nb::arg("simplex"), nb::arg("dimension"),
               R"pbdoc(TODO)pbdoc")
          .def("expansion",
               &Gudhi::Simplex_tree_interface::expansion,
               nb::arg("max_dim"),
               R"pbdoc(TODO)pbdoc")
          .def("remove_maximal_simplex",
               &Gudhi::Simplex_tree_interface::remove_maximal_simplex,
               nb::arg("simplex"),
               R"pbdoc(TODO)pbdoc")
          .def("prune_above_filtration",
               &Gudhi::Simplex_tree_interface::prune_above_filtration,
               nb::arg("filtration"),
               R"pbdoc(TODO)pbdoc")
          .def("prune_above_dimension",
               &Gudhi::Simplex_tree_interface::prune_above_dimension,
               nb::arg("dimension"),
               R"pbdoc(TODO)pbdoc")
          .def("make_filtration_non_decreasing",
               &Gudhi::Simplex_tree_interface::make_filtration_non_decreasing,
               R"pbdoc(TODO)pbdoc")
          .def("compute_extended_filtration",
               &Gudhi::Simplex_tree_interface::compute_extended_filtration,
               R"pbdoc(TODO)pbdoc")
          .def("collapse_edges",
               &Gudhi::Simplex_tree_interface::collapse_edges,
               nb::arg("nb_collapse_iteration"),
               R"pbdoc(TODO)pbdoc")
          .def("reset_filtration",
               &Gudhi::Simplex_tree_interface::reset_filtration,
               nb::arg("filtration"), nb::arg("dimension")=0,
               R"pbdoc(TODO)pbdoc")
          .def("get_simplex_and_filtration",
               &Gudhi::Simplex_tree_interface::get_simplex_and_filtration,
               nb::arg("f_simplex"),
               R"pbdoc(TODO)pbdoc")
          .def("get_simplices_iterator_begin",
               &Gudhi::Simplex_tree_interface::get_simplices_iterator_begin,
               R"pbdoc(TODO)pbdoc")
          .def("get_simplices_iterator_end",
               &Gudhi::Simplex_tree_interface::get_simplices_iterator_end,
               R"pbdoc(TODO)pbdoc")
          .def("get_filtration_iterator_begin",
               &Gudhi::Simplex_tree_interface::get_filtration_iterator_begin,
               R"pbdoc(TODO)pbdoc")
          .def("get_filtration_iterator_end",
               &Gudhi::Simplex_tree_interface::get_filtration_iterator_end,
               R"pbdoc(TODO)pbdoc")
          .def("get_skeleton_iterator_begin",
               &Gudhi::Simplex_tree_interface::get_skeleton_iterator_begin,
               nb::arg("dimension"),
               R"pbdoc(TODO)pbdoc")
          .def("get_skeleton_iterator_end",
               &Gudhi::Simplex_tree_interface::get_skeleton_iterator_end,
               nb::arg("dimension"),
               R"pbdoc(TODO)pbdoc")
          .def("get_boundary_iterators",
               &Gudhi::Simplex_tree_interface::get_boundary_iterators,
               nb::arg("simplex"),
               R"pbdoc(TODO)pbdoc")
     ;
}

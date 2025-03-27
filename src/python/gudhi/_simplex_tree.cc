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

NB_MODULE(_simplex_tree_ext, m) {
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
                    R"pbdoc("""This function returns the number of vertices of the simplicial
        complex.

        :returns:  The simplicial complex number of vertices.
        :rtype:  int
        """)pbdoc")
            .def("num_simplices",
                    nb::overload_cast<>(&gsti::num_simplices, nb::const_),
                    R"pbdoc("""This function returns the number of simplices of the simplicial
        complex.

        :returns:  the simplicial complex number of simplices.
        :rtype:  int
        """)pbdoc")
            .def("is_empty",
                    &gsti::is_empty,
                    R"pbdoc("""This function returns whether the simplicial complex is empty.

        :returns:  True if the simplicial complex is empty.
        :rtype:  bool
        """)pbdoc")
            .def("set_dimension",
                    &gsti::set_dimension,
                    nb::arg("dimension"), nb::arg("exact")=true,
                    R"pbdoc("""This function sets the dimension of the simplicial complex.

        :param dimension: The new dimension value.
        :type dimension: int

        .. note::

            This function must be used with caution because it disables
            dimension recomputation when required
            (this recomputation can be triggered by
            :func:`remove_maximal_simplex`
            or
            :func:`prune_above_filtration`
            ).
        """)pbdoc")
            .def("dimension",
                    nb::overload_cast<>(&gsti::dimension, nb::const_),
                    R"pbdoc("""This function returns the dimension of the simplicial complex.

        :returns:  the simplicial complex dimension.
        :rtype:  int

        .. note::

            This function is not constant time because it can recompute
            dimension if required (can be triggered by
            :func:`remove_maximal_simplex`
            or
            :func:`prune_above_filtration`
            methods).
        """)pbdoc")
            .def("upper_bound_dimension",
                    &gsti::upper_bound_dimension,
                    R"pbdoc("""This function returns a valid dimension upper bound of the
        simplicial complex.

        :returns:  an upper bound on the dimension of the simplicial complex.
        :rtype:  int
        """)pbdoc")
            .def("find_simplex",
                    &gsti::find_simplex,
                    nb::arg("simplex"),
                    R"pbdoc(TODO)pbdoc")
            .def("insert",
                    &gsti::insert,
                    nb::arg("simplex"), nb::arg("filtration"),
                    R"pbdoc("""This function inserts the given N-simplex and its subfaces with the
        given filtration value (default value is '0.0'). If some of those
        simplices are already present with a higher filtration value, their
        filtration value is lowered.

        :param simplex: The N-simplex to insert, represented by a list of
            vertex.
        :type simplex: list of int
        :param filtration: The filtration value of the simplex.
        :type filtration: float
        :returns:  true if the simplex was not yet in the complex, false
            otherwise (whatever its original filtration value).
        :rtype:  bool
        """)pbdoc")
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
    
    nb::class_<gsti>(m, (m, "_Simplex_tree_persistence_interface")
            .def(nb::init<CC&, bool>())
            .def("compute_persistence",
                &gsti::compute_persistence,
                nb::arg("homology_coeff_field"),
                nb::arg("double min_persistence"))
            .def("get_persistence", &gsti::get_persistence)
            .def("cofaces_of_cubical_persistence_pairs", &gsti::cofaces_of_cubical_persistence_pairs)
            .def("vertices_of_cubical_persistence_pairs", &gsti::vertices_of_cubical_persistence_pairs)
            .def("betti_numbers", &gsti::betti_numbers)
            .def("persistent_betti_numbers", &gsti::persistent_betti_numbers, nb::arg("from_value"), nb::arg("to_value"))
            .def("intervals_in_dimension", &gsti::intervals_in_dimension, nb::arg("dimension"));
}

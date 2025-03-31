/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - 2023/02 Vincent Rouvreau: Add serialize/deserialize for pickle feature
 *      - 2025/03 Alexis Gobé & Hannah Schreiber: Use nanobind instead of Cython for python bindings.
 *      - YYYY/MM Author: Description of the modification
 */

#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/operators.h>
#include <nanobind/ndarray.h>

#include <Persistent_cohomology_interface.h>
#include <Simplex_tree_interface.h>

namespace nb = nanobind;

using gsti = Gudhi::Simplex_tree_interface;
using gpers = Gudhi::Persistent_cohomology_interface<gsti>;

gsti deserialize_from_python(const nb::ndarray<char, nb::ndim<1>, nb::numpy> &state)
{
  gsti st;
  st.deserialize(state.data(), state.shape(0));
  return st;
}

NB_MODULE(_simplex_tree_ext, m)
{
  m.attr("__license__") = "MIT";

  nb::class_<gsti::Simplex_handle>(m, "_Simplex_handle").def(nb::init<>());

  nb::class_<gsti>(m, "_Simplex_tree_python_interface")
      .def(nb::init<>())
      .def(nb::init<gsti &>())
      .def("filtration", &gsti::simplex_filtration, nb::arg("simplex"), R"pbdoc(
"""This function returns the filtration value for a given N-simplex in
this simplicial complex, or +infinity if it is not in the complex.

:param simplex: The N-simplex, represented by a list of vertex.
:type simplex: list of int
:returns:  The simplicial complex filtration value.
:rtype:  float
"""
           )pbdoc")
      .def("assign_filtration",
           &gsti::assign_simplex_filtration,
           nb::arg("simplex"),
           nb::arg("filtration"),
           R"pbdoc(
"""This function assigns a new filtration value to a
given N-simplex.

:param simplex: The N-simplex, represented by a list of vertex.
:type simplex: list of int
:param filtration:  The new filtration value.
:type filtration:  float

.. note::
        Beware that after this operation, the structure may not be a valid
        filtration anymore, a simplex could have a lower filtration value
        than one of its faces. Callers are responsible for fixing this
        (with more :meth:`assign_filtration` or
        :meth:`make_filtration_non_decreasing` for instance) before calling
        any function that relies on the filtration property, like
        :meth:`persistence`.
"""
           )pbdoc")
      .def("initialize_filtration", &gsti::initialize_filtration)
      .def("num_vertices",
           &gsti::num_vertices,
           R"pbdoc(
"""This function returns the number of vertices of the simplicial
complex.

:returns:  The simplicial complex number of vertices.
:rtype:  int
"""
           )pbdoc")
      .def("num_simplices",
           nb::overload_cast<>(&gsti::num_simplices, nb::const_),
           R"pbdoc(
"""This function returns the number of simplices of the simplicial
complex.

:returns:  the simplicial complex number of simplices.
:rtype:  int
"""
           )pbdoc")
      .def("is_empty",
           &gsti::is_empty,
           R"pbdoc("""This function returns whether the simplicial complex is empty.

        :returns:  True if the simplicial complex is empty.
        :rtype:  bool
        """)pbdoc")
      .def("set_dimension",
           &gsti::set_dimension,
           nb::arg("dimension"),
           nb::arg("exact") = true,
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
      .def("find", &gsti::find_simplex, nb::arg("simplex"), R"pbdoc(
        """This function returns if the N-simplex was found in the simplicial
        complex or not.

        :param simplex: The N-simplex to find, represented by a list of vertex.
        :type simplex: list of int
        :returns:  true if the simplex was found, false otherwise.
        :rtype:  bool
        """
        )pbdoc")
      .def("insert",
           &gsti::insert,
           nb::arg("simplex"),
           nb::arg("filtration") = 0.0,
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
      .def("insert_matrix", &gsti::insert_matrix, nb::arg("filtrations"), nb::arg("max_filtration"))
      .def("insert_batch_vertices", &gsti::insert_batch_vertices<std::vector<int>>)
      .def("get_star", &gsti::get_star, nb::arg("simplex"), R"pbdoc(
        """This function returns the star of a given N-simplex.

        :param simplex: The N-simplex, represented by a list of vertex.
        :type simplex: list of int
        :returns:  The (simplices of the) star of a simplex.
        :rtype:  list of tuples(simplex, filtration)
        """
        )pbdoc")
      .def("get_cofaces", &gsti::get_cofaces, nb::arg("simplex"), nb::arg("dimension"), R"pbdoc(
        """This function returns the cofaces of a given N-simplex with a
        given codimension.

        :param simplex: The N-simplex, represented by a list of vertex.
        :type simplex: list of int
        :param codimension: The codimension. If codimension = 0, all cofaces
            are returned (equivalent of get_star function)
        :type codimension: int
        :returns:  The (simplices of the) cofaces of a simplex
        :rtype:  list of tuples(simplex, filtration)
        """
        )pbdoc")
      .def("expansion",
           &gsti::expansion,
           nb::arg("max_dimension"),
           R"pbdoc("""Expands the simplex tree containing only its one skeleton
        until dimension max_dim.

        The expanded simplicial complex until dimension :math:`d`
        attached to a graph :math:`G` is the maximal simplicial complex of
        dimension at most :math:`d` admitting the graph :math:`G` as
        :math:`1`-skeleton.
        The filtration value assigned to a simplex is the maximal filtration
        value of one of its edges.

        The simplex tree must contain no simplex of dimension bigger than
        1 when calling the method.

        :param max_dimension: The maximal dimension.
        :type max_dimension: int
        """)pbdoc")
      .def("remove_maximal_simplex",
           &gsti::remove_maximal_simplex,
           nb::arg("simplex"),
           R"pbdoc("""This function removes a given maximal N-simplex from the simplicial
        complex.

        :param simplex: The N-simplex, represented by a list of vertex.
        :type simplex: list of int

        .. note::

            The dimension of the simplicial complex may be lower after calling
            remove_maximal_simplex than it was before. However,
            :func:`upper_bound_dimension`
            method will return the old value, which
            remains a valid upper bound. If you care, you can call
            :func:`dimension`
            to recompute the exact dimension.
        """)pbdoc")
      .def("prune_above_filtration",
           &gsti::prune_above_filtration,
           nb::arg("filtration"),
           R"pbdoc(        """Prune above filtration value given as parameter.

        :param filtration: Maximum threshold value.
        :type filtration: float
        :returns: The filtration modification information.
        :rtype: bool


        .. note::

            Note that the dimension of the simplicial complex may be lower
            after calling
            :func:`prune_above_filtration`
            than it was before. However,
            :func:`upper_bound_dimension`
            will return the old value, which remains a
            valid upper bound. If you care, you can call
            :func:`dimension`
            method to recompute the exact dimension.
        """)pbdoc")
      .def("prune_above_dimension",
           &gsti::prune_above_dimension,
           nb::arg("dimension"),
           R"pbdoc("""Remove all simplices of dimension greater than a given value.

        :param dimension: Maximum dimension value.
        :type dimension: int
        :returns: The modification information.
        :rtype: bool
        """)pbdoc")
      .def("make_filtration_non_decreasing",
           &gsti::make_filtration_non_decreasing,
           R"pbdoc(        """This function ensures that each simplex has a higher filtration
        value than its faces by increasing the filtration values.

        :returns: True if any filtration value was modified,
            False if the filtration was already non-decreasing.
        :rtype: bool
        """)pbdoc")
      .def(
          "extend_filtration",
          &gsti::compute_extended_filtration,
          R"pbdoc(""" Extend filtration for computing extended persistence. This function only uses the filtration values at the
        0-dimensional simplices, and computes the extended persistence diagram induced by the lower-star filtration
        computed with these values.

        .. note::

            Note that after calling this function, the filtration values are actually modified within the simplex tree.
            The function :func:`extended_persistence` retrieves the original values.

        .. note::

            Note that this code creates an extra vertex internally, so you should make sure that the simplex tree does
            not contain a vertex with the largest possible value (i.e., 4294967295).

        This `notebook <https://github.com/GUDHI/TDA-tutorial/blob/master/Tuto-GUDHI-extended-persistence.ipynb>`_
        explains how to compute an extension of persistence called extended persistence.
        """)pbdoc")
      .def("collapse_edges", &gsti::collapse_edges, nb::arg("nb_collapse_iteration"))
      .def(
          "reset_filtration",
          &gsti::reset_filtration,
          nb::arg("filtration"),
          nb::arg("min_dim") = 0,
          R"pbdoc("""This function resets the filtration value of all the simplices of dimension at least min_dim. Resets all the
        simplex tree when `min_dim = 0`.
        `reset_filtration` may break the filtration property with `min_dim > 0`, and it is the user's responsibility to
        make it a valid filtration (using a large enough `filt_value`, or calling `make_filtration_non_decreasing`
        afterwards for instance).

        :param filtration: New threshold value.
        :type filtration: float.
        :param min_dim: The minimal dimension. Default value is 0.
        :type min_dim: int.
        """)pbdoc")
      .def(nb::self == nb::self, R"pbdoc(
        """:returns: True if the 2 complexes have the same simplices with the same filtration values, False otherwise.
        :rtype: bool
        """
        )pbdoc")
      .def("get_simplex_and_filtration", &gsti::get_simplex_and_filtration, nb::arg("f_simplex"))
      .def("simplex_iter", &gsti::get_simplex_python_iterator)
      .def("filtration_iter", &gsti::get_filtration_python_iterator)
      .def("skeleton_iter", &gsti::get_skeleton_python_iterator)
      .def("boundary_iter", &gsti::get_boundary_python_iterator)
      .def("get_boundary_iterators", &gsti::get_boundary_iterators, nb::arg("simplex"))
      .def("expansion_with_blocker",
           &gsti::expansion_with_blockers_callback,
           nb::arg("max_dim"),
           nb::arg("blocker_func"),
           R"pbdoc(
        """Expands the Simplex_tree containing only a graph. Simplices corresponding to cliques in the graph are added
        incrementally, faces before cofaces, unless the simplex has dimension larger than `max_dim` or `blocker_func`
        returns `True` for this simplex.

        The function identifies a candidate simplex whose faces are all already in the complex, inserts it with a
        filtration value corresponding to the maximum of the filtration values of the faces, then calls `blocker_func`
        with this new simplex (represented as a list of int). If `blocker_func` returns `True`, the simplex is removed,
        otherwise it is kept. The algorithm then proceeds with the next candidate.

        .. warning::
            Several candidates of the same dimension may be inserted simultaneously before calling `blocker_func`, so
            if you examine the complex in `blocker_func`, you may hit a few simplices of the same dimension that have
            not been vetted by `blocker_func` yet, or have already been rejected but not yet removed.

        :param max_dim: Expansion maximal dimension value.
        :type max_dim: int
        :param blocker_func: Blocker oracle.
        :type blocker_func: Callable[[List[int]], bool]
        """
        )pbdoc")
      .def("clear", &gsti::clear)
      .def("__getstate__",
           [](const gsti &st) -> nb::ndarray<char, nb::ndim<1>, nb::numpy> {
             auto buffer_size = st.get_serialization_size();
             char *buffer = new char[buffer_size];
             st.serialize(buffer, buffer_size);
             nb::ndarray<char, nb::ndim<1>, nb::numpy> np_buffer(
                 buffer, {buffer_size}, nb::capsule(buffer, [](void *p) noexcept {
                   delete reinterpret_cast<char *>(p);
                 }));
             return np_buffer;
           })
      .def("__setstate__", [](gsti &st, const nb::ndarray<char, nb::ndim<1>, nb::numpy> &state) -> void {
        new (&st) gsti(deserialize_from_python(state));
      });

  nb::class_<gpers>(m, "_Simplex_tree_persistence_interface")
      .def(nb::init<gsti &, bool>())
      .def("compute_persistence",
           &gpers::compute_persistence,
           nb::arg("homology_coeff_field"),
           nb::arg("double min_persistence"))
      .def("get_persistence", &gpers::get_persistence)
      .def("betti_numbers", &gpers::betti_numbers)
      .def("persistent_betti_numbers", &gpers::persistent_betti_numbers, nb::arg("from_value"), nb::arg("to_value"))
      .def("intervals_in_dimension", &gpers::intervals_in_dimension, nb::arg("dimension"))
      .def("write_output_diagram", &gpers::write_output_diagram, nb::arg("diagram_file_name"))
      .def("persistence_pairs", &gpers::persistence_pairs)
      .def("lower_star_generators", &gpers::lower_star_generators)
      .def("flag_generators", &gpers::flag_generators)
      .def("compute_extended_persistence_subdiagrams",
           &gpers::compute_extended_persistence_subdiagrams,
           nb::arg("min_persistence"));
}
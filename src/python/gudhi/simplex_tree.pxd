# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2016 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

from cython cimport numeric
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp cimport bool
from libcpp.string cimport string

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "MIT"

cdef extern from "Simplex_tree_interface.h" namespace "Gudhi":
    cdef cppclass Simplex_tree_options_full_featured:
        pass

    cdef cppclass Simplex_tree_simplex_handle "Gudhi::Simplex_tree_interface<Gudhi::Simplex_tree_options_full_featured>::Simplex_handle":
        pass

    cdef cppclass Simplex_tree_simplices_iterator "Gudhi::Simplex_tree_interface<Gudhi::Simplex_tree_options_full_featured>::Complex_simplex_iterator":
        Simplex_tree_simplices_iterator()
        Simplex_tree_simplex_handle& operator*()
        Simplex_tree_simplices_iterator operator++()
        bint operator!=(Simplex_tree_simplices_iterator)

    cdef cppclass Simplex_tree_skeleton_iterator "Gudhi::Simplex_tree_interface<Gudhi::Simplex_tree_options_full_featured>::Skeleton_simplex_iterator":
        Simplex_tree_skeleton_iterator()
        Simplex_tree_simplex_handle& operator*()
        Simplex_tree_skeleton_iterator operator++()
        bint operator!=(Simplex_tree_skeleton_iterator)


    cdef cppclass Simplex_tree_interface_full_featured "Gudhi::Simplex_tree_interface<Gudhi::Simplex_tree_options_full_featured>":
        Simplex_tree()
        double simplex_filtration(vector[int] simplex)
        void assign_simplex_filtration(vector[int] simplex, double filtration)
        void initialize_filtration()
        int num_vertices()
        int num_simplices()
        void set_dimension(int dimension)
        int dimension()
        int upper_bound_dimension()
        bool find_simplex(vector[int] simplex)
        bool insert_simplex_and_subfaces(vector[int] simplex,
                                         double filtration)
        vector[pair[vector[int], double]] get_star(vector[int] simplex)
        vector[pair[vector[int], double]] get_cofaces(vector[int] simplex,
                                                          int dimension)
        void expansion(int max_dim) except +
        void remove_maximal_simplex(vector[int] simplex)
        bool prune_above_filtration(double filtration)
        bool make_filtration_non_decreasing()
        # Iterators over Simplex tree
        pair[vector[int], double] get_simplex_and_filtration(Simplex_tree_simplex_handle f_simplex)
        Simplex_tree_simplices_iterator get_simplices_iterator_begin()
        Simplex_tree_simplices_iterator get_simplices_iterator_end()
        vector[Simplex_tree_simplex_handle].const_iterator get_filtration_iterator_begin()
        vector[Simplex_tree_simplex_handle].const_iterator get_filtration_iterator_end()
        Simplex_tree_skeleton_iterator get_skeleton_iterator_begin(int dimension)
        Simplex_tree_skeleton_iterator get_skeleton_iterator_end(int dimension)

cdef extern from "Persistent_cohomology_interface.h" namespace "Gudhi":
    cdef cppclass Simplex_tree_persistence_interface "Gudhi::Persistent_cohomology_interface<Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_full_featured>>":
        Simplex_tree_persistence_interface(Simplex_tree_interface_full_featured * st, bool persistence_dim_max)
        vector[pair[int, pair[double, double]]] get_persistence(int homology_coeff_field, double min_persistence)
        vector[int] betti_numbers()
        vector[int] persistent_betti_numbers(double from_value, double to_value)
        vector[pair[double,double]] intervals_in_dimension(int dimension)
        void write_output_diagram(string diagram_file_name)
        vector[pair[vector[int], vector[int]]] persistence_pairs()
        pair[vector[vector[int]], vector[vector[int]]] lower_star_generators(double)
        pair[vector[vector[int]], vector[vector[int]]] flag_generators(double)

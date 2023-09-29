# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2016 Inria
#
# Modification(s):
#   - 2023/02 Vincent Rouvreau: Add serialize/deserialize for pickle feature
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
    cdef cppclass Simplex_tree_simplex_handle "Gudhi::Simplex_tree_interface::Simplex_handle":
        pass

    cdef cppclass Simplex_tree_simplices_iterator "Gudhi::Simplex_tree_interface::Complex_simplex_iterator":
        Simplex_tree_simplices_iterator() nogil
        Simplex_tree_simplex_handle& operator*() nogil
        Simplex_tree_simplices_iterator operator++() nogil
        bint operator!=(Simplex_tree_simplices_iterator) nogil

    cdef cppclass Simplex_tree_skeleton_iterator "Gudhi::Simplex_tree_interface::Skeleton_simplex_iterator":
        Simplex_tree_skeleton_iterator() nogil
        Simplex_tree_simplex_handle& operator*() nogil
        Simplex_tree_skeleton_iterator operator++() nogil
        bint operator!=(Simplex_tree_skeleton_iterator) nogil

    cdef cppclass Simplex_tree_boundary_iterator "Gudhi::Simplex_tree_interface::Boundary_simplex_iterator":
        Simplex_tree_boundary_iterator() nogil
        Simplex_tree_simplex_handle& operator*() nogil
        Simplex_tree_boundary_iterator operator++() nogil
        bint operator!=(Simplex_tree_boundary_iterator) nogil


    cdef cppclass Simplex_tree_python_interface "Gudhi::Simplex_tree_interface":
        Simplex_tree_python_interface() nogil
        Simplex_tree_python_interface(Simplex_tree_python_interface&) nogil
        double simplex_filtration(vector[int] simplex) nogil
        void assign_simplex_filtration(vector[int] simplex, double filtration) nogil except +
        void initialize_filtration() nogil
        int num_vertices() nogil
        int num_simplices() nogil
        bool is_empty() nogil
        void set_dimension(int dimension) nogil
        int dimension() nogil
        int upper_bound_dimension() nogil
        bool find_simplex(vector[int] simplex) nogil
        bool insert(vector[int] simplex, double filtration) nogil
        void insert_matrix(double* filtrations, int n, int stride0, int stride1, double max_filtration) nogil except +
        void insert_batch_vertices(vector[int] v, double f) nogil except +
        vector[pair[vector[int], double]] get_star(vector[int] simplex) nogil
        vector[pair[vector[int], double]] get_cofaces(vector[int] simplex, int dimension) nogil
        void expansion(int max_dim) nogil except +
        void remove_maximal_simplex(vector[int] simplex) nogil
        bool prune_above_filtration(double filtration) nogil
        bool prune_above_dimension(int dimension) nogil
        bool make_filtration_non_decreasing() nogil
        void compute_extended_filtration() nogil
        Simplex_tree_python_interface* collapse_edges(int nb_collapse_iteration) nogil except +
        void reset_filtration(double filtration, int dimension) nogil
        bint operator==(Simplex_tree_python_interface) nogil
        # Iterators over Simplex tree
        pair[vector[int], double] get_simplex_and_filtration(Simplex_tree_simplex_handle f_simplex) nogil
        Simplex_tree_simplices_iterator get_simplices_iterator_begin() nogil
        Simplex_tree_simplices_iterator get_simplices_iterator_end() nogil
        vector[Simplex_tree_simplex_handle].const_iterator get_filtration_iterator_begin() nogil
        vector[Simplex_tree_simplex_handle].const_iterator get_filtration_iterator_end() nogil
        Simplex_tree_skeleton_iterator get_skeleton_iterator_begin(int dimension) nogil
        Simplex_tree_skeleton_iterator get_skeleton_iterator_end(int dimension) nogil
        pair[Simplex_tree_boundary_iterator, Simplex_tree_boundary_iterator] get_boundary_iterators(vector[int] simplex) nogil except +
        # Expansion with blockers
        ctypedef bool (*blocker_func_t)(vector[int], void *user_data) except +
        void expansion_with_blockers_callback(int dimension, blocker_func_t user_func, void *user_data) except +
        void serialize(char* buffer, const size_t buffer_size) nogil except +
        void deserialize(const char* buffer, const size_t buffer_size) nogil except +
        size_t get_serialization_size() nogil

cdef extern from "Persistent_cohomology_interface.h" namespace "Gudhi":
    cdef cppclass Simplex_tree_persistence_interface "Gudhi::Persistent_cohomology_interface<Gudhi::Simplex_tree_interface>":
        Simplex_tree_persistence_interface(Simplex_tree_python_interface * st, bool persistence_dim_max) nogil
        void compute_persistence(int homology_coeff_field, double min_persistence) nogil except +
        vector[pair[int, pair[double, double]]] get_persistence() nogil
        vector[int] betti_numbers() nogil
        vector[int] persistent_betti_numbers(double from_value, double to_value) nogil
        vector[pair[double,double]] intervals_in_dimension(int dimension) nogil
        void write_output_diagram(string diagram_file_name) nogil except +
        vector[pair[vector[int], vector[int]]] persistence_pairs() nogil
        pair[vector[vector[int]], vector[vector[int]]] lower_star_generators() nogil
        pair[vector[vector[int]], vector[vector[int]]] flag_generators() nogil
        vector[vector[pair[int, pair[double, double]]]] compute_extended_persistence_subdiagrams(double min_persistence) nogil

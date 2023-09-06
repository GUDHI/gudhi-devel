# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2016 Inria
#
# Modification(s):
#	- 2022 David Loiseaux, Hannah Schreiber: adapt for multipersistence. 
#   - YYYY/MM Author: Description of the modification

from cython cimport numeric
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp cimport bool
from libcpp.string cimport string

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "MIT"

ctypedef int dimension_type
ctypedef float value_type
ctypedef Finitely_critical_multi_filtration filtration_type
ctypedef vector[value_type] python_filtration_type ## TODO: move constructor for C++ filtration type ?
ctypedef vector[int] simplex_type
ctypedef vector[simplex_type] simplex_list
ctypedef vector[pair[pair[int,int], pair[double, double]]] edge_list
ctypedef vector[int] euler_char_list
ctypedef pair[simplex_type, value_type*] simplex_filtration_type




cdef extern from "multi_filtrations/finitely_critical_filtrations.h" namespace "Gudhi::multi_filtrations":
	cdef cppclass Finitely_critical_multi_filtration "Gudhi::multi_filtrations::Finitely_critical_multi_filtration<Gudhi::Simplex_tree_options_multidimensional_filtration::value_type>":
		Finitely_critical_multi_filtration()  except + nogil
		Finitely_critical_multi_filtration(vector[value_type]) except +
		Finitely_critical_multi_filtration& operator=(const Finitely_critical_multi_filtration&)
		filtration_type& get_vector()  nogil
		int size() nogil
		void push_back(const value_type&) nogil
		void clear() nogil
		value_type& operator[](int)
		


cdef extern from "Simplex_tree_interface_multi.h" namespace "Gudhi":
	cdef cppclass Simplex_tree_options_multidimensional_filtration:
		pass

	cdef cppclass Simplex_tree_multi_simplex_handle "Gudhi::Simplex_tree_interface_multi<Gudhi::Simplex_tree_options_multidimensional_filtration>::Simplex_handle":
		pass

	cdef cppclass Simplex_tree_multi_simplices_iterator "Gudhi::Simplex_tree_interface_multi<Gudhi::Simplex_tree_options_multidimensional_filtration>::Complex_simplex_iterator":
		Simplex_tree_multi_simplices_iterator() nogil
		Simplex_tree_multi_simplex_handle& operator*() nogil
		Simplex_tree_multi_simplices_iterator operator++() nogil
		bint operator!=(Simplex_tree_multi_simplices_iterator) nogil

	cdef cppclass Simplex_tree_multi_skeleton_iterator "Gudhi::Simplex_tree_interface_multi<Gudhi::Simplex_tree_options_multidimensional_filtration>::Skeleton_simplex_iterator":
		Simplex_tree_multi_skeleton_iterator() nogil
		Simplex_tree_multi_simplex_handle& operator*() nogil
		Simplex_tree_multi_skeleton_iterator operator++() nogil
		bint operator!=(Simplex_tree_multi_skeleton_iterator) nogil

	cdef cppclass Simplex_tree_multi_boundary_iterator "Gudhi::Simplex_tree_interface_multi<Gudhi::Simplex_tree_options_multidimensional_filtration>::Boundary_simplex_iterator":
		Simplex_tree_multi_boundary_iterator() nogil
		Simplex_tree_multi_simplex_handle& operator*() nogil
		Simplex_tree_multi_boundary_iterator operator++() nogil
		bint operator!=(Simplex_tree_multi_boundary_iterator) nogil


	cdef cppclass Simplex_tree_multi_interface "Gudhi::Simplex_tree_interface_multi<Gudhi::Simplex_tree_options_multidimensional_filtration>":
		Simplex_tree_multi_interface() nogil
		Simplex_tree_multi_interface(Simplex_tree_multi_interface&) nogil
		value_type* simplex_filtration(const vector[int]& simplex) nogil
		void assign_simplex_filtration(vector[int]& simplex, const filtration_type& filtration) nogil
		void initialize_filtration() nogil
		int num_vertices() nogil
		int num_simplices() nogil
		void set_dimension(int dimension) nogil
		dimension_type dimension() nogil
		dimension_type upper_bound_dimension() nogil
		bool find_simplex(vector[int]& simplex) nogil
		bool insert(vector[int]& simplex, filtration_type& filtration) nogil
		vector[simplex_filtration_type] get_star(const vector[int]& simplex) nogil
		vector[simplex_filtration_type] get_cofaces(const vector[int]& simplex, int dimension) nogil
		void expansion(int max_dim)  except + nogil
		void remove_maximal_simplex(simplex_type simplex) nogil
		# bool prune_above_filtration(filtration_type filtration) nogil
		bool make_filtration_non_decreasing()  except + nogil
		# void compute_extended_filtration() nogil
		# Simplex_tree_multi_interface* collapse_edges(int nb_collapse_iteration)  except + nogil
		void reset_filtration(const filtration_type& filtration, int dimension) nogil
		bint operator==(Simplex_tree_multi_interface) nogil
		# Iterators over Simplex tree
		simplex_filtration_type get_simplex_and_filtration(Simplex_tree_multi_simplex_handle f_simplex) nogil
		Simplex_tree_multi_simplices_iterator get_simplices_iterator_begin() nogil
		Simplex_tree_multi_simplices_iterator get_simplices_iterator_end() nogil
		vector[Simplex_tree_multi_simplex_handle].const_iterator get_filtration_iterator_begin() nogil
		vector[Simplex_tree_multi_simplex_handle].const_iterator get_filtration_iterator_end() nogil
		Simplex_tree_multi_skeleton_iterator get_skeleton_iterator_begin(int dimension) nogil
		Simplex_tree_multi_skeleton_iterator get_skeleton_iterator_end(int dimension) nogil
		pair[Simplex_tree_multi_boundary_iterator, Simplex_tree_multi_boundary_iterator] get_boundary_iterators(vector[int] simplex)  except + nogil
		# Expansion with blockers
		ctypedef bool (*blocker_func_t)(vector[int], void *user_data)
		void expansion_with_blockers_callback(int dimension, blocker_func_t user_func, void *user_data)

		## MULTIPERS STUFF
		void set_keys_to_enumerate() nogil const
		int get_key(const simplex_type) nogil
		void set_key(simplex_type, int) nogil
		void fill_lowerstar(const vector[value_type]&, int) nogil
		simplex_list get_simplices_of_dimension(int) nogil
		edge_list get_edge_list() nogil
		# euler_char_list euler_char(const vector[filtration_type]&) nogil
		void resize_all_filtrations(int) nogil
		void set_number_of_parameters(int) nogil
		int get_number_of_parameters() nogil



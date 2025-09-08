# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2016 Inria
#
# Modification(s):
#	- 2023 David Loiseaux : Conversions with standard simplextree, scc2020 format, edge collapses, euler characteristic, grid filtrations.
#	- 2022/11 Hannah Schreiber / David Loiseaux : adapt for multipersistence. 
#   - YYYY/MM Author: Description of the modification



__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "MIT"

from libc.stdint cimport intptr_t, int32_t, int64_t
from cython.operator import dereference, preincrement
from libc.stdint cimport intptr_t
from libc.stdint cimport uintptr_t
from libcpp.map cimport map


ctypedef fused some_int:
	int32_t
	int64_t
	int

ctypedef fused some_float:
	float
	double

ctypedef vector[pair[pair[int,int],pair[value_type,value_type]]] edge_list_type

from typing import Any

cimport numpy as cnp
import numpy as np
cnp.import_array()

from gudhi.simplex_tree_multi cimport *
cimport cython
from gudhi.simplex_tree import SimplexTree ## Small hack for typing
from typing import Iterable
from tqdm import tqdm



from warnings import warn


cdef extern from "Simplex_tree_multi.h" namespace "Gudhi::multiparameter":
	void multify_from_ptr(const uintptr_t, const uintptr_t, const unsigned int, const vector[value_type]&)  except + nogil
	void flatten_from_ptr(const uintptr_t, const uintptr_t, const unsigned int) nogil
	void linear_projection_from_ptr(const uintptr_t, const uintptr_t, const vector[value_type]&) nogil
	void flatten_diag_from_ptr(const uintptr_t, const uintptr_t, const vector[value_type], int) nogil
	void squeeze_filtration(uintptr_t, const vector[vector[value_type]]&, bool)  except + nogil
	vector[vector[vector[value_type]]] get_filtration_values(uintptr_t, const vector[int]&)  except + nogil


# cdef bool callback(vector[int] simplex, void *blocker_func):
	# 	return (<object>blocker_func)(simplex)

# SimplexTree python interface
cdef class SimplexTreeMulti:
	"""The simplex tree is an efficient and flexible data structure for
	representing general (filtered) simplicial complexes. The data structure
	is described in Jean-Daniel Boissonnat and Clément Maria. The Simplex
	Tree: An Efficient Data Structure for General Simplicial Complexes.
	Algorithmica, pages 1–22, 2014.

	This class is a multi-filtered, with keys, and non contiguous vertices version
	of the simplex tree. 
	"""
	# unfortunately 'cdef public Simplex_tree_multi_interface* thisptr' is not possible
	# Use intptr_t instead to cast the pointer
	cdef public intptr_t thisptr
	cdef public vector[vector[value_type]] filtration_grid

	# Get the pointer casted as it should be
	cdef Simplex_tree_multi_interface* get_ptr(self) nogil:
		return <Simplex_tree_multi_interface*>(self.thisptr)

	# cdef Simplex_tree_persistence_interface * pcohptr
	# Fake constructor that does nothing but documenting the constructor
	def __init__(self, other = None, num_parameters:int=2,default_values=[], safe_conversion=False):
		"""SimplexTreeMulti constructor.
		
		:param other: If `other` is `None` (default value), an empty `SimplexTreeMulti` is created.
			If `other` is a `SimplexTree`, the `SimplexTreeMulti` is constructed from a deep copy of `other`.
			If `other` is a `SimplexTreeMulti`, the `SimplexTreeMulti` is constructed from a deep copy of `other`.
		:type other: SimplexTree or SimplexTreeMulti (Optional)
		:param num_parameters: The number of parameter of the multi-parameter filtration.
		:type num_parameters: int
		:returns: An empty or a copy simplex tree.
		:rtype: SimplexTreeMulti

		:raises TypeError: In case `other` is neither `None`, nor a `SimplexTree`, nor a `SimplexTreeMulti`.
		"""
		
	# The real cython constructor
	def __cinit__(self, other = None, int num_parameters=2, 
		default_values=[-np.inf], # I'm not sure why `[]` does not work. Cython bug ? 
bool safe_conversion=False,
		): #TODO doc
		cdef vector[value_type] c_default_values=default_values
		if other is not None:
			if isinstance(other, SimplexTreeMulti):
				self.thisptr = _get_copy_intptr(other)
				num_parameters = other.num_parameters
				self.filtration_grid = other.filtration_grid
			elif isinstance(other, SimplexTree): # Constructs a SimplexTreeMulti from a SimplexTree
				self.thisptr = <intptr_t>(new Simplex_tree_multi_interface())
				if safe_conversion:
					new_st_multi  = _safe_simplextree_multify(other, num_parameters = num_parameters)
					self.thisptr, new_st_multi.thisptr = new_st_multi.thisptr, self.thisptr
				else:
					multify_from_ptr(other.thisptr, self.thisptr, num_parameters, c_default_values)
			else:
				raise TypeError("`other` argument requires to be of type `SimplexTree`, `SimplexTreeMulti`, or `None`.")
		else:
			self.thisptr = <intptr_t>(new Simplex_tree_multi_interface())
		self.get_ptr().set_number_of_parameters(num_parameters)
		self.filtration_grid=[[]*num_parameters] 

	def __dealloc__(self):
		cdef Simplex_tree_multi_interface* ptr = self.get_ptr()
		if ptr != NULL:
			del ptr
		# if self.pcohptr != NULL:
		#     del self.pcohptr 

	def __is_defined(self):
		"""Returns true if SimplexTree pointer is not NULL.
			"""
		return self.get_ptr() != NULL

	# def __is_persistence_defined(self):
	#     """Returns true if Persistence pointer is not NULL.
	#      """
	#     return self.pcohptr != NULL

	def copy(self)->SimplexTreeMulti:
		"""
		:returns: A simplex tree that is a deep copy of itself.
		:rtype: SimplexTreeMulti

		:note: The persistence information is not copied. If you need it in the clone, you have to call
			:func:`compute_persistence` on it even if you had already computed it in the original.
		"""
		stree = SimplexTreeMulti(self,num_parameters=self.num_parameters)
		return stree

	def __deepcopy__(self):
		return self.copy()

	def filtration(self, simplex:list|np.ndarray)->np.ndarray:
		"""This function returns the filtration value for a given N-simplex in
		this simplicial complex, or +infinity if it is not in the complex.

		:param simplex: The N-simplex, represented by a list of vertex.
		:type simplex: list of int
		:returns:  The simplicial complex multi-critical filtration value.
		:rtype:  numpy array of shape (-1, num_parameters)
		"""
		return self[simplex]

	def assign_filtration(self, simplex:list|np.ndarray, filtration:list|np.ndarray)->None:
		"""This function assigns a new multi-critical filtration value to a
		given N-simplex.

		:param simplex: The N-simplex, represented by a list of vertex.
		:type simplex: list of int
		:param filtration:  The new filtration(s) value(s), concatenated.
		:type filtration:  list[float] or np.ndarray[float, ndim=1]

		.. note::
			Beware that after this operation, the structure may not be a valid
			filtration anymore, a simplex could have a lower filtration value
			than one of its faces. Callers are responsible for fixing this
			(with more :meth:`assign_filtration` or
			:meth:`make_filtration_non_decreasing` for instance) before calling
			any function that relies on the filtration property, like
			:meth:`persistence`.
		"""
		assert len(filtration)>0 and len(filtration) % self.get_ptr().get_number_of_parameters() == 0
		self.get_ptr().assign_simplex_filtration(simplex, Finitely_critical_multi_filtration(<python_filtration_type>filtration))

	def __getitem__(self, simplex):
		cdef vector[int] csimplex = simplex 
		cdef value_type[:] filtration_view = <value_type[:self.get_ptr().get_number_of_parameters()]> self.get_ptr().simplex_filtration(csimplex)
		return np.asarray(filtration_view)
	

	@property
	def num_vertices(self)->int:
		"""This function returns the number of vertices of the simplicial
		complex.

		:returns:  The simplicial complex number of vertices.
		:rtype:  int
		"""
		return self.get_ptr().num_vertices()
	
	@property
	def num_simplices(self)->int:
		"""This function returns the number of simplices of the simplicial
		complex.

		:returns:  the simplicial complex number of simplices.
		:rtype:  int
		"""
		return self.get_ptr().num_simplices()

	@property
	def dimension(self)->int:
		"""This function returns the dimension of the simplicial complex.

		:returns:  the simplicial complex dimension.
		:rtype:  int

		.. note::

			This function is not constant time because it can recompute
			dimension if required (can be triggered by
			:func:`remove_maximal_simplex`
			or
			:func:`prune_above_filtration`
			methods).
		"""
		return self.get_ptr().dimension()
	def upper_bound_dimension(self)->int:
		"""This function returns a valid dimension upper bound of the
		simplicial complex.

		:returns:  an upper bound on the dimension of the simplicial complex.
		:rtype:  int
		"""
		return self.get_ptr().upper_bound_dimension()

	def set_dimension(self, dimension)->None:
		"""This function sets the dimension of the simplicial complex.

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
		"""
		self.get_ptr().set_dimension(<int>dimension)

	# def find(self, simplex)->bool:
	# 	"""This function returns if the N-simplex was found in the simplicial
	# 	complex or not.

	# 	:param simplex: The N-simplex to find, represented by a list of vertex.
	# 	:type simplex: list of int
	# 	:returns:  true if the simplex was found, false otherwise.
	# 	:rtype:  bool
	# 	"""
	# 	return self.get_ptr().find_simplex(simplex)
	def __contains__(self, simplex)->bool:
		"""This function returns if the N-simplex was found in the simplicial
		complex or not.

		:param simplex: The N-simplex to find, represented by a list of vertex.
		:type simplex: list of int
		:returns:  true if the simplex was found, false otherwise.
		:rtype:  bool
		"""
		return self.get_ptr().find_simplex(simplex)

	def insert(self, simplex, filtration:list|np.ndarray|None=None)->bool:
		"""This function inserts the given N-simplex and its subfaces with the
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
		"""
		# TODO C++, to be compatible with insert_batch and multicritical filtrations
		num_parameters = self.get_ptr().get_number_of_parameters()
		assert filtration is None or len(filtration) % num_parameters == 0
		if filtration is None:	
			filtration = np.array([-np.inf]*num_parameters, dtype = float)
		return self.get_ptr().insert(simplex, Finitely_critical_multi_filtration(<python_filtration_type>filtration))
		
	@cython.boundscheck(False)
	@cython.wraparound(False)
	def insert_batch(self, some_int[:,:] vertex_array, some_float[:,:]  filtrations)->SimplexTreeMulti:
		"""Inserts k-simplices given by a sparse array in a format similar
		to `torch.sparse <https://pytorch.org/docs/stable/sparse.html>`_.
		The n-th simplex has vertices `vertex_array[0,n]`, ...,
		`vertex_array[k,n]` and filtration value `filtrations[n,num_parameters]`.
		/!\ Only compatible with 1-critical filtrations. If a simplex is repeated, 
		only one filtration value will be taken into account.

		:param vertex_array: the k-simplices to insert.
		:type vertex_array: numpy.array of shape (k+1,n)
		:param filtrations: the filtration values.
		:type filtrations: numpy.array of shape (n,num_parameters)
		"""
		# TODO : multi-critical
		# cdef vector[int] vertices = np.unique(vertex_array)
		cdef Py_ssize_t k = vertex_array.shape[0]
		cdef Py_ssize_t n = vertex_array.shape[1]
		assert filtrations.shape[0] == n, 'inconsistent sizes for vertex_array and filtrations'
		assert filtrations.shape[1] == self.num_parameters, "wrong number of parameters"
		cdef Py_ssize_t i
		cdef Py_ssize_t j
		cdef vector[int] v
		cdef Finitely_critical_multi_filtration w
		cdef int n_parameters = self.num_parameters
		with nogil:
			for i in range(n):
				for j in range(k):
					v.push_back(vertex_array[j, i])
				for j in range(n_parameters):
					w.push_back(filtrations[i,j])
				self.get_ptr().insert(v, w)
				v.clear()
				w.clear()
		return self


	@cython.boundscheck(False)
	@cython.wraparound(False)
	def assign_batch_filtration(self, some_int[:,:] vertex_array, some_float[:,:]  filtrations, bool propagate=True)->SimplexTreeMulti:
		"""Assign k-simplices given by a sparse array in a format similar
		to `torch.sparse <https://pytorch.org/docs/stable/sparse.html>`_.
		The n-th simplex has vertices `vertex_array[0,n]`, ...,
		`vertex_array[k,n]` and filtration value `filtrations[n,num_parameters]`.
		/!\ Only compatible with 1-critical filtrations. If a simplex is repeated, 
		only one filtration value will be taken into account.

		:param vertex_array: the k-simplices to assign.
		:type vertex_array: numpy.array of shape (k+1,n)
		:param filtrations: the filtration values.
		:type filtrations: numpy.array of shape (n,num_parameters)
		"""
		cdef Py_ssize_t k = vertex_array.shape[0]
		cdef Py_ssize_t n = vertex_array.shape[1]
		assert filtrations.shape[0] == n, 'inconsistent sizes for vertex_array and filtrations'
		assert filtrations.shape[1] == self.num_parameters, "wrong number of parameters"
		cdef Py_ssize_t i
		cdef Py_ssize_t j
		cdef vector[int] v
		cdef Finitely_critical_multi_filtration w
		cdef int n_parameters = self.num_parameters
		with nogil:
			for i in range(n):
				for j in range(k):
					v.push_back(vertex_array[j, i])
				for j in range(n_parameters):
					w.push_back(filtrations[i,j])
				self.get_ptr().assign_simplex_filtration(v, w)
				v.clear()
				w.clear()
		if propagate: self.make_filtration_non_decreasing()
		return self



	def get_simplices(self):
		"""This function returns a generator with simplices and their given
		filtration values.

		:returns:  The simplices.
		:rtype:  generator with tuples(simplex, filtration)
		"""
		cdef Simplex_tree_multi_simplices_iterator it = self.get_ptr().get_simplices_iterator_begin()
		cdef Simplex_tree_multi_simplices_iterator end = self.get_ptr().get_simplices_iterator_end()
		cdef Simplex_tree_multi_simplex_handle sh = dereference(it)
		# cdef pair[simplex_type,Finitely_critical_multi_filtration] out_
		# while it != end:
		# 	out_ = self.get_ptr().get_simplex_and_filtration(dereference(it))
		# 	out = (out_.first,out_.second.get_vector())
		# 	yield out
		# 	preincrement(it)
		# cdef pair[simplex_type,filtration_type] out
		cdef int num_parameters = self.get_ptr().get_number_of_parameters()
		while it != end:
			pair = self.get_ptr().get_simplex_and_filtration(dereference(it))

			yield (np.asarray(pair.first, dtype=int),np.asarray(<value_type[:num_parameters]> pair.second))
			# yield SimplexTreeMulti._pair_simplex_filtration_to_python(out)
			preincrement(it)

	def get_filtration(self):
		"""This function returns a generator with simplices and their given
		filtration values sorted by increasing filtration values.

		:returns:  The simplices sorted by increasing filtration values.
		:rtype:  generator with tuples(simplex, filtration)
		"""
		cdef vector[Simplex_tree_multi_simplex_handle].const_iterator it = self.get_ptr().get_filtration_iterator_begin()
		cdef vector[Simplex_tree_multi_simplex_handle].const_iterator end = self.get_ptr().get_filtration_iterator_end()
		cdef int num_parameters = self.get_ptr().get_number_of_parameters()
		while it != end:
			# yield self.get_ptr().get_simplex_and_filtration(dereference(it))
			pair = self.get_ptr().get_simplex_and_filtration(dereference(it))
			yield (np.asarray(pair.first, dtype=int),np.asarray(<value_type[:num_parameters]> pair.second))
			preincrement(it)

	def get_skeleton(self, dimension):
		"""This function returns a generator with the (simplices of the) skeleton of a maximum given dimension.

		:param dimension: The skeleton dimension value.
		:type dimension: int
		:returns:  The (simplices of the) skeleton of a maximum dimension.
		:rtype:  generator with tuples(simplex, filtration)
		"""
		cdef Simplex_tree_multi_skeleton_iterator it = self.get_ptr().get_skeleton_iterator_begin(dimension)
		cdef Simplex_tree_multi_skeleton_iterator end = self.get_ptr().get_skeleton_iterator_end(dimension)
		cdef int num_parameters = self.get_ptr().get_number_of_parameters()
		while it != end:
			# yield self.get_ptr().get_simplex_and_filtration(dereference(it))
			pair = self.get_ptr().get_simplex_and_filtration(dereference(it))
			yield (np.asarray(pair.first, dtype=int),np.asarray(<value_type[:num_parameters]> pair.second))
			preincrement(it)

	def get_star(self, simplex):
		"""This function returns the star of a given N-simplex.

		:param simplex: The N-simplex, represented by a list of vertex.
		:type simplex: list of int
		:returns:  The (simplices of the) star of a simplex.
		:rtype:  list of tuples(simplex, filtration)
		"""
		cdef simplex_type csimplex = simplex
		cdef int num_parameters = self.num_parameters
		# for i in simplex:
		# 	csimplex.push_back(i)
		cdef vector[simplex_filtration_type] star \
			= self.get_ptr().get_star(csimplex)
		ct = []

		for filtered_simplex in star:
			v = []
			for vertex in filtered_simplex.first:
				v.append(vertex)
			ct.append((v, np.asarray(<value_type[:num_parameters]>filtered_simplex.second)))
		return ct

	def get_cofaces(self, simplex, codimension):
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
		cdef vector[int] csimplex = simplex
		cdef int num_parameters = self.num_parameters
		# for i in simplex:
		# 	csimplex.push_back(i)
		cdef vector[simplex_filtration_type] cofaces \
			= self.get_ptr().get_cofaces(csimplex, <int>codimension)
		ct = []
		for filtered_simplex in cofaces:
			v = []
			for vertex in filtered_simplex.first:
				v.append(vertex)
			ct.append((v, np.asarray(<value_type[:num_parameters]>filtered_simplex.second)))
		return ct

	def get_boundaries(self, simplex):
		"""This function returns a generator with the boundaries of a given N-simplex.
		If you do not need the filtration values, the boundary can also be obtained as
		:code:`itertools.combinations(simplex,len(simplex)-1)`.

		:param simplex: The N-simplex, represented by a list of vertex.
		:type simplex: list of int.
		:returns:  The (simplices of the) boundary of a simplex
		:rtype:  generator with tuples(simplex, filtration)
		"""
		cdef pair[Simplex_tree_multi_boundary_iterator, Simplex_tree_multi_boundary_iterator] it =  self.get_ptr().get_boundary_iterators(simplex)

		# while it.first != it.second:
		# 	yield self.get_ptr().get_simplex_and_filtration(dereference(it.first))
		# 	preincrement(it.first)
		cdef int num_parameters = self.get_ptr().get_number_of_parameters()
		while it.first != it.second:
			# yield self.get_ptr().get_simplex_and_filtration(dereference(it))
			pair = self.get_ptr().get_simplex_and_filtration(dereference(it.first))
			yield (np.asarray(pair.first, dtype=int),np.asarray(<value_type[:num_parameters]> pair.second))
			preincrement(it.first)
	def remove_maximal_simplex(self, simplex):
		"""This function removes a given maximal N-simplex from the simplicial
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
		"""
		self.get_ptr().remove_maximal_simplex(simplex)

	# def prune_above_filtration(self, filtration)->bool:
	# 	"""Prune above filtration value given as parameter.

	# 	:param filtration: Maximum threshold value.
	# 	:type filtration: float
	# 	:returns: The filtration modification information.
	# 	:rtype: bool


	# 	.. note::

	# 		Note that the dimension of the simplicial complex may be lower
	# 		after calling
	# 		:func:`prune_above_filtration`
	# 		than it was before. However,
	# 		:func:`upper_bound_dimension`
	# 		will return the old value, which remains a
	# 		valid upper bound. If you care, you can call
	# 		:func:`dimension`
	# 		method to recompute the exact dimension.
	# 	"""
	# 	return self.get_ptr().prune_above_filtration(filtration)
	def prune_above_dimension(self, int dimension):
		"""Remove all simplices of dimension greater than a given value.

		:param dimension: Maximum dimension value.
		:type dimension: int
		:returns: The modification information.
		:rtype: bool
		"""
		return self.get_ptr().prune_above_dimension(dimension)
	def expansion(self, int max_dim)->SimplexTreeMulti:
		"""Expands the simplex tree containing only its one skeleton
		until dimension max_dim.

		The expanded simplicial complex until dimension :math:`d`
		attached to a graph :math:`G` is the maximal simplicial complex of
		dimension at most :math:`d` admitting the graph :math:`G` as
		:math:`1`-skeleton.
		The filtration value assigned to a simplex is the maximal filtration
		value of one of its edges.

		The simplex tree must contain no simplex of dimension bigger than
		1 when calling the method.

		:param max_dim: The maximal dimension.
		:type max_dim: int
		"""
		with nogil:
			self.get_ptr().expansion(max_dim)
			# This is a fix for multipersistence. FIXME expansion in c++
			self.get_ptr().make_filtration_non_decreasing()
		return self

	def make_filtration_non_decreasing(self)->bool: 
		"""This function ensures that each simplex has a higher filtration
		value than its faces by increasing the filtration values.

		:returns: True if any filtration value was modified,
			False if the filtration was already non-decreasing.
		:rtype: bool
		"""
		cdef bool out
		with nogil:
			out = self.get_ptr().make_filtration_non_decreasing()
		return out

	def reset_filtration(self, filtration, min_dim = 0):
		"""This function resets the filtration value of all the simplices of dimension at least min_dim. Resets all the
		simplex tree when `min_dim = 0`.
		`reset_filtration` may break the filtration property with `min_dim > 0`, and it is the user's responsibility to
		make it a valid filtration (using a large enough `filt_value`, or calling `make_filtration_non_decreasing`
		afterwards for instance).

		:param filtration: New threshold value.
		:type filtration: float.
		:param min_dim: The minimal dimension. Default value is 0.
		:type min_dim: int.
		"""
		self.get_ptr().reset_filtration(Finitely_critical_multi_filtration(<python_filtration_type>filtration), min_dim)

	

	# def extend_filtration(self):
	#     """ Extend filtration for computing extended persistence. This function only uses the filtration values at the
	#     0-dimensional simplices, and computes the extended persistence diagram induced by the lower-star filtration
	#     computed with these values.
	#
	#     .. note::
	#
	#         Note that after calling this function, the filtration values are actually modified within the simplex tree.
	#         The function :func:`extended_persistence` retrieves the original values.
	#
	#     .. note::
	#
	#         Note that this code creates an extra vertex internally, so you should make sure that the simplex tree does
	#         not contain a vertex with the largest possible value (i.e., 4294967295).
	#
	#     This `notebook <https://github.com/GUDHI/TDA-tutorial/blob/master/Tuto-GUDHI-extended-persistence.ipynb>`_
	#     explains how to compute an extension of persistence called extended persistence.
	#     """
	#     self.get_ptr().compute_extended_filtration()

	# def extended_persistence(self, homology_coeff_field=11, min_persistence=0):
	#     """This function retrieves good values for extended persistence, and separate the diagrams into the Ordinary,
	#     Relative, Extended+ and Extended- subdiagrams.
	#
	#     :param homology_coeff_field: The homology coefficient field. Must be a prime number. Default value is 11. Max is 46337.
	#     :type homology_coeff_field: int
	#     :param min_persistence: The minimum persistence value (i.e., the absolute value of the difference between the
	#         persistence diagram point coordinates) to take into account (strictly greater than min_persistence).
	#         Default value is 0.0. Sets min_persistence to -1.0 to see all values.
	#     :type min_persistence: float
	#     :returns: A list of four persistence diagrams in the format described in :func:`persistence`. The first one is
	#         Ordinary, the second one is Relative, the third one is Extended+ and the fourth one is Extended-.
	#         See https://link.springer.com/article/10.1007/s10208-008-9027-z and/or section 2.2 in
	#         https://link.springer.com/article/10.1007/s10208-017-9370-z for a description of these subtypes.
	#
	#     .. note::
	#
	#         This function should be called only if :func:`extend_filtration` has been called first!
	#
	#     .. note::
	#
	#         The coordinates of the persistence diagram points might be a little different than the
	#         original filtration values due to the internal transformation (scaling to [-2,-1]) that is
	#         performed on these values during the computation of extended persistence.
	#
	#     This `notebook <https://github.com/GUDHI/TDA-tutorial/blob/master/Tuto-GUDHI-extended-persistence.ipynb>`_
	#     explains how to compute an extension of persistence called extended persistence.
	#     """
	#     cdef vector[pair[int, pair[value_type, value_type]]] persistence_result
	#     if self.pcohptr != NULL:
	#         del self.pcohptr
	#     self.pcohptr = new Simplex_tree_persistence_interface(self.get_ptr(), False)
	#     self.pcohptr.compute_persistence(homology_coeff_field, -1.)
	#     return self.pcohptr.compute_extended_persistence_subdiagrams(min_persistence)

	# TODO : cython3
	# def expansion_with_blocker(self, max_dim, blocker_func):
	# 	"""Expands the Simplex_tree containing only a graph. Simplices corresponding to cliques in the graph are added
	# 	incrementally, faces before cofaces, unless the simplex has dimension larger than `max_dim` or `blocker_func`
	# 	returns `True` for this simplex.

	# 	The function identifies a candidate simplex whose faces are all already in the complex, inserts it with a
	# 	filtration value corresponding to the maximum of the filtration values of the faces, then calls `blocker_func`
	# 	with this new simplex (represented as a list of int). If `blocker_func` returns `True`, the simplex is removed,
	# 	otherwise it is kept. The algorithm then proceeds with the next candidate.

	# 	.. warning::
	# 		Several candidates of the same dimension may be inserted simultaneously before calling `blocker_func`, so
	# 		if you examine the complex in `blocker_func`, you may hit a few simplices of the same dimension that have
	# 		not been vetted by `blocker_func` yet, or have already been rejected but not yet removed.

	# 	:param max_dim: Expansion maximal dimension value.
	# 	:type max_dim: int
	# 	:param blocker_func: Blocker oracle.
	# 	:type blocker_func: Callable[[List[int]], bool]
	# 	"""
	# 	self.get_ptr().expansion_with_blockers_callback(max_dim, callback, <void*>blocker_func)

	# def persistence(self, homology_coeff_field=11, min_persistence=0, persistence_dim_max = False):
	#     """This function computes and returns the persistence of the simplicial complex.
	#
	#     :param homology_coeff_field: The homology coefficient field. Must be a
	#         prime number. Default value is 11. Max is 46337.
	#     :type homology_coeff_field: int
	#     :param min_persistence: The minimum persistence value to take into
	#         account (strictly greater than min_persistence). Default value is
	#         0.0.
	#         Set min_persistence to -1.0 to see all values.
	#     :type min_persistence: float
	#     :param persistence_dim_max: If true, the persistent homology for the
	#         maximal dimension in the complex is computed. If false, it is
	#         ignored. Default is false.
	#     :type persistence_dim_max: bool
	#     :returns: The persistence of the simplicial complex.
	#     :rtype:  list of pairs(dimension, pair(birth, death))
	#     """
	#     self.compute_persistence(homology_coeff_field, min_persistence, persistence_dim_max)
	#     return self.pcohptr.get_persistence()
				
		
## This function is only meant for the edge collapse interface.
	def get_edge_list(self):
		cdef edge_list out
		with nogil:
			out = self.get_ptr().get_edge_list()
		return out
	
	def collapse_edges(self, max_dimension:int=None, num:int=1, progress:bool=False, strong:bool=True, full:bool=False, ignore_warning:bool=False)->SimplexTreeMulti:
		"""Edge collapse for 1-critical 2-parameter clique complex (see https://arxiv.org/abs/2211.05574).
		It uses the code from the github repository https://github.com/aj-alonso/filtration_domination .

		Parameters
		----------
		max_dimension:int
			Max simplicial dimension of the complex. Unless specified, keeps the same dimension.
		num:int
			The number of collapses to do.
		strong:bool
			Whether to use strong collapses or standard collapses (slower, but may remove more edges)
		full:bool
			Collapses the maximum number of edges if true, i.e., will do (at most) 100 strong collapses and (at most) 100 non-strong collapses afterward.
		progress:bool
			If true, shows the progress of the number of collapses.

		WARNING
		-------
			- This will destroy all of the k-simplices, with k>=2. Be sure to use this with a clique complex, if you want to preserve the homology >= dimension 1.
			- This is for 1 critical simplices, with 2 parameter persistence.
		Returns
		-------
		self:SimplexTreeMulti
			A (smaller) simplex tree that has the same homology over this bifiltration.

		"""
		# TODO : find a way to do multiple edge collapses without python conversions.
		if num <= 0:
			return self
		assert self.num_parameters == 2, "Number of parameters has to be 2 to use edge collapses ! This is a limitation of Filtration-domination"
		if self.dimension > 1 and not ignore_warning: warn("This method ignores simplices of dimension > 1 !")
		
		max_dimension = self.dimension if max_dimension is None else max_dimension

		# Retrieves the edge list, and send it to filration_domination
		edges = self.get_edge_list()
		from gudhi.multiparameter_edge_collapse import _collapse_edge_list
		edges = _collapse_edge_list(edges, num=num, full=full, strong=strong, progress=progress)
		# Retrieves the collapsed simplicial complex
		self._reconstruct_from_edge_list(edges, swap=True, expand_dimension=max_dimension)
		return self
	def _reconstruct_from_edge_list(self, edges, swap:bool=True, expand_dimension:int=None)->SimplexTreeMulti:
		"""
		Generates a 1-dimensional copy of self, with the edges given as input. Useful for edge collapses

		Input
		-----
		 - edges : Iterable[(int,int),(float,float)] ## This is the format of the rust library filtration-domination
		 - swap : bool
		 	If true, will swap self and the collapsed simplextrees.
		 - expand_dim : int
		 	expands back the simplextree to this dimension
		Ouput
		-----
		The reduced SimplexTreeMulti having only these edges.
		"""
		reduced_tree = SimplexTreeMulti(num_parameters=self.num_parameters)

		## Adds vertices back, with good filtration
		if self.num_vertices > 0:
			vertices = np.asarray([splx for splx, f in self.get_skeleton(0)], dtype=int).T
			vertices_filtration = np.asarray([f for splx, f in self.get_skeleton(0)], dtype=np.float32)
			reduced_tree.insert_batch(vertices, vertices_filtration)
		
		## Adds edges again
		if self.num_simplices - self.num_vertices > 0:
			edges_filtration = np.asarray([f for e,f in edges], dtype=np.float32)
			edges = np.asarray([e for e, _ in edges], dtype=int).T
			reduced_tree.insert_batch(edges, edges_filtration)
		if swap:
			# Swaps the simplextrees pointers
			self.thisptr, reduced_tree.thisptr = reduced_tree.thisptr, self.thisptr # Swaps self and reduced tree (self is a local variable)
		if expand_dimension is not None:
			self.expansion(expand_dimension) # Expands back the simplextree to the original dimension.
		return self if swap else reduced_tree
	
	@property
	def num_parameters(self)->int:
		return self.get_ptr().get_number_of_parameters()
	def get_simplices_of_dimension(self, dim:int)->np.ndarray:
		return np.asarray(self.get_ptr().get_simplices_of_dimension(dim), dtype=int)
	def key(self, simplex:list|np.ndarray):
		return self.get_ptr().get_key(simplex)
	def set_keys_to_enumerate(self)->None:
		self.get_ptr().set_keys_to_enumerate()
		return
	def set_key(self,simplex:list|np.ndarray, key:int)->None:
		self.get_ptr().set_key(simplex, key)
		return
	
	
	def to_scc(self, path="scc_dataset.txt", progress:bool=True, overwrite:bool=False, ignore_last_generators:bool=True, strip_comments:bool=False, reverse_block:bool=True, rivet_compatible=False)->None:
		""" Create a file with the scc2020 standard, representing the n-filtration of the simplextree.
		Link : https://bitbucket.org/mkerber/chain_complex_format/src/master/

		Parameters
		----------
		path:str
			path of the file.
		ignore_last_generators:bool = True
			If false, will include the filtration values of the last free persistence module.
		progress:bool = True
			Shows the progress bar.
		overwrite:bool = False
			If true, will overwrite the previous file if it already exists.
		ignore_last_generators:bool=True
			If true, does not write the final generators to the file. Rivet ignores them.
		reverse_block:bool=True
			Some obscure programs reverse the inside-block order.
		rivet_compatible:bool=False
			Returns a firep (old scc2020) format instead. Only Rivet uses this.

		Returns
		-------
		Nothing
		"""
		### initialize keys
		self.set_keys_to_enumerate()
		### File 
		from os.path import exists
		from os import remove
		if exists(path):
			if not(overwrite):
				raise Exception(f"The file {path} already exists. Use the `overwrite` flag if you want to overwrite.")
			remove(path)
		file = open(path, "a")
		file.write("scc2020\n") if not rivet_compatible else file.write("firep\n")
		if not strip_comments and not rivet_compatible: file.write("# Number of parameters\n")
		num_parameters = self.get_ptr().get_number_of_parameters()
		if rivet_compatible:
			assert num_parameters == 2
			file.write("Filtration 1\n")
			file.write("Filtration 2\n")
		else:
			file.write(f"{num_parameters}\n") 
		if not strip_comments: file.write("# Sizes of generating sets\n")
		## WRITES TSR VARIABLES
		tsr:list= [0]*(self.dimension+1) # dimension --- 0
		for splx,f in self.get_simplices():
			dim = len(splx)-1
			tsr[dim] += (int)(len(f) // num_parameters)
		tsr.reverse()
		file.write(" ".join([str(n) for n in tsr])+"\n")

		## Adds the boundaries to the dictionnary + tsr
		dict_splx_to_firep_number = {}
		tsr:list = [[] for _ in range(len(tsr))] # tsr stores simplices vertices, according to dimension, and the dictionnary
		for dim in range(self.dimension,-1 , -1): # range(2,-1,-1):
			for splx,F in self.get_skeleton(dim):
				if len(splx) != dim+1:	continue
				for b,_ in self.get_boundaries(splx):
					if not self.key(b) in dict_splx_to_firep_number:
						dict_splx_to_firep_number[self.key(b)] = len(tsr[dim-1])
						tsr[dim-1].append(b)
		
		## Adds simplices that are not borders to tsr, i.e., simplices not in the dictionnary 
		for splx,_ in self.get_simplices():
			if not self.key(splx) in dict_splx_to_firep_number:
				tsr[len(splx)-1].append(splx)
		## Writes simplices of tsr to file
		dim_range = range(self.dimension,0,-1) if ignore_last_generators else range(self.dimension,-1,-1)
		for dim in dim_range: # writes block by block
			if not strip_comments: file.write(f"# Block of dimension {dim}\n")
			if reverse_block:	tsr[dim].reverse()
			for splx in tsr[dim]: # for simplices of dimension
				F = np.concatenate(self.filtration(splx), axis=0)
				nbirth = (int)(len(F)//num_parameters)
				for i in range(nbirth):
					simplex_filtration = F[i*num_parameters:(i+1)*num_parameters]
					file.write(" ".join([str(f) for f in simplex_filtration]))
					file.write(" ;")
					for b,_ in self.get_boundaries(splx):
						file.write(f" {dict_splx_to_firep_number[self.key(b)]}")
					file.write("\n")
		file.close()
		return
	
	def to_rivet(self, path="rivet_dataset.txt", degree:int|None = None, progress:bool=False, overwrite:bool=False, xbins:int|None=None, ybins:int|None=None)->None:
		""" Create a file that can be imported by rivet, representing the filtration of the simplextree.

		Parameters
		----------
		path:str
			path of the file.
		degree:int
			The homological degree to ask rivet to compute.
		progress:bool = True
			Shows the progress bar.
		overwrite:bool = False
			If true, will overwrite the previous file if it already exists.
		Returns
		-------
		Nothing
		"""
		...
		from os.path import exists
		from os import remove
		if exists(path):
			if not(overwrite):
				print(f"The file {path} already exists. Use the `overwrite` flag if you want to overwrite.")
				return
			remove(path)
		file = open(path, "a")
		file.write("--datatype bifiltration\n")
		file.write(f"--homology {degree}\n") 	if degree is not None else None
		file.write(f"-x {xbins}\n") 			if xbins is not None else None
		file.write(f"-y {ybins}\n") 			if ybins is not None else None
		file.write("--xlabel time of appearance\n")
		file.write("--ylabel density\n\n")
		from tqdm import tqdm
		with tqdm(total=self.num_simplices, position=0, disable = not(progress), desc="Writing simplex to file") as bar:
			for dim in range(0,self.dimension+1): # Not sure if dimension sort is necessary for rivet. Check ?
				file.write(f"# block of dimension {dim}\n")
				for s,F in self.get_skeleton(dim):
					if len(s) != dim+1:	continue
					for i in s:
						file.write(str(i) + " ")
					file.write("; ")
					for f in F:
						file.write(str(f) + " ")
					file.write("\n")
					bar.update(1)
		file.close()
		return



	def _get_filtration_values(self, vector[int] degrees, bool inf_to_nan:bool=False)->Iterable[np.ndarray]:
		# cdef vector[int] c_degrees = degrees
		cdef intptr_t ptr = self.thisptr
		cdef vector[vector[vector[value_type]]] out
		with nogil:
			out = get_filtration_values(ptr, degrees)
		filtrations_values =  [np.asarray(filtration) for filtration in out]
		# Removes infs
		if inf_to_nan:
			for i,f in enumerate(filtrations_values):
				filtrations_values[i][f == np.inf] = np.nan
				filtrations_values[i][f == - np.inf] = np.nan
		return filtrations_values
	
	@staticmethod
	def _reduce_grid(filtrations_values,resolutions=None, strategy:str="exact", bool unique=True, some_float _q_factor=1., drop_quantiles=[0,0]):
		num_parameters = len(filtrations_values)
		if resolutions is None and strategy not in ["exact", "precomputed"]:
			raise ValueError("Resolutions must be provided for this strategy.")
		elif resolutions is not None:
			try:
				int(resolutions)
				resolutions = [resolutions]*num_parameters
			except:
				pass
		try:
			a,b=drop_quantiles
		except:
			a,b=drop_quantiles,drop_quantiles
		
		if a != 0 or b != 0:
			boxes = np.asarray([np.quantile(filtration, [a, b], axis=1, method='closest_observation') for filtration in filtrations_values])
			min_filtration, max_filtration = np.min(boxes, axis=(0,1)), np.max(boxes, axis=(0,1)) # box, birth/death, filtration
			filtrations_values = [
				filtration[(m<filtration) * (filtration <M)] 
				for filtration, m,M in zip(filtrations_values, min_filtration, max_filtration)
			]

		## match doesn't work with cython BUG
		if strategy == "exact":
			F=[np.unique(f) for f in filtrations_values]
		elif strategy == "quantile":
			F = [np.unique(f) for f in filtrations_values]
			max_resolution = [min(len(f),r) for f,r in zip(F,resolutions)]
			F = [np.quantile(f, q=np.linspace(0,1,num=int(r*_q_factor)), axis=0, method='closest_observation') for f,r in zip(F, resolutions)]
			if unique:
				F = [np.unique(f) for f in F]
				if np.all(np.asarray(max_resolution) > np.asarray([len(f) for f in F])):
					return SimplexTreeMulti._reduce_grid(filtrations_values=filtrations_values, resolutions=resolutions, strategy="quantile",_q_factor=1.5*_q_factor)
		elif strategy == "regular":
			F = [np.linspace(f.min(),f.max(),num=r) for f,r in zip(filtrations_values, resolutions)]
		elif strategy == "regular_closest":
			
			F = [_todo_regular_closest(f,r, unique) for f,r in zip(filtrations_values, resolutions)]
		elif strategy == "partition":
			F = [_todo_partition(f,r, unique) for f,r in zip(filtrations_values, resolutions)]
		elif strategy == "precomputed":
			F=filtrations_values
		else:
			raise Exception("Invalid strategy. Pick either regular, regular_closest, partition, quantile, precomputed or exact.")

		return F
	
	def get_filtration_grid(self, resolution:Iterable[int]|None=None, degrees:Iterable[int]|None=None, drop_quantiles:float|tuple=0, grid_strategy:str="exact")->Iterable[np.ndarray]:
		"""
		Returns a grid over the n-filtration, from the simplextree. Usefull for grid_squeeze. TODO : multicritical

		Parameters
		----------
			resolution: list[int]
				resolution of the grid, for each parameter
			box=None : pair[list[float]]
				Grid bounds. format : [low bound, high bound]
				If None is given, will use the filtration bounds of the simplextree.
			grid_strategy="regular" : string
				Either "regular", "quantile", or "exact".
		Returns
		-------
			List of filtration values, for each parameter, defining the grid.
		"""
		if degrees is None:
			degrees = range(self.dimension+1)
		

		## preprocesses the filtration values:
		filtrations_values = np.concatenate(self._get_filtration_values(degrees, inf_to_nan=True), axis=1)
		# removes duplicate + sort (nan at the end)
		filtrations_values = [np.unique(filtration) for filtration in filtrations_values]
		# removes nan
		filtrations_values = [filtration[:-1] if np.isnan(filtration[-1]) else filtration for filtration in filtrations_values]
		
		return self._reduce_grid(filtrations_values=filtrations_values, resolutions=resolution,strategy=grid_strategy,drop_quantiles=drop_quantiles)
	
	

	def grid_squeeze(self, filtration_grid:np.ndarray|list|None=None, bool coordinate_values=True, force=False, **filtration_grid_kwargs)->SimplexTreeMulti:
		"""
		Fit the filtration of the simplextree to a grid.
		
		:param filtration_grid: The grid on which to squeeze. An example of grid can be given by the `get_filtration_grid` method.
		:type filtration_grid: list[list[float]]
		:param coordinate_values: If true, the filtrations values of the simplices will be set to the coordinate of the filtration grid.
		:type coordinate_values: bool
		"""
		if not force and self._is_squeezed:
			raise Exception("SimplexTree already squeezed, use `force=True` if that's really what you want to do.") 
		#TODO : multi-critical
		if filtration_grid is None:	filtration_grid = self.get_filtration_grid(**filtration_grid_kwargs)
		cdef vector[vector[value_type]] c_filtration_grid = filtration_grid
		cdef intptr_t ptr = self.thisptr
		if coordinate_values:
			self.filtration_grid = c_filtration_grid
		with nogil:
			squeeze_filtration(ptr, c_filtration_grid, coordinate_values)
		return self

	@property
	def _is_squeezed(self)->bool:
		return self.num_vertices > 0 and self.filtration_grid[0].size() > 0

	def filtration_bounds(self, degrees:Iterable[int]|None=None, q:float|tuple=0, split_dimension:bool=False)->np.ndarray:
		"""
		Returns the filtrations bounds of the finite filtration values.
		"""
		try:
			a,b =q
		except:
			a,b,=q,q
		degrees = range(self.dimension+1) if degrees is None else degrees
		filtrations_values = self._get_filtration_values(degrees, inf_to_nan=True) ## degree, parameter, pt
		boxes = np.array([np.nanquantile(filtration, [a, 1-b], axis=1) for filtration in filtrations_values],dtype=float)
		if split_dimension: return boxes
		return np.asarray([np.nanmin(boxes, axis=(0,1)), np.nanmax(boxes, axis=(0,1))]) # box, birth/death, filtration


	

	def fill_lowerstar(self, F, parameter:int)->SimplexTreeMulti:
		""" Fills the `dimension`th filtration by the lower-star filtration defined by F.

		Parameters
		----------
		F:1d array
			The density over the vertices, that induces a lowerstar filtration.
		parameter:int
			Which filtration parameter to fill. /!\ python starts at 0.

		Returns
		-------
		self:SimplexTreeMulti
		"""
		# for s, sf in self.get_simplices():
		# 	self.assign_filtration(s, [f if i != dimension else np.max(np.array(F)[s]) for i,f in enumerate(sf)])
		cdef int c_parameter = parameter
		cdef vector[value_type] c_F = F
		with nogil:
			self.get_ptr().fill_lowerstar(c_F, c_parameter)
		return self


	def project_on_line(self, parameter:int=0, basepoint:None|list|np.ndarray= None)->SimplexTree:
		"""Converts an multi simplextree to a gudhi simplextree.
		Parameters
		----------
			parameter:int = 0
				The parameter to keep. WARNING will crash if the multi simplextree is not well filled.
			basepoint:None
				Instead of keeping a single parameter, will consider the filtration defined by the diagonal line crossing the basepoint.
		WARNING 
		-------
			There are no safeguard yet, it WILL crash if asking for a parameter that is not filled.
		Returns
		-------
			A SimplexTree with chosen 1D filtration.
		"""
		# FIXME : deal with multicritical filtrations
		import gudhi as gd
		new_simplextree = gd.SimplexTree()
		assert parameter < self.get_ptr().get_number_of_parameters()
		cdef int c_parameter = parameter
		cdef intptr_t old_ptr = self.thisptr
		cdef intptr_t new_ptr = new_simplextree.thisptr
		cdef vector[value_type] c_basepoint = [] if basepoint is None else basepoint
		if basepoint is None:
			with nogil:
				flatten_from_ptr(old_ptr, new_ptr, c_parameter)
		else:
			with nogil:
				flatten_diag_from_ptr(old_ptr, new_ptr, c_basepoint, c_parameter)
		return new_simplextree

	def linear_projections(self, linear_forms:np.ndarray)->Iterable[SimplexTree]:
		"""
		Compute the 1-parameter projections, w.r.t. given the linear forms, of this simplextree.

		Input
		-----
		 - Array of shape (num_linear_forms, num_parameters)
		
		Output
		------
		 - List of projected (gudhi) simplextrees.
		"""
		cdef Py_ssize_t num_projections = linear_forms.shape[0]
		cdef Py_ssize_t num_parameters = linear_forms.shape[1]
		if num_projections == 0:	return []
		cdef vector[vector[value_type]] c_linear_forms = linear_forms
		assert num_parameters==self.num_parameters, f"The linear forms has to have the same number of parameter as the simplextree ({self.num_parameters})."
		
		# Gudhi copies are faster than inserting simplices one by one
		import gudhi as gd
		flattened_simplextree = gd.SimplexTree()
		cdef intptr_t multi_prt = self.thisptr
		cdef intptr_t flattened_ptr = flattened_simplextree.thisptr
		with nogil:
			flatten_from_ptr(multi_prt, flattened_ptr, num_parameters)
		out = [flattened_simplextree] + [gd.SimplexTree(flattened_simplextree) for _ in range(num_projections-1)]

		# Fills the 1-parameter simplextrees.
		cdef vector[intptr_t] out_ptrs = [st.thisptr for st in out]
		with nogil:
			for i in range(num_projections):
				linear_projection_from_ptr(out_ptrs[i], multi_prt, c_linear_forms[i])
		return out


	def set_num_parameter(self, num:int):
		"""
		Sets the numbers of parameters. 
		WARNING : it will resize all the filtrations to this size. 
		"""
		self.get_ptr().resize_all_filtrations(num)
		self.get_ptr().set_number_of_parameters(num)
		return

	def __eq__(self, other:SimplexTreeMulti):
		"""Test for structural equality
		:returns: True if the 2 simplex trees are equal, False otherwise.
		:rtype: bool
		"""
		return dereference(self.get_ptr()) == dereference(other.get_ptr())
	
cdef intptr_t _get_copy_intptr(SimplexTreeMulti stree) nogil:
	return <intptr_t>(new Simplex_tree_multi_interface(dereference(stree.get_ptr())))


def _todo_regular_closest(cnp.ndarray[some_float,ndim=1] f, int r, bool unique):
	f_regular = np.linspace(np.min(f),np.max(f),num=r)
	f_regular_closest = np.asarray([f[np.argmin(np.abs(f-x))] for x in f_regular])
	if unique: f_regular_closest = np.unique(f_regular_closest)
	return f_regular_closest

def _todo_partition(cnp.ndarray[some_float,ndim=1] data,int resolution, bool unique):
	if data.shape[0] < resolution: resolution=data.shape[0]
	k = data.shape[0] // resolution
	partitions = np.partition(data, k)
	f = partitions[[i*k for i in range(resolution)]]
	if unique: f= np.unique(f)
	return f



def _simplextree_multify(simplextree:SimplexTree, num_parameters:int=2, default_values=[])->SimplexTreeMulti:
	"""Converts a gudhi simplextree to a multi simplextree.
	Parameters
	----------
		parameters:int = 2
			The number of filtrations
	Returns
	-------
		A multi simplextree, with first filtration value being the one from the original simplextree.
	"""
	if isinstance(simplextree, SimplexTreeMulti):
		return simplextree
	st = SimplexTreeMulti(num_parameters=num_parameters)
	cdef int c_num_parameters = num_parameters
	cdef intptr_t old_ptr = simplextree.thisptr
	cdef intptr_t new_ptr = st.thisptr
	cdef vector[value_type] c_default_values=default_values
	with nogil:
		multify_from_ptr(old_ptr, new_ptr, c_num_parameters, c_default_values)
	return st

def _safe_simplextree_multify(simplextree:SimplexTree,int num_parameters=2)->SimplexTreeMulti:
	if isinstance(simplextree, SimplexTreeMulti):
		return simplextree
	simplices = [[] for _ in range(simplextree.dimension()+1)]
	filtration_values = [[] for _ in range(simplextree.dimension()+1)]
	for s,f in simplextree.get_simplices():
		filtration_values[len(s)-1].append(f)
		simplices[len(s)-1].append(s)
	st_multi = SimplexTreeMulti(num_parameters=1)
	for batch_simplices, batch_filtrations in zip(simplices,filtration_values):
		st_multi.insert_batch(np.asarray(batch_simplices).T, np.asarray(batch_filtrations)[:,None])
	if num_parameters > 1:
		st_multi.set_num_parameter(num_parameters)
	return st_multi 
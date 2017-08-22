from cython cimport numeric
from libcpp.vector cimport vector #use similar line if you need string, pair, etc....
from libcpp.utility cimport pair
from libcpp cimport bool  
import os

"""This file is part of the Gudhi Library. The Gudhi library
   (Geometric Understanding in Higher Dimensions) is a generic C++
   library for computational topology.

   Author(s):       Pawel Dlotko

   Copyright (C) 2017 Swansea University

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

__author__ = "Pawel Dlotko"
__copyright__ = "Copyright (C) 2017 Swansea University"
__license__ = "GPL v3"






"""
This is a promisse that there will be a class in this file with the following function signature. Something like C++ predeclaration.
According to Vincent, most of the tutorials in cython suggest to separate pre-declaration below with the definition of the method. Hovewer it seems to create problems, that is why we keep them both here. 
"""
    
cdef extern from "persistence_representations_intervals.h" namespace "Gudhi::Persistence_representations":
	cdef cppclass Persistence_intervals_interface "Gudhi::Persistence_representations::Persistence_intervals_interface":        
		Persistence_intervals_interface()
		Persistence_intervals_interface(const char* , unsigned )
		Persistence_intervals_interface(const vector[pair[double, double]] intervals)
		pair[double, double] get_x_range_interface() const 
		pair[double, double] get_y_range_interface() const 	
		vector[double] length_of_dominant_intervals_interface(size_t where_to_cut) const		
		vector[pair[double, double] ] dominant_intervals_interface(size_t where_to_cut) const
		vector[size_t] histogram_of_lengths_interface(size_t number_of_bins ) const
		vector[size_t] cumulative_histogram_of_lengths_interface(size_t number_of_bins) const
		vector[double] characteristic_function_of_diagram_interface(double x_min, double x_max, size_t number_of_bins) const
		vector[double] cumulative_characteristic_function_of_diagram_interface(double x_min, double x_max, size_t number_of_bins) const                                                              
		vector[pair[double, size_t] ] compute_persistent_betti_numbers_interface() const
		double project_to_R_interface(int number_of_function) const
		size_t number_of_projections_to_R_interface() const
		vector[double] vectorize_interface(int number_of_function) const
		size_t number_of_vectorize_functions_interface() const       

"""                 
make sure that here we call the functions from the intermediate .h file, with dummy names, so that later below we can use the same names of the functions as in C++ version. 
Over here I need to list all the functions that will be used in the file. So there should be a list of constructors, methors, etc.
to separate the function, use newline. Put there only C++ signature
"""

 
 
#convention for python class is PersistenceIntervals instead of Persistence_intervals
#for methods it is def num_simplices(self).
cdef class PersistenceIntervals:
	"""
	Persistence intrvals is a standard representation of persistent homology. This file provide implementation of a number of operations on persistence diagrams.
	"""
        
	cdef Persistence_intervals_interface * thisptr
    
	#do we need a fake constructor here, as in case of bitmaps??
	#We do need it so that we have a doc for python becaus ethe documentation only read from __init__, it do not read from __cinit__
	#__ means private memeber
	def __init__(self, vector_of_intervals=None, dimension = None, file_with_intervals=''):
		"""Persistence interals is a standard representation of persistent homology. This file provide implementation of a number of operations on persistence diagrams. 

		:param dimensions: A vector of birth-death pairs.	

		Or

		:param Gudhi style file togethr with a dimension of birth-death pairs to consider.	
		"""


	# The real cython constructor
	def __cinit__(self, vector_of_intervals=None, dimension = None, file_with_intervals=''):
		if (vector_of_intervals is None) and (file_with_intervals is not ''):
			self.thisptr = new Persistence_intervals_interface(file_with_intervals, dimension) 
		elif (file_with_intervals is not '') and (vector_of_intervals is not None):
			if os.path.isfile(file_with_intervals):
				self.thisptr = new Persistence_intervals_interface(str.encode(file_with_intervals))
			else:
				print("file " + file_with_intervals + " not found.")
		else:
			print("Persistence interals can be constructed from vector of birth-death pairs,  vector_of_intervals or a Gudhi-style file.")
			
			
			
	def __dealloc__(self):
		if self.thisptr != NULL:
			del self.thisptr      
		
		
	#from here on this is my try. Do we need to specify the returned type?? 
	#no, we do not. 
	def get_x_range(self):
		if self.thisptr != NULL:
			return self.thisptr.get_x_range_interface()
								
	def get_y_range(self):
		if self.thisptr != NULL:
			return self.thisptr.get_y_range_interface()

	def length_of_dominant_intervals(self,where_to_cut):
		if (self.thisptr != NULL) and (where_to_cut is not None):
			return self.thisptr.length_of_dominant_intervals_interface(where_to_cut)   
		else:
			if (self.thisptr != NULL):
				return self.thisptr.dominant_intervals_interface(100)#default argument	 		     

	def dominant_intervals(self,where_to_cut):
		if (self.thisptr != NULL) and (where_to_cut is not None):
			return self.thisptr.dominant_intervals_interface(where_to_cut)   	 
		else:
			if (self.thisptr != NULL):
				return self.thisptr.dominant_intervals_interface(100)#default argument 	 
		 
	def histogram_of_lengths(self,number_of_bins):
		if (self.thisptr != NULL) and (number_of_bins is not None):
			return self.thisptr.histogram_of_lengths_interface(number_of_bins)	
		else:
			if (self.thisptr != NULL):
				return self.thisptr.dominant_intervals_interface(100) #default argument 	 		 				
		 
	def cumulative_histogram_of_lengths(self,number_of_bins):
		if (self.thisptr != NULL) and (number_of_bins is not None):
			return self.thisptr.cumulative_histogram_of_lengths_interface(number_of_bins)
		else:
			if (self.thisptr != NULL):
				return self.thisptr.cumulative_histogram_of_lengths_interface(10)#default argument 	
		 
	def characteristic_function_of_diagram(self,x_min,x_max,number_of_bins):
		if (self.thisptr != NULL) and ( x_min is not None ) and ( x_max is not None ) and ( number_of_bins is not None ):
			return self.thisptr.characteristic_function_of_diagram_interface( x_min , x_max , number_of_bins )
		else:
			if (self.thisptr != NULL) and ( x_min is not None ) and ( x_max is not None ):
				return self.thisptr.characteristic_function_of_diagram_interface( x_min , x_max ,10 )	#default argument 	
		 
	def cumulative_characteristic_function_of_diagram(self,x_min,x_max,number_of_bins):
		if (self.thisptr != NULL) and ( x_min is not None ) and ( x_max is not None ) and ( number_of_bins is not None ):
			return self.thisptr.cumulative_characteristic_function_of_diagram_interface( x_min , x_max , number_of_bins )
		else:
			if (self.thisptr != NULL ) and ( x_min is not None ) and ( x_max is not None ):
				return self.thisptr.cumulative_characteristic_function_of_diagram_interface( x_min , x_max , 10)	#default argument  				
		 
	def compute_persistent_betti_numbers(self):
		if self.thisptr != NULL:
			return self.thisptr.compute_persistent_betti_numbers_interface()
					
	def project_to_R(self,number_of_function):
		if (self.thisptr != NULL) and (number_of_function is not None):
			return self.thisptr.project_to_R_interface(number_of_function)
		 
	def number_of_projections_to_R(self):
		if self.thisptr != NULL:
			return self.thisptr.number_of_projections_to_R_interface()
		 
	def vectorize(self,number_of_function):
		if (self.thisptr != NULL) and ( number_of_function is not None ):
			return self.thisptr.vectorize_interface(number_of_function)
		 
	def number_of_vectorize_functions(self):
		if (self.thisptr != NULL):
			return self.thisptr.number_of_vectorize_functions_interface()			 			 			 
			



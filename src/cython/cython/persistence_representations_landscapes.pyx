from cython cimport numeric
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp cimport bool  
from cython.operator cimport dereference as deref
import os
import sys

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



cdef extern from "Persistence_landscape_interface.h" namespace "Gudhi::Persistence_representations":
    cdef cppclass Persistence_landscape_interface "Gudhi::Persistence_representations::Persistence_landscape_interface":
        Persistence_landscape_interface()
        #Persistence_landscape_interface(vector[pair[double, double]], bool, size_t)
        #Persistence_landscape_interface(const char*, size_t , size_t)
        void load_landscape_from_file(const char*)
        void print_to_file(const char*)const
        double compute_integral_of_landscape()const
        double compute_integral_of_a_level_of_a_landscape(size_t)const
        double compute_integral_of_landscape(double)const
        double compute_value_at_a_given_point(unsigned , double)const
        double compute_maximum()const
        double compute_minimum()const
        double compute_norm_of_landscape(double)
        Persistence_landscape_interface* new_abs_interface()
        size_t size()const
        double find_max(unsigned)const
        double project_to_R(int)const
        size_t number_of_projections_to_R()const
        vector[double] vectorize(int)const
        size_t number_of_vectorize_functions()const
        void compute_average(const vector[Persistence_landscape_interface*]&)
        void new_compute_average(const vector[Persistence_landscape_interface*]&) 
        double distance(const Persistence_landscape_interface&, double)
        double compute_scalar_product(const Persistence_landscape_interface&)const
        pair[double, double] get_y_range(size_t)const
        #**************
        #static methods
        @staticmethod
        Persistence_landscape_interface* construct_from_file( const char*, size_t, size_t)
        @staticmethod
        Persistence_landscape_interface* construct_from_vector_of_pairs( const vector[pair[double, double]], size_t)
        #***************




#convention for python class is PersistenceIntervals instead of Persistence_intervals
#for methods it is def num_simplices(self).
cdef class PersistenceLandscapes:

    cdef Persistence_landscape_interface* thisptr



#Can we have only one constructor, or can we have more
    def __init__(self, vector_of_intervals=None, dimension=None, file_with_intervals='',number_of_levels=sys.maxsize):
        """
        This is a class implementing persistence landscapes data structures.
        For theoretical description, please consult <i>Statistical topological
        data analysis using persistence landscapes</i>\cite bubenik_landscapes_2015,
        and for details of algorithms, A persistence landscapes toolbox for
        topological statistics</i>\cite bubenik_dlotko_landscapes_2016.

        Persistence landscapes allow vectorization, computations of distances,
        computations of projections to Real, computations of averages and
        scalar products. Therefore they implement suitable interfaces. It
        implements the following concepts: Vectorized_topological_data,
        Topological_data_with_distances, Real_valued_topological_data,
        Topological_data_with_averages, Topological_data_with_scalar_product

        Note that at the moment, due to rounding errors during the construction
        of persistence landscapes, elements which are different by 0.000005 are
        considered the same. If the scale in your persistence diagrams is
        comparable to this value, please rescale them before use this code.
        """
   


    def __cinit__(self, vector_of_intervals=None, dimension=None, file_with_intervals='',number_of_levels=sys.maxsize):
        """
        This is a constructor of a class PersistenceLandscapes.
        It either take text file and a positive integer, or a vector
        of pairs. The last optional parameter is the nunmer of levels of
        the landscapes to be generated. If not set, all the levels will
        be generated. In case of file, each line of the input file is,
        supposed to contain two numbers of a type double
        (or convertible to double) representing the birth and the death
        of the persistence interval. If the pairs are not sorted so that
        birth <= death, then the constructor will sort then that way.
        In case of vector of pairs, it simply accept vector of pair of
        doubles.
        :param vector_of_intervals -- vector of pairs of doubles with
               birth-death pairs. None if we construct it from file.
        :type vector of pairs of doubles or None
        :param dimension -- diension of intervals to be extracted from file
        :type nonnegative integer or None
        :param file_with_intervals - a path to Gudhi style file with
               persistence interfals.
        :type string of None.
        :param number_of_levels - number of levels of landscape to be
        generated (if not set, all of the are generated).
        :type positive integer
        """
        if (vector_of_intervals is None) and (file_with_intervals is not ''):
            if (dimension is not None):
                if os.path.isfile(file_with_intervals):
                    #self.thisptr = new Persistence_landscape_interface(file_with_intervals, dimension, number_of_levels)
                    self.thisptr = Persistence_landscape_interface.construct_from_file(file_with_intervals, dimension, number_of_levels)
                else:
                    print("file " + file_with_intervals + " not found.")
            else:
                    #self.thisptr = new Persistence_landscape_interface(file_with_intervals, number_of_levels)
                    self.thisptr = Persistence_landscape_interface.construct_from_file(file_with_intervals,0, number_of_levels)
        elif (file_with_intervals is '') and (vector_of_intervals is not None):
            #self.thisptr = new Persistence_landscape_interface(vector_of_intervals, true, number_of_levels)
            self.thisptr = Persistence_landscape_interface.construct_from_vector_of_pairs(vector_of_intervals, number_of_levels)
        else:
            print("Persistence interals can be constructed from vector of birth-death pairs,  vector_of_intervals or a Gudhi-style file.")                     
            self.thisptr = new Persistence_landscape_interface()
            
    def __dealloc__(self):
        """
        destructor
        """
        if self.thisptr != NULL:
            del self.thisptr

    def load_landscape_from_file(self,filename):
        """
        This procedure loads a landscape from file. It erase all the data
        that was previously stored in this landscape.
        :param Name of the file.
        :type String
        """
        if ( self.thisptr != NULL ) and ( filename is not None ):
            self.thisptr.load_landscape_from_file(filename)

    def print_to_file(self,filename) :
        """
        The procedure stores a landscape to a file. The file can be later
        used by a procedure load_landscape_from_file.
        :param Name of the file.
        :type String
        """
        if ( self.thisptr != NULL ) and ( filename is not None ):
            self.thisptr.print_to_file(filename)

    def compute_integral_of_landscape(self):
        """
        This function compute integral of the landscape (defined formally as
        sum of integrals on R of all landscape functions)
        """
        if ( self.thisptr != NULL ):
            return self.thisptr.compute_integral_of_landscape()

    def compute_integral_of_a_level_of_a_landscape(self,level):
        """
        This function compute integral of the 'level'-level of a landscape.
        :param Level of the landscape, n, so that the integral of lambda_n is
        computed
        :type nonnegative integer.
        """
        if ( self.thisptr != NULL ) and ( level is not None ):
            return self.thisptr.compute_integral_of_landscape(level)

    def compute_integral_of_landscape(self,p):
        """
        This function compute integral of the landscape p-th power of a
        landscape (defined formally as sum of integrals on R of p-th powers
        of all landscape functions)
        :param An positive real p such that the integral of p-th power of
        landscape is computed.
        :type Real value
        """
        if ( self.thisptr != NULL ) and ( p is not None ):
            return self.thisptr.compute_integral_of_landscape(p)

    def compute_value_at_a_given_point(self, level,  x):
        """
        A function that computes the value of a landscape at a given point.
        The parameters of the function are: unsigned
        level and double x.
        The procedure will compute the value of the level-landscape at the
        point x.
        param: level n of lanscape (positive integer) and real number x.
        The value \lambda_n(x) is
        type: nonnegative integer
        :param A real number x. The value \lambda_n(x) is computed.
        :type real
        """
        if ( self.thisptr != NULL ) and ( level is not None ) and ( x is not None ):
            return self.thisptr.compute_value_at_a_given_point(level,x)

    def compute_maximum( self ):
        """
        Computations of maximum (y) value of landscape.
        """
        if ( self.thisptr != NULL ):
            return self.thisptr.compute_maximum()

    def compute_minimum( self ):
        """
        Computations of minimum (y) value of landscape.
        """
        if ( self.thisptr != NULL ):
            return self.thisptr.compute_minimum()

    def compute_norm_of_landscape(self,i):
        """
        Computations of a \f$L^i\f$ norm of landscape, where i is the input parameter.
        :type integer.
        :param i
        """
        if ( self.thisptr != NULL ) and ( i is not None ):
            return self.thisptr.compute_norm_of_landscape(i)          

    def abs( self ):
        """
        Function to compute absolute value of a PL function. The representation
        of persistence landscapes allow to store
        general PL-function. When computing distance between two landscapes,
        we compute difference between
        them. In this case, a general PL-function with negative value can
        appear as a result. Then in order to compute
        distance, we need to take its absolute value. This is the purpose of
        this procedure.
        """
        if ( self.thisptr != NULL ):
            abs_ = PersistenceLandscapes()
            abs_.thisptr = self.thisptr.new_abs_interface()
            return abs_


    def size( self ):
        """
        Computes the number of landscape functions.
        """
        if ( self.thisptr != NULL ):
            return self.thisptr.size()

    def find_max(self, lambda_):
        """
        Compute maximal value of lambda-level landscape.
        :param level of landscape
        :type nonnegative integer
        """
        if ( self.thisptr != NULL ) and ( lambda_ is not None ):
            return self.thisptr.find_max(lambda_)

    def project_to_R(self, number_of_function):
        """
        The number of projections to R is defined to the number of nonzero
        landscape functions. I-th projection is an
        integral of i-th landscape function over whole R.
        This function is required by the Real_valued_topological_data concept.
        At the moment this function is not tested, since it is quite likely
        to be changed in the future. Given this, when
        using it, keep in mind that it
        will be most likely changed in the next versions.
        :param number of function
        :type nonnegative integer
        """
        if ( self.thisptr != NULL ) and ( number_of_function is not None ):
            return self.thisptr.project_to_R(number_of_function)

    def number_of_projections_to_R(self):
        """
        The function gives the number of possible projections to R. This
        function is required by the
        Real_valued_topological_data concept
        """
        if ( self.thisptr != NULL ):
            return self.thisptr.number_of_projections_to_R()

    def vectorize(self, number_of_function):
        """
        This function produce a vector of doubles based on a landscape. It
        is required in a concept
        Vectorized_topological_data
        :param number of function
        :type nonnegative intege
        """
        if ( self.thisptr != NULL ) and ( number_of_function is not None ):
            return self.thisptr.vectorize(number_of_function)

    def number_of_vectorize_functions(self):
        """
        The number of projections to R is defined to the number of nonzero
        landscape functions. I-th projection is an
        integral of i-th landscape function over whole R.
        This function is required by the Real_valued_topological_data concept.
        At the moment this function is not tested, since it is quite likely
        to be changed in the future. Given this, when
        using it, keep in mind that it
        will be most likely changed in the next versions
        """
        if ( self.thisptr != NULL ):
            return self.thisptr.number_of_vectorize_functions()

    def compute_average( self,to_average=[] ):
        """
        A function to compute averaged persistence landscape, based on vector
        of persistence landscapes.
        This function is required by Topological_data_with_averages concept.
        :param vector of persistence landscapes to average
        :type vectors of references to persistence landscapes
        """            
        #TODO -- add a check if all objects in the to_average are of the same type. 
        cdef vector[Persistence_landscape_interface*] cpp_list    
        if ( self.thisptr != NULL ) and ( to_average is not None ):	
            for elt in to_average: 
                cpp_list.push_back((<PersistenceLandscapes>elt).thisptr)
            self.thisptr.new_compute_average( cpp_list )                
            


    def distance(self, PersistenceLandscapes second, power):
        """
        A function to compute distance between persistence landscape.
        The parameter of this function is a Persistence_landscape.
        This function is required in Topological_data_with_distances concept.
        For max norm distance, set power to std::numeric_limits<double>::max()
        :param the landascape to compute distance to
        :type PersistenceLandscape
        """
        if ( self.thisptr != NULL ) and ( second is not None ) and ( power is not None ):
            return self.thisptr.distance( deref(second.thisptr), power)
 
    def compute_scalar_product(self, PersistenceLandscapes second):
        """
        A function to compute scalar product of persistence landscapes.
        The parameter of this function is a Persistence_landscape.
        This function is required in Topological_data_with_scalar_product concept.
        :param the landascape to compute scalar product with
        :type PersistenceLandscape
        """
        if ( self.thisptr != NULL ) and ( second is not None ):
            return self.thisptr.compute_scalar_product( deref(second.thisptr) )

    def get_y_range(self, level):
        """
        This procedure returns y-range of a given level persistence landscape.
        If a default value is used, the y-range
        of 0th level landscape is given (and this range contains the ranges
        of all other landscapes).
        :param The level of lrandscape
        :type nonnegative integer
        """
        if ( self.thisptr != NULL ):
            return self.thisptr.get_y_range(level)


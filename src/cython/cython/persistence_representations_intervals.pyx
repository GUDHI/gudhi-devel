from cython cimport numeric
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp cimport bool
import os

"""
This file is part of the Gudhi Library. The Gudhi library
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
This is a promisse that there will be a class in this file with the following
function signature. Something like C++ predeclaration.
According to Vincent, most of the tutorials in cython suggest to
separate pre-declaration below with the definition of the method.
Hovewer it seems to create problems, that is why we keep them both here.
"""




cdef extern from "Persistence_intervals_interface.h" namespace "Gudhi::Persistence_representations":
    cdef cppclass Persistence_intervals_interface "Gudhi::Persistence_representations::Persistence_intervals_interface":
        Persistence_intervals_interface(const char*, unsigned)
        Persistence_intervals_interface(const vector[pair[double, double]] intervals)
        pair[double, double] get_x_range() const
        pair[double, double] get_y_range() const
        vector[double] length_of_dominant_intervals(size_t where_to_cut)const
        vector[pair[double, double]] dominant_intervals(size_t where_to_cut) const
        vector[size_t] histogram_of_lengths(size_t number_of_bins) const
        vector[size_t] cumulative_histogram_of_lengths(size_t number_of_bins) const
        vector[double] characteristic_function_of_diagram(double x_min, double x_max, size_t number_of_bins)const
        vector[double] cumulative_characteristic_function_of_diagram(double x_min, double x_max, size_t number_of_bins)const
        vector[pair[double, size_t]] compute_persistent_betti_numbers()const
        double project_to_R(int number_of_function) const
        size_t number_of_projections_to_R() const
        vector[double] vectorize(int number_of_function) const
        size_t number_of_vectorize_functions() const

"""
make sure that here we call the functions from the intermediate .h file,
with dummy names, so that later below we can use the same names of the
functions as in C++ version.
Over here I need to list all the functions that will be used in the file.
So there should be a list of constructors, methors, etc.
to separate the function, use newline. Put there only C++ signature
"""


#convention for python class is PersistenceIntervals instead of
#Persistence_intervals for methods it is def num_simplices(self).
cdef class PersistenceIntervals:
    """
    Persistence intrvals is a standard representation of persistent homology. This file provide implementation of a number of operations on persistence diagrams.
    """
    cdef Persistence_intervals_interface * thisptr

    #do we need a fake constructor here, as in case of bitmaps??
    #We do need it so that we have a doc for python because the
    #documentation only read from __init__, it do not read from
    #__cinit__, where __ means private memeber
    def __init__(self, vector_of_intervals=None, dimension=None,
                 file_with_intervals=''):
        """Persistence interals is a standard representation of
           persistent homology. This file provide implementation of a
           number of operations on persistence diagrams.

          :param dimensions: A vector of birth-death pairs.

        Or

        :param Gudhi style file togethr with a dimension of birth-death
               pairs to consider.
        """

    #The real cython constructor
    def __cinit__(self, vector_of_intervals=None, dimension=None,
                  file_with_intervals=''):
        """
        This is a constructor of a class Persistence_intervals.
        It either take text file and a positive integer, or a vector
        of pairs. In case of file, each line of the input file is
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
        """
        if (vector_of_intervals is None) and (file_with_intervals is not ''):
            if (dimension is not None):
                if os.path.isfile(file_with_intervals):
                    self.thisptr = new Persistence_intervals_interface(file_with_intervals, dimension)
                else:
                    print("file " + file_with_intervals + " not found.")
            else:
                    self.thisptr = new Persistence_intervals_interface(file_with_intervals)
        elif (file_with_intervals is '') and (vector_of_intervals is not None):            
            self.thisptr = new Persistence_intervals_interface(vector_of_intervals)            
        else:
            print("Persistence interals can be constructed from vector of birth-death pairs,  vector_of_intervals or a Gudhi-style file.")

    def __dealloc__(self):
        """
        destructor
        """
        if self.thisptr != NULL:
            del self.thisptr

    #from here on this is my try. Do we need to specify the returned type??
    #no, we do not.
    def get_x_range(self):
        """
        This procedure returns x-range of a given persistence diagram.
        """
        if self.thisptr != NULL:
                        return self.thisptr.get_x_range()

    def get_y_range(self):
        """
        This procedure returns y-range of a given persistence diagram.
        """
        if self.thisptr != NULL:
            return self.thisptr.get_y_range()

    def length_of_dominant_intervals(self, where_to_cut):
        """
        Procedure that compute the vector of lengths of the dominant
        (i.e. the longest) persistence intervals. The list is
        truncated at the parameter of the call where_to_cut
        (set by default to 100).
        :param where_to_cut -- number of domiannt intervals to be returned.
        :type positive integer.
        """
        if (self.thisptr != NULL) and (where_to_cut is not None):
            return self.thisptr.length_of_dominant_intervals(where_to_cut)
        else:
            if (self.thisptr != NULL):
                return self.thisptr.dominant_intervals(100)

    def dominant_intervals(self, where_to_cut):
        """
        Procedure that compute the vector of the dominant (i.e. the longest)
        persistence intervals. The parameter of the procedure (set by default
        to 100) is the number of dominant intervals returned by the procedure.
        :param where_to_cut -- number of lengths of domiannt intervals to
        be returned.
       :type positive integer.
        """
        if (self.thisptr != NULL) and (where_to_cut is not None):
            return self.thisptr.dominant_intervals(where_to_cut)
        else:
            if (self.thisptr != NULL):
                return self.thisptr.dominant_intervals(100)

    def histogram_of_lengths(self, number_of_bins):
        """
        Procedure to compute a histogram of interval's length.
        A histogram is a block plot. The number of blocks is
        determined by the first parameter of the function
        (set by default to 10).
        For the sake of argument let us assume that the length of the
        longest interval is 1 and the number of bins is
        10. In this case the i-th block correspond to a range between
        i-1/10 and i10. The vale of a block supported at the interval is
        the number of persistence intervals of a length between x_0
        and x_1.
        :param where_to_cut -- number of bins in the histogram.
        :type positive integer.
        """
        if (self.thisptr != NULL) and (number_of_bins is not None):
            return self.thisptr.histogram_of_lengths(number_of_bins)
        else:
            if (self.thisptr != NULL):
                return self.thisptr.dominant_intervals(100)

    def cumulative_histogram_of_lengths(self, number_of_bins):
        """
        Based on a histogram of intervals lengths computed by the
        function histogram_of_lengths H the procedure below
        computes the cumulative histogram. The i-th position
        of the resulting histogram
        is the sum of values of H for the positions from 0 to i.
        :param where_to_cut -- number of bins in the histogram.
        :type positive integer.
        """
        if (self.thisptr != NULL) and (number_of_bins is not None):
            return self.thisptr.cumulative_histogram_of_lengths(number_of_bins)
        else:
            if (self.thisptr != NULL):
                return self.thisptr.cumulative_histogram_of_lengths(10)

    def characteristic_function_of_diagram(self, x_min, x_max, number_of_bins):
        """
        In this procedure we assume that each barcode is a characteristic
        function of a hight equal to its length. The persistence diagram
        is a sum of such a functions. The procedure below construct a
        function being a sum of the characteristic functions of
        persistence intervals. The first two parameters are the range in
        which the function is to be computed and the last parameter is
        the number of bins in the discretization of the interval
        [_min,_max]
        :param x_min -- Begin of range of function.
        :type real number
        :param x_max -- End of range of function.
        :type real number
        :param number_of_bins -- Number of bins in characteristic function.
        :type positive integer
        """
        if (self.thisptr != NULL) and (x_min is not None) and (x_max is not None) and (number_of_bins is not None):
            return self.thisptr.characteristic_function_of_diagram(x_min, x_max	, number_of_bins)
        else:
            if (self.thisptr != NULL) and (x_min is not None) and (x_max is not None):
                return self.thisptr.characteristic_function_of_diagram(x_min, x_max, 10)

    def cumulative_characteristic_function_of_diagram(self, x_min, x_max, number_of_bins):
        """
        Cumulative version of the function characteristic_function_of_diagram.
        :param x_min -- Begin of range of function.
        :type real number
        :param x_max -- End of range of function.
        :type real number
        :param number_of_bins -- Number of bins in characteristic function.
        :type positive integer
        """
        if (self.thisptr != NULL) and (x_min is not None) and (x_max is not None) and (number_of_bins is not None):
            return self.thisptr.cumulative_characteristic_function_of_diagram(x_min, x_max, number_of_bins)
        else:
            if (self.thisptr != NULL) and (x_min is not None) and (x_max is not None):
                return
                self.thisptr.cumulative_characteristic_function_of_diagram(x_min, x_max, 10)

    def compute_persistent_betti_numbers(self):
        """
        Compute the function of persistence Betti numbers. The returned
        value is a vector of pair. First element of each
        pair is a place where persistence Betti numbers change.
        Second element of each pair is the value of Persistence Betti
        numbers at that point.
        """
        if self.thisptr != NULL:
            return self.thisptr.compute_persistent_betti_numbers()

    def project_to_R(self, number_of_function):
        """
        This is a simple function projecting the persistence intervals
        to a real number. The function we use here is a sum
        of squared lengths of intervals. It can be naturally interpreted as
        sum of step function, where the step hight it equal to the length
        of the interval. At the moment this function is not tested, since
        it is quite likely to be changed in the future. Given this, when
        using it, keep in mind that it
        will be most likely changed in the next versions.
        :param number of projection
        :type positive integer.
        """
        if (self.thisptr != NULL) and (number_of_function is not None):
            return self.thisptr.project_to_R(number_of_function)

    def number_of_projections_to_R(self):
        """
        The function gives the number of possible projections to R.
        This function is required by the
        Real_valued_topological_data concept.
        """
        if self.thisptr != NULL:
            return self.thisptr.number_of_projections_to_R()

    def vectorize(self, number_of_function):
        """
        Return a family of vectors obtained from the persistence diagram.
        The i-th vector consist of the length of i
        dominant persistence intervals.
        :param number of function to vectorizes
        :type positive integer.
        """
        if (self.thisptr != NULL) and (number_of_function is not None):
            return self.thisptr.vectorize(number_of_function)

    def number_of_vectorize_functions(self):
        """
        This function return the number of functions that allows
        vectorization of a persistence diagram. It is required
        in a concept Vectorized_topological_data
        """
        if (self.thisptr != NULL):
            return self.thisptr.number_of_vectorize_functions()

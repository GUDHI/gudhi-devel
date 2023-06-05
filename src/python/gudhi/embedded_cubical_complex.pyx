from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "Embedded_cubical_complex_interface.h" namespace "Gudhi":
    cdef cppclass Embedded_cubical_complex_base_interface "Gudhi::Cubical_complex::Embedded_cubical_complex_interface<>":
        Embedded_cubical_complex_interface(vector[unsigned] dimensions, vector[double] top_dimensional_cells) nogil
        vector[double] compute_hybrid_transform(string kernel, vector[vector[double]] directions_list, int num_jobs) nogil

cdef extern from "hybrid_transform_kernels.h" namespace "Gudhi":
    cdef cppclass Hybrid_transform_kernels "Gudhi::Cubical_complex::Hybrid_transform_kernels":
        Hybrid_transform_kernels() nogil
        double kernel_exp(double x) nogil
        double kernel_sin(double x) nogil
        double kernel_cos(double x) nogil

cdef class EmbeddedComplex:

    cdef Embedded_cubical_complex_interface * this_ptr
    cdef Hybrid_transform_kernels * kernel_obj
    cdef dict kernels_dict
    
    # Fake constructor that does nothing but documenting the constructor
    def __init__(self, sizes, data):
        """EmbeddedComplex constructor.
        :param sizes: The sizes of the embedded complex.
        :type sizes: list of int
        :param data: The filtration values of the top dimensional cells of the cubical complex.
        :type data: list of double
        """
    
    # The real cython constructor
    def __cinit__(self, dimensions=None, top_dimensional_cells=None):
        if((dimension is not None) and (top_dimensional_cells is not None)):
            with nogil:
                self.thisptr = new Embedded_cubical_complex_base_interface(dimensions, top_dimensional_cells)
                self.kernel_obj = new Hybrid_transform_kernels()

    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr

    def compute_hybrid_transform(self, string kernel = None, vector[double] direction = None, int num = -1):
        if((kernel is not None) and (direction is not None)):
            with gil:
                #Finding the kernel function
            with nogil:
                #Computing the transform
# distutils: language = c++
from libcpp.vector cimport vector

cdef extern from "embedded_complex_release.hpp" namespace "Gudhi":
    cdef cppclass Embedded_cubical_complex "Embedded_cubical_complex<Gudhi::cubical_complex::Bitmap_cubical_complex_base<double>>":
        Embedded_cubical_complex(vector[int] dimensions, vector[double] top_dimensional_cells)
        #double compute_hybrid_transform(double (*kernel)(double), vector[double] e)

cdef class EmbeddedComplex:

    cdef Embedded_cubical_complex * this_ptr
    
    def __init__(self, sizes, data):
        """EmbeddedComplex :
        :The constructors takes two arguments :
        :sizes: the sizes of the embedded complex
        :data: the filtration values of the top dimensional cells of the cubical complex
        """
    def __cinit__(self, sizes, data):
        cdef vector[int] v_sizes
        cdef vector[double] v_data
        self.this_ptr = new Embedded_cubical_complex(sizes,data)
from cython cimport numeric
from libcpp.vector cimport vector
from libcpp.utility cimport pair

cimport cgudhi

# SimplexTree python interface
cdef class SimplexTree:
    cdef cgudhi.Simplex_tree_interface[cgudhi.Simplex_tree_options_full_featured] *thisptr
    def __cinit__(self, points=None, max_dimension=3, max_edge_length=float('inf')):
        self.thisptr = new cgudhi.Simplex_tree_interface[cgudhi.Simplex_tree_options_full_featured]()
        # Constructor from graph expansion
        if points is not None:
            self.thisptr.graph_expansion(points,max_dimension,max_edge_length)
    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr
    def get_filtration(self):
        return self.thisptr.filtration()
    def filtration(self, simplex):
        return self.thisptr.simplex_filtration(simplex)
    def set_filtration(self, filtration):
        self.thisptr.set_filtration(<double>filtration)
    def initialize_filtration(self):
        self.thisptr.initialize_filtration()
    def num_vertices(self):
        return self.thisptr.num_vertices()
    def num_simplices(self):
        return self.thisptr.num_simplices()
    def dimension(self):
        return self.thisptr.dimension()
    def set_dimension(self, dim):
        self.thisptr.set_dimension(<int>dim)
    def find(self, simplex):
        cdef vector[int] complex
        for i in simplex:
          complex.push_back(i)
        return self.thisptr.find_simplex(complex)
    def insert(self, simplex, filtration = 0.0):
        return self.thisptr.insert_simplex(simplex, <double>filtration)
    def insert_with_subfaces(self, simplex, filtration = 0.0):
        return self.thisptr.insert_simplex_and_subfaces(simplex, <double>filtration)
    def get_filtered_tree(self):
        return self.thisptr.get_filtered_tree()
    def get_skeleton_tree(self, dim):
        return self.thisptr.get_skeleton_tree(<int>dim)
    def get_star_tree(self, simplex):
        return self.thisptr.get_star_tree(simplex)
    def get_coface_tree(self, simplex, dim):
        return self.thisptr.get_coface_tree(simplex, <int>dim)


cdef class MiniSimplexTree:
    cdef cgudhi.Simplex_tree_interface[cgudhi.Simplex_tree_options_mini] *thisptr
    def __cinit__(self):
        self.thisptr = new cgudhi.Simplex_tree_interface[cgudhi.Simplex_tree_options_mini]()
    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr
    def num_vertices(self):
        return self.thisptr.num_vertices()
    def num_simplices(self):
        return self.thisptr.num_simplices()
    def dimension(self):
        return self.thisptr.dimension()
    def set_dimension(self, dim):
        self.thisptr.set_dimension(<int>dim)
    def find(self, simplex):
        cdef vector[int] complex
        for i in simplex:
          complex.push_back(i)
        return self.thisptr.find_simplex(complex)

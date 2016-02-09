from cython cimport numeric
cimport cgudhi

cdef class SimplexTree:
    cdef cgudhi.Simplex_tree[cgudhi.Simplex_tree_options_full_featured] *thisptr
    def __cinit__(self):
        self.thisptr = new cgudhi.Simplex_tree[cgudhi.Simplex_tree_options_full_featured]()
    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr

cdef class MiniSimplexTree:
    cdef cgudhi.Simplex_tree[cgudhi.Simplex_tree_options_fast_persistence] *thisptr
    def __cinit__(self):
        self.thisptr = new cgudhi.Simplex_tree[cgudhi.Simplex_tree_options_fast_persistence]()
    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr

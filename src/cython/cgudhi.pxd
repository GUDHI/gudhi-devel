
cdef extern from "include/gudhi/Simplex_tree.h" namespace "Gudhi":
    cdef cppclass Simplex_tree_options_full_featured:
        pass
    cdef Simplex_tree_options_fast_persistence:
        pass
    cdef cppclass Simplex_tree[T]:
        Simplex_tree()


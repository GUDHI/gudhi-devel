from libcpp.vector cimport vector
from libcpp.utility cimport pair

cdef extern from "Simplex_tree_interface.h" namespace "Gudhi":
    cdef cppclass Simplex_tree_options_full_featured:
        pass
    cdef cppclass Simplex_tree_options_fast_persistence:
        pass
    cdef cppclass Simplex_tree_options_mini:
        pass
    cdef cppclass Simplex_tree_interface[T]:
        Simplex_tree()
        double filtration()
        double simplex_filtration(vector[int] simplex)
        void set_filtration(double filtration)
        void initialize_filtration()
        int num_vertices()
        int num_simplices()
        void set_dimension(int dimension)
        int dimension()
        bint find_simplex(vector[int] simplex)
        bint insert_simplex(vector[int] simplex, double filtration)
        bint insert_simplex_and_subfaces(vector[int] simplex, double filtration)
        vector[pair[vector[int], double]] get_filtered_tree()
        vector[pair[vector[int], double]] get_skeleton_tree(int dimension)
        vector[pair[vector[int], double]] get_star_tree(vector[int] simplex)
        vector[pair[vector[int], double]] get_coface_tree(vector[int] simplex, int dimension)
        void graph_expansion(vector[vector[double]] points,int max_dimension,double max_edge_length)



#!/usr/bin/env python

import gudhi

print("#####################################################################")
print("MiniSimplexTree creation from insertion")

""" Complex to build.
     1   3
     o---o
    /X\ /
   o---o   o
   2   0   4
"""

triangle012 = [0, 1, 2]
edge03 = [0, 3]
edge13 = [1, 3]
vertex4 = [4]
mini_st = gudhi.MiniSimplexTree()
mini_st.insert(triangle012)
mini_st.insert(edge03)
mini_st.insert(edge13)
mini_st.insert(vertex4)

# FIXME: Remove this line
mini_st.set_dimension(2)

# initialize_filtration required before plain_homology
mini_st.initialize_filtration()

print("plain_homology(2)=", mini_st.plain_homology(2))

edge02 = [0, 2]
if mini_st.find(edge02):
    # Only coface is 012
    print("coface(edge02,1)=", mini_st.get_coface_tree(edge02, 1))

if mini_st.get_coface_tree(triangle012, 1) == []:
    # Precondition: Check the simplex has no coface before removing it.
    mini_st.remove_maximal_simplex(triangle012)

# initialize_filtration required after removing
mini_st.initialize_filtration()

print("filtered_tree after triangle012 removal =", mini_st.get_filtered_tree())

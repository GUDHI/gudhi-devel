#!/usr/bin/env python

import gudhi

st = gudhi.SimplexTree()

print("#####################################################################")
print("SimplexTree creation from insertion")
if st.insert([0, 1]):
    print("Inserted !!")
else:
    print("Not inserted...")

if st.find([0, 1]):
    print("Found !!")
else:
    print("Not found...")

if st.insert([0, 1, 2], filtration=4.0):
    print("Inserted !!")
else:
    print("Not inserted...")

# FIXME: Remove this line
st.set_dimension(3)
print("dimension=", st.dimension())

st.set_filtration(4.0)
st.initialize_filtration()
print("filtration=", st.get_filtration())
print("filtration[1, 2]=", st.filtration([1, 2]))
print("filtration[4, 2]=", st.filtration([4, 2]))

print("num_simplices=", st.num_simplices())
print("num_vertices=", st.num_vertices())

print("skeleton_tree[2]=", st.get_skeleton_tree(2))
print("skeleton_tree[1]=", st.get_skeleton_tree(1))
print("skeleton_tree[0]=", st.get_skeleton_tree(0))

print("persistence(2)=", st.persistence(homology_coeff_field=2,
                                        min_persistence=0))

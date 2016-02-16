#!/usr/bin/env python

import gudhi

st = gudhi.SimplexTree()

print("#######################################################################")
print("SimplexTree creation from insertion")
if st.insert([0,1]):
  print("Inserted !!")
else:
  print("Not inserted...")

if st.find([0,1]):
  print("Found !!")
else:
  print("Not found...")

if st.insert_with_subfaces([0,1,2], filtration=4.0):
  print("Inserted !!")
else:
  print("Not inserted...")

# FIXME: Remove this line
st.set_dimension(3)
print("dimension=", st.dimension())

st.set_filtration(4.0)
st.initialize_filtration()
print("filtration=", st.get_filtration())
print("filtration[1,2]=", st.filtration([1,2]))
print("filtration[4,2]=", st.filtration([4,2]))

print("num_simplices=", st.num_simplices())
print("num_vertices=", st.num_vertices())

print("skeleton_tree[2]=", st.get_skeleton_tree(2))
print("skeleton_tree[1]=", st.get_skeleton_tree(1))
print("skeleton_tree[0]=", st.get_skeleton_tree(0))

print("#######################################################################")
print("SimplexTree creation from graph expansion")
st_from_graph_expansion = gudhi.SimplexTree(points=[[0,0],[1,0],[0,1],[1,1]],max_dimension=1,max_edge_length=42)

print("filtered_tree=", st_from_graph_expansion.get_filtered_tree())
print("star([0])=", st_from_graph_expansion.get_star_tree([0]))
print("coface([0],1)=", st_from_graph_expansion.get_coface_tree([0], 1))


print("#######################################################################")
print("MiniSimplexTree creation from insertion")
triangle012 = [0, 1, 2]
edge03 = [0, 3]
mini_st = gudhi.MiniSimplexTree()
mini_st.insert_with_subfaces(triangle012)
mini_st.insert_with_subfaces(edge03)
# FIXME: Remove this line
mini_st.set_dimension(2);

edge02 = [0, 2]
if mini_st.find(edge02):
  # Only coface is 012
  print("coface(edge02,1)=", mini_st.get_coface_tree(edge02, 1))


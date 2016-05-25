#!/usr/bin/env python

import gudhi

print("#####################################################################")
print("AlphaComplex creation from points")
alpha_complex = gudhi.AlphaComplex(points=[[0, 0], [1, 0], [0, 1], [1, 1]],
                                   max_alpha_square=60.0)

if alpha_complex.find([0, 1]):
    print("[0, 1] Found !!")
else:
    print("[0, 1] Not found...")

if alpha_complex.find([4]):
    print("[4] Found !!")
else:
    print("[4] Not found...")

if alpha_complex.insert([0, 1, 2], filtration=4.0):
    print("[0, 1, 2] Inserted !!")
else:
    print("[0, 1, 2] Not inserted...")

if alpha_complex.insert([0, 1, 4], filtration=4.0):
    print("[0, 1, 4] Inserted !!")
else:
    print("[0, 1, 4] Not inserted...")

if alpha_complex.find([4]):
    print("[4] Found !!")
else:
    print("[4] Not found...")

print("dimension=", alpha_complex.dimension())
print("filtered_tree=", alpha_complex.get_filtered_tree())
print("star([0])=", alpha_complex.get_star_tree([0]))
print("coface([0], 1)=", alpha_complex.get_coface_tree([0], 1))

print("point[0]=", alpha_complex.get_point(0))
print("point[5]=", alpha_complex.get_point(5))

alpha_complex.initialize_filtration()
print("persistence(2)=", alpha_complex.persistence(homology_coeff_field=2,
                                                   min_persistence=0))

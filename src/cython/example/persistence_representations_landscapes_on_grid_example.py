#!/usr/bin/env python

import gudhi


"""This file is part of the Gudhi Library. The Gudhi library
   (Geometric Understanding in Higher Dimensions) is a generic C++
   library for computational topology.

   Author(s):       Pawel Dlotko

   Copyright (C) 2017 Swansea University

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

__author__ = "Pawel Dlotko"
__copyright__ = "Copyright (C) 2017 Swansea University"
__license__ = "GPL v3"

print("#####################################################################")
print("Persistence representations landscapes on a grid example")

persistence1 = [(1, 2),(6, 8),(0, 4),(3, 8)]
persistence2 = [(2, 9),(1, 6),(3, 5),(6, 10)]

#create two persistence landscapes based on persistence1 and persistence2:
l1 = PersistenceLandscapeOnGrid(persistence1, 0, 11, 20)
l2 = PersistenceLandscapeOnGrid(persistence2, 0, 11, 20)

#This is how to compute integral of landscapes:
print "Integral of the first landscape : " , l1.compute_integral_of_landscape() 
print "Integral of the second landscape : " , l2.compute_integral_of_landscape() 

#here are the maxima of the functions:
print "Maximum of l1 : " , l1.compute_maximum() 
print "Maximum of l2 : " , l2.compute_maximum() 

#here are the norms of landscapes:
print "L^1 Norm of l1 : " , l1.compute_norm_of_landscape(1.) 
print "L^1 Norm of l2 : " , l2.compute_norm_of_landscape(1.) 

#here is the average of landscapes:
average = PersistenceLandscapeOnGrid();
average.compute_average(to_average=[l1, l2]);


#here is the distance of landscapes:
print "Distance : " , l1.distance(l2) 

#here is the scalar product of landscapes:
print "Scalar product : " , l1.compute_scalar_product(l2) 


















persistence1 = [(1,2),(6,8),(0,4),(3,8)]
persistence2 = [(2,9),(1,6),(3,5),(6,10)]


#create two persistence landscapes based on persistence1 and persistence2:
l1 = gudhi.PersistenceLandscapes(vector_of_intervals=persistence1, dimension=3)
l2 = gudhi.PersistenceLandscapes(vector_of_intervals=persistence2)

#This is how to compute integral of landscapes:
print "Integral of the first landscape : ", l1.compute_integral_of_landscape()
print "Integral of the second landscape : ", l2.compute_integral_of_landscape()

#here are the maxima of the functions:
print "Maximum of l1 : ", l1.compute_maximum()
print "Maximum of l2 : ", l2.compute_maximum()

#here are the norms of landscapes:
print "L^1 Norm of l1 : ", l1.compute_norm_of_landscape(1.)
print "L^1 Norm of l2 : ", l2.compute_norm_of_landscape(1.)

#here is the average of landscapes:
average = gudhi.PersistenceLandscapes()
average.compute_average(to_average=[l1, l2])

#here is the distance of landscapes:
print "Distance : ", l1.distance(average,1)

#here is the scalar product of landscapes:
print "Scalar product : ", l1.compute_scalar_product(l2)

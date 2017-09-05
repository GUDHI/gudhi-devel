#!/usr/bin/env python

import gudhi
from libcpp.vector cimport vector
from libcpp.utility cimport pair


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
print("Persistence representations landscapes example")

persistence1 = vector[pair[double, double]]
persistence2 = vector[pair[double, double]]

persistence1.push_back(make_pair(1, 2))
persistence1.push_back(make_pair(6, 8))
persistence1.push_back(make_pair(0, 4))
persistence1.push_back(make_pair(3, 8))

persistence2.push_back(make_pair(2, 9))
persistence2.push_back(make_pair(1, 6))
persistence2.push_back(make_pair(3, 5))
persistence2.push_back(make_pair(6, 10))

#create two persistence landscapes based on persistence1 and persistence2:
l1 = gudhi.Persistence_landscape(persistence1)
l2 = gudhi.Persistence_landscape(persistence2)

#This is how to compute integral of landscapes:
print "Integral of the first landscape : ", l1.compute_integral_of_landscape()
print "Integral of the second landscape : ", l2.compute_integral_of_landscape()

#Arithmetic operations on landscapes:
sum_ = l1 + l2;  

#here are the maxima of the functions:
print "Maximum of l1 : ", l1.compute_maximum()
print "Maximum of l2 : ", l2.compute_maximum()

#here are the norms of landscapes:
print "L^1 Norm of l1 : ", l1.compute_norm_of_landscape(1.)
print "L^1 Norm of l2 : ", l2.compute_norm_of_landscape(1.)

#here is the average of landscapes:
average = Persistence_landscape
average.compute_average({l1, l2});  

#here is the distance of landscapes:
print "Distance : ", l1.distance(l2)

#here is the scalar product of landscapes:
print "Scalar product : ", l1.compute_scalar_product(l2)

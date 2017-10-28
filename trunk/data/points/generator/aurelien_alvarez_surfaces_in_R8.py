#   This file is part of the Gudhi Library. The Gudhi library
#   (Geometric Understanding in Higher Dimensions) is a generic C++
#   library for computational topology.
#  
#   Author(s):       Aurélien Alvarez
#  
#   Copyright (C) 2016  Université d'Orléans (France)
#  
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#  
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#  
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import random
from math import factorial

I = complex(0,1)

#################################################
#################################################

#Surface réelle d'équation x.conj(y)^d + y.conj(z)^d + z.conj(x)^d = 0 dans P2(C)
#Équation affine (z=1) multipliée par sa conjuguée (d = 2) : x.conj(x)^2.y^4 + 2x^3.conj(x).y^2 + y + conj(x)^2 + x^5 = 0
def equationAffineSurfaceReelle(x):
    polynome = [0]*(degre**2+1)
    for k in range(degre+1):
        polynome[k*degre] = (-1)**degre*x*factorial(degre)/(factorial(k)*factorial(degre-k))*x**(k*degre)*np.conjugate(x)**(degre-k)
    polynome[-2] += 1
    polynome[-1] += np.conjugate(x)**degre
    return polynome

#################################################
#################################################

def calculRacines(equation,nombrePoints,module_x):
    racines = [[1,0,0],[0,1,0],[0,0,1]]
    for _ in range(nombrePoints):
        x = module_x*(2*random.random()-1+I*(2*random.random()-1))
        fool = [[[x,y,1],[y,1,x],[1,x,y]] for y in np.roots(equation(x)) if abs(x*np.conjugate(y)**degre+y+np.conjugate(x)**degre) < 0.0001]
        for bar in fool:
            racines += bar
    return racines

#################################################
#################################################

def plongementDansR8(pointDansCP2):
    z0 = pointDansCP2[0]
    z1 = pointDansCP2[1]
    z2 = pointDansCP2[2]
    a = z0*np.conjugate(z0)
    b = z1*np.conjugate(z1)
    c = z2*np.conjugate(z2)
    normeCarree = a+b+c
    a = a/normeCarree
    b = b/normeCarree
    u = z0*np.conjugate(z1)/normeCarree
    v = z0*np.conjugate(z2)/normeCarree
    w = z1*np.conjugate(z2)/normeCarree
    return [a.real,b.real,u.real,u.imag,v.real,v.imag,w.real,w.imag]

def plongementListeDansR8(listePointsDansCP2):
    listePointsDansR8 = []
    for point in listePointsDansCP2:
        listePointsDansR8 += [plongementDansR8(point)]
    return listePointsDansR8

#################################################
#################################################

degre = 3
nombrePoints = 10**4
module_x = 10

with open("surface.txt","w") as fichier: 
    bar = calculRacines(equationAffineSurfaceReelle,nombrePoints,module_x)
    listePoints = plongementListeDansR8(bar)
    fichier.write(str(len(bar)) + "\n")
    for point in listePoints:
        fichier.write(str(point[0]) + " " + str(point[1]) + " " + str(point[2]) + " " + str(point[3]) + " " + str(point[4]) + " " + str(point[5]) + " " + str(point[6]) + " " + str(point[7]) + "\n")


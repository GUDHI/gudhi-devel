#!/usr/bin/env python

from os import listdir

"""This file is part of the Gudhi Library. The Gudhi library
   (Geometric Understanding in Higher Dimensions) is a generic C++
   library for computational topology.

   Author(s):       Vincent Rouvreau

   Copyright (C) 2017 INRIA

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

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2017 INRIA"
__license__ = "GPL v3"

"""
generate_examples.py generates examples.inc to be included in examples.rst.
Refer to Makefile and make.bat to see if it is correctly launched.
"""

output_file = open('examples.inc','w')

output_file.write('.. only:: builder_html\n\n')

for file in listdir('../example/'):
    output_file.write("    * :download:`" + file + " <../example/" + file + ">`\n")

output_file.close()

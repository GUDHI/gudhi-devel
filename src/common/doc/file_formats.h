/*    This file is part of the Gudhi Library. The Gudhi library
*    (Geometric Understanding in Higher Dimensions) is a generic C++
*    library for computational topology.
*
*    Author(s):       Clément Jamin
*
*    Copyright (C) 2017 Inria
*
*    This program is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    This program is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef DOC_COMMON_FILE_FORMAT_H_
#define DOC_COMMON_FILE_FORMAT_H_

namespace Gudhi {

/*! \page fileformats File formats

 \tableofcontents

 \section FileFormatsPers Persistence Diagram

 Such a file, whose extension is usually `.pers`, contains a list of persistence intervals.<br>
 Lines starting with `#` are ignored (comments).<br>
 Other lines might contain 2, 3 or 4 values (the number of values on each line must be the same for all lines):
 \verbatim
   [[field] dimension] birth death
 \endverbatim

 Here is a simple sample file:
 \verbatim
   # Persistence diagram example
   2 2.7 3.7
   2 9.6 14.
   # Some comments
   3 34.2 34.974
   4 3. inf
 \endverbatim

 Other sample files can be found in the `data/persistence_diagram` folder.

 Such files can be generated with `Gudhi::persistent_cohomology::Persistent_cohomology::output_diagram()` and read with
 `Gudhi::read_persistence_intervals_and_dimension()`, `Gudhi::read_persistence_intervals_grouped_by_dimension()` or
 `Gudhi::read_persistence_intervals_in_dimension()`.
 

 \section FileFormatsIsoCuboid Iso-cuboid

 Such a file describes an iso-oriented cuboid with diagonal opposite vertices (min_x, min_y, min_z,...) and (max_x, max_y, max_z, ...). The format is:<br>
 \verbatim
   min_x min_y [min_z ...]
   max_x max_y [max_z ...]
 \endverbatim

 Here is a simple sample file in the 3D case:
 \verbatim
   -1. -1. -1.
   1. 1. 1.
 \endverbatim


 \section FileFormatsPerseus Perseus

 This file format is the format used by the Perseus software
 (http://www.sas.upenn.edu/~vnanda/perseus/) by Vidit Nanda.
 The first line contains a number d begin the dimension of the
 bitmap (2 in the example below). Next d lines are the numbers of top dimensional cubes in each dimensions (3 and 3
 in the example below). Next, in lexicographical order, the filtration of top dimensional cubes is given (1 4 6 8
 20 4 7 6 5 in the example below).
 
 \image html "exampleBitmap.png" "Example of a input data."
 
 The input file for the following complex is:
 \verbatim
 2
 3
 3
 1
 4
 6
 8
 20
 4
 7
 6
 5
 \endverbatim

 To indicate periodic boundary conditions in a
 given direction, then number of top dimensional cells in this direction have to be multiplied by -1. For instance:

 \verbatim
 2
 -3
 3
 1
 4
 6
 8
 20
 4
 7
 6
 5
 \endverbatim

 Indicate that we have imposed periodic boundary conditions in the direction x, but not in the direction y.

 Other sample files can be found in the `data/bitmap` folder.
 
 Please note that unlike in Perseus format the filtration on the maximal cubes can be any double precision number. 
 Consequently one cannot mark the cubes that are not present with -1's. To do that please set their filtration value to +Inf.

*/
}  // namespace Gudhi

#endif  // DOC_COMMON_FILE_FORMAT_H_

/*    This file is part of the Gudhi Library. The Gudhi library
*    (Geometric Understanding in Higher Dimensions) is a generic C++
*    library for computational topology.
*
*    Author(s):       Clément Jamin
*
*    Copyright (C) 2017  INRIA
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
 \code{.unparsed}
   [[field] dimension] birth death
 \endcode

 Here is a simple sample file:
 \code{.unparsed}
   # Beautiful persistence diagram
   2 2.7 3.7
   2 9.6 14.
   3 34.2 34.974
   4 3. inf
 \endcode

 Other sample files can be found in the `data/persistence_diagram` folder.
*/
}  // namespace Gudhi

#endif  // DOC_COMMON_FILE_FORMAT_H_

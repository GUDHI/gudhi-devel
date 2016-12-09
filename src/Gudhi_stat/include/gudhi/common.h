/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2015  INRIA Sophia-Saclay (France)
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

#ifndef HAUSDORFF_DISTANCES_H
#define HAUSDORFF_DISTANCES_H
 //this file contain an implementation of some common procedures used in Gudhi_stat.




 //double epsi = std::numeric_limits<double>::epsilon();
double epsi = 0.000005;


/**
 *  A procedure used to compare doubles. Typically gien two doubles A and B, comparing A == B is not good idea. In this case, we use the procedure almostEqual with the epsi defined at
 *  the top of the file. Setting up the epsi give the user a tolerance on what should be consider equal.
**/
inline bool almost_equal( double a , double b )
{
    if ( fabs(a-b) < epsi )
        return true;
    return false;
}

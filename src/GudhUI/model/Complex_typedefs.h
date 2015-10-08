/* This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014  INRIA Sophia Antipolis-Mediterranee (France)
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
 * 
 */

#ifndef MODEL_COMPLEX_TYPEDEFS_H_
#define MODEL_COMPLEX_TYPEDEFS_H_

#include <gudhi/Skeleton_blocker/Skeleton_blocker_simple_geometric_traits.h>
#include <gudhi/Skeleton_blocker_geometric_complex.h>

#include <CGAL/Kernel_d/Point_d.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Cartesian_d.h>

struct Geometry_trait : public CGAL::Cartesian_d<double> {
  typedef CGAL::Cartesian<double>::Point_3 Point_3;
  typedef CGAL::Cartesian<double>::Vector_3 Vector_3;
  typedef CGAL::Point_d<Cartesian_d<double>> Point;
  typedef CGAL::Vector_d<Cartesian_d<double>> Vector;
};

typedef Geometry_trait::Point Point;

using namespace Gudhi;
using namespace Gudhi::skbl;

typedef Skeleton_blocker_simple_geometric_traits<Geometry_trait> Complex_geometric_traits;
typedef Skeleton_blocker_geometric_complex< Complex_geometric_traits > Complex;

#endif  // MODEL_COMPLEX_TYPEDEFS_H_

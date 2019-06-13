/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
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
using namespace Gudhi::skeleton_blocker;

typedef Skeleton_blocker_simple_geometric_traits<Geometry_trait> Complex_geometric_traits;
typedef Skeleton_blocker_geometric_complex< Complex_geometric_traits > Complex;

#endif  // MODEL_COMPLEX_TYPEDEFS_H_

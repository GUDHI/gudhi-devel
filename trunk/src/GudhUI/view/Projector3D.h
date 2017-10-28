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

#ifndef VIEW_PROJECTOR3D_H_
#define VIEW_PROJECTOR3D_H_

#include "model/Complex_typedefs.h"

class Projector3D {
 public:
  typedef Geometry_trait::Point Point;
  typedef Geometry_trait::Point_3 Point_3;

  virtual Point_3 operator()(const Point&) const = 0;

  virtual ~Projector3D() { }
};

#endif  // VIEW_PROJECTOR3D_H_

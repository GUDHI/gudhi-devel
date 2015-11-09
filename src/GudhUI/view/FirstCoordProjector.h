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

#ifndef VIEW_FIRSTCOORDPROJECTOR_H_
#define VIEW_FIRSTCOORDPROJECTOR_H_

#include "utils/UI_utils.h"
#include "Projector3D.h"

class FirstCoordProjector3D : public Projector3D {
  typedef Projector3D::Point Point;
  typedef Projector3D::Point_3 Point_3;

  Point_3 operator()(const Point& p) const {
    if (p.dimension() >= 3)
    return Point_3(p.x(), p.y(), p.z());
    else if  (p.dimension() >= 2)
    return Point_3(p.x(), p.y(), 0.0);

  }
};

#endif  // VIEW_FIRSTCOORDPROJECTOR_H_

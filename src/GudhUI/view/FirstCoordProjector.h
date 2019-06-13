/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
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
    else
      return Point_3(0.0, 0.0, 0.0);
  }
};

#endif  // VIEW_FIRSTCOORDPROJECTOR_H_

/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
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

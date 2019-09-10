/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef VIEW_COLOR_H_
#define VIEW_COLOR_H_

struct Color {
  double r;
  double g;
  double b;

  Color(double r_, double g_, double b_) : r(r_), g(g_), b(b_) { }
};

#endif  // VIEW_COLOR_H_

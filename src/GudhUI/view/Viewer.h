/* This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
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

#ifndef VIEW_VIEWER_H_
#define VIEW_VIEWER_H_

#include <QGLViewer/qglviewer.h>

#include <vector>

#include "View_parameter.h"
#include "model/Complex_typedefs.h"
#include "Color.h"
#include "Viewer_instructor.h"

class Viewer_instructor;

class Viewer : public QGLViewer {
  Q_OBJECT

  Viewer_instructor * instructor;

  /**
   * light angles
   */
  double theta, phi;
  typedef Complex::GT Gudhi_kernel;
  typedef Gudhi_kernel::Point_3 Point_3;

 public:
  Viewer(QWidget* parent);

  void set_instructor(Viewer_instructor* instructor_);

  void show_entire_scene();

  void draw();

  void set_bounding_box(const Point_3 & lower_left, const Point_3 & upper_right);

  void update_GL();

  void init_scene();

  void init_light();

  void set_light();

  void set_light_direction(double theta, double phi);

  /**
   * set the light in the direction of the observer
   */
  void set_light_direction();

 protected:
  virtual void postSelection(const QPoint& point);

 public:
  ////////////////////////
  // draw
  ////////////////////////
  void set_size_point(double size_points);

  void set_color(const Color& color);

  void draw_point(const Point_3& p, const Color& color, double size_points);

  void begin_draw_points(double size, bool light = false);

  void draw_points(const Point_3 & point);

  void end_draw_points();

  void draw_edge(const Point_3 &a, const Point_3 &b, const Color& color, double size);

  void begin_draw_edges(double size, bool light = false);

  void draw_edges(const Point_3 &a, const Point_3 &b);

  void end_draw_edges();

  void begin_draw_triangles(double size, bool light, bool transparent = false);

  void draw_triangles(const Point_3& p1, const Point_3& p2, const Point_3& p3);

  // todo remove
  void draw_triangles(const std::vector<Point_3*>& points);

  void end_draw_triangles();


 signals:
  void click(const Point_3& position);
};

#endif  // VIEW_VIEWER_H_

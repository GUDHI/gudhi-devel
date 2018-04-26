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

#include "Viewer.h"
#include "utils/UI_utils.h"

Viewer::Viewer(QWidget* parent) : QGLViewer(QGLFormat(QGL::SampleBuffers), parent), instructor(0), theta(0), phi(0) { }

void Viewer::set_instructor(Viewer_instructor* instructor_) {
  instructor = instructor_;
}

void Viewer::show_entire_scene() {
  this->showEntireScene();
}

void Viewer::draw() {
  instructor->give_instructions();
}

void Viewer::set_bounding_box(const Point_3 & lower_left, const Point_3 & upper_right) {
  this->camera()->setSceneBoundingBox(qglviewer::Vec(lower_left[0], lower_left[1], lower_left[2]),
                                      qglviewer::Vec(upper_right[0], upper_right[1], upper_right[2]));
}

void Viewer::update_GL() {
  this->updateGL();
}

void Viewer::init_scene() {
  this->setBackgroundColor(Qt::white);
  ::glEnable(GL_LINE_SMOOTH);
  init_light();
}

void Viewer::init_light() {
  ::glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
}

void Viewer::set_light() {
  if (theta >= 0 && phi >= 0) {
    const GLfloat pos[4] = {static_cast<float> (sin(phi) * cos(theta)),
                            static_cast<float> (sin(phi) * sin(theta)),
                            static_cast<float> (cos(phi)), 0.};
    glLightfv(GL_LIGHT0, GL_POSITION, pos);
  }
}

void Viewer::set_light_direction(double theta_, double phi_) {
  theta = theta_;
  phi = phi_;
}

/**
 * set the light in the direction of the observer
 */
void Viewer::set_light_direction() {
  theta = -1;
  phi = -1;
}

void Viewer::postSelection(const QPoint& point) {
  bool found;

  auto vec = this->camera()->pointUnderPixel(point, found);

  if (found) {
    Point_3 position(vec[0], vec[1], vec[2]);
    emit(click(position));
  }
}

////////////////////////
// draw
////////////////////////

void Viewer::set_size_point(double size_points) {
  ::glPointSize(size_points);
}

void Viewer::draw_point(const Point_3& p, const Color& color, double size_points) {
  ::glColor3f(color.r, color.g, color.b);
  ::glDisable(GL_LIGHTING);
  ::glEnable(GL_POINT_SMOOTH);
  ::glPointSize(size_points);
  ::glBegin(GL_POINTS);
  ::glVertex3d(p.x(), p.y(), p.z());
  ::glEnd();
  ::glDisable(GL_POINT_SMOOTH);
}

void Viewer::begin_draw_points(double size, bool light) {
  light ? glEnable(GL_LIGHTING) : glDisable(GL_LIGHTING);
  ::glEnable(GL_POINT_SMOOTH);
  ::glPointSize(size);
  ::glBegin(GL_POINTS);
}

void Viewer::set_color(const Color& color) {
  ::glColor3f(color.r, color.g, color.b);
}

void Viewer::draw_points(const Point_3 & point) {
  ::glVertex3d(point.x(), point.y(), point.z());
}

void Viewer::end_draw_points() {
  ::glEnd();
  ::glDisable(GL_POINT_SMOOTH);
}

void Viewer::draw_edge(const Point_3 &a, const Point_3 &b, const Color& color, double size) {
  ::glColor3f(color.r, color.g, color.b);
  ::glPointSize(3.0);
  ::glLineWidth(size);
  ::glBegin(GL_LINES);
  ::glVertex3f(a.x(), a.y(), a.z());
  ::glVertex3f(b.x(), b.y(), b.z());
  ::glEnd();
}

void Viewer::begin_draw_edges(double size, bool light) {
  ::glLineWidth(size);
  ::glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  ::glEnable(GL_POLYGON_OFFSET_LINE);
  ::glPolygonOffset(3.0f, -3.0f);
  light ? glEnable(GL_LIGHTING) : glDisable(GL_LIGHTING);
  ::glBegin(GL_LINES);
}

void Viewer::draw_edges(const Point_3 &a, const Point_3 &b) {
  ::glVertex3f(a.x(), a.y(), a.z());
  ::glVertex3f(b.x(), b.y(), b.z());
}

void Viewer::end_draw_edges() {
  ::glEnd();
}

void Viewer::begin_draw_triangles(double size, bool light, bool transparent) {
  if (transparent) {
    ::glEnable(GL_BLEND);
    ::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  }
  ::glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  ::glEnable(GL_POLYGON_OFFSET_FILL);
  ::glPolygonOffset(3.0f, -3.0f);
  light ? glEnable(GL_LIGHTING) : glDisable(GL_LIGHTING);
  ::glBegin(GL_TRIANGLES);
}

void Viewer::draw_triangles(const Point_3& p1, const Point_3& p2, const Point_3& p3) {
  if (!CGAL::collinear(p1, p2, p3)) {
    auto triangle_normal = CGAL::unit_normal(p1, p2, p3);
    ::glNormal3d(triangle_normal.x(), triangle_normal.y(), triangle_normal.z());
    ::glVertex3d(p1.x(), p1.y(), p1.z());
    ::glVertex3d(p2.x(), p2.y(), p2.z());
    ::glVertex3d(p3.x(), p3.y(), p3.z());
  }
}

void Viewer::end_draw_triangles() {
  ::glEnd();
}

#include "Viewer.moc"

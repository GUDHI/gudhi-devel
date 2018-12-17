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

#include <utility>

#include "Viewer_instructor.h"
#include "utils/UI_utils.h"
#include "FirstCoordProjector.h"

Viewer_instructor::Viewer_instructor(QWidget* parent,
                                     Viewer* viewer,
                                     const Complex& mesh
                                     )
    : viewer_(viewer), mesh_(mesh), projector_(new FirstCoordProjector3D()) {
  viewer_->set_instructor(this);
}

void Viewer_instructor::initialize_bounding_box() {
  auto pair_bounding_box = compute_bounding_box_corners();
  viewer_->set_bounding_box(proj(pair_bounding_box.first), proj(pair_bounding_box.second));
  viewer_->init_scene();
}

std::pair<Complex::Point, Complex::Point> Viewer_instructor::compute_bounding_box_corners() {
  if (mesh_.empty()) {
    return std::make_pair(Point(-1, -1, -1, 1), Point(1, 1, 1, 1));
  } else {
    double x_min = 1e10;
    double y_min = 1e10;
    double z_min = 1e10;
    double x_max = -1e10;
    double y_max = -1e10;
    double z_max = -1e10;
    for (auto vi : mesh_.vertex_range()) {
      auto pt = proj(mesh_.point(vi));
      x_min = (std::min)(x_min, pt.x());
      y_min = (std::min)(y_min, pt.y());
      z_min = (std::min)(z_min, pt.z());

      x_max = (std::max)(x_max, pt.x());
      y_max = (std::max)(y_max, pt.y());
      z_max = (std::max)(z_max, pt.z());
    }
    return std::make_pair(Point(x_min, y_min, z_min, 1.),
                          Point(x_max, y_max, z_max, 1.));
  }
}

void Viewer_instructor::show_entire_scene() {
  viewer_->show_entire_scene();
}

const qglviewer::Camera* Viewer_instructor::camera() const {
  return viewer_->camera();
}

int
Viewer_instructor::width() const {
  return viewer_->width();
}

int
Viewer_instructor::height() const {
  return viewer_->height();
}

/**
 * to change display parameters
 */
View_parameter& Viewer_instructor::view_params() {
  return view_params_;
}

void
Viewer_instructor::give_instructions() {
  if (view_params_.relative_light)
    viewer_->set_light_direction();
  else
    viewer_->set_light_direction(view_params_.theta, view_params_.phi);
  viewer_->set_light();

  if (view_params_.edge_mode) draw_edges();
  if (view_params_.triangle_mode) draw_triangles();
  if (view_params_.vertex_mode) draw_points();
}

void Viewer_instructor::draw_edges() {
  viewer_->begin_draw_edges(view_params_.size_edges, false);

  for (auto edge : mesh_.edge_range()) {
    set_color_edge(edge);
    const Point& a = mesh_.point(mesh_.first_vertex(edge));
    const Point& b = mesh_.point(mesh_.second_vertex(edge));
    viewer_->draw_edges(proj(a), proj(b));
  }

  viewer_->end_draw_edges();
}

void Viewer_instructor::draw_triangles() {
  const double size_triangles = 1.0;
  viewer_->begin_draw_triangles(size_triangles, view_params_.light);

  for (const auto& fit : mesh_.triangle_range()) {
    set_color_triangle(fit);
    if (view_params_.triangle_mode) {
      auto fit_it = fit.begin();
      const Point& p1 = mesh_.point(*fit_it);
      const Point& p2 = mesh_.point(*(++fit_it));
      const Point& p3 = mesh_.point(*(++fit_it));
      viewer_->draw_triangles(proj(p1), proj(p2), proj(p3));
    }
  }
  viewer_->end_draw_triangles();
}

void Viewer_instructor::draw_points() {
  viewer_->begin_draw_points(view_params_.size_vertices);
  for (auto vi : mesh_.vertex_range()) {
    viewer_->set_size_point(view_params_.size_vertices);
    set_color_vertex(vi);
    viewer_->draw_points(proj(mesh_.point(vi)));
  }
  viewer_->end_draw_points();
}

void Viewer_instructor::draw_edge(const Point&, const Point&) { }

void Viewer_instructor::draw_point(const Point&) { }

/**
 * set the right color of vertex/edge/triangle considering the view_params choice
 */
void Viewer_instructor::set_color_vertex(Vertex_handle vh) {
  viewer_->set_color(Color(view_params_.light_edges, view_params_.light_edges, view_params_.light_edges));
}

void Viewer_instructor::set_color_edge(Edge_handle eh) {
  viewer_->set_color(Color(view_params_.light_edges, view_params_.light_edges, view_params_.light_edges));
}

void Viewer_instructor::set_color_triangle(const Simplex& triangle) {
  viewer_->set_color(Color(view_params_.light_triangles, view_params_.light_triangles, view_params_.light_triangles));
}

Viewer_instructor::Point_3
Viewer_instructor::proj(const Point& p) const {
  return (*projector_)(p);
}

void Viewer_instructor::sceneChanged() {
  UIDBG("sceneChanged");
  viewer_->update_GL();
}

void Viewer_instructor::change_draw_vertices() {
  view_params_.change_vertex_mode();
}

void Viewer_instructor::change_draw_edges() {
  view_params_.change_edge_mode();
}

void Viewer_instructor::change_draw_triangles() {
  view_params_.change_triangle_mode();
}

void Viewer_instructor::change_light() {
  view_params_.light = !view_params_.light;
}

#include "Viewer_instructor.moc"

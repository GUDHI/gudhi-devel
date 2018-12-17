/*    This file is part of the Gudhi Library. The Gudhi library
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
 */

#include "MainWindow.h"

#include <QWidget>
#include <QInputDialog>
#include <QFileDialog>

#include <fenv.h>

#include "gui/Menu_k_nearest_neighbors.h"
#include "gui/Menu_uniform_neighbors.h"
#include "gui/Menu_edge_contraction.h"
#include "gui/Menu_persistence.h"

MainWindow::MainWindow(QWidget* parent) :
    menu_k_nearest_neighbors_(new Menu_k_nearest_neighbors(this)),
    menu_uniform_neighbors_(new Menu_uniform_neighbors(this)),
    menu_edge_contraction_(new Menu_edge_contraction(this, model_)),
    menu_persistence_(new Menu_persistence(this)) {
  // #ifndef NDEBUG // catch nan
  // feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  // #endif

  setupUi(this);

  viewer_instructor_ = new Viewer_instructor(
                                             this,
                                             this->viewer,
                                             model_.complex_);

  connectActions();

  update_view();
}

void
MainWindow::closeEvent(QCloseEvent *event) {
  exit(0);
}

void
MainWindow::connectActions() {
  QObject::connect(this->actionLoad_complex, SIGNAL(triggered()), this,
                   SLOT(off_file_open()));
  QObject::connect(this->actionLoad_points, SIGNAL(triggered()), this,
                   SLOT(off_points_open()));
  QObject::connect(this->actionSave_complex, SIGNAL(triggered()), this,
                   SLOT(off_file_save()));
  QObject::connect(this->actionSave_points, SIGNAL(triggered()), this,
                   SLOT(off_points_save()));

  QObject::connect(this->actionShow_graph_stats, SIGNAL(triggered()), this,
                   SLOT(show_graph_stats()));
  QObject::connect(this->actionShow_complex_stats, SIGNAL(triggered()), this,
                   SLOT(show_complex_stats()));
  QObject::connect(this->actionShow_complex_dimension, SIGNAL(triggered()), this,
                   SLOT(show_complex_dimension()));

  QObject::connect(this->actionUniform_proximity_graph, SIGNAL(triggered()), this,
                   SLOT(build_rips_menu()));
  QObject::connect(this->actionK_nearest_neighbors_graph, SIGNAL(triggered()), this,
                   SLOT(build_k_nearest_neighbors_menu()));


  QObject::connect(this->actionContract_edges, SIGNAL(triggered()), this,
                   SLOT(contract_edge_menu()));

  QObject::connect(this->actionCollapse_vertices, SIGNAL(triggered()), this,
                   SLOT(collapse_vertices()));

  QObject::connect(this->actionCollapse_edges, SIGNAL(triggered()), this,
                   SLOT(collapse_edges()));

  QObject::connect(this->actionNoise, SIGNAL(triggered()), this,
                   SLOT(uniform_noise()));
  QObject::connect(this->actionLloyd, SIGNAL(triggered()), this,
                   SLOT(lloyd()));


  // view
  QObject::connect(this->actionPoints, SIGNAL(triggered()), this->viewer_instructor_,
                   SLOT(change_draw_vertices()));
  QObject::connect(this->actionEdges, SIGNAL(triggered()), this->viewer_instructor_,
                   SLOT(change_draw_edges()));
  QObject::connect(this->actionTriangles, SIGNAL(triggered()), this->viewer_instructor_,
                   SLOT(change_draw_triangles()));

  // topology
  QObject::connect(this->actionShow_homology_group, SIGNAL(triggered()), this,
                   SLOT(show_homology_group()));
  QObject::connect(this->actionEuler_characteristic, SIGNAL(triggered()), this,
                   SLOT(show_euler_characteristic()));
  QObject::connect(this->actionPersistence, SIGNAL(triggered()), this,
                   SLOT(persistence_menu()));
  QObject::connect(this->actionEstimate_topological_changes, SIGNAL(triggered()), this,
                   SLOT(critical_points_menu()));
  QObject::connect(this->actionIs_manifold, SIGNAL(triggered()), this,
                   SLOT(is_manifold_menu()));


  QObject::connect(this, SIGNAL(sceneChanged()), this->viewer_instructor_,
                   SLOT(sceneChanged()));
}

void
MainWindow::init_view() const {
  viewer_instructor_->initialize_bounding_box();
  viewer_instructor_->show_entire_scene();
  update_view();
}

void
MainWindow::update_view() const {
  emit(sceneChanged());
}

/**
 * open a file chooser to choose an off to load
 */
void
MainWindow::off_file_open() {
  QString fileName = QFileDialog::getOpenFileName(this, tr("Open off File"),
                                                  "~/", tr("off files (*.off)"));
  if (!fileName.isEmpty()) {
    model_.off_file_open(fileName.toStdString());
    init_view();
  }
}

void
MainWindow::off_points_open() {
  QString fileName = QFileDialog::getOpenFileName(this, tr("Open points in a off file"),
                                                  "~/", tr("off files (*.off)"));
  if (!fileName.isEmpty()) {
    model_.off_points_open(fileName.toStdString());
    init_view();
  }
}

/**
 * open a file chooser to choose an off to save
 */
void
MainWindow::off_file_save() {
  QString fileName = QFileDialog::getOpenFileName(this, tr("Save to off File"),
                                                  "~/", tr("off files (*.off)"));
  if (!fileName.isEmpty()) {
    model_.off_file_save(fileName.toStdString());
  }
}

/**
 * open a file chooser to choose an off to save
 */
void
MainWindow::off_points_save() {
  QString fileName = QFileDialog::getOpenFileName(this, tr("Save to off File"),
                                                  "~/", tr("off files (*.off)"));
  if (!fileName.isEmpty()) {
    model_.off_points_save(fileName.toStdString());
  }
}

void
MainWindow::show_graph_stats() {
  model_.show_graph_stats();
}

void
MainWindow::show_complex_stats() {
  model_.show_complex_stats();
}

void
MainWindow::show_complex_dimension() {
  model_.show_complex_dimension();
}

void
MainWindow::build_rips_menu() {
  menu_uniform_neighbors_->show();
}

void
MainWindow::build_rips(double alpha) {
  model_.build_rips(alpha);
  update_view();
}

void
MainWindow::build_k_nearest_neighbors_menu() {
  menu_k_nearest_neighbors_->show();
}

void
MainWindow::build_k_nearest_neighbors(unsigned k) {
  model_.build_k_nearest_neighbors(k);
  update_view();
}

void
MainWindow::contract_edge_menu() {
  menu_edge_contraction_->show();
}

void
MainWindow::contract_edges(unsigned num_collapses) {
  std::cerr << "Collapse " << num_collapses << " vertices\n";
  model_.contract_edges(num_collapses);
  update_view();
}

void
MainWindow::collapse_edges() {
  model_.collapse_edges(model_.num_edges());
  update_view();
}

void
MainWindow::collapse_vertices() {
  std::cerr << "Collapse vertices edges\n";
  model_.collapse_vertices(0);
  update_view();
}

void
MainWindow::uniform_noise() {
  bool ok;
  double amplitude = QInputDialog::getDouble(this, tr("Uniform noise"),
                                             tr("Amplitude:"), 0, 0, 10000, 3, &ok);
  srand(time(NULL));
  if (ok)
    model_.uniform_noise(amplitude);
}

void
MainWindow::lloyd() {
  // todo 1 ask lloyd parameters
  model_.lloyd(0, 0);
  update_view();
}

void
MainWindow::show_homology_group() {
  model_.show_homology_group();
}

void
MainWindow::show_euler_characteristic() {
  model_.show_euler_characteristic();
}

void
MainWindow::persistence_menu() {
  menu_persistence_->show();
}

void
MainWindow::compute_persistence(int p, double threshold, int max_dim, double min_pers) {
  model_.show_persistence(p, threshold, max_dim, min_pers);
}

void
MainWindow::critical_points_menu() {
  bool ok;
  double max_length = QInputDialog::getDouble(this, tr("Maximal edge length for the Rips"),
                                              tr("Maximal edge length:"), 0, 0, 10000, 3, &ok);
  if (ok)
    model_.show_critical_points(max_length);
}

void
MainWindow::is_manifold_menu() {
  model_.show_is_manifold();
}


#include "MainWindow.moc"

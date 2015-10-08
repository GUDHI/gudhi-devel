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

#include "Menu_uniform_neighbors.h"

Menu_uniform_neighbors::Menu_uniform_neighbors(QMainWindow* parent_) :
    parent(parent_) {
  setupUi(this);
  connectActions(parent_);
}

void Menu_uniform_neighbors::connectActions(QMainWindow* parent) {
  QObject::connect(this->pushButtonCompute,
                   SIGNAL(clicked()),
                   this,
                   SLOT(send_compute_uniform_neighbors()));
  QObject::connect(this->doubleSpinBoxAlpha,
                   SIGNAL(valueChanged(double)),
                   this,
                   SLOT(update_alpha(double)));
  QObject::connect(this,
                   SIGNAL(compute_uniform_neighbors(double)),
                   parent,
                   SLOT(build_rips(double)));
}

void Menu_uniform_neighbors::send_compute_uniform_neighbors() {
  emit(compute_uniform_neighbors(doubleSpinBoxAlpha->value()));
}

void Menu_uniform_neighbors::accept() {
  send_compute_uniform_neighbors();
}

void Menu_uniform_neighbors::update_alpha(double alpha) {
  if (checkBoxAutoUpdate->isChecked())
    emit(compute_uniform_neighbors(doubleSpinBoxAlpha->value()));
}

#include "Menu_uniform_neighbors.moc"

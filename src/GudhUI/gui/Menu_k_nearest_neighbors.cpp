/*    This file is part of the Gudhi Library. The Gudhi library
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
 */

#include "Menu_k_nearest_neighbors.h"

Menu_k_nearest_neighbors::Menu_k_nearest_neighbors(QMainWindow* parent_)
    : parent(parent_) {
  setupUi(this);
  connectActions(parent_);
}

void Menu_k_nearest_neighbors::connectActions(QMainWindow* parent) {
  QObject::connect(this->pushButtonCompute,
                   SIGNAL(clicked()),
                   this,
                   SLOT(send_compute_k_nearest_neighbors()));
  QObject::connect(this->spinBoxK,
                   SIGNAL(valueChanged(int)),
                   this,
                   SLOT(update_k(int)));
  QObject::connect(this,
                   SIGNAL(compute_k_nearest_neighbors(unsigned)),
                   parent,
                   SLOT(build_k_nearest_neighbors(unsigned)));
}

void Menu_k_nearest_neighbors::send_compute_k_nearest_neighbors() {
  emit(compute_k_nearest_neighbors((unsigned) spinBoxK->value()));
}

void Menu_k_nearest_neighbors::accept() {
  send_compute_k_nearest_neighbors();
}

void Menu_k_nearest_neighbors::update_k(int new_k_value) {
  if (checkBoxAutoUpdate->isChecked())
    emit(compute_k_nearest_neighbors((unsigned) spinBoxK->value()));
}

#include "Menu_k_nearest_neighbors.moc"

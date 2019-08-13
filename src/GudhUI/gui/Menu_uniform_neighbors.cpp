/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
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

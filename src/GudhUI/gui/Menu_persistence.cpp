/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */


#include "Menu_persistence.h"

Menu_persistence::Menu_persistence(QMainWindow* parent_) : parent(parent_) {
  setupUi(this);
  connectActions(parent_);
}

void Menu_persistence::connectActions(QMainWindow* parent) {
  QObject::connect(this,
                   SIGNAL(compute_persistence(int, double, int, double)),
                   parent,
                   SLOT(compute_persistence(int, double, int, double)));
}

void Menu_persistence::send_compute_persistence() {
  emit(compute_persistence(p_spinBox->value(), threshold_doubleSpinBox->value(),
                           maxdimension_spinBox->value(), minpersistence_doubleSpinBox->value()));
}

void Menu_persistence::accept() {
  send_compute_persistence();
}

// void Menu_persistence::compute_persistence(int p_fied,double threshold,int dim_max,double min_persistence) {
//  if(checkBoxAutoUpdate->isChecked())
//    emit(compute_k_nearest_neighbors((unsigned)spinBoxK->value()));
// }

#include "Menu_persistence.moc"

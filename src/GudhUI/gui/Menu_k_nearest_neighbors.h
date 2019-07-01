/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef GUI_MENU_K_NEAREST_NEIGHBORS_H_
#define GUI_MENU_K_NEAREST_NEIGHBORS_H_

#include <QMainWindow>
#include "ui_KNearestNeighborsMenu.h"

class QWidget;

class Menu_k_nearest_neighbors : public QDialog, public Ui::KNearestNeighborsMenu {
  Q_OBJECT

 private:
  QMainWindow* parent;

 public:
  Menu_k_nearest_neighbors(QMainWindow* parent_);

  void connectActions(QMainWindow* parent);

 public slots:
  void send_compute_k_nearest_neighbors();
  void update_k(int k);
  void accept();

 signals:
  void compute_k_nearest_neighbors(unsigned k);
};

#endif  // GUI_MENU_K_NEAREST_NEIGHBORS_H_

/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef GUI_MENU_UNIFORM_NEIGHBORS_H_
#define GUI_MENU_UNIFORM_NEIGHBORS_H_

#include <QMainWindow>
#include "ui_UniformNeighborsMenu.h"

class Menu_uniform_neighbors : public QDialog, public Ui::UniformMenu {
  Q_OBJECT

 private:
  QMainWindow* parent;

 public:
  Menu_uniform_neighbors(QMainWindow* parent_);

  void connectActions(QMainWindow* parent);

 public slots:
  void send_compute_uniform_neighbors();
  void update_alpha(double alpha);
  void accept();

 signals:
  void compute_uniform_neighbors(double alpha);
};


#endif  // GUI_MENU_UNIFORM_NEIGHBORS_H_

/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef GUI_MENU_PERSISTENCE_H_
#define GUI_MENU_PERSISTENCE_H_

#include <QMainWindow>
#include "ui_PersistenceMenu.h"

class QWidget;

class Menu_persistence : public QDialog, public Ui::PersistenceMenu {
  Q_OBJECT

 private:
  QMainWindow* parent;

 public:
  Menu_persistence(QMainWindow* parent_);

  void connectActions(QMainWindow* parent);

 public slots:
  void send_compute_persistence();
  void accept();

 signals:
  void compute_persistence(int p_fied, double threshold, int dim_max, double min_persistence);
};

#endif  // GUI_MENU_PERSISTENCE_H_

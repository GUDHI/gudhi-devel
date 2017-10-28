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

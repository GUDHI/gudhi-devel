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

#ifndef GUI_MENU_EDGE_CONTRACTION_H_
#define GUI_MENU_EDGE_CONTRACTION_H_

#include "gui/MainWindow.h"
#include "ui_MenuEdgeContraction.h"

#include "model/Model.h"

class Menu_edge_contraction : public QDialog, public Ui::MenuEdgeContraction {
  Q_OBJECT

 private:
  MainWindow* parent_;
  const Model& model_;

  void update_slider_value();

 public:
  Menu_edge_contraction(MainWindow* parent, const Model& model);

  void connectActions(MainWindow* parent);

 private:
  unsigned num_vertices();
  unsigned num_collapses();

 public slots:
  void slider_value_changed(int new_slider_value);
  void update_gui_numbers();
  void update_gui_numbers(int gui_numbers);

  void send_contract_edges();

 signals:
  void contract_edges(unsigned num_collapses);
};

#endif  // GUI_MENU_EDGE_CONTRACTION_H_

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

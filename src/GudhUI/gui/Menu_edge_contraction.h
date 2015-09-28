/*
 * Menu_edge_contraction.h
 *
 *  Created on: Sep 11, 2014
 *      Author: dsalinas
 */

#ifndef MENU_EDGE_CONTRACTION_H_
#define MENU_EDGE_CONTRACTION_H_

// Workaround for moc-qt4 not parsing boost headers
#include <CGAL/config.h>

#include "gui/MainWindow.h"
#include "gui/ui_MenuEdgeContraction.h"

#include "model/Model.h"


class Menu_edge_contraction : public QDialog,public Ui::MenuEdgeContraction{
	Q_OBJECT
private:
	MainWindow* parent_;
	const Model& model_;


	void update_slider_value();
public:

	Menu_edge_contraction(MainWindow* parent,const Model& model);

	void connectActions(MainWindow* parent);


private:
	unsigned num_vertices();
	unsigned num_collapses();

	public slots:

	void slider_value_changed(int new_slider_value);
	void update_gui_numbers();
	void update_gui_numbers(int);

	void send_contract_edges();
	signals:

	void contract_edges(unsigned num_collapses);

};



#endif /* MENU_EDGE_CONTRACTION_H_ */

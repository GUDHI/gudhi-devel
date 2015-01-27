/*
 * Menu_uniform_neighbors.h
 *
 *  Created on: Sep 11, 2014
 *      Author: dsalinas
 */

#ifndef MENU_UNIFORM_NEIGHBORS_H_
#define MENU_UNIFORM_NEIGHBORS_H_

#include <QMainWindow>
#include "gui/ui_UniformNeighborsMenu.h"

class Menu_uniform_neighbors : public QDialog,public Ui::UniformMenu {
	Q_OBJECT
private:
	QMainWindow* parent;

public:

	Menu_uniform_neighbors(QMainWindow* parent_);

	void connectActions(QMainWindow* parent);

	public slots:
	void send_compute_uniform_neighbors();
	void update_alpha(double);
	void accept();

	signals:
	void compute_uniform_neighbors(double alpha);

};


#endif /* MENU_UNIFORM_NEIGHBORS_H_ */

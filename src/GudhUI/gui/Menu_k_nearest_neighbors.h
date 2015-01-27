/*
 * Menu_k_nearest_neighbors.h
 *
 *  Created on: Sep 11, 2014
 *      Author: dsalinas
 */

#ifndef MENU_K_NEAREST_NEIGHBORS_H_
#define MENU_K_NEAREST_NEIGHBORS_H_

#include <QMainWindow>
#include "gui/ui_KNearestNeighborsMenu.h"

class QWidget;


class Menu_k_nearest_neighbors : public QDialog,public Ui::KNearestNeighborsMenu{
	Q_OBJECT
private:
	QMainWindow* parent;

public:

	Menu_k_nearest_neighbors(QMainWindow* parent_);

	void connectActions(QMainWindow* parent);

	public slots:
	void send_compute_k_nearest_neighbors();
	void update_k(int);
	void accept();

	signals:
	void compute_k_nearest_neighbors(unsigned k);

};


#endif /* MENU_K_NEAREST_NEIGHBORS_H_ */

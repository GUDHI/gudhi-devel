/*
 * Menu_persistence.cpp
 *  Created on: Jan 27, 2015
 * This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014  INRIA Sophia Antipolis-Méditerranée (France)
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
 * 
 */


#include "Menu_persistence.h"

Menu_persistence::Menu_persistence(QMainWindow* parent_):parent(parent_)
{
	setupUi(this);
	connectActions(parent_);
}

void Menu_persistence::connectActions(QMainWindow* parent){
	QObject::connect(
			this,
			SIGNAL(compute_persistence(int,double,int,double)),
			parent,
			SLOT(compute_persistence(int,double,int,double))
	);
}

void Menu_persistence::send_compute_persistence(){
	emit(compute_persistence(p_spinBox->value(),threshold_doubleSpinBox->value(),
			maxdimension_spinBox->value(),minpersistence_doubleSpinBox->value()));
}

void Menu_persistence::accept(){
	send_compute_persistence();
}

//void Menu_persistence::compute_persistence(int p_fied,double threshold,int dim_max,double min_persistence){
////	if(checkBoxAutoUpdate->isChecked())
////		emit(compute_k_nearest_neighbors((unsigned)spinBoxK->value()));
//}


#include "Menu_persistence.moc"

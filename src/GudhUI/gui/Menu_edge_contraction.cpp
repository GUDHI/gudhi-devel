/*
 * Menu_edge_contraction.cpp
 *
 *  Created on: Sep 11, 2014
 *      Author: dsalinas
 */

#ifndef MENU_EDGE_CONTRACTION_CPP_
#define MENU_EDGE_CONTRACTION_CPP_


#include "Menu_edge_contraction.h"

Menu_edge_contraction::Menu_edge_contraction(MainWindow* parent,const Model& model):
parent_(parent),model_(model)
{
	setupUi(this);
	connectActions(parent_);
}

void Menu_edge_contraction::connectActions(MainWindow* parent)
{
	QObject::connect(
			this->horizontalSlider,
			SIGNAL(valueChanged(int)),
			this,
			SLOT(slider_value_changed(int))
	);


	QObject::connect(this, SIGNAL(contract_edges(unsigned)), parent, SLOT(contract_edges(unsigned)));

	QObject::connect(this->pushButton_collapse, SIGNAL(clicked()), this, SLOT(send_contract_edges()));

}

void Menu_edge_contraction::slider_value_changed(int new_slider_value){
	int num_collapses =
			(horizontalSlider->value()==1)? 1 : horizontalSlider->value() * model_.num_vertices() / 100 ;
	this->txt_nb_vertices->setNum((int)model_.num_vertices());
	this->txt_nb_collapses->setNum(num_collapses);
	this->spinBox_nb_remaining_vertices->setValue(model_.num_vertices()-num_collapses);
}


void
Menu_edge_contraction::update_slider_value(){
	int num_vertices = model_.num_vertices();
	int num_collapses = (horizontalSlider->value()==1)? 1 : horizontalSlider->value() * num_vertices / 100 ;
	int horizontal_slider_position = num_vertices>0?  num_collapses/(double)num_vertices * 100 : 1  ;
	horizontalSlider->setValue(horizontal_slider_position);
}


void
Menu_edge_contraction::update_gui_numbers(){
	update_slider_value();
	bool ok;
	int num_collapses = this->txt_nb_collapses->text().toInt(&ok,10);
	if(!ok) return;
	this->txt_nb_vertices->setNum((int)model_.num_vertices());
	this->txt_nb_collapses->setNum(num_collapses);
	this->spinBox_nb_remaining_vertices->setValue(model_.num_vertices()-num_collapses);
}

void
Menu_edge_contraction::update_gui_numbers(int new_value){
	update_gui_numbers();
}


void
Menu_edge_contraction::send_contract_edges(){
	emit(contract_edges(num_collapses()));
	update_gui_numbers();
}

unsigned
Menu_edge_contraction::num_vertices(){
	return model_.num_vertices();
}

unsigned
Menu_edge_contraction::num_collapses(){
	return (horizontalSlider->value()==1)? 1 : horizontalSlider->value() * num_vertices() / 100 ;
}


#include "Menu_edge_contraction.moc"

#endif /* MENU_EDGE_CONTRACTION_CPP_ */

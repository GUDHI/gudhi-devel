/*    This file is part of the MMA Library - https://gitlab.inria.fr/dloiseau/multipers - which is released under MIT.
 *    See file LICENSE for full license details.
 *    Author(s):       David Loiseaux
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 */
/**
 * @file line_filtration_translation.h
 * @author David Loiseaux
 * @brief
 */


#ifndef LINE_FILTRATION_TRANSLATION_H_INCLUDED
#define LINE_FILTRATION_TRANSLATION_H_INCLUDED

#include "utilities.h"
#include "box.h"

namespace utils{
	using point_type = std::vector<value_type>;
	template<typename T=value_type>
	class Line
	{
	
	public:
		Line();
		Line(point_type x);
		Line(point_type x, point_type v);
		point_type push_forward(point_type x) const;
		point_type push_back(point_type x) const;
		dimension_type get_dim() const;
		std::pair<point_type, point_type> get_bounds(const Box<T> &box) const;


	private:
		point_type basepoint_; // any point on the line
		point_type direction_; // direction of the line

	};
	template<typename T>
	Line<T>::Line(){}

	template<typename T>
	Line<T>::Line(point_type x){
		this->basepoint_.swap(x);
		// this->direction_ = {}; // diagonal line
	}
	template<typename T>
	Line<T>::Line(point_type x, point_type v){
		this->basepoint_.swap(x);
		this->direction_.swap(v);
	}
	template<typename T>
	point_type Line<T>::push_forward(point_type x) const{ //TODO remove copy
		x -= basepoint_;
		value_type t=negInf;
		for (unsigned int i = 0; i<x.size(); i++){
			value_type dir  = this->direction_.size() > i ? direction_[i] : 1;
			t = std::max(t, x[i]/dir);
		}
		point_type out(basepoint_.size());
		for (unsigned int i = 0; i < out.size(); i++)
			out[i] = basepoint_[i] + t * (this->direction_.size() > i ? direction_[i] : 1) ;
		return out;
	}
	template<typename T>
	point_type Line<T>::push_back(point_type x) const{
		x -= basepoint_;
		value_type t=inf;
		for (unsigned int i = 0; i<x.size(); i++){
			value_type dir  = this->direction_.size() > i ? direction_[i] : 1;
			t = std::min(t, x[i]/dir);
		}
		point_type out(basepoint_.size());
		for (unsigned int i = 0; i < out.size(); i++)
			out[i] = basepoint_[i] + t * (this->direction_.size() > i ? direction_[i] : 1) ;
		return out;
	}
	template<typename T>
	dimension_type Line<T>::get_dim() const{
		return basepoint_.size();
	}
	template<typename T>
	std::pair<point_type, point_type> Line<T>::get_bounds(const Box<T> &box) const{
		return {this->push_forward(box.get_bottom_corner()), this->push_back(box.get_upper_corner())};
	}


}

#endif // LINE_FILTRATION_TRANSLATION_H_INCLUDED

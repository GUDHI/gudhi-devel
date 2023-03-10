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

#include "box.h"
#include "finitely_critical_filtrations.h"

namespace Gudhi::multi_filtrations{

	
	template<typename T>
	class Line
	{
		
	public:
		using point_type = Finitely_critical_multi_filtration<T>;
		Line();
		Line(point_type x);
		Line(point_type x, point_type v);
		point_type push_forward(point_type x) const;
		point_type push_back(point_type x) const;
		int get_dim() const;
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
	typename Line<T>::point_type Line<T>::push_forward(point_type x) const { //TODO remove copy
		x -= basepoint_;
		T t = - std::numeric_limits<T>::infinity();;
		for (std::size_t i = 0; i<x.size(); i++){
			T dir  = this->direction_.size() > i ? direction_[i] : 1;
			t = std::max(t, x[i]/dir);
		}
		point_type out(basepoint_.size());
		for (unsigned int i = 0; i < out.size(); i++)
			out[i] = basepoint_[i] + t * (this->direction_.size() > i ? direction_[i] : 1) ;
		return out;
	}
	template<typename T>
	typename Line<T>::point_type Line<T>::push_back(point_type x) const{
		x -= basepoint_;
		T t = std::numeric_limits<T>::infinity();
		for (unsigned int i = 0; i<x.size(); i++){
			T dir  = this->direction_.size() > i ? direction_[i] : 1;
			t = std::min(t, x[i]/dir);
		}
		point_type out(basepoint_.size());
		for (unsigned int i = 0; i < out.size(); i++)
			out[i] = basepoint_[i] + t * (this->direction_.size() > i ? direction_[i] : 1) ;
		return out;
	}
	template<typename T>
	int Line<T>::get_dim() const{
		return basepoint_.size();
	}
	template<typename T>
	std::pair<typename Line<T>::point_type, typename Line<T>::point_type> Line<T>::get_bounds(const Box<T> &box) const{
		return {this->push_forward(box.get_bottom_corner()), this->push_back(box.get_upper_corner())};
	}
}

#endif // LINE_FILTRATION_TRANSLATION_H_INCLUDED

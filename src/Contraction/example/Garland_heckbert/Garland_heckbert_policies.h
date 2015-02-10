/*
 * Garland_heckbert_policies.h
 *  Created on: Feb 10, 2015
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


#ifndef GARLAND_HECKBERT_POLICIES_H_
#define GARLAND_HECKBERT_POLICIES_H_

#include "gudhi/Edge_contraction.h"
#include "Error_quadric.h"

template<typename EdgeProfile>
class GH_visitor: public Gudhi::contraction::Contraction_visitor<EdgeProfile> {
	typedef typename EdgeProfile::Complex Complex;
	typedef typename Complex::Point Point;
	Complex& complex_;
public:
	GH_visitor(Complex& complex):complex_(complex){}


	void on_started(Complex & complex) override{
		//Compute quadrics for every vertex v
		//The quadric of v consists in the sum of quadric
		//of every triangles passing through v weighted by its area
		for(auto v : complex.vertex_range()){
			auto & quadric_v(complex[v].quadric);
			for(auto t : complex.triangle_range(v)){
				auto t_it = t.begin();
				const auto& p0(complex.point(*t_it++));
				const auto& p1(complex.point(*t_it++));
				const auto& p2(complex.point(*t_it++));
				quadric_v+=Error_quadric<Point>(p0,p1,p2);
			}
		}
	}

	/**
	 * @brief Called when an edge is about to be contracted and replaced by a vertex whose position is *placement.
	 */
	void on_contracting(EdgeProfile const &profile, boost::optional< Point > placement)
	override{
		profile.v0().quadric += profile.v1().quadric;
	}
};


template<typename EdgeProfile>
class GH_placement : public Gudhi::contraction::Placement_policy<EdgeProfile>{
	typedef typename EdgeProfile::Complex Complex;
	Complex& complex_;
public:
	typedef typename EdgeProfile::Point Point;
	typedef typename Gudhi::contraction::Placement_policy<EdgeProfile>::Placement_type Placement_type;

	GH_placement(Complex& complex):complex_(complex){}

	Placement_type operator()(const EdgeProfile& profile) const override{
		auto sum_quad(profile.v0().quadric);
		sum_quad += profile.v1().quadric;

		boost::optional<Point> min_quadric_pt(sum_quad.min_cost());
		if (min_quadric_pt)
			return Placement_type(*min_quadric_pt);
		else
			return profile.p0();
	}
};

template<typename EdgeProfile>
class GH_cost : public Gudhi::contraction::Cost_policy<EdgeProfile>{
	typedef typename EdgeProfile::Complex Complex;
	Complex& complex_;
public:

	typedef typename Gudhi::contraction::Cost_policy<EdgeProfile>::Cost_type Cost_type;
	typedef typename Complex::Point Point;

	GH_cost(Complex& complex):complex_(complex){}

	Cost_type operator()( EdgeProfile const& profile, boost::optional<Point> const& new_point ) const override {
		Cost_type res;
		if (new_point){
			auto sum_quad(profile.v0().quadric);
			sum_quad += profile.v1().quadric;
			res = sum_quad.cost(*new_point);
		}
		return res;
	}
};



#endif /* GARLAND_HECKBERT_POLICIES_H_ */

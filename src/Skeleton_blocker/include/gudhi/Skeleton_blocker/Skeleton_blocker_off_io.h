/*
 * Skeleton_blocker_off_io.h
 *  Created on: Nov 28, 2014
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


#ifndef SKELETON_BLOCKER_OFF_IO_H_
#define SKELETON_BLOCKER_OFF_IO_H_

#include "gudhi/Off_reader.h"

namespace Gudhi {

namespace skbl {

template<typename Complex>
class Skeleton_blocker_off_visitor_reader{
	Complex& complex_;
	typedef typename Complex::Vertex_handle Vertex_handle;
	typedef typename Complex::Point Point;

	const bool load_only_points_;

public:
	Skeleton_blocker_off_visitor_reader(Complex& complex,bool load_only_points = false):
		complex_(complex),
		load_only_points_(load_only_points){}


	void init(int dim,int num_vertices,int num_faces,int num_edges){
		//todo do an assert to check that this number are correctly read
		//todo reserve size for vector points
	}


	void point(const std::vector<double>& point){
		complex_.add_vertex(point);
	}

	void maximal_face(const std::vector<int>& face){
		if (!load_only_points_){
			for(size_t i = 0 ; i < face.size();++i)
				for(size_t j = i+1 ; j < face.size();++j){
					complex_.add_edge(Vertex_handle(face[i]),Vertex_handle(face[j]));
				}
			}
	}

	void done(){
	}
};

template<typename Complex>
class Skeleton_blocker_off_reader{
public:
	/**
	 * name_file : file to read
	 * read_complex : complex that will receive the file content
	 * read_only_points : specify true if only the points must be read
	 */
	Skeleton_blocker_off_reader(const std::string & name_file,Complex& read_complex,bool read_only_points = false):valid_(false){
		std::ifstream stream(name_file);
		if(stream.is_open()){
			Skeleton_blocker_off_visitor_reader<Complex> off_visitor(read_complex,read_only_points);
			Off_reader off_reader(stream);
			valid_ = off_reader.read(off_visitor);
		}
	}

	/**
	 * return true iff reading did not meet problems.
	 */
	bool is_valid() const{
		return valid_;
	}

private:
	bool valid_;
};

}  // namespace skbl


}  // namespace Gudhi


#endif /* SKELETON_BLOCKER_OFF_IO_H_ */

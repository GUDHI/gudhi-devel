/*    This file is a prototype for the Gudhi Library.
 *    Author(s):       Clément Maria
 *    Copyright (C) 2021 Inria
 *    This version is under developement, please do not redistribute this software. 
 *    This program is for academic research use only. 
 */

#include <iostream>
#include <fstream>
#include "gudhi/reader_utils.h"
#include <gudhi/distance_functions.h>
#include <gudhi/Points_off_io.h>
#include <gudhi/Point_cloud.h>
#include <CGAL/Epick_d.h>

// Types definition
using K = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
using Point_d = typename K::Point_d;
using Points_off_reader = Gudhi::Points_off_reader<Point_d>;
using Coordinate_type = double;

int main(int argc, char* argv[])
{
	if(argc != 2) {
		std::cout << "./exe filename\n\n"; 
		std::cout << "Compute the number of duplicates in a point cloud, then perturb the points with a small noise, then check the distance genericity.\n\nfilename is the name of a file containing a set of points, in off format.\n";
		return 0;
	}

	std::string off_file_points = argv[1];
	//ensure there is no duplicate points
	Points_off_reader off_reader(off_file_points); //read points

	if(off_reader.get_point_cloud().empty()) { std::cout << "Empty point cloud.\n"; return 0;}

	int dim = off_reader.get_point_cloud().begin()->size();

	std::cout << "Number of input points: " << off_reader.get_point_cloud().size() << "\n";
	std::cout << "There are " << duplicates(off_reader.get_point_cloud()) << " duplicated points in file " << off_file_points << "\n"; 
	//remove all duplicates
  auto points_unique(off_reader.get_point_cloud());

	remove_duplicates(points_unique);
	std::cout <<"Number of points after removals of duplicates: " << points_unique.size() << "\n";
	//compute diameter and min distance between two distinct points
	K k_d;   
	auto dist = k_d.squared_distance_d_object();

	auto spread = max_min_distances<Coordinate_type>(points_unique,dist);

	std::cout << "Diameter: " << std::sqrt(spread.first) << "   closest points: " << std::sqrt(spread.second) << "    Points spread: " << std::sqrt((spread.first/spread.second)) << "\n";

	//perturb the points by 0.001 of the diameter
	random_perturbation(points_unique, dim, 0.001 * spread.first);

	//check whether the point cloud is generic w.r.t. pairwise distances
	auto min_delta = generic_distances<Coordinate_type>(points_unique, dist);

	std::cout << "Generic distances: " << (min_delta > 0.) << "\n";
	return 0;
}




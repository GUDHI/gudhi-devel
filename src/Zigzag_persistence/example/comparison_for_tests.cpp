/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <algorithm>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>  // for pair
#include <vector>
#include <unordered_map>

#include <boost/range/adaptors.hpp>

#include <gudhi/Zigzag_persistence/oscillating_rips_iterators.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/distance_functions.h>

using ST = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_wide_indexation>;
using Filtration_value = ST::Filtration_value;
using OR = Gudhi::zigzag_persistence::Oscillating_rips_edge_range<Filtration_value>;
using ZE = Gudhi::zigzag_persistence::Zigzag_edge<Filtration_value>;
using Point = std::vector<double>;

void print_points(const std::vector<Point>& points){
	std::cout << "Number of points: " << points.size() << "\n";
	for (const Point& p : points){
		std::cout << "(" << p[0] << ", " << p[1] << ")\n";
	}
	std::cout << "\n";
}

std::vector<Point> build_point_cloud(unsigned int numberOfPoints, int seed){
	std::vector<Point> finalPoints;
	std::set<Point> points;
	std::random_device dev;
	std::mt19937 rng(dev());
	if (seed > -1) rng.seed(seed);
	std::uniform_real_distribution<double> dist(0,10);

	for (unsigned int i = 0; i < numberOfPoints; ++i){
		auto res = points.insert({dist(rng), dist(rng)});
		while(!res.second){
			res = points.insert({dist(rng), dist(rng)});
		}
		finalPoints.push_back(*res.first);
	}

	// print_points(finalPoints);

	return finalPoints;
}

int main(int argc, char* const argv[]) {
	if (argc != 4 && argc != 5) {
		std::cout << "Usage: ./comp nu mu nomberOfPoints [seed]\n";
		return 0;
	}

	double nu = std::stod(argv[1]);
	double mu = std::stod(argv[2]);
	unsigned int numberOfPoints = std::stoi(argv[3]);
	int seed = -1;

	if (argc == 5)
		seed = std::stoi(argv[4]);

	std::cout << "nu, mu: " << nu << ", " << mu << "\n";
	std::cout << "number of points: " << numberOfPoints << "\n";
	std::cout << "seed: " << seed << "\n";

	std::vector<Point> points = build_point_cloud(numberOfPoints, seed);

	std::vector<ZE> edges_v1 = OR::compute_oscillating_rips_edges(nu, mu, points, Gudhi::Euclidean_distance(), OR::Order_policy::FARTHEST_POINT_ORDERING);
	std::cout << edges_v1.size() << "\n";

	// unsigned int i = 0;
	// for (const auto& e : OR::compute_oscillating_rips_edges_as_iterator(nu, mu, points, Gudhi::Euclidean_distance(), OR::Order_policy::FARTHEST_POINT_ORDERING)){
	// 	++i;
	// }
	// std::cout << i << "\n";

	// unsigned int i = 0;
	// for (const auto& e : OR::compute_oscillating_rips_edges_as_iterator(nu, mu, points, Gudhi::Euclidean_distance(), OR::Order_policy::FARTHEST_POINT_ORDERING)){
	// 	if (i < edges_v1.size()){
	// 		if (!(edges_v1[i] == e)){
	// 		std::cout << "[" << i << "] different:\n";
	// 		std::cout << edges_v1[i].get_smallest_vertex() << ", " << edges_v1[i].get_biggest_vertex() << ", " << edges_v1[i].get_filtration_value() << ", " << edges_v1[i].get_direction() << "\n";
	// 		std::cout << e.get_smallest_vertex() << ", " << e.get_biggest_vertex() << ", " << e.get_filtration_value() << ", " << e.get_direction() << "\n";
	// 	}/*  else {
	// 		std::cout << "[" << i << "] same:\n";
	// 		std::cout << edges_v1[i].get_smallest_vertex() << ", " << edges_v1[i].get_biggest_vertex() << ", " << edges_v1[i].get_filtration_value() << ", " << edges_v1[i].get_direction() << "\n";
	// 		std::cout << e.get_smallest_vertex() << ", " << e.get_biggest_vertex() << ", " << e.get_filtration_value() << ", " << e.get_direction() << "\n";
	// 	} */
	// 	} else {
	// 		std::cout << "[" << i << "] too long:\n";
	// 		std::cout << e.get_smallest_vertex() << ", " << e.get_biggest_vertex() << ", " << e.get_filtration_value() << ", " << e.get_direction() << "\n";
	// 	}
	// 	++i;
	// }

	return 0;
}

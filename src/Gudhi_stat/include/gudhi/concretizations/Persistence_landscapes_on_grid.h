/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2015  INRIA Sophia-Saclay (France)
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
 */

#pragma once
#ifndef Persistence_landscapes_on_grid_H
#define Persistence_landscapes_on_grid_H

#include <cmath>
#include <iostream>
#include <vector>
#include <limits>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <unistd.h>





#include "../abstract_classes/Abs_Vectorized_topological_data.h"
#include "../abstract_classes/Abs_Topological_data_with_averages.h"
#include "../abstract_classes/Abs_Topological_data_with_distances.h"
#include "../abstract_classes/Abs_Real_valued_topological_data.h"
#include "../abstract_classes/Abs_Topological_data_with_scalar_product.h"
using namespace std;

class Persistence_landscapes_on_grid :  public Abs_Vectorized_topological_data , 
										public Abs_Topological_data_with_distances, 
										public Abs_Real_valued_topological_data, 
										public Abs_Topological_data_with_averages, 
										public Abs_Topological_data_with_scalar_product
{
public:
	aaa
private:
	double grid_min;
	double grid_max;
	size_t number_of_points;
	std::vector< std::vector< double > > values_of_landscapes;
};

#endif

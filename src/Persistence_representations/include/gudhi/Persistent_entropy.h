/*    
 *    Author(s):       Alejandro Luque Cerpa
 */

#ifndef PERSISTENT_ENTROPY_H_
#define PERSISTENT_ENTROPY_H_

// gudhi include
#include <gudhi/read_persistence_from_file.h>
#include <gudhi/Persistence_intervals.h>

// standard include
#include <limits>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <functional>
#include <utility>
#include <string>

namespace Gudhi {
namespace Persistent_entropy{

double persistent_entropy_from_file(const char* filename, bool keep_inf=false, double val_inf=-1.0, bool normalize=false){
    std::vector<std::pair<double, double>> intervals;
	if(keep_inf){
		intervals = Gudhi::Persistence_representations::read_persistence_intervals_in_one_dimension_from_file(filename, -1, val_inf);
	}else{
		intervals = Gudhi::Persistence_representations::read_persistence_intervals_in_one_dimension_from_file(filename);
	}	
	double res = 0;
	double SL = 0;
	double aux;
	double pi;
	int N = intervals.size();
	std::vector<double> props;
	// This function optimizes the memory access
	props.reserve(N);
	for (int i = 0; i < N; i++){
		aux = intervals.at(i).second - intervals.at(i).first;
		SL += aux;
		props.push_back(aux);
	}
	for (int i = 0; i < N; i++){
		pi = props.at(i) / SL;
		res -= pi * log(pi);
	}	
        if (normalize){
		res = res / log(N);
	}
	return res;
} // double persistent_entropy_from_file


std::vector<double> persistent_entropy_from_files_list(const char* filename, bool keep_inf=false, double val_inf=-1.0, bool normalize=false){
	std::vector<double> res;
	std::ifstream in(filename);
	if(!in.is_open()){
		std::string error_str("Persistent entropy from files list - Unable to open file ");
		error_str.append(filename);
		std::cerr << error_str << std::endl;
		throw std::invalid_argument(error_str);
	}

	while (!in.eof()){
		std::string line;
		getline(in, line);
		char cstr[line.size() + 1];
		line.copy(cstr, line.size() + 1);
		cstr[line.size()] = '\0';
		if(line.length() != 0 && line[0] != '#'){
			res.push_back(persistent_entropy_from_file(cstr, keep_inf, val_inf, normalize));
		}
	}
	return res;
} // double persistent_entropy_from_files_list

} // namespace Persistent_entropy
} // namespace Gudhi

#endif  // PERSISTENT_ENTROPY_H_

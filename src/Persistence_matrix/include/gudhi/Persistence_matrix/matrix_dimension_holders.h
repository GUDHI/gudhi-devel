/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022-23 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_MATRIX_DIM_HOLDER_H
#define PM_MATRIX_DIM_HOLDER_H

#include <utility>	//std::swap, std::move & std::exchange
#include <vector>

namespace Gudhi {
namespace persistence_matrix {

struct Dummy_matrix_dimension_holder{
	template<typename dimension_type>
	Dummy_matrix_dimension_holder([[maybe_unused]] dimension_type maximalDimension){}

	friend void swap([[maybe_unused]] Dummy_matrix_dimension_holder& d1, [[maybe_unused]] Dummy_matrix_dimension_holder& d2){}
};

template<typename dimension_type>
class Matrix_max_dimension_holder
{
public:
	Matrix_max_dimension_holder(dimension_type maximalDimension = -1) : maxDim_(maximalDimension){};
	Matrix_max_dimension_holder(const Matrix_max_dimension_holder& toCopy) : maxDim_(toCopy.maxDim_){};
	Matrix_max_dimension_holder(Matrix_max_dimension_holder&& other) noexcept 
		: maxDim_(std::exchange(other.maxDim_, -1)){};

	dimension_type get_max_dimension() const{ return maxDim_; };

	Matrix_max_dimension_holder& operator=(const Matrix_max_dimension_holder& other){
		std::swap(maxDim_, other.maxDim_); 
		return *this; 
	};

	friend void swap(Matrix_max_dimension_holder& matrix1, Matrix_max_dimension_holder& matrix2){
		std::swap(matrix1.maxDim_, matrix2.maxDim_);
	}

protected:
	dimension_type maxDim_;

	void update_up(dimension_type dimension){
		if (maxDim_ < dimension) maxDim_ = dimension;
	};
};

template<typename dimension_type>
class Matrix_all_dimension_holder
{
public:
	Matrix_all_dimension_holder(dimension_type maximalDimension = -1) : maxDim_(maximalDimension){};
	Matrix_all_dimension_holder(const Matrix_all_dimension_holder& toCopy) 
		: dimensions_(toCopy.dimensions_), maxDim_(toCopy.maxDim_){};
	Matrix_all_dimension_holder(Matrix_all_dimension_holder&& other) noexcept 
		: dimensions_(std::move(other.dimensions_)), maxDim_(std::exchange(other.maxDim_, -1)){};

	dimension_type get_max_dimension() const{ return maxDim_; };

	Matrix_all_dimension_holder& operator=(Matrix_all_dimension_holder other){
		std::swap(maxDim_, other.maxDim_); 
		dimensions_.swap(other.dimensions_);
		return *this; 
	};

	friend void swap(Matrix_all_dimension_holder& matrix1, Matrix_all_dimension_holder& matrix2){
		std::swap(matrix1.maxDim_, matrix2.maxDim_);
		matrix1.dimensions_.swap(matrix2.dimensions_);
	}

protected:
	std::vector<unsigned int> dimensions_;
	dimension_type maxDim_;

	void update_up(unsigned int dimension){
		if (dimensions_.size() <= dimension) dimensions_.resize(dimension + 1, 0);
		++(dimensions_[dimension]);
		maxDim_ = dimensions_.size() - 1;
	};

	void update_down(unsigned int dimension){
		--(dimensions_[dimension]);		//assumes dimension already exists and is not 0
		while (dimensions_.back() == 0)	//assumes that there is always at least one cell in the complex
			dimensions_.pop_back();
		maxDim_ = dimensions_.size() - 1;
	};
};

} //namespace persistence_matrix
} //namespace Gudhi

#endif // PM_MATRIX_DIM_HOLDER_H

/*    This file is part of the MMA Library - https://gitlab.inria.fr/dloiseau/multipers - which is released under MIT.
 *    See file LICENSE for full license details.
 *    Author(s):       David Loiseaux
 *
 *    Copyright (C) 2021 Inria
 *
 *    Modification(s):
 * 
 */
/**
 * @file box.h
 * @author David Loiseaux
 * @brief BOX.
 */

#ifndef BOX_H_INCLUDED
#define BOX_H_INCLUDED

#include <vector>
#include <ostream>
// #include <iomanip>
#include <cmath>
#include <limits>
#include <cassert>
// #include "gudhi/Simplex_tree.h"






/**
 * @brief Holds the square box on which to compute.
 */

namespace utils{

template<typename T>
class Box
{
	using corner_type = std::vector<T>;
	using point_type = std::vector<T>;
public:
	Box();
	Box(const corner_type& bottomCorner, const corner_type& upperCorner);
	Box(const std::pair<corner_type, corner_type>& box);

	void inflate(T delta);
	const corner_type& get_bottom_corner() const;
	const corner_type& get_upper_corner() const;
	bool contains(const point_type& point) const;
	void infer_from_filters(const std::vector<std::vector<T>> &Filters_list);
    bool is_trivial() const ;

private:
	corner_type bottomCorner_;
	corner_type upperCorner_;
};

template<typename T>
inline Box<T>::Box()
{}

template<typename T>
inline Box<T>::Box(const corner_type &bottomCorner, const corner_type &upperCorner)
	: bottomCorner_(bottomCorner),
	  upperCorner_(upperCorner)
{
	assert(bottomCorner.size() == upperCorner.size()
		   && is_smaller(bottomCorner, upperCorner)
		   && "This box is trivial !");
}

template<typename T>
inline Box<T>::Box(const std::pair<corner_type, corner_type> &box)
	: bottomCorner_(box.first),
	  upperCorner_(box.second)
{}


template<typename T>
inline void Box<T>::inflate(T delta)
{
#pragma omp simd
	for (unsigned int i = 0; i < bottomCorner_.size(); i++){
		bottomCorner_[i] -= delta;
		upperCorner_[i] += delta;
	}
}

template<typename T>
inline void Box<T>::infer_from_filters(const std::vector<std::vector<T>> &Filters_list){
	unsigned int dimension = Filters_list.size();
	unsigned int nsplx = Filters_list[0].size();
	std::vector<T> lower(dimension);
	std::vector<T> upper(dimension);
	for (unsigned int i =0; i < dimension; i++){
		T min = Filters_list[i][0];
		T max = Filters_list[i][0];
		for (unsigned int j=1; j<nsplx; j++){
			min = std::min(min, Filters_list[i][j]);
			max = std::max(max, Filters_list[i][j]);
		}
		lower[i] = min;
		upper[i] = max;
	}
	bottomCorner_.swap(lower);
	upperCorner_.swap(upper);
}
template<typename T>
inline bool Box<T>::is_trivial() const {
    return bottomCorner_.empty() || upperCorner_.empty() || bottomCorner_.size() != upperCorner_.size();
}

template<typename T>
inline const std::vector<T> &Box<T>::get_bottom_corner() const
{
	return bottomCorner_;
}

template<typename T>
inline const std::vector<T> &Box<T>::get_upper_corner() const
{
	return upperCorner_;
}

template<typename T>
inline bool Box<T>::contains(const std::vector<T> &point) const
{
	if (point.size() != bottomCorner_.size()) return false;

	for (unsigned int i = 0; i < point.size(); i++){
		if (point[i] < bottomCorner_[i]) return false;
		if (point[i] > upperCorner_[i]) return false;
	}

	return true;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const Box<T>& box)
{
    os << "Box -- Bottom corner : ";
    os << box.get_bottom_corner();
    os << ", Top corner : ";
    os << box.get_upper_corner();
    return os;
}

} // namespace utils


#endif // BOX_H_INCLUDED

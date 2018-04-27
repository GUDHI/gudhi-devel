/*    This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
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

#ifndef POINT_H_
#define POINT_H_

#include <cmath>
#include <vector>
#include <cassert>
#include <cstddef>
#include <initializer_list>

class Point_d {
 public:
  Point_d(size_t dim = 3) : coords_(dim, 0) { }

  Point_d(const Point_d& other) : coords_(other.coords_) { }

  Point_d(const std::initializer_list<double>& list) : coords_(list) { }

  template<typename CoordsIt>
  Point_d(CoordsIt begin, CoordsIt end) : coords_(begin, end) { }

  size_t dimension() const {
    return coords_.size();
  }

  double x() const {
    return coords_[0];
  }

  double y() const {
    return coords_[1];
  }

  double z() const {
    return coords_[2];
  }

  double& x() {
    return coords_[0];
  }

  double& y() {
    return coords_[1];
  }

  double& z() {
    return coords_[2];
  }

  std::vector<double>::const_iterator begin() const {
    return coords_.begin();
  }

  std::vector<double>::const_iterator end() const {
    return coords_.end();
  }

  double& operator[](unsigned i) {
    return coords_[i];
  }

  const double& operator[](unsigned i) const {
    return coords_[i];
  }

  double squared_norm() const {
    double res = 0;
    for (auto x : coords_)
      res += x * x;
    return res;
  }

  friend double squared_dist(const Point_d& p1, const Point_d& p2) {
    assert(p1.dimension() == p2.dimension());
    double res = 0;
    for (unsigned i = 0; i < p1.coords_.size(); ++i)
      res += (p1[i] - p2[i])*(p1[i] - p2[i]);
    return res;
  }

  /**
   * dot product
   */
  double operator*(const Point_d& other) const {
    assert(dimension() == other.dimension());
    double res = 0;
    for (unsigned i = 0; i < coords_.size(); ++i)
      res += coords_[i] * other[i];
    return res;
  }

  /**
   * only if points have dimension 3
   */
  Point_d cross_product(const Point_d& other) {
    assert(dimension() == 3 && other.dimension() == 3);
    Point_d res(3);
    res[0] = (*this)[1] * other[2] - (*this)[2] * other[1];
    res[1] = (*this)[2] * other[0] - (*this)[0] * other[2];
    res[2] = (*this)[0] * other[1] - (*this)[1] * other[0];
    return res;
  }

  Point_d operator+(const Point_d& other) const {
    assert(dimension() == other.dimension());
    Point_d res(dimension());
    for (unsigned i = 0; i < coords_.size(); ++i)
      res[i] = (*this)[i] + other[i];
    return res;
  }

  Point_d operator*(double lambda) const {
    Point_d res(dimension());
    for (unsigned i = 0; i < coords_.size(); ++i)
      res[i] = (*this)[i] * lambda;
    return res;
  }

  Point_d operator/(double lambda) const {
    Point_d res(dimension());
    for (unsigned i = 0; i < coords_.size(); ++i)
      res[i] = (*this)[i] / lambda;
    return res;
  }

  Point_d operator-(const Point_d& other) const {
    assert(dimension() == other.dimension());
    Point_d res(dimension());
    for (unsigned i = 0; i < coords_.size(); ++i)
      res[i] = (*this)[i] - other[i];
    return res;
  }

  friend Point_d unit_normal(const Point_d& p1, const Point_d& p2, const Point_d& p3) {
    assert(p1.dimension() == 3);
    assert(p2.dimension() == 3);
    assert(p3.dimension() == 3);
    Point_d p1p2 = p2 - p1;
    Point_d p1p3 = p3 - p1;
    Point_d res(p1p2.cross_product(p1p3));
    return res / std::sqrt(res.squared_norm());
  }

 private:
  std::vector<double> coords_;
};

#endif  // POINT_H_

/*
 * Error_quadric.h
 *
 *  Created on: 24 janv. 2014
 *      Author: dsalinas
 */

#ifndef ERROR_QUADRIC_H_
#define ERROR_QUADRIC_H_

#include <vector>
#include <utility>
#include <boost/optional/optional.hpp>


template <typename Point> class Error_quadric{
private :
	double coeff[10];

public :
	Error_quadric(){
		clear();
	}

	/**
	 * Quadric corresponding to the L2 distance to the plane.
	 *
	 * According to the notation of Garland Heckbert, they
	 * denote a quadric symetric matrix as :
	 * Q = [ q11 q12 q13 q14]
	 *     [ q12 q22 q23 q24]
	 *     [ q13 q23 q33 q34]
	 *     [ q14 q24 q34 q44]
	 *
	 * which is represented by a vector with 10 elts that
	 * are denoted ci for clarity with :
	 * Q = [ c0  c1  c2  c3 ]
	 *     [ c1  c4  c5  c6 ]
	 *     [ c2  c5  c7  c8 ]
	 *     [ c3  c6  c8  c9 ]
	 *
	 * The constructor return the quadrics that represents
	 * the squared distance to the plane defined by triangle p0,p1,p2
	 * times the area of triangle p0,p1,p2.
	 */
	Error_quadric(const Point & p0,const Point & p1,const Point & p2){

		Point normal(unit_normal(p0,p1,p2));
		double a=normal[0];
		double b=normal[1];
		double c=normal[2];
		double d= -a*p0[0]-b*p0[1]-c*p0[2];
		coeff[0] = a*a ;
		coeff[1] = a*b ;
		coeff[2] = a*c ;
		coeff[3] = a*d ;
		coeff[4] = b*b ;
		coeff[5] = b*c ;
		coeff[6] = b*d ;
		coeff[7] = c*c ;
		coeff[8] = c*d ;
		coeff[9] = d*d ;

		double area_p0p1p2 = std::sqrt(squared_area(p0,p1,p2));
		for(auto& x : coeff)
			x*= area_p0p1p2;
	}


	inline double squared_area(const Point& p0,const Point& p1,const Point& p2) {
		//if (x1,x2,x3) = p1-p0 and (y1,y2,y3) = p2-p0
		//then the squared area is = (u^2+v^2+w^2)/4
		//with:		u = x2 * y3 - x3 * y2;
		//		    v = x3 * y1 - x1 * y3;
		//	     	w = x1 * y2 - x2 * y1;
		Point p0p1(p1-p0);
		Point p0p2(p2-p0);
		double A = p0p1[1] * p0p2[2] - p0p1[2] * p0p2[1];
		double B = p0p1[2] * p0p2[0] - p0p1[0] * p0p2[2];
		double C = p0p1[0] * p0p2[1] - p0p1[1] * p0p2[0];
		return 1./4. * (A*A+B*B+C*C);
	}


	void clear(){
		for(auto& x:coeff)
			x=0;
	}

	Error_quadric& operator+=(const Error_quadric& other){
		if(this!=&other)
			for(int i = 0 ; i < 10; ++i)
				coeff[i] += other.coeff[i];
		return *this;
	}

	/**
	 * @return The quadric quost defined by the scalar product v^T Q v where Q is the quadratic form of Garland/Heckbert
	 */
	inline double cost(const Point& point) const{
		double cost =
				coeff[0]*point.x()*point.x()+coeff[4]*point.y()*point.y()+coeff[7]*point.z()*point.z()
				+2*(coeff[1]*point.x()*point.y()+coeff[5]*point.y()*point.z()+coeff[2]*point.z()*point.x())
				+2*(coeff[3]*point.x()+coeff[6]*point.y()+coeff[8]*point.z())
				+coeff[9];
		if(cost<0) return 0;
		else {
			return cost;
		}
	}

	inline double grad_determinant() const{
		return
				coeff[0] * coeff[4] * coeff[7]
				                            -  coeff[0] * coeff[5] * coeff[5]
				                                                           -  coeff[1] * coeff[1] * coeff[7]
				                                                                                          +2*coeff[1] * coeff[5] * coeff[2]
				                                                                                                                         -  coeff[4] * coeff[2] * coeff[2];
	}

	/**
	 * Return the point such that it minimizes the gradient of the quadric.
	 * Det must be passed with the determinant value of the gradient (should be non zero).
	 */
	inline Point solve_linear_gradient(double det) const{
		return Point({
				(-coeff[1]*coeff[5]*coeff[8]+coeff[1]*coeff[7]*coeff[6]+coeff[2]*coeff[8]*coeff[4]-coeff[2]*coeff[5]*coeff[6]-coeff[3]*coeff[4]*coeff[7]+coeff[3]*coeff[5]*coeff[5])/ det,
				(coeff[0]*coeff[5]*coeff[8]-coeff[0]*coeff[7]*coeff[6]-coeff[5]*coeff[2]*coeff[3]-coeff[1]*coeff[2]*coeff[8]+coeff[6]*coeff[2]*coeff[2]+coeff[1]*coeff[3]*coeff[7])/det,
				(-coeff[8]*coeff[0]*coeff[4]+coeff[8]*coeff[1]*coeff[1]+coeff[2]*coeff[3]*coeff[4]+coeff[5]*coeff[0]*coeff[6]-coeff[5]*coeff[1]*coeff[3]-coeff[1]*coeff[2]*coeff[6])/det
		});
	}


	/**
	 * returns the point that minimizes the quadric.
	 * It inverses the quadric if its determinant is higher that a given threshold .
	 * If the determinant is lower than this value the returned value is uninitialized.
	 */
	boost::optional<Point> min_cost(double scale=1) const{
		//		const double min_determinant = 1e-4 * scale*scale;
		const double min_determinant = 1e-5;
		boost::optional<Point> pt_res;
		double det = grad_determinant();
		if (std::abs(det)>min_determinant)
			pt_res = solve_linear_gradient(det);
		return pt_res;
	}

	friend std::ostream& operator<< (std::ostream& stream, const Error_quadric& quadric) {
		stream << "\n[ "<<quadric.coeff[0]<<","<<quadric.coeff[1]<<","<<quadric.coeff[2]<<","<<quadric.coeff[3]<<";\n";
		stream <<    " "<<quadric.coeff[1]<<","<<quadric.coeff[4]<<","<<quadric.coeff[5]<<","<<quadric.coeff[6]<<";\n";
		stream <<    " "<<quadric.coeff[2]<<","<<quadric.coeff[5]<<","<<quadric.coeff[7]<<","<<quadric.coeff[8]<<";\n";
		stream <<    " "<<quadric.coeff[3]<<","<<quadric.coeff[6]<<","<<quadric.coeff[8]<<","<<quadric.coeff[9]<<"]";
		return stream;
	}


};




#endif /* ERROR_QUADRIC_H_ */


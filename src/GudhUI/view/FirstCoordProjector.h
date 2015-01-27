/*
 * FirstCoordProjector.h
 *
 *  Created on: Aug 27, 2014
 *      Author: dsalinas
 */

#ifndef FIRSTCOORDPROJECTOR_H_
#define FIRSTCOORDPROJECTOR_H_

#include "utils/UI_utils.h"
#include "Projector3D.h"

class FirstCoordProjector3D : public Projector3D{
	typedef Projector3D::Point Point;
	typedef Projector3D::Point_3 Point_3;

	Point_3 operator()(const Point& p) const{
		assert(p.dimension()>=3);
		return Point_3(p.x(),p.y(),p.z());
	}
};

#endif /* FIRSTCOORDPROJECTOR_H_ */

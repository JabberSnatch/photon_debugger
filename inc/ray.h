#ifndef __YS_RAY_HPP__
#define __YS_RAY_HPP__

#include "vector.h"
#include "point.h"


namespace maths
{


struct Ray
{
	Ray(Point3f _o, Vec3f _d, Decimal _tMax, Decimal _time) :
		origin{ _o }, direction{ _d }, tMax{ _tMax }, time{ _time }
	{}
	inline Point3f operator()(Decimal _t) const { return origin + direction * _t; }
	inline bool HasNaNs() const { return origin.HasNaNs() || direction.HasNaNs() || std::isnan(tMax); }

	Point3f			origin;
	Vec3f			direction;
	mutable Decimal	tMax;
	Decimal			time;
};


}


#endif // __YS_RAY_HPP__

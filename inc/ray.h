#ifndef __YS_RAY_HPP__
#define __YS_RAY_HPP__

#include "maths.h"
#include "vector.h"
#include "point.h"
#include "bounds.h"


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

	template <typename T>
	bool IntersectP(Bounds<T, 3> const &_bounds, Decimal &_hit0, Decimal &_hit1) const
	{
		// r(t) = o + d * t
		// r(t).x = _bounds.min.x	<=> (o + d * t).x = _bounds.min.x
		//							<=> o.x + d.x * t = _bounds.min.x
		//							<=> t = (_bounds.min.x - o.x) / d.x
		// r(t) = _bounds.min		<=> o + d * t = _bounds.min
		//							<=> t = (_bounds.min - o) * (vec3(1) / d)
		Decimal t0{ 0._d }, t1{ tMax };
		for (int i = 0; i < 3; ++i)
		{
			Decimal	inv_dir = 1._d / direction[i];
			Decimal	tNear = (_bounds.min[i] - origin[i]) * inv_dir;
			Decimal tFar = (_bounds.max[i] - origin[i]) * inv_dir;

			if (tNear > tFar) std::swap(tNear, tFar);

			t0 = maths::Max(t0, tNear);
			t1 = maths::Min(t1, tFar);
			if (t0 > t1) return false;
		}
		_hit0 = t0; _hit1 = t1;
		return true;
	}
	template <typename T>
	bool IntersectP(Bounds<T, 3> const &_bounds) const
	{
		Decimal t0{ 0._d }, t1{ tMax };
		for (int i = 0; i < 3; ++i)
		{
			Decimal	inv_dir = 1._d / direction[i];
			Decimal	tNear = (_bounds.min[i] - origin[i]) * inv_dir;
			Decimal tFar = (_bounds.max[i] - origin[i]) * inv_dir;

			if (tNear > tFar) std::swap(tNear, tFar);

			t0 = maths::Max(t0, tNear);
			t1 = maths::Min(t1, tFar);
			if (t0 > t1) return false;
		}
		return true;
	}
	// NOTE: pbrt book provides an optimized IntersectP function for BVH tree traversal (page 128).
	//		 It will surely be useful later.
};


}


#endif // __YS_RAY_HPP__

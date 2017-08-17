#ifndef __YS_RAY_HPP__
#define __YS_RAY_HPP__

#include "maths/maths.h"
#include "maths/vector.h"
#include "maths/point.h"
#include "maths/bounds.h"


namespace maths {


struct Ray
{
	explicit Ray() = default;
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
	bool DoesIntersect(Bounds<T, 3> const &_bounds, Decimal &_hit0, Decimal &_hit1) const
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
	bool DoesIntersect(Bounds<T, 3> const &_bounds) const
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
	// Thanks, past self, you're the best.
	template <typename T>
	bool DoesIntersect(Bounds<T, 3> const &_bounds, Vec3f const &_direction_inverse,
					   Vector<int, 3> const &_is_negative) const
	{
		Decimal t_min = (_bounds[_is_negative.x].x - origin.x) * _direction_inverse.x;
		Decimal t_max = (_bounds[1 - _is_negative.x].x - origin.x) * _direction_inverse.x;
		Decimal const y_min = (_bounds[_is_negative.y].y - origin.y) * _direction_inverse.y;
		Decimal const y_max = (_bounds[1 - _is_negative.y].y - origin.y) * _direction_inverse.y;
		if (t_min > y_max || y_min > t_max)
			return false;
		t_min = Max(y_min, t_min);
		t_max = Min(y_max, t_max);

		Decimal const z_min = (_bounds[_is_negative.z].z - origin.z) * _direction_inverse.z;
		Decimal const z_max = (_bounds[1 - _is_negative.z].z - origin.z) * _direction_inverse.z;
		if (t_min > z_max || z_min > t_max)
			return false;
		t_min = Max(z_min, t_min);
		t_max = Min(z_max, t_max);

		return (t_min < tMax) && (t_max > 0);// && (t_min <= t_max);
	}
};


} // namespace maths


#endif // __YS_RAY_HPP__

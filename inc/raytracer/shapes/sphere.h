#pragma once
#ifndef __YS_SPHERE_HPP__
#define __YS_SPHERE_HPP__

#include <vector>

#include "api/api.h"
#include "raytracer/raytracer.h"
#include "raytracer/shape.h"


namespace raytracer
{

class Sphere final : public Shape
{
public:
	Sphere(maths::Transform const &_world_transform, bool _flip_normals,
		   maths::Decimal _radius, maths::Decimal _z_min, maths::Decimal _z_max,
		   maths::Decimal _phi_max);

	virtual bool Intersect(maths::Ray const &_ray,
						   maths::Decimal &_tHit,
						   SurfaceInteraction &_hit_info) const override;
	virtual bool DoesIntersect(maths::Ray const &_ray) const override;

	virtual maths::Decimal	Area() const override;

	virtual maths::Bounds3f	ObjectBounds() const override;

	maths::Decimal const	radius;
	maths::Decimal const	z_min, z_max;
	maths::Decimal const	theta_min, theta_max, phi_max;
};

} // namespace raytracer


#endif // __YS_SPHERE_HPP__

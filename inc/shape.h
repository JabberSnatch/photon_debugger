#pragma once
#ifndef __YS_SHAPE_HPP__
#define __YS_SHAPE_HPP__

#include "maths.h"
#include "raytracer.h"


namespace raytracer
{

class Shape
{
public:
	Shape(maths::Transform const &_world_transform, bool _flip_normals = false);
	virtual ~Shape() = default;

	virtual bool Intersect(maths::Ray const &_ray, 
						   maths::Decimal &_tHit, 
						   SurfaceInteraction &_hit_info) const = 0;
	virtual bool DoesIntersect(maths::Ray const &_ray) const;

	virtual maths::Decimal	Area() const = 0;

	virtual maths::Bounds3f	ObjectBounds() const = 0;
	virtual maths::Bounds3f	WorldBounds() const;

	maths::Transform const	&world_transform;
	bool const				swaps_handedness;
	bool const				flip_normals;
};

} // namespace raytracer


#endif // __YS_SHAPE_HPP__

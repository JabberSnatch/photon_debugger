#pragma once
#ifndef __YS_SHAPE_HPP__
#define __YS_SHAPE_HPP__

#include "maths/maths.h"
#include "raytracer/raytracer.h"

#include "maths/point.h"
#include "maths/vector.h"

namespace raytracer
{

class Shape
{
public:
	struct SurfacePoint
	{
		maths::Point3f const position;
		maths::Norm3f const normal;
		maths::Vec3f const position_error;
	};
public:
	Shape(maths::Transform const &_world_transform, bool _flip_normals = false);
	virtual ~Shape() = default;

	virtual bool Intersect(maths::Ray const &_ray, 
						   maths::Decimal &_tHit, 
						   SurfaceInteraction &_hit_info) const = 0;
	virtual bool DoesIntersect(maths::Ray const &_ray) const;
	virtual maths::Decimal	Area() const = 0;
	virtual SurfacePoint	SampleSurface(maths::Vec2f const &_ksi) const = 0;
	virtual maths::Decimal	SurfacePdf(SurfaceInteraction const &_origin) const;
	virtual SurfacePoint	SampleVisibleSurface(SurfaceInteraction const &_origin,
												 maths::Vec2f const &_ksi) const;
	virtual maths::Decimal	VisibleSurfacePdf(SurfaceInteraction const &_origin,
											  maths::Vec3f const &_wi) const;
	virtual maths::Bounds3f	ObjectBounds() const = 0;
	virtual maths::Bounds3f	WorldBounds() const;
public:
	// REFACTOR: Maybe Primitive should have a transform instead of shape
	maths::Transform const	&world_transform;
	bool const				swaps_handedness;
	bool const				flip_normals;
};

} // namespace raytracer


#endif // __YS_SHAPE_HPP__

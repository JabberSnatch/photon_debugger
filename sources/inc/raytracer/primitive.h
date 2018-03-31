#pragma once
#ifndef __YS_PRIMITIVE_HPP__
#define __YS_PRIMITIVE_HPP__

#include "maths/maths.h"
#include "raytracer/raytracer.h"


namespace raytracer {

class Primitive
{
public:
	virtual ~Primitive() = default;
	virtual bool	Intersect(maths::Ray &_ray, SurfaceInteraction &_hit_info) const = 0;
	virtual bool	DoesIntersect(maths::Ray const &_ray) const = 0;
	virtual maths::Bounds3f	WorldBounds() const = 0;
};


class GeometryPrimitive final :
	public Primitive
{
public:
	GeometryPrimitive(Shape const &_shape);
	bool	Intersect(maths::Ray &_ray, SurfaceInteraction &_hit_info) const override;
	bool	DoesIntersect(maths::Ray const &_ray) const override;
	maths::Bounds3f	WorldBounds() const override;
private:
	Shape const	&shape_;
};


} // namespace raytracer


#endif // __YS_PRIMITIVE_HPP__

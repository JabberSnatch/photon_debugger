#include "raytracer/primitive.h"

#include "maths/ray.h"
#include "raytracer/shape.h"
#include "raytracer/surface_interaction.h"

namespace raytracer {


GeometryPrimitive::GeometryPrimitive(Shape const &_shape) :
	shape_{ _shape }
{}

bool
GeometryPrimitive::Intersect(maths::Ray &_ray, SurfaceInteraction &_hit_info) const
{
	maths::Decimal	t_hit = maths::infinity<maths::Decimal>;

	if (!shape_.Intersect(_ray, t_hit, _hit_info))
		return false;

	_ray.tMax = t_hit;
	_hit_info.primitive = this;

	return true;
}

bool
GeometryPrimitive::DoesIntersect(maths::Ray const &_ray) const
{
	return shape_.DoesIntersect(_ray);
}

maths::Bounds3f
GeometryPrimitive::WorldBounds() const
{
	return shape_.WorldBounds();
}


} // namespace raytracer


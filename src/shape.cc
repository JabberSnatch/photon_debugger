#include "shape.h"

#include "transform.h"
#include "bounds.h"
#include "surface_interaction.h"

namespace raytracer
{


Shape::Shape(maths::Transform const &_world_transform, bool _flip_normals) :
	world_transform{ _world_transform },
	swaps_handedness{ _world_transform.SwapsHandedness() },
	flip_normals{ _flip_normals }
{}

maths::Bounds3f
Shape::WorldBounds() const
{
	return world_transform(ObjectBounds());
}

bool
Shape::DoesIntersect(maths::Ray const &_ray) const
{
	maths::Decimal		t;
	SurfaceInteraction	hit_info;
	return Intersect(_ray, t, hit_info);
}


} // namespace raytracer

#include "raytracer/shape.h"

#include "maths/transform.h"
#include "maths/bounds.h"
#include "raytracer/surface_interaction.h"

namespace raytracer
{


Shape::Shape(maths::Transform const &_world_transform, bool _flip_normals) :
	world_transform{ _world_transform },
	swaps_handedness{ _world_transform.SwapsHandedness() },
	flip_normals{ _flip_normals }
{}


bool
Shape::DoesIntersect(maths::Ray const &_ray) const
{
	maths::Decimal		t;
	SurfaceInteraction	hit_info;
	return Intersect(_ray, t, hit_info);
}


maths::Decimal
Shape::SurfacePdf(SurfaceInteraction const &_origin) const
{
	return 1._d / Area();
}


Shape::SurfacePoint
Shape::SampleVisibleSurface(SurfaceInteraction const &_origin, maths::Vec2f const &_ksi) const
{
	return SampleSurface(_ksi);
}


maths::Decimal
Shape::VisibleSurfacePdf(SurfaceInteraction const &_origin, maths::Vec3f const &_wi) const
{
	maths::Decimal result = 0._d;
	maths::Ray const visibility_ray = _origin.SpawnRay(_wi);
	maths::Decimal _t;
	SurfaceInteraction hit_info;
	if (Intersect(visibility_ray, _t, hit_info))
	{
		maths::Decimal const r_sqr = maths::SqrDistance(_origin.position, hit_info.position);
		maths::Decimal const cos_theta = maths::Abs(maths::Dot(hit_info.geometry.normal(), -_wi));
		maths::Decimal const dA = Area();
		result = r_sqr / (cos_theta * dA);
	}
	return result;
}


maths::Bounds3f
Shape::WorldBounds() const
{
	return world_transform(ObjectBounds());
}


} // namespace raytracer

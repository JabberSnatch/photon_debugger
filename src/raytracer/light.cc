#include "raytracer/light.h"


#include "raytracer/shape.h"
#include "raytracer/surface_interaction.h"


namespace raytracer {


Light::Light(maths::Transform const &_world_transform) :
	world_transform_{ _world_transform }
{}


maths::Transform const &
Light::world_transform() const
{
	return world_transform_;
}


AreaLight::AreaLight(maths::Transform const &_world_transform, maths::Vec3f const &_le, raytracer::Shape const &_shape) :
	Light(_world_transform),
	le_{ _le }, shape_{ _shape }
{}


maths::Vec3f
AreaLight::Le() const
{
	return le_;
}


Light::LiSample
AreaLight::Sample(SurfaceInteraction const &_hit_info, maths::Vec2f const &_ksi) const
{
	Shape::SurfacePoint const shape_point =  shape_.SampleVisibleSurface(_hit_info, _ksi);
	maths::Vec3f const wi = maths::Normalized(shape_point.position - _hit_info.position);
	maths::Decimal const sample_probability = shape_.VisibleSurfacePdf(_hit_info, wi);
	return LiSample{ _hit_info.position, shape_point.position, le_, sample_probability };
}


maths::Decimal
AreaLight::Pdf(SurfaceInteraction const &_hit_info, maths::Vec3f const &_wi) const
{
	return shape_.VisibleSurfacePdf(_hit_info, _wi);
}


} // namespace raytracer

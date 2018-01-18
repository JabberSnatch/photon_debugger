#pragma once
#ifndef __YS_LIGHT_HPP__
#define __YS_LIGHT_HPP__


#include "maths/maths.h"
#include "maths/point.h"
#include "maths/vector.h"
#include "raytracer/shape.h"


namespace raytracer
{


class SurfaceInteraction;


class Light
{
public:
	struct LiSample
	{
		maths::Vec3f wi() const
		{
			return maths::Normalized(position - origin);
		}
		maths::Point3f const origin;
		maths::Point3f const position;
		maths::Vec3f const li;
		maths::Decimal const probability;
	};
public:
	Light(maths::Transform const & _world_transform);
	virtual ~Light() = default;
	// _ksi is a 2D [0,1) point used to produce the sample
	virtual maths::Vec3f Le() const = 0;
	virtual LiSample Sample(SurfaceInteraction const &_hit_info, maths::Vec2f const &_ksi) const = 0;
	virtual maths::Decimal Pdf(SurfaceInteraction const &_hit_info, maths::Vec3f const &_wi) const = 0;
	// TODO: remove world_transform
	maths::Transform const &world_transform() const;
private:
	maths::Transform const &world_transform_;
};


class AreaLight final : public Light
{
public:
	AreaLight(maths::Transform const &_world_transform, maths::Vec3f const &_le, raytracer::Shape const &_shape);
	maths::Vec3f Le() const override;
	LiSample Sample(SurfaceInteraction const &_hit_info, maths::Vec2f const &_ksi) const override;
	maths::Decimal Pdf(SurfaceInteraction const &_hit_info, maths::Vec3f const &_wi) const override;
private:
	maths::Vec3f const le_;
	raytracer::Shape const &shape_;
};


} // namespace raytracer


#endif // __YS_LIGHT_HPP__

#pragma once
#ifndef __YS_TRIANGLE_HPP__
#define __YS_TRIANGLE_HPP__

#include "raytracer/shape.h"


namespace raytracer {


class TriangleMeshRawData;


class Triangle : public Shape
{
public:
	Triangle(maths::Transform const &_world_transform, bool _flip_normals,
			 TriangleMeshRawData const &_mesh_data, int32_t _face_index);
	virtual bool Intersect(maths::Ray const &_ray,
						   maths::Decimal &_tHit,
						   SurfaceInteraction &_hit_info) const override;
	virtual bool DoesIntersect(maths::Ray const &_ray) const override;
	virtual maths::Decimal	Area() const override;
	virtual SurfacePoint	SampleSurface(maths::Vec2f const &_ksi) const override;
	virtual maths::Bounds3f	ObjectBounds() const override;
	virtual maths::Bounds3f	WorldBounds() const override;
	maths::Point2f	uv(uint32_t _index) const;
private:
	TriangleMeshRawData const	&mesh_data_;
	int32_t const				*vertex_index_;
};

} // namespace raytracer


#endif // __YS_TRIANGLE_HPP__

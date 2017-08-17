#pragma once
#ifndef __YS_TRIANGLE_HPP__
#define __YS_TRIANGLE_HPP__

#include "raytracer/shape.h"
#include "raytracer/raytracer.h"


namespace raytracer {

class Triangle : public Shape
{
public:
	Triangle(maths::Transform const &_world_transform, bool _flip_normals,
			 TriangleMesh const &_mesh, int32_t _face_index);

	virtual bool Intersect(maths::Ray const &_ray,
						   maths::Decimal &_tHit,
						   SurfaceInteraction &_hit_info) const override;
	virtual bool DoesIntersect(maths::Ray const &_ray) const override;

	virtual maths::Decimal	Area() const override;

	virtual maths::Bounds3f	ObjectBounds() const override;
	virtual maths::Bounds3f	WorldBounds() const override;

	maths::Point2f	uv(uint32_t _index) const;

private:
	TriangleMesh const	&mesh_;
	int32_t const		*vertex_index_;
};

} // namespace raytracer


#endif // __YS_TRIANGLE_HPP__

#pragma once
#ifndef __YS_TRIANGLE_MESH_HPP__
#define __YS_TRIANGLE_MESH_HPP__

#include <string>
#include <vector>

#include "maths/maths.h"
#include "maths/bounds.h"
#include "raytracer/shape.h"
#include "raytracer/triangle_mesh_data.h"


namespace core {
class MemoryRegion;
}


namespace raytracer {


template <typename InstancingPolicy>
class TriangleMesh : public Shape
{
private:
	using TriangleContainer_t = std::vector<Triangle const*>;
public:
	TriangleMesh(maths::Transform const &_world_transform,
				 bool _flip_normals,
				 TriangleMeshData const &_data);
	TriangleMesh(maths::Transform const &_world_transform,
				 bool _flip_normals,
				 TriangleMesh const &_sibling_instance);
	virtual bool Intersect(maths::Ray const &_ray,
						   maths::Decimal &_tHit,
						   SurfaceInteraction &_hit_info) const override;
	virtual bool DoesIntersect(maths::Ray const &_ray) const override;
	virtual maths::Decimal	Area() const override;
	virtual SurfacePoint	SampleSurface(maths::Vec2f const &_ksi) const override;
	virtual maths::Bounds3f	ObjectBounds() const override;
	virtual maths::Bounds3f	WorldBounds() const override;
private:
	TriangleMeshData const	&data_;
};


} // namespace raytracer


#endif // __YS_TRIANGLE_MESH_HPP__

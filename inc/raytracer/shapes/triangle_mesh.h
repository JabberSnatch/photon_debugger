#pragma once
#ifndef __YS_TRIANGLE_MESH_HPP__
#define __YS_TRIANGLE_MESH_HPP__

#include <string>
#include <vector>

#include "core/memory_region.h"
#include "maths/maths.h"
#include "maths/bounds.h"
#include "raytracer/shape.h"
#include "raytracer/bvh_accelerator.h"


namespace core {
class MemoryRegion;
}


namespace raytracer {

class TriangleMeshData;


class TriangleMesh : public Shape
{
private:
	using TriangleContainer_t = std::vector<Triangle const*>;
public:
	TriangleMesh(maths::Transform const &_world_transform,
				 bool _flip_normals,
				 TriangleMeshData const &_data,
				 core::MemoryRegion &_mem_region);
	virtual bool Intersect(maths::Ray const &_ray,
						   maths::Decimal &_tHit,
						   SurfaceInteraction &_hit_info) const override;
	virtual bool DoesIntersect(maths::Ray const &_ray) const override;
	virtual maths::Decimal	Area() const override;
	virtual SurfacePoint	SampleSurface(maths::Vec2f const &_ksi) const override;
	virtual maths::Bounds3f	ObjectBounds() const override;
	virtual maths::Bounds3f	WorldBounds() const override;
private:
	static constexpr uint32_t kBvhNodeSize = 10u;
	static TriangleContainer_t MakeTriangles_(
		maths::Transform const &_world_transform,
		bool _flip_normals,
		TriangleMeshData const &_data,
		core::MemoryRegion &_mem_region);
	static BvhAccelerator::PrimitiveArray_t MakePrimitives_(TriangleContainer_t const &_triangles,
															core::MemoryRegion &_mem_region);
private:
	maths::Bounds3f		world_bounds_;
	TriangleContainer_t	triangles_;
	BvhAccelerator		bvh_;
	core::MemoryRegion	mem_region_;
};

TriangleMesh	*ReadTriangleMeshFromFile(std::string const &_path,
										  maths::Transform const &_world_transform,
										  bool _flip_normals,
										  core::MemoryRegion &_mem_region);

} // namespace raytracer


#endif // __YS_TRIANGLE_MESH_HPP__

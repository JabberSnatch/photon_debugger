#pragma once
#ifndef __YS_TRIANGLE_MESH_DATA_HPP__
#define __YS_TRIANGLE_MESH_DATA_HPP__


#include <vector>
#include "maths/bounds.h"
#include "maths/point.h"
#include "maths/vector.h"
#include "core/memory_region.h"
#include "raytracer/bvh_accelerator.h"

namespace raytracer {

class Triangle;

class TriangleMeshRawData
{
public:
	using IndicesContainer_t = std::vector<int32_t>;
	using VerticesContainer_t = std::vector<maths::Point3f>;
	using NormalsContainer_t = std::vector<maths::Norm3f>;
	using TangentsContainer_t = std::vector<maths::Vec3f>;
	using UvsContainer_t = std::vector<maths::Point2f>;
public:
	static int32_t	IndexOffset(int32_t _face_index) { return 3 * _face_index; }
public:
	TriangleMeshRawData(int32_t _triangle_count,
						IndicesContainer_t const &_indices,
						VerticesContainer_t const &_vertices,
						NormalsContainer_t const *_normals = nullptr,
						TangentsContainer_t const *_tangents = nullptr,
						UvsContainer_t const *_uvs = nullptr);
	bool		has_normals() const { return normals.size() > 0; }
	bool		has_tangents() const { return tangents.size() > 0; }
	bool		has_uvs() const { return uvs.size() > 0; }
public:
	int32_t					triangle_count;
	IndicesContainer_t		indices;
	VerticesContainer_t		vertices;
	NormalsContainer_t		normals;
	TangentsContainer_t		tangents;
	UvsContainer_t			uvs;
	maths::Bounds3f 		bounds;
};


// TODO: move back to triangle_mesh.h from here to the end of namespace
struct InstancingPolicyClass
{
	struct SharedSource
	{};
	struct Transformed
	{
		static TriangleMeshRawData const &GetRawData(TriangleMeshRawData const &_base,
													 maths::Transform const &_world_transform,
													 core::MemoryRegion &_mem_region);
	};
};

class TriangleMeshData
{
public:
	using TriangleContainer_t = std::vector<Triangle const*>;
public:
	TriangleMeshData(maths::Transform const &_world_transform,
					 bool _flip_normals,
					 TriangleMeshRawData const &_mesh_raw_data,
					 InstancingPolicyClass::SharedSource const &);
	TriangleMeshData(maths::Transform const &_world_transform,
					 bool _flip_normals,
					 TriangleMeshRawData const &_mesh_raw_data,
					 InstancingPolicyClass::Transformed const &);
	maths::Bounds3f const &bounds() const { return data_source_.bounds; }
	TriangleContainer_t const &triangles() const { return triangles_; }
	BvhAccelerator const &bvh() const { return bvh_; }
private:
	static constexpr uint32_t kBvhNodeSize = 10u;
	static TriangleContainer_t MakeTriangles_(
		maths::Transform const &_world_transform,
		bool _flip_normals,
		TriangleMeshRawData const &_raw_data,
		core::MemoryRegion &_mem_region);
	static BvhAccelerator::PrimitiveArray_t MakePrimitives_(TriangleContainer_t const &_triangles,
															core::MemoryRegion &_mem_region);
private:
	core::MemoryRegion			mem_region_;
	TriangleMeshRawData const	&data_source_;
	TriangleContainer_t 		triangles_;
	BvhAccelerator				bvh_;
};


} // namespace raytracer


#endif // __YS_TRIANGLE_MESH_DATA_HPP__

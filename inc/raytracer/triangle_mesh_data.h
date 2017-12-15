#pragma once
#ifndef __YS_TRIANGLE_MESH_DATA_HPP__
#define __YS_TRIANGLE_MESH_DATA_HPP__


#include <vector>
#include "maths/bounds.h"
#include "maths/point.h"
#include "maths/vector.h"


namespace raytracer {


class TriangleMeshData
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
	TriangleMeshData(maths::Transform const &_world_transform,
					 int32_t _triangle_count,
					 IndicesContainer_t const &_indices,
					 VerticesContainer_t const &_vertices,
					 NormalsContainer_t const *_normals = nullptr,
					 TangentsContainer_t const *_tangents = nullptr,
					 UvsContainer_t const *_uvs = nullptr);
	bool		has_normals() const { return normals.size() > 0; }
	bool		has_tangents() const { return tangents.size() > 0; }
	bool		has_uvs() const { return uvs.size() > 0; }
public:
	maths::Bounds3f			world_bounds;
	int32_t					triangle_count;
	IndicesContainer_t		indices;
	VerticesContainer_t		vertices;
	NormalsContainer_t		normals;
	TangentsContainer_t		tangents;
	UvsContainer_t			uvs;
};


} // namespace raytracer


#endif // __YS_TRIANGLE_MESH_DATA_HPP__

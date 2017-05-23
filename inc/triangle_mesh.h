#pragma once
#ifndef __YS_TRIANGLE_MESH_HPP__
#define __YS_TRIANGLE_MESH_HPP__

#include <vector>
#include "maths.h"


namespace raytracer {

struct TriangleMesh
{
	TriangleMesh(maths::Transform const &_world_transform, int32_t _triangle_count,
				 std::vector<int32_t> const &_indices,
				 std::vector<maths::Point3f> const &_vertices,
				 std::vector<maths::Norm3f> const *_normals = nullptr,
				 std::vector<maths::Vec3f> const *_tangents = nullptr,
				 std::vector<maths::Point2f> const *_uv = nullptr);

	static int32_t	IndexOffset(int32_t _face_index) { return 3 * _face_index; }

	inline bool		has_normals() const { return normals.size() > 0; }
	inline bool		has_tangents() const { return tangents.size() > 0; }
	inline bool		has_uv() const { return uv.size() > 0; }

	int32_t							triangle_count;
	std::vector<int32_t>			indices;
	std::vector<maths::Point3f>		vertices;
	std::vector<maths::Norm3f>		normals;
	std::vector<maths::Vec3f>		tangents;
	std::vector<maths::Point2f>		uv;
};

} // namespace raytracer


#endif // __YS_TRIANGLE_MESH_HPP__

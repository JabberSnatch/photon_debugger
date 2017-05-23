#include "triangle_mesh.h"

#include "globals.h"
#include "common_macros.h"
#include "profiler.h"
#include "transform.h"
#include "vector.h"
#include "point.h"

namespace raytracer {

TriangleMesh::TriangleMesh(maths::Transform const &_world_transform, int32_t _triangle_count,
			 std::vector<int32_t> const &_indices,
			 std::vector<maths::Point3f> const &_vertices,
			 std::vector<maths::Norm3f> const *_normals,
			 std::vector<maths::Vec3f> const *_tangents,
			 std::vector<maths::Point2f> const *_uv) :
	triangle_count{ _triangle_count }, indices{ _indices }
{
	TIMED_SCOPE(TriangleMesh_ctor);

	YS_ASSERT(3 * triangle_count == indices.size());

	vertices.reserve(_vertices.size());
	for (auto &&vertex : _vertices)
		vertices.emplace_back(_world_transform(vertex));
	
	if (_normals != nullptr)
	{
		YS_ASSERT(_normals->size() == vertices.size());
		normals.insert(normals.begin(), _normals->begin(), _normals->end());
	}

	if (_tangents != nullptr)
	{
		YS_ASSERT(_tangents->size() == vertices.size());
		tangents.insert(tangents.begin(), _tangents->begin(), _tangents->end());
	}

	if (_uv != nullptr)
	{
		YS_ASSERT(_uv->size() == vertices.size());
		uv.insert(uv.begin(), _uv->begin(), _uv->end());
	}
}

} // namespace raytracer

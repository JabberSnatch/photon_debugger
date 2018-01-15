#include "raytracer/triangle_mesh_data.h"

#include <functional>
#include <sstream>

#include "maths/transform.h"


namespace raytracer {


TriangleMeshData::TriangleMeshData(maths::Transform const &_world_transform,
								   int32_t _triangle_count,
								   IndicesContainer_t const &_indices,
								   VerticesContainer_t const &_vertices,
								   NormalsContainer_t const *_normals,
								   TangentsContainer_t const *_tangents,
								   UvsContainer_t const *_uvs) :
	triangle_count{ _triangle_count }, indices{ _indices }, world_bounds{}
{
	TIMED_SCOPE(TriangleMeshData_ctor);

	YS_ASSERT(triangle_count > 0);
	YS_ASSERT(static_cast<size_t>(3 * triangle_count) == indices.size());

	vertices.reserve(_vertices.size());
	std::transform(_vertices.cbegin(), _vertices.cend(),
				   std::back_inserter(vertices), [this, &_world_transform] 
				   (maths::Point3f const &_vertex) {
		maths::Point3f const world_vertex{ _world_transform(_vertex) };
		world_bounds = maths::Union(world_bounds, world_vertex);
		return world_vertex;
	});
	if (_normals != nullptr)
	{
		YS_ASSERT(_normals->size() == vertices.size());
		std::transform(_normals->cbegin(), _normals->cend(),
					   std::back_inserter(normals), [this, &_world_transform]
					   (maths::Norm3f const &_normal) {
			return _world_transform(_normal);
		});
	}
	if (_tangents != nullptr)
	{
		YS_ASSERT(_tangents->size() == vertices.size());
		std::transform(_tangents->cbegin(), _tangents->cend(),
					   std::back_inserter(tangents), [this, &_world_transform]
					   (maths::Vec3f const &_tangent) {
			return _world_transform(_tangent);
		});
	}
	if (_uvs != nullptr)
	{
		YS_ASSERT(_uvs->size() == vertices.size());
		uvs.insert(uvs.begin(), _uvs->begin(), _uvs->end());
	}

	std::stringstream message_stream{};
	message_stream << "WORLD_BOUNDS : " << std::endl <<
		"	" << world_bounds.min.x << "; " << world_bounds.min.y << "; " << world_bounds.min.z <<
		std::endl <<
		"	" << world_bounds.max.x << "; " << world_bounds.max.y << "; " << world_bounds.max.z <<
		std::endl;
	LOG_INFO(tools::kChannelGeneral, message_stream.str());
}

} // namespace raytracer

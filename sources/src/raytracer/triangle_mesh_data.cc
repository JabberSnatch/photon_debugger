#include "raytracer/triangle_mesh_data.h"

#include <functional>
#include <numeric>
#include <sstream>

#include "maths/transform.h"
#include "raytracer/primitive.h"
#include "raytracer/shapes/triangle.h"


namespace raytracer {


TriangleMeshRawData::TriangleMeshRawData(int32_t _triangle_count,
										 IndicesContainer_t const &_indices,
										 VerticesContainer_t const &_vertices,
										 NormalsContainer_t const *_normals,
										 TangentsContainer_t const *_tangents,
										 UvsContainer_t const *_uvs) :
	triangle_count{ _triangle_count },
	indices{ _indices },
	vertices{ _vertices },
	normals{ _normals ? *_normals : NormalsContainer_t{} },
	tangents{ _tangents ? *_tangents : TangentsContainer_t{} },
	uvs{ _uvs ? *_uvs : UvsContainer_t{} },
	bounds{ std::accumulate(vertices.cbegin(), vertices.cend(), maths::Bounds3f{},
							[](maths::Bounds3f const &_acc, maths::Point3f const &_vertex) {
								return maths::Union(_acc, _vertex);
							}) }
{
	YS_ASSERT(triangle_count > 0);
	YS_ASSERT(static_cast<size_t>(3 * triangle_count) == indices.size());
	if(!vertices.empty())
	{
		YS_ASSERT(normals.empty() || vertices.size() == normals.size());
		YS_ASSERT(tangents.empty() || vertices.size() == tangents.size());
		YS_ASSERT(uvs.empty() || vertices.size() == uvs.size());
	}
	else
	{
		YS_ASSERT(false);
	}
}


TriangleMeshRawData const &
InstancingPolicyClass::Transformed::GetRawData(TriangleMeshRawData const &_base,
											   maths::Transform const &_world_transform,
											   core::MemoryRegion &_mem_region)
{
	TriangleMeshRawData::VerticesContainer_t transformed_vertices{};
	TriangleMeshRawData::NormalsContainer_t transformed_normals{};
	TriangleMeshRawData::TangentsContainer_t transformed_tangents{};
	TriangleMeshRawData::UvsContainer_t transformed_uvs{};

	transformed_vertices.reserve(_base.vertices.size());
	std::transform(_base.vertices.cbegin(), _base.vertices.cend(),
				   std::back_inserter(transformed_vertices), [&_world_transform]
				   (maths::Point3f const &_vertex) {
					   return _world_transform(_vertex);
				   });
	std::transform(_base.normals.cbegin(), _base.normals.cend(),
				   std::back_inserter(transformed_normals), [&_world_transform]
				   (maths::Norm3f const &_normal) {
					   return _world_transform(_normal);
				   });
	std::transform(_base.tangents.cbegin(), _base.tangents.cend(),
				   std::back_inserter(transformed_tangents), [&_world_transform]
				   (maths::Vec3f const &_tangent) {
					   return _world_transform(_tangent);
				   });

	// TODO: uvs are not transformed by this operation, they shouldn't be copied but shared
	TriangleMeshRawData const* const transformed_data =
		new (_mem_region) TriangleMeshRawData{ _base.triangle_count,
											   _base.indices,
											   transformed_vertices,
											   &transformed_normals,
											   &transformed_tangents,
											   &_base.uvs };

	return *transformed_data;
}


TriangleMeshData::TriangleMeshData(maths::Transform const &_world_transform,
								   bool _flip_normals,
								   TriangleMeshRawData const &_mesh_raw_data,
								   InstancingPolicyClass::SharedSource const &) :
	mem_region_{},
	data_source_{ _mesh_raw_data },
	triangles_{ MakeTriangles_(_world_transform, _flip_normals, data_source_, mem_region_) },
	bvh_{ MakePrimitives_(triangles_, mem_region_), kBvhNodeSize }
{}


TriangleMeshData::TriangleMeshData(maths::Transform const &_world_transform,
								   bool _flip_normals,
								   TriangleMeshRawData const &_mesh_raw_data,
								   InstancingPolicyClass::Transformed const &) :
	mem_region_{},
	data_source_{ InstancingPolicyClass::Transformed::GetRawData(_mesh_raw_data,
																 _world_transform,
																 mem_region_) },
	triangles_{ MakeTriangles_(_world_transform, _flip_normals, data_source_, mem_region_) },
	bvh_{ MakePrimitives_(triangles_, mem_region_), kBvhNodeSize }
{}


TriangleMeshData::TriangleContainer_t
TriangleMeshData::MakeTriangles_(maths::Transform const &_world_transform,
								 bool _flip_normals,
								 TriangleMeshRawData const &_raw_data,
								 core::MemoryRegion &_mem_region)
{
	TriangleContainer_t result{};
	result.reserve(_raw_data.triangle_count);
	for (int32_t face_index = 0; face_index < _raw_data.triangle_count; ++face_index)
	{
		result.emplace_back(new (_mem_region)
							raytracer::Triangle(_world_transform,
												_flip_normals,
												_raw_data,
												face_index));
	}
	return result;
}


BvhAccelerator::PrimitiveArray_t
TriangleMeshData::MakePrimitives_(TriangleContainer_t const &_triangles,
								  core::MemoryRegion &_mem_region)
{
	BvhAccelerator::PrimitiveArray_t result{};
	result.reserve(_triangles.size());
	std::transform(_triangles.cbegin(), _triangles.cend(),
				   std::back_inserter(result), [&_mem_region] (Triangle const* const _triangle) {
		return new (_mem_region) GeometryPrimitive(*static_cast<Shape const*>(_triangle));
	});
	return result;
}


} // namespace raytracer

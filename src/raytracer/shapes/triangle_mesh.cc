#include "raytracer/shapes/triangle_mesh.h"

#include <numeric>

#include "assimp/Importer.hpp"
#include "assimp/postprocess.h"
#include "assimp/scene.h"

#include "maths/ray.h"
#include "maths/transform.h"
#include "raytracer/primitive.h"
#include "raytracer/shapes/triangle.h"
#include "raytracer/triangle_mesh_data.h"


namespace raytracer {


TriangleMesh::TriangleMesh(maths::Transform const &_world_transform,
						   bool _flip_normals,
						   TriangleMeshData const &_data,
						   core::MemoryRegion &_mem_region) :
	Shape(_world_transform, _flip_normals),
	world_bounds_(_data.world_bounds),
	triangles_(MakeTriangles_(_world_transform, _flip_normals, _data, _mem_region)),
	bvh_(MakePrimitives_(triangles_, _mem_region), kBvhNodeSize)
{}


bool
TriangleMesh::Intersect(maths::Ray const &_ray,
						maths::Decimal &_tHit,
						SurfaceInteraction &_hit_info) const
{
	maths::Ray bvh_ray{ _ray };
	bool const result = bvh_.Intersect(bvh_ray, _hit_info);
	_tHit = (result) ? bvh_ray.tMax : _tHit;
	return result;
}


bool
TriangleMesh::DoesIntersect(maths::Ray const &_ray) const
{
	return bvh_.DoesIntersect(_ray);
}


maths::Decimal
TriangleMesh::Area() const
{
	return std::accumulate(triangles_.cbegin(), triangles_.cend(), 0._d,
						   [](maths::Decimal const &_acc, Triangle const* const _triangle) {
		return _acc + _triangle->Area();
	});
}


Shape::SurfacePoint
TriangleMesh::SampleSurface(maths::Vec2f const &_ksi) const
{
	YS_ASSERT(false);
	return Shape::SurfacePoint();
}


maths::Bounds3f
TriangleMesh::ObjectBounds() const
{
	return world_transform(world_bounds_, maths::Transform::kInverse);
}


maths::Bounds3f
TriangleMesh::WorldBounds() const
{
	return world_bounds_;
}


TriangleMesh::TriangleContainer_t
TriangleMesh::MakeTriangles_(maths::Transform const &_world_transform,
							 bool _flip_normals,
							 TriangleMeshData const &_data,
							 core::MemoryRegion &_mem_region)
{
	TriangleContainer_t result{};
	result.reserve(_data.triangle_count);
	for (int32_t face_index = 0; face_index < _data.triangle_count; ++face_index)
	{
		result.emplace_back(new (_mem_region)
							raytracer::Triangle(_world_transform,
												_flip_normals,
												_data,
												face_index));
	}
	return result;
}


BvhAccelerator::PrimitiveArray_t
TriangleMesh::MakePrimitives_(TriangleContainer_t const &_triangles,
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

// TODO: implement instancing
TriangleMesh*
ReadTriangleMeshFromFile(std::string const &_path,
						 maths::Transform const &_world_transform,
						 bool _flip_normals,
						 core::MemoryRegion &_mem_region)
{
	TriangleMesh *result = nullptr;
	uint32_t const		load_flags =
		aiProcess_Triangulate |
		aiProcess_PreTransformVertices |
		aiProcess_JoinIdenticalVertices;
	Assimp::Importer	importer;
	aiScene const		*scene = importer.ReadFile(_path.c_str(), load_flags);
	if (scene != nullptr)
	{
		uint32_t const		mesh_count = scene->mNumMeshes;

		YS_ASSERT(mesh_count > 0u);

		int32_t						out_triangle_count{ 0 };
		std::vector<int32_t>		out_indices;
		std::vector<maths::Point3f>	out_vertices;
		std::vector<maths::Norm3f>	out_normals;

		uint32_t					indices_capacity{ 0u };
		uint32_t					vertices_capacity{ 0u };
		bool						load_normals = true;

		aiMesh *const				*in_meshes = scene->mMeshes;

		for (uint32_t mesh_index = 0; mesh_index < mesh_count; ++mesh_index)
			load_normals &= in_meshes[mesh_index]->HasNormals();

		for (uint32_t mesh_index = 0; mesh_index < mesh_count; ++mesh_index)
		{
			aiMesh const		*in_mesh = in_meshes[mesh_index];
			aiFace const		*in_faces = in_mesh->mFaces;
			aiVector3D const	*in_vertices = in_mesh->mVertices;
			aiVector3D const	*in_normals = in_mesh->mNormals;

			out_triangle_count += in_mesh->mNumFaces;
			indices_capacity += in_mesh->mNumFaces * 3;
			vertices_capacity += in_mesh->mNumVertices;
			out_indices.reserve(indices_capacity);
			out_vertices.reserve(vertices_capacity);

			if (load_normals)
				out_normals.reserve(vertices_capacity);

			for (uint32_t i = 0; i < in_mesh->mNumFaces; ++i)
			{
				aiFace const	&face = in_faces[i];
				YS_ASSERT(face.mNumIndices == 3);
				for (uint32_t j = 0; j < face.mNumIndices; ++j)
				{
					uint32_t const	base_index = static_cast<uint32_t>(face.mIndices[j]);
					uint32_t const	offset_index =
						base_index + static_cast<uint32_t>(out_vertices.size());
					YS_ASSERT(offset_index >= base_index);
					out_indices.emplace_back(offset_index);
				}
			}

			for (uint32_t i = 0; i < in_mesh->mNumVertices; ++i)
			{
				aiVector3D const	&vertex = in_vertices[i];
				out_vertices.emplace_back(vertex.x, vertex.y, vertex.z);
				if (load_normals)
				{
					aiVector3D const	&normal = in_normals[i];
					maths::Decimal const x = std::isnan(normal.x) ? 0._d : normal.x;
					maths::Decimal const y = std::isnan(normal.y) ? 0._d : normal.y;
					maths::Decimal const z = std::isnan(normal.z) ? 0._d : normal.z;
					out_normals.emplace_back(normal.x, normal.y, normal.z);
				}
			}
		}

#ifndef YS_NO_LOGS
		std::string	const	info_message =
			"Loading " + std::to_string(mesh_count) + " submesh from " + _path + ". " +
			std::to_string(out_triangle_count) + " total triangles." +
			((load_normals) ? " Normals provided by file." : "");
		LOG(tools::kChannelGeneral, tools::kLevelInfo, info_message);
#endif
		
		TriangleMeshData const* const mesh_data = new (_mem_region)
			TriangleMeshData(_world_transform,
							 out_triangle_count,
							 out_indices,
							 out_vertices,
							 (load_normals ? &out_normals : nullptr));
		result = new (_mem_region) TriangleMesh(_world_transform,
												_flip_normals,
												*mesh_data,
												_mem_region);
	}
	else
	{
		LOG_WARNING(tools::kChannelGeneral, "could not load file " + _path);
	}
	return result;
}


} // namespace raytracer

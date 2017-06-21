#include "triangle_mesh.h"

#include "assimp/Importer.hpp"
#include "assimp/postprocess.h"
#include "assimp/scene.h"

#include "globals.h"
#include "common_macros.h"
#include "profiler.h"
#include "logger.h"
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
	triangle_count{ _triangle_count }, indices{ _indices }, world_bounds{}
{
	TIMED_SCOPE(TriangleMesh_ctor);

	YS_ASSERT(3 * triangle_count == indices.size());

	vertices.reserve(_vertices.size());
	for (auto &&vertex : _vertices)
	{
		maths::Point3f const	world_vertex{ _world_transform(vertex) };
		vertices.emplace_back(world_vertex);
		world_bounds = maths::Union(world_bounds, world_vertex);
	}
	
	if (_normals != nullptr)
	{
		YS_ASSERT(_normals->size() == vertices.size());
		for (auto &&normal : *_normals)
			normals.emplace_back(_world_transform(normal));
		//normals.insert(normals.begin(), _normals->begin(), _normals->end());
	}

	if (_tangents != nullptr)
	{
		YS_ASSERT(_tangents->size() == vertices.size());
		for (auto &&tangent : *_tangents)
			tangents.emplace_back(_world_transform(tangent));
		//tangents.insert(tangents.begin(), _tangents->begin(), _tangents->end());
	}

	if (_uv != nullptr)
	{
		YS_ASSERT(_uv->size() == vertices.size());
		uv.insert(uv.begin(), _uv->begin(), _uv->end());
	}
}

TriangleMesh
ReadTriangleMeshFromFile(std::string const &_path,
						 maths::Transform const &_world_transform)
{
	uint32_t const		load_flags =
		aiProcess_Triangulate |
		aiProcess_PreTransformVertices | 
		aiProcess_JoinIdenticalVertices;
	Assimp::Importer	importer;
	aiScene const		*scene = importer.ReadFile(_path.c_str(), load_flags);
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
				YS_ASSERT(face.mIndices[j] >= 0);
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

	return TriangleMesh(_world_transform, out_triangle_count, out_indices, out_vertices,
						(load_normals ? &out_normals : nullptr));
}

} // namespace raytracer

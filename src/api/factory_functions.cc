#include "api/factory_functions.h"

#include <random>

#include "assimp/Importer.hpp"
#include "assimp/postprocess.h"
#include "assimp/scene.h"
#include "boost/filesystem.hpp"
#include "boost/numeric/conversion/cast.hpp"

#include "api/param_set.h"
#include "api/resource_context.h"
#include "common_macros.h"
#include "core/logger.h"
#include "maths/transform.h"
#include "raytracer/camera.h"
#include "raytracer/film.h"
#include "raytracer/integrator.h"
#include "raytracer/integrators/direct_lighting_integrator.h"
#include "raytracer/light.h"
#include "raytracer/sampler.h"
#include "raytracer/samplers/random_sampler.h"
#include "raytracer/samplers/halton_sampler.h"
#include "raytracer/shape.h"
#include "raytracer/shapes/sphere.h"
#include "raytracer/shapes/triangle.h"
#include "raytracer/shapes/triangle_mesh.h"
#include "raytracer/triangle_mesh_data.h"


namespace api {


raytracer::Shape* MakeSphere(api::ResourceContext &_context, api::ParamSet const &_params);
raytracer::Shape* MakeTriangleMesh(api::ResourceContext &_context, api::ParamSet const &_params);

raytracer::Sampler* MakeRandomSampler(api::ResourceContext &_context,
									  api::ParamSet const &_params);
raytracer::Sampler* MakeHaltonSampler(api::ResourceContext &_context,
									  api::ParamSet const &_params);

raytracer::Integrator* MakeNormalIntegrator(api::ResourceContext &_context,
											api::ParamSet const &_params);
raytracer::Integrator* MakeAOIntegrator(api::ResourceContext &_context,
										api::ParamSet const &_params);
raytracer::Integrator* MakeDirectLightingIntegrator(api::ResourceContext &_context,
													api::ParamSet const &_params);

raytracer::Light* MakeAreaLight(api::ResourceContext &_context,
								api::ParamSet const &_params);


raytracer::Camera*
MakeCamera(api::ResourceContext &_context, api::ParamSet const &_params)
{
	maths::Point3f const	position = static_cast<maths::Point3f>(
		_params.FindFloat<3>("position", { 0._d, 0._d, 0._d })
		);
	maths::Point3f const	lookat = static_cast<maths::Point3f>(
		_params.FindFloat<3>("lookat", maths::Vec3f(position) + maths::Vec3f{ 0._d, 1._d, 0._d })
		);
	maths::Vec3f const		up = _params.FindFloat<3>("up", { 0._d, 0._d, 1._d });
	maths::Decimal const	fov = _params.FindFloat("fov", 60._d);
	//
	ResourceContext::ObjectDescriptor const &film_desc =
		_context.GetAnyDescOfType(ResourceContext::ObjectType::kFilm);
	return new (_context.mem_region()) raytracer::Camera{
		position, lookat, up, fov, _context.Fetch<raytracer::Film>(film_desc.unique_id) };
}


raytracer::Film*
MakeFilm(api::ResourceContext &_context, api::ParamSet const &_params)
{
	maths::Vec2i const		resolution = _params.FindInt<2>("resolution", { 800, 600 });
	maths::Decimal const	side = _params.FindFloat("side", .036_d);

	return new (_context.mem_region()) raytracer::Film{ resolution.x, resolution.y, side };
}


raytracer::TriangleMeshRawData*
MakeTriangleMeshRawData(api::ResourceContext &_context, api::ParamSet const &_params)
{
	// TODO: compare this process execution time agains bvh construction
	std::string const path_string = _params.FindString("path", "");

	raytracer::TriangleMeshRawData* result = nullptr;
	uint32_t const		load_flags =
		aiProcess_Triangulate |
		aiProcess_PreTransformVertices |
		aiProcess_JoinIdenticalVertices;
	Assimp::Importer	importer;
	aiScene const		*scene = importer.ReadFile(path_string.c_str(), load_flags);
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
			"Loading " + std::to_string(mesh_count) + " submesh from " + path_string + ". " +
			std::to_string(out_triangle_count) + " total triangles." +
			((load_normals) ? " Normals provided by file." : "");
		LOG(tools::kChannelGeneral, tools::kLevelInfo, info_message);
#endif

		result = new (_context.mem_region()) raytracer::TriangleMeshRawData{
			out_triangle_count, out_indices, out_vertices,
			(load_normals ? &out_normals : nullptr) };
	}
	else
	{
		LOG_WARNING(tools::kChannelGeneral, "could not load file " + path_string);
	}
	return result;
}


template <typename CallbackType>
CallbackType const &
LookupObjectFuncImpl_(std::string const &_id, NamedCallbackContainer_t<CallbackType> const &_container)
{
	using ContainerType = NamedCallbackContainer_t<CallbackType>;
	typename ContainerType::const_iterator const msccit = _container.find(_id);
	YS_ASSERT(msccit != _container.cend());
	return msccit->second;
}

ShapeCallbackContainer_t const &
shape_callbacks()
{
	static ShapeCallbackContainer_t const callbacks{
		{ "sphere", &MakeSphere },
		{ "triangle_mesh", &MakeTriangleMesh },
	};
	return callbacks;
}
MakeShapeCallback_t const &
LookupShapeFunc(std::string const &_id)
{
	return LookupObjectFuncImpl_(_id, shape_callbacks());
}


SamplerCallbackContainer_t const &
sampler_callbacks()
{
	static SamplerCallbackContainer_t const callbacks{
		{ "random", &MakeRandomSampler },
		{ "halton", &MakeHaltonSampler }
	};
	return callbacks;
}
MakeSamplerCallback_t const &
LookupSamplerFunc(std::string const &_id)
{
	return LookupObjectFuncImpl_(_id, sampler_callbacks());
}


IntegratorCallbackContainer_t const &
integrator_callbacks()
{
	static IntegratorCallbackContainer_t const callbacks{
		{ "ambient_occlusion", &MakeAOIntegrator },
		{ "normal", &MakeNormalIntegrator },
		{ "direct_lighting", &MakeDirectLightingIntegrator },
	};
	return callbacks;
}
MakeIntegratorCallback_t const &
LookupIntegratorFunc(std::string const &_id)
{
	return LookupObjectFuncImpl_(_id, integrator_callbacks());
}


LightCallbackContainer_t const &
light_callbacks()
{
	static LightCallbackContainer_t const callbacks{
		{ "area_light", &MakeAreaLight },
	};
	return callbacks;
}
MakeLightCallback_t const &
LookupLightFunc(std::string const &_id)
{
	return LookupObjectFuncImpl_(_id, light_callbacks());
}


raytracer::Shape*
MakeSphere(api::ResourceContext &_context, api::ParamSet const &_params)
{
	maths::Transform const	&world_transform = _params.FindTransform("world_transform", maths::Transform::Identity());
	bool const				flip_normals = _params.FindBool("flip_normals", false);
	maths::Decimal const	radius = _params.FindFloat("radius", 1._d);
	maths::Decimal const	z_min = _params.FindFloat("z_min", -radius);
	maths::Decimal const	z_max = _params.FindFloat("z_max", radius);
	maths::Decimal const	phi_max = _params.FindFloat("phi_max", 360._d);
	raytracer::Shape *const sphere_shape = new (_context.mem_region())
		raytracer::Sphere{ world_transform, flip_normals, radius, z_min, z_max, phi_max };
	return sphere_shape;
}
raytracer::Shape*
MakeTriangleMesh(api::ResourceContext &_context, api::ParamSet const &_params)
{
	raytracer::Shape* result = nullptr;
	maths::Transform const	&world_transform = _params.FindTransform("world_transform", maths::Transform::Identity());
	bool const				flip_normals = _params.FindBool("flip_normals", false);
	std::string const		path_string = _params.FindString("path", "");
	if (path_string != "")
	{
		boost::filesystem::path path(path_string);
		if (path.is_relative())
		{
			boost::filesystem::path workdir(_context.workdir());
			path = workdir / path;
			YS_ASSERT(path.is_absolute());
		}

		if (boost::filesystem::exists(path))
		{
			// if rawdata descriptor is missing, push it
			if (!_context.IsUniqueIdFree(path_string))
			{
				YS_ASSERT(_context.GetDesc(path_string).type_id ==
						  ResourceContext::ObjectType::kTriangleMeshRawData);
			}
			else
			{
				ParamSet *const params = new (_context.mem_region()) ParamSet();
				params->PushString("path", path_string);
				_context.PushDescriptor(path_string,
										ResourceContext::ObjectType::kTriangleMeshRawData,
										*params);
			}
			// pick instancing policy
				// count trianglemeshes sharing the same path_string
				// 1 => transform, n => instancing
			ResourceContext::ObjectDescriptorContainer_t const shape_descs =
				_context.GetAllDescsOfType(ResourceContext::ObjectType::kShape);
			uint32_t const instance_count = boost::numeric_cast<uint32_t>(std::count_if(
				shape_descs.cbegin(), shape_descs.cend(),
				[&path_string](ResourceContext::ObjectDescriptor const *_desc) {
					return _desc->param_set.FindString("path", "") == path_string;
				}));
			YS_ASSERT(instance_count != 0u);
			// switch on instancing policy class
			// case transform
				// fetch rawdata
				// build TriangleMeshData, and then TriangleMesh
			// case sharedsource
				// try to find a preexisting TriangleMesh that shares the same rawdata
				// if it exists
					// build new TriangleMesh using found "sibling"'s TriangleMeshData
				// else
					// fallback on fetching rawdata, then building TriangleMeshData and TriangleMesh
			if (instance_count == 1u) // transform
			{
				using InstancingPolicy = raytracer::InstancingPolicyClass::Transformed;
				using LocalTriangleMesh = raytracer::TriangleMesh<InstancingPolicy>;
				raytracer::TriangleMeshRawData const &raw_data =
					_context.Fetch<raytracer::TriangleMeshRawData>(path_string);
				raytracer::TriangleMeshData const *const mesh_data =
					new (_context.mem_region()) raytracer::TriangleMeshData{
					world_transform,
					flip_normals,
					raw_data,
					InstancingPolicy{} };
				result = new (_context.mem_region()) LocalTriangleMesh{ world_transform,
																		flip_normals,
																		*mesh_data };
			}
			else // sharedsource
			{
				using InstancingPolicy = raytracer::InstancingPolicyClass::SharedSource;
				using LocalTriangleMesh = raytracer::TriangleMesh<InstancingPolicy>;
				using DescConstIt = ResourceContext::ObjectDescriptorContainer_t::const_iterator;
				DescConstIt const odcit = std::find_if(
					shape_descs.cbegin(), shape_descs.cend(),
					[&path_string, &_context](ResourceContext::ObjectDescriptor const *_desc) {
						std::string const path_param = _desc->param_set.FindString("path", "");
						return (path_param == path_string)
							&& (_context.HasInstance(_desc->unique_id));
					});
				if (odcit != shape_descs.cend())
				{
					ResourceContext::ObjectDescriptor const &desc = **odcit;
					LocalTriangleMesh const &sibling = dynamic_cast<LocalTriangleMesh const &>(
						_context.GetInstance<raytracer::Shape>(desc.unique_id));
					result = new (_context.mem_region()) LocalTriangleMesh{ world_transform,
																			flip_normals,
																			sibling };
				}
				else
				{
					raytracer::TriangleMeshRawData const &raw_data =
						_context.Fetch<raytracer::TriangleMeshRawData>(path_string);
					raytracer::TriangleMeshData const *const mesh_data =
						new (_context.mem_region()) raytracer::TriangleMeshData{
						world_transform,
						flip_normals,
						raw_data,
						InstancingPolicy{} };
					result = new (_context.mem_region()) LocalTriangleMesh{ world_transform,
																			flip_normals,
																			*mesh_data };
				}
			}
			if (!result)
			{
				LOG_ERROR(tools::kChannelGeneral,
						  "Failed to load a triangle mesh at path : " + path.generic_string());
			}
		}
		else
		{
			LOG_ERROR(tools::kChannelGeneral,
					  "No file found to create a triangle mesh at path : " + path.generic_string());
		}
	}
	else
	{
		LOG_WARNING(tools::kChannelGeneral, "Declared a triangle mesh without a path");
	}
	return result;
}


raytracer::Sampler*
MakeRandomSampler(api::ResourceContext &_context, api::ParamSet const &_params)
{
	uint64_t const	seed = _params.FindUint("seed", std::random_device()());
	uint64_t const	samples_per_pixel = _params.FindUint("sample_count", 1u);
	uint64_t const	dimensions_per_sample = _params.FindUint("dimension_count", 5u);
	raytracer::Sampler *const random_sampler = new (_context.mem_region())
		raytracer::RandomSampler{ seed, samples_per_pixel, dimensions_per_sample };
	return random_sampler;
}
raytracer::Sampler*
MakeHaltonSampler(api::ResourceContext &_context, api::ParamSet const &_params)
{
	uint64_t const	seed = _params.FindUint("seed", std::random_device()());
	uint64_t const	samples_per_pixel = _params.FindUint("sample_count", 1u);
	uint64_t const	dimensions_per_sample = _params.FindUint("dimension_count", 5u);
	maths::Vec2u const	tile_resolution = _params.FindUint<2>("tile_resolution", { 128u, 128u });
	raytracer::Sampler *const halton_sampler = new (_context.mem_region())
		raytracer::HaltonSampler{ seed, samples_per_pixel, dimensions_per_sample, tile_resolution };
	return halton_sampler;
}


namespace {
template <typename T>
T&
FetchForIDOrAny(std::string const &_object_id, api::ResourceContext &_context)
{
	T* result = nullptr;
	if (_object_id.empty())
	{
		api::ResourceContext::ObjectDescriptor const &object_desc =
			_context.GetAnyDescOfType(api::ResourceContext::GetType<T>());
		result = &(_context.Fetch<T>(object_desc.unique_id));
	}
	else
	{
		result = &(_context.Fetch<T>(_object_id));
	}
	return *result;
}

std::tuple<raytracer::Camera*, raytracer::Film*, raytracer::Sampler*>
IntegratorCommonMake(api::ResourceContext &_context, api::ParamSet const &_params)
{
	std::string const camera_id = _params.FindString("camera_id", "");
	std::string const film_id = _params.FindString("film_id", "");
	std::string const sampler_id = _params.FindString("sampler_id", "");
	raytracer::Camera &camera =
		FetchForIDOrAny<raytracer::Camera>(camera_id, _context);
	raytracer::Film &film =
		FetchForIDOrAny<raytracer::Film>(film_id, _context);
	raytracer::Sampler &sampler =
		FetchForIDOrAny<raytracer::Sampler>(sampler_id, _context);
	return std::make_tuple(&camera, &film, &sampler);
}
}

raytracer::Integrator*
MakeNormalIntegrator(api::ResourceContext &_context, api::ParamSet const &_params)
{
	auto integrator_base_params = IntegratorCommonMake(_context, _params);
	bool const	remap = _params.FindBool("remap", false);
	bool const	absolute = _params.FindBool("absolute", false);
	raytracer::Integrator *const normal_integrator = new (_context.mem_region())
		raytracer::NormalIntegrator{ *std::get<0>(integrator_base_params),
									 *std::get<1>(integrator_base_params),
									 *std::get<2>(integrator_base_params),
									 remap, absolute };
	return normal_integrator;
}

raytracer::Integrator*
MakeAOIntegrator(api::ResourceContext &_context, api::ParamSet const &_params)
{
	auto integrator_base_params = IntegratorCommonMake(_context, _params);
	uint64_t const	sample_count = _params.FindUint("sample_count", 1u);
	bool const	shading_geometry = _params.FindBool("shading_geometry", false);
	raytracer::Integrator *const ao_integrator = new (_context.mem_region())
		raytracer::AOIntegrator{ *std::get<0>(integrator_base_params),
								 *std::get<1>(integrator_base_params),
								 *std::get<2>(integrator_base_params),
								 sample_count, shading_geometry };
	return ao_integrator;
}

raytracer::Integrator*
MakeDirectLightingIntegrator(api::ResourceContext &_context, api::ParamSet const &_params)
{
	auto integrator_base_params = IntegratorCommonMake(_context, _params);
	uint64_t const			shadow_ray_count = _params.FindUint("shadow_rays", 1u);
	raytracer::Integrator *const direct_lighting_integrator = new (_context.mem_region())
		raytracer::DirectLightingIntegrator{ *std::get<0>(integrator_base_params),
											 *std::get<1>(integrator_base_params),
											 *std::get<2>(integrator_base_params),
											 shadow_ray_count };
	return direct_lighting_integrator;
}


raytracer::Light*
MakeAreaLight(api::ResourceContext &_context, api::ParamSet const &_params)
{
	maths::Transform const	&world_transform = _params.FindTransform("world_transform", maths::Transform::Identity());
	maths::Vec3f const		&emission_color =
		_params.FindFloat<3>("emission_color", maths::Vec3f{1._d, 0._d, 1._d});
	maths::Decimal const	intensity =
		_params.FindFloat("intensity", 1._d);
	std::string const 		&shape_id = _params.FindString("shape_id", "");
	if (shape_id == "")
	{
		LOG_ERROR(tools::kChannelParsing, "An area light descriptor fails to provide a shape_id");
	}
	raytracer::Shape const	&shape = _context.Fetch<raytracer::Shape>(shape_id);
	raytracer::Light *const area_light = new (_context.mem_region())
		raytracer::AreaLight{ world_transform, emission_color * intensity, shape };
	_context.FlagLightShape(shape);
	return area_light;
}

} // namespace api

#include "api/factory_functions.h"

#include <random>

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


namespace api {

raytracer::Camera *MakeCamera(api::ResourceContext &_context, api::ParamSet const &_params);

raytracer::Film *MakeFilm(api::ResourceContext &_context, api::ParamSet const &_params);

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


raytracer::Shape*
MakeSphere(api::ResourceContext &_context, api::ParamSet const &_params)
{
	// TODO: the identity matrix shouldn't be looked up,
	// maths::Transform should provide a static instance
	maths::Transform const	&identity = _context.transform_cache().Lookup(maths::Transform());
	maths::Transform const	&world_transform = _params.FindTransform("world_transform", identity);
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
	maths::Transform const	&identity = _context.transform_cache().Lookup(maths::Transform());
	maths::Transform const	&world_transform = _params.FindTransform("world_transform", identity);
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
		//
		if (boost::filesystem::exists(path))
		{
			result = raytracer::ReadTriangleMeshFromFile(path.generic_string(),
														 world_transform,
														 flip_normals,
														 _context.mem_region());
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


raytracer::Integrator*
MakeNormalIntegrator(api::ResourceContext &_context, api::ParamSet const &_params)
{
	// REFACTOR: Integrator should lookup its dependecies upon construction
	bool const	remap = _params.FindBool("remap", false);
	bool const	absolute = _params.FindBool("absolute", false);
	raytracer::Integrator *const normal_integrator = new (_context.mem_region())
		raytracer::NormalIntegrator{ remap, absolute };
	return normal_integrator;
}
raytracer::Integrator*
MakeAOIntegrator(api::ResourceContext &_context, api::ParamSet const &_params)
{
	uint64_t const	sample_count = _params.FindUint("sample_count", 1u);
	bool const	shading_geometry = _params.FindBool("shading_geometry", false);
	raytracer::Integrator *const ao_integrator = new (_context.mem_region())
		raytracer::AOIntegrator{ sample_count, shading_geometry };
	return ao_integrator;
}
raytracer::Integrator*
MakeDirectLightingIntegrator(api::ResourceContext &_context, api::ParamSet const &_params)
{
	raytracer::Integrator *const direct_lighting_integrator = new (_context.mem_region())
		raytracer::DirectLightingIntegrator{};
	return direct_lighting_integrator;
}


raytracer::Light*
MakeAreaLight(api::ResourceContext &_context, api::ParamSet const &_params)
{
	maths::Transform const	&identity = _context.transform_cache().Lookup(maths::Transform());
	maths::Transform const	&world_transform = _params.FindTransform("world_transform", identity);
	maths::Vec3f const		&emission_color =
		_params.FindFloat<3>("emission_color", maths::Vec3f{1._d, 0._d, 1._d});
	std::string const &shape_id = _params.FindString("shape_id", "");
	if (shape_id == "")
	{
		LOG_ERROR(tools::kChannelParsing, "An area light descriptor fails to provide a shape_id");
	}
	raytracer::Shape const	&shape = _context.Fetch<raytracer::Shape>(shape_id);
	raytracer::Light *const area_light = new (_context.mem_region())
		raytracer::AreaLight{ world_transform, emission_color, shape };
	_context.FlagLightShape(shape);
	return area_light;
}
	
} // namespace api


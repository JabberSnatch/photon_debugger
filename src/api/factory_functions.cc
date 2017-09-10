#include "api/factory_functions.h"

#include <random>

#include "boost/numeric/conversion/cast.hpp"

#include "api/param_set.h"
#include "api/render_context.h"
#include "maths/transform.h"
#include "raytracer/camera.h"
#include "raytracer/film.h"
#include "raytracer/sampler.h"
#include "raytracer/samplers/random_sampler.h"
#include "raytracer/shape.h"
#include "raytracer/shapes/sphere.h"


namespace api {


raytracer::Camera*
MakeCamera(api::RenderContext &_context, api::ParamSet const &_params)
{
	maths::Point3f const	position = static_cast<maths::Point3f>(
		_params.FindFloat<3>("position", { 0._d, 0._d, 0._d })
		);
	maths::Point3f const	lookat = static_cast<maths::Point3f>(
		_params.FindFloat<3>("lookat", maths::Vec3f(position) + maths::Vec3f{ 0._d, 1._d, 0._d })
		);
	maths::Vec3f const		up = _params.FindFloat<3>("up", { 0._d, 0._d, 1._d });
	maths::Decimal const	fov = _params.FindFloat("fov", 60._d);

	return new (_context.mem_region()) raytracer::Camera{
		position, lookat, up, fov, *_context.film() };
}

raytracer::Film*
MakeFilm(api::RenderContext &_context, api::ParamSet const &_params)
{
	maths::Vec2i const		resolution = _params.FindInt<2>("resolution", { 800, 600 });
	maths::Decimal const	side = _params.FindFloat("side", .036_d);

	return new (_context.mem_region()) raytracer::Film{ resolution.x, resolution.y, side };
}


std::vector<raytracer::Shape*>
MakeSphere(api::RenderContext &_context, api::ParamSet const &_params)
{
	maths::Transform const	&identity = _context.transform_cache().Lookup(maths::Transform());
	maths::Transform const	&world_transform = _params.FindTransform("world_transform", identity);
	bool const				flip_normals = _params.FindBool("flip_normals", false);
	maths::Decimal const	radius = _params.FindFloat("radius", 1._d);
	maths::Decimal const	z_min = _params.FindFloat("z_min", -radius);
	maths::Decimal const	z_max = _params.FindFloat("z_max", radius);
	maths::Decimal const	phi_max = _params.FindFloat("phi_max", 360._d);
	raytracer::Shape *const sphere_shape = new (_context.mem_region())
		raytracer::Sphere{ world_transform, flip_normals, radius, z_min, z_max, phi_max };
	return std::vector<raytracer::Shape*>{ sphere_shape };
}


raytracer::Sampler*
MakeRandomSampler(api::RenderContext &_context, api::ParamSet const &_params)
{
	uint64_t const	seed = _params.FindUint("seed", std::random_device()());
	uint64_t const	samples_per_pixel = _params.FindUint("sample_count", 1);
	uint64_t const	dimensions_per_sample = _params.FindUint("dimension_count", 5);
	raytracer::Sampler *const random_sampler = new (_context.mem_region())
		raytracer::RandomSampler{ seed, samples_per_pixel, dimensions_per_sample };
	return random_sampler;
}


ShapeCallbackContainer_t const &
shape_callbacks()
{
	static ShapeCallbackContainer_t const callbacks {
			{ "sphere", &MakeSphere },
	};
	return callbacks;
}

MakeShapeCallback_t const &
LookupShapeFunc(std::string const &_id)
{
	ShapeCallbackContainer_t const &callbacks = shape_callbacks();
	ShapeCallbackContainer_t::const_iterator const msccit = callbacks.find(_id);
	YS_ASSERT(msccit != callbacks.cend());
	return msccit->second;
}


SamplerCallbackContainer_t const &
sampler_callbacks()
{
	static SamplerCallbackContainer_t const callbacks{
		{ "random_sampler", &MakeRandomSampler },
	};
	return callbacks;
}

MakeSamplerCallback_t const &
LookupSamplerFunc(std::string const &_id)
{
	SamplerCallbackContainer_t const &callbacks = sampler_callbacks();
	SamplerCallbackContainer_t::const_iterator const msccit = callbacks.find(_id);
	YS_ASSERT(msccit != callbacks.cend());
	return msccit->second;
}


} // namespace api


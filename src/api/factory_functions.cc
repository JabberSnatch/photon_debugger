#include "api/factory_functions.h"

#include "api/param_set.h"
#include "raytracer/camera.h"
#include "raytracer/film.h"
#include "raytracer/shape.h"
#include "raytracer/shapes/sphere.h"
#include "raytracer/render_context.h"


namespace api {


raytracer::Camera*
MakeCamera(raytracer::RenderContext &_context, api::ParamSet const &_params)
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
MakeFilm(raytracer::RenderContext &_context, api::ParamSet const &_params)
{
	maths::Vec2i32 const	resolution = _params.FindInt<2>("resolution", { 800, 600 });
	maths::Decimal const	side = _params.FindFloat("side", .036_d);

	return new (_context.mem_region()) raytracer::Film{ resolution.x, resolution.y, side };
}


std::vector<raytracer::Shape*>
MakeSphere(raytracer::RenderContext &_context,
		   maths::Transform const &_t, bool _flip_normals,
		   api::ParamSet const &_params)
{
	bool const				flip_normals = _params.FindBool("flip_normals", false);
	maths::Decimal const	radius = _params.FindFloat("radius", 1._d);
	maths::Decimal const	z_min = _params.FindFloat("z_min", -radius);
	maths::Decimal const	z_max = _params.FindFloat("z_max", radius);
	maths::Decimal const	phi_max = _params.FindFloat("phi_max", 360._d);

	raytracer::Shape *const sphere_shape = new (_context.mem_region())
		raytracer::Sphere{ _t, _flip_normals, radius, z_min, z_max, phi_max };

	return std::vector<raytracer::Shape*>{ sphere_shape };
}


} // namespace api


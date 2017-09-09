#pragma once
#ifndef __YS_FACTORY_FUNCTIONS_HPP__
#define __YS_FACTORY_FUNCTIONS_HPP__

#include <vector>


namespace maths {
class Transform;
} // namespace maths

namespace raytracer {

class RenderContext;
class Camera;
class Film;
class Shape;

} // namespace raytracer


namespace api {

class ParamSet;

raytracer::Camera *MakeCamera(
	raytracer::RenderContext &_context, api::ParamSet const &_params);
raytracer::Film *MakeFilm(
	raytracer::RenderContext &_context, api::ParamSet const &_params);
std::vector<raytracer::Shape*>
MakeSphere(raytracer::RenderContext &_context,
		   maths::Transform const &_t, bool _flip_normals,
		   api::ParamSet const &_params);


} // namespace api


#endif // __YS_FACTORY_FUNCTIONS_HPP__

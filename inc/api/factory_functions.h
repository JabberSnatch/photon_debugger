#pragma once
#ifndef __YS_FACTORY_FUNCTIONS_HPP__
#define __YS_FACTORY_FUNCTIONS_HPP__

#include <vector>
#include <unordered_map>
#include <functional>


namespace maths {
class Transform;
} // namespace maths

namespace raytracer {

class Camera;
class Film;
class Shape;
class Sampler;

} // namespace raytracer


namespace api {

class ParamSet;
class RenderContext;

using MakeShapeCallback_t = std::function<
	std::vector<raytracer::Shape*>(api::RenderContext &_context, ParamSet const &_params)
>;
using ShapeCallbackContainer_t = std::unordered_map<std::string, MakeShapeCallback_t>;

using MakeSamplerCallback_t = std::function<
	raytracer::Sampler*(api::RenderContext &_context, ParamSet const &_params)
>;
using SamplerCallbackContainer_t = std::unordered_map<std::string, MakeSamplerCallback_t>;


raytracer::Camera *MakeCamera(api::RenderContext &_context, api::ParamSet const &_params);
raytracer::Film *MakeFilm(api::RenderContext &_context, api::ParamSet const &_params);

std::vector<raytracer::Shape*> MakeSphere(api::RenderContext &_context,
										  api::ParamSet const &_params);

raytracer::Sampler* MakeRandomSampler(api::RenderContext &_context,
									  api::ParamSet const &_params);


ShapeCallbackContainer_t const &shape_callbacks();
MakeShapeCallback_t const &LookupShapeFunc(std::string const &_id);

SamplerCallbackContainer_t const &sampler_callbacks();
MakeSamplerCallback_t const &LookupSamplerFunc(std::string const &_id);

} // namespace api


#endif // __YS_FACTORY_FUNCTIONS_HPP__

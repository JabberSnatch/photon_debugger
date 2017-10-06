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
class Integrator;
class Shape;
class Sampler;

} // namespace raytracer


namespace api {

class ParamSet;
class RenderContext;


template <typename ObjectType>
using MakeObjectCallback_t =
	std::function<ObjectType*(api::RenderContext &_context, ParamSet const &_params)>;
template <typename CallbackType>
using NamedCallbackContainer_t = std::unordered_map<std::string, CallbackType>;


using MakeShapeCallback_t = std::function<
	std::vector<raytracer::Shape*>(api::RenderContext &_context, ParamSet const &_params)
>;
using MakeSamplerCallback_t = MakeObjectCallback_t<raytracer::Sampler>;
using MakeIntegratorCallback_t = MakeObjectCallback_t<raytracer::Integrator>;
using ShapeCallbackContainer_t = NamedCallbackContainer_t<MakeShapeCallback_t>;
using SamplerCallbackContainer_t = NamedCallbackContainer_t<MakeSamplerCallback_t>;
using IntegratorCallbackContainer_t = NamedCallbackContainer_t<MakeIntegratorCallback_t>;


ShapeCallbackContainer_t const &shape_callbacks();
SamplerCallbackContainer_t const &sampler_callbacks();
IntegratorCallbackContainer_t const &integrator_callbacks();
MakeShapeCallback_t const &LookupShapeFunc(std::string const &_id);
MakeSamplerCallback_t const &LookupSamplerFunc(std::string const &_id);
MakeIntegratorCallback_t const &LookupIntegratorFunc(std::string const &_id);


raytracer::Camera* MakeCamera(api::RenderContext &_context, api::ParamSet const &_params);
raytracer::Film* MakeFilm(api::RenderContext &_context, api::ParamSet const &_params);


} // namespace api


#endif // __YS_FACTORY_FUNCTIONS_HPP__

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
class Light;
class Shape;
class Sampler;
class TriangleMeshRawData;

} // namespace raytracer


namespace api {

class ParamSet;
class ResourceContext;


template <typename ObjectType>
using MakeObjectCallback_t =
	std::function<ObjectType*(api::ResourceContext &_context, ParamSet const &_params)>;
template <typename CallbackType>
using NamedCallbackContainer_t = std::unordered_map<std::string, CallbackType>;


using MakeShapeCallback_t = MakeObjectCallback_t<raytracer::Shape>;
using MakeSamplerCallback_t = MakeObjectCallback_t<raytracer::Sampler>;
using MakeIntegratorCallback_t = MakeObjectCallback_t<raytracer::Integrator>;
using MakeLightCallback_t = MakeObjectCallback_t<raytracer::Light>;
using ShapeCallbackContainer_t = NamedCallbackContainer_t<MakeShapeCallback_t>;
using SamplerCallbackContainer_t = NamedCallbackContainer_t<MakeSamplerCallback_t>;
using IntegratorCallbackContainer_t = NamedCallbackContainer_t<MakeIntegratorCallback_t>;
using LightCallbackContainer_t = NamedCallbackContainer_t<MakeLightCallback_t>;


raytracer::Camera* MakeCamera(api::ResourceContext &_context, api::ParamSet const &_params);
raytracer::Film* MakeFilm(api::ResourceContext &_context, api::ParamSet const &_params);

raytracer::TriangleMeshRawData* MakeTriangleMeshRawData(api::ResourceContext &_context,
														api::ParamSet const &_params);


ShapeCallbackContainer_t const &shape_callbacks();
SamplerCallbackContainer_t const &sampler_callbacks();
IntegratorCallbackContainer_t const &integrator_callbacks();
LightCallbackContainer_t const &light_callbacks();
MakeShapeCallback_t const &LookupShapeFunc(std::string const &_id);
MakeSamplerCallback_t const &LookupSamplerFunc(std::string const &_id);
MakeIntegratorCallback_t const &LookupIntegratorFunc(std::string const &_id);
MakeLightCallback_t const &LookupLightFunc(std::string const &_id);

} // namespace api


#endif // __YS_FACTORY_FUNCTIONS_HPP__

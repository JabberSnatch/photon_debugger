#pragma once
#ifndef __YS_RENDER_CONTEXT_HPP__
#define __YS_RENDER_CONTEXT_HPP__

#include <vector>

#include "api/transform_cache.h"
#include "core/memory_region.h"
#include "core/noncopyable.h"
#include "core/nonmovable.h"
#include "raytracer/camera.h"
#include "raytracer/light.h"
#include "raytracer/primitive.h"
#include "raytracer/sampler.h"
#include "raytracer/integrator.h"
#include "raytracer/bvh_accelerator.h"

namespace api {


class RenderContext final
{
public:
	using PrimitiveContainer_t = std::vector<raytracer::Primitive const*>;
	using LightContainer_t = std::vector<raytracer::Light const*>;
public:
	RenderContext();
	RenderContext(raytracer::Integrator &_integrator,
				  PrimitiveContainer_t &_primitives,
				  LightContainer_t &_lights);
	void	Clear();
	void	AddPrimitive(raytracer::Primitive *_prim);
public:
	bool	GoodForRender() const;
	void	RenderAndWrite(std::string const &_path);
private:
	raytracer::Integrator	*integrator_ = nullptr;
	PrimitiveContainer_t	primitives_{};
	LightContainer_t		lights_{};
};


} // namespace raytracer

#endif // __YS_RENDER_CONTEXT_HPP__

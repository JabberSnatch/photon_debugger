#pragma once
#ifndef __YS_RENDER_CONTEXT_HPP__
#define __YS_RENDER_CONTEXT_HPP__

#include <vector>

#include "api/transform_cache.h"
#include "core/memory_region.h"
#include "core/noncopyable.h"
#include "core/nonmovable.h"
#include "raytracer/camera.h"
#include "raytracer/primitive.h"
#include "raytracer/sampler.h"
#include "raytracer/integrator.h"
#include "raytracer/bvh_accelerator.h"

namespace api {


class RenderContext final :
	private core::noncopyable,
	private core::nonmovable
{
public:
	RenderContext();
	void	Clear();
	void	SetFilm(raytracer::Film *_f);
	void	SetCamera(raytracer::Camera *_c);
	void	SetSampler(raytracer::Sampler *_s);
	void	SetIntegrator(raytracer::Integrator *_i);
	void	AddPrimitive(raytracer::Primitive *_prim);
public:
	bool	GoodForRender() const;
	void	RenderAndWrite(std::string const &_path);
public:
	raytracer::Film		&film() { return *film_; }
	TransformCache		&transform_cache() { return transform_cache_; }
	core::MemoryRegion	&mem_region() { return mem_region_; }
private:
	core::MemoryRegion					mem_region_;
	TransformCache						transform_cache_;
	raytracer::Camera					*camera_ = nullptr;
	raytracer::Film						*film_ = nullptr;
	raytracer::Sampler					*sampler_ = nullptr;
	raytracer::Integrator				*integrator_ = nullptr;
	std::vector<raytracer::Primitive*>	primitives_;
	//
public:
	std::string workdir; // To be removed, at some point in the future
};


} // namespace raytracer

#endif // __YS_RENDER_CONTEXT_HPP__

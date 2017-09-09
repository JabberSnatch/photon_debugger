#pragma once
#ifndef __YS_RENDER_CONTEXT_HPP__
#define __YS_RENDER_CONTEXT_HPP__

#include <vector>

#include "raytracer/raytracer.h"
#include "core/memory_region.h"

#include "raytracer/camera.h"
#include "raytracer/primitive.h"
#include "raytracer/integrator.h"
#include "raytracer/samplers/random_sampler.h"

namespace raytracer {


class RenderContext
{
public:
	RenderContext() = default;
	~RenderContext() = default;
	RenderContext(RenderContext const &) = delete;
	RenderContext(RenderContext &&) = delete;
	RenderContext &operator=(RenderContext const &) = delete;
	RenderContext &operator=(RenderContext &&) = delete;

	void	Clear()
	{
		film_ = nullptr;
		camera_ = nullptr;
		primitives_.clear();
		mem_region_.Clear();
	}
	void	SetFilm(Film *_f)
	{
		film_ = _f;
		integrator_.SetFilm(film_);
	}
	void	SetCamera(Camera *_c)
	{
		camera_ = _c;
		integrator_.SetCamera(camera_);
	}
	bool	GoodForRender() const
	{
		return (camera_ && film_);
	}
	void	RenderAndWrite(std::string const &_path)
	{
		core::RNG sampler_rng(14587126349871134129u);
		Sampler *sampler = new (mem_region_) RandomSampler(sampler_rng, 16u, 2u);
		integrator_.SetSampler(sampler);

		integrator_.Integrate(primitives_, 0._d);
		//camera_->Expose(primitives_, 0._d);
		camera_->WriteToFile(_path);
	}

	void AddPrimitive(Primitive *_prim)
	{
		primitives_.push_back(_prim);
	}

	raytracer::Film		*film() { return film_; }
	core::MemoryRegion<> &mem_region() { return mem_region_; }

private:
	raytracer::Camera			*camera_ = nullptr;
	raytracer::Film				*film_ = nullptr;
	raytracer::Integrator		integrator_;
	std::vector<Primitive*>		primitives_;

	core::MemoryRegion<>		mem_region_;
};


} // namespace raytracer

#endif // __YS_RENDER_CONTEXT_HPP__

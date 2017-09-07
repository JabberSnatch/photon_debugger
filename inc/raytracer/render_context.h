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
		//if (camera_)
		//	camera_->SetFilm(film_);
	}
	void	SetCamera(Camera *_c)
	{
		camera_ = _c;
		integrator_.SetCamera(camera_);
		//if (film_)
		//	camera_->SetFilm(film_);
	}
	bool	GoodForRender() const
	{
		return (camera_ && film_);
	}
	void	RenderAndWrite(std::string const &_path)
	{
		Sampler *sampler = mem_region_.Alloc<RandomSampler>();
		core::RNG sampler_rng(14587126349871134129u);
		new (sampler) RandomSampler(sampler_rng, 16u, 2u);
		integrator_.SetSampler(sampler);

		integrator_.Integrate(primitives_, 0._d);
		//camera_->Expose(primitives_, 0._d);
		camera_->WriteToFile(_path);
	}

	template <typename T>
	Shape	*AllocShapes(size_t _count = 1)
	{
		return static_cast<Shape*>(mem_region_.Alloc<T>(_count));
	}
	template <typename T>
	Primitive	*AllocPrimitives(size_t _count = 1)
	{
		Primitive	*first = static_cast<Primitive*>(mem_region_.Alloc<T>(_count));
		//primitives_.insert(primitives_.cend(), first, first + _count);
		for (int i = 0; i < _count; ++i)
			primitives_.emplace_back(first + i);
		return	first;
	}
	Film	*AllocFilm()
	{
		return mem_region_.Alloc<Film>();
	}
	Camera	*AllocCamera()
	{
		return mem_region_.Alloc<Camera>();
	}

	raytracer::Film		*film() { return film_; }

private:
	raytracer::Camera			*camera_ = nullptr;
	raytracer::Film				*film_ = nullptr;
	raytracer::Integrator		integrator_;
	std::vector<Primitive*>		primitives_;

	core::MemoryRegion<>		mem_region_;
};


} // namespace raytracer

#endif // __YS_RENDER_CONTEXT_HPP__

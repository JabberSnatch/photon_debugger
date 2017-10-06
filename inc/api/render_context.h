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
	RenderContext() :
		mem_region_{},
		transform_cache_{},
		camera_{ nullptr },
		film_{ nullptr },
		sampler_{ nullptr },
		integrator_{ nullptr },
		primitives_{}
	{}
	void	Clear()
	{
		mem_region_.Clear();
		film_ = nullptr;
		camera_ = nullptr;
		sampler_ = nullptr;
		integrator_ = nullptr;
		primitives_.clear();
	}
	void	SetFilm(raytracer::Film *_f)
	{
		film_ = _f;
		integrator_->SetFilm(film_);
	}
	void	SetCamera(raytracer::Camera *_c)
	{
		camera_ = _c;
		integrator_->SetCamera(camera_);
	}
	void	SetSampler(raytracer::Sampler *_s)
	{
		sampler_ = _s;
		integrator_->SetSampler(sampler_);
	}
	void	SetIntegrator(raytracer::Integrator *_i)
	{
		integrator_ = _i;
	}
	bool	GoodForRender() const
	{
		return (camera_ && film_ && sampler_);
	}
	void	RenderAndWrite(std::string const &_path)
	{
		constexpr uint32_t bvh_node_max_size = 20;
		integrator_->Prepare();
		if (primitives_.size() > bvh_node_max_size)
		{
			LOG_INFO(tools::kChannelGeneral, "Primitive count exceeded threshold, building BVH");
			raytracer::Primitive *const bvh =
				new (mem_region()) raytracer::BvhAccelerator(primitives_, bvh_node_max_size);
			integrator_->Integrate({ bvh }, 0._d);
		}
		else
		{
			integrator_->Integrate(primitives_, 0._d);
		}
		camera_->WriteToFile(_path);
	}

	void AddPrimitive(raytracer::Primitive *_prim)
	{
		primitives_.push_back(_prim);
	}

	raytracer::Film		*film() { return film_; }
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
};


} // namespace raytracer

#endif // __YS_RENDER_CONTEXT_HPP__

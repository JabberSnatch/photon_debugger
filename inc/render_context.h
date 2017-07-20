#pragma once
#ifndef __YS_RENDER_CONTEXT_HPP__
#define __YS_RENDER_CONTEXT_HPP__

#include <vector>

#include "raytracer.h"
#include "memory_region.h"

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

private:
	raytracer::Camera			*camera_ = nullptr;
	raytracer::Film				*film_ = nullptr;
	std::vector<Primitive *>	primitives_;

	core::MemoryRegion<>		mem_region_;
};


} // namespace raytracer

#endif // __YS_RENDER_CONTEXT_HPP__

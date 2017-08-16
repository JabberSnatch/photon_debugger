#pragma once
#ifndef __YS_RENDER_CONTEXT_HPP__
#define __YS_RENDER_CONTEXT_HPP__

#include <vector>

#include "raytracer.h"
#include "memory_region.h"

#include "camera.h"
#include "primitive.h"

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
		//if (camera_)
		//	camera_->SetFilm(film_);
	}
	void	SetCamera(Camera *_c)
	{
		camera_ = _c;
		//if (film_)
		//	camera_->SetFilm(film_);
	}
	bool	GoodForRender() const
	{
		return (camera_ && film_);
	}
	void	RenderAndWrite(std::string const &_path) const
	{
		camera_->Expose(primitives_, 0._d);
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
	std::vector<Primitive*>		primitives_;

	core::MemoryRegion<>		mem_region_;
};


} // namespace raytracer

#endif // __YS_RENDER_CONTEXT_HPP__

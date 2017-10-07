#include "api/render_context.h"

namespace api
{


RenderContext::RenderContext() :
	mem_region_{},
	transform_cache_{},
	camera_{ nullptr },
	film_{ nullptr },
	sampler_{ nullptr },
	integrator_{ nullptr },
	primitives_{}
{}


void
RenderContext::Clear()
{
	mem_region_.Clear();
	film_ = nullptr;
	camera_ = nullptr;
	sampler_ = nullptr;
	integrator_ = nullptr;
	primitives_.clear();
}


void
RenderContext::SetFilm(raytracer::Film *_f)
{
	film_ = _f;
}


void
RenderContext::SetCamera(raytracer::Camera *_c)
{
	camera_ = _c;
}


void
RenderContext::SetSampler(raytracer::Sampler *_s)
{
	sampler_ = _s;
}


void
RenderContext::SetIntegrator(raytracer::Integrator *_i)
{
	integrator_ = _i;
}


void
RenderContext::AddPrimitive(raytracer::Primitive *_prim)
{
	primitives_.push_back(_prim);
}


bool
RenderContext::GoodForRender() const
{
	return (camera_ && film_ && sampler_ && integrator_);
}


void
RenderContext::RenderAndWrite(std::string const &_path)
{
	integrator_->SetFilm(film_);
	integrator_->SetCamera(camera_);
	integrator_->SetSampler(sampler_);
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


} // namespace api

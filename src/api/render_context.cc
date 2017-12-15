#include "api/render_context.h"

namespace api
{


RenderContext::RenderContext() :
	camera_{ nullptr },
	sampler_{ nullptr },
	integrator_{ nullptr },
	primitives_{}
{}


RenderContext::RenderContext(raytracer::Camera &_camera,
							 raytracer::Sampler &_sampler,
							 raytracer::Integrator &_integrator,
							 PrimitiveContainer_t &_primitives) :
	camera_{ &_camera },
	sampler_{ &_sampler },
	integrator_{ &_integrator },
	primitives_{ _primitives }
{}


void
RenderContext::Clear()
{
	camera_ = nullptr;
	sampler_ = nullptr;
	integrator_ = nullptr;
	primitives_.clear();
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
	return (camera_ && sampler_ && integrator_);
}


void
RenderContext::RenderAndWrite(std::string const &_path)
{
	// REFACTOR: Still some work to do here, as soon as Integrator is constructed from its dependencies
	integrator_->SetFilm(&(camera_->film()));
	integrator_->SetCamera(camera_);
	integrator_->SetSampler(sampler_);
	integrator_->Prepare();
	integrator_->Integrate({ primitives_, {} }, 0._d);
	camera_->WriteToFile(_path);
}


} // namespace api

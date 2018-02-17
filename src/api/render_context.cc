#include "api/render_context.h"

namespace api
{


RenderContext::RenderContext() :
	integrator_{ nullptr },
	primitives_{},
	lights_{}
{}


RenderContext::RenderContext(raytracer::Integrator &_integrator,
							 PrimitiveContainer_t &_primitives,
							 LightContainer_t &_lights) :
	integrator_{ &_integrator },
	primitives_{ _primitives },
	lights_{ _lights }
{}


void
RenderContext::Clear()
{
	integrator_ = nullptr;
	primitives_.clear();
	lights_.clear();
}


void
RenderContext::AddPrimitive(raytracer::Primitive *_prim)
{
	primitives_.push_back(_prim);
}


bool
RenderContext::GoodForRender() const
{
	return (integrator_);
}


void
RenderContext::RenderAndWrite(std::string const &_path)
{
	integrator_->Prepare(primitives_, lights_);
	integrator_->Integrate({ primitives_, lights_ }, 0._d);
	integrator_->camera().WriteToFile(_path);
}


} // namespace api

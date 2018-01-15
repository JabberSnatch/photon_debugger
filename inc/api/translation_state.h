#pragma once
#ifndef __YS_TRANSLATION_STATE_HPP__
#define __YS_TRANSLATION_STATE_HPP__

#include <functional>
#include <string>

#include "api/param_set.h"
#include "api/render_context.h"
#include "api/resource_context.h"
#include "api/transform_cache.h"
#include "core/noncopyable.h"
#include "core/nonmovable.h"
#include "maths/maths.h"

namespace raytracer {
class Shape;
} // namespace api


namespace api {

// NOTE: Implementation doesn't enforce syntax compliance, be very careful about using it on its own.
class TranslationState final :
	private core::noncopyable,
	private core::nonmovable
{
private:
	using TransformStack_t = std::vector<maths::Transform>;
	using ResourceCounterContainer_t =
		std::array<uint32_t, static_cast<size_t>(ResourceContext::ObjectType::kCount)>;
public:
	TranslationState();
	void	SceneBegin();
	void	SceneEnd();
	void	ScopeBegin();
	void	ScopeEnd();
	void	Workdir(std::string const &_absolute_path);
	void	Output(std::string const &_relative_path);
	void	ObjectId(std::string const &_object_id);
	void	Film();
	void	Camera();
	void	Shape(std::string const &_type);
	void	Light(std::string const &_type);
	void	Sampler(std::string const &_type);
	void	Integrator(std::string const &_type);
	void	Identity();
	void	Translate(maths::Vec3f const &_t);
	void	Rotate(maths::Decimal _angle, maths::Vec3f const &_axis);
	void	Scale(maths::Decimal _x, maths::Decimal _y, maths::Decimal _z);
	ParamSet	&param_set() { return *parameters_; }
public:
	api::RenderContext	&render_context() { return render_context_; }
	std::string			&output_path() { return output_path_; }
private:
	void	SceneSetup_();
	void	PushObjectDesc_(ResourceContext::ObjectType const _type,
							std::string const &_subtype_id);
private:
	std::string				workdir_;
	std::string				output_path_;
	api::RenderContext		render_context_;
	api::ResourceContext	resource_context_;
	uint32_t				scope_depth_;
	TransformStack_t		transform_stack_;
	std::string				cached_object_id_;
	ParamSet				*parameters_;
public:
	void	ResetResourceCounters();
private:
	std::string					MakeUniqueResourceID_(ResourceContext::ObjectType const _type) const;
	ResourceCounterContainer_t	resource_counters_{};
};


} // namespace api


#endif // __YS_TRANSLATION_STATE_HPP__

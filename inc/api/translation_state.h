#pragma once
#ifndef __YS_TRANSLATION_STATE_HPP__
#define __YS_TRANSLATION_STATE_HPP__

#include <functional>
#include <string>

#include "api/param_set.h"
#include "api/render_context.h"
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
public:
	TranslationState();
	void	Workdir(std::string const &_absolute_path);
	void	Output(std::string const &_relative_path);
	void	Film();
	void	Camera();
	void	Shape(std::string const &_type);
	void	Identity();
	void	Translate(maths::Vec3f const &_t);
	void	Rotate(maths::Decimal _angle, maths::Vec3f const &_axis);
	void	Scale(maths::Decimal _x, maths::Decimal _y, maths::Decimal _z);
	void	ScopeBegin();
	void	ScopeEnd();
	void	SceneBegin();
	void	SceneEnd();
	ParamSet	&param_set() { return parameters_; }
private:
	std::string			workdir_;
	std::string			output_path_;
	api::RenderContext	render_context_;
	ParamSet			parameters_;
	uint32_t			scope_depth_;
	TransformStack_t	transform_stack_;
};


} // namespace api


#endif // __YS_TRANSLATION_STATE_HPP__
